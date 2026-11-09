#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 efic_hpge.py - Curva de eficiencia de pico absoluta para detector HPGe ORTEC
                a partir de espectros ASCII .Spe (Maestro / GammaVision)
===============================================================================

O que o codigo faz
------------------
1.  Le um ou mais arquivos .Spe (formato ASCII da ORTEC) extraindo o espectro,
    o tempo vivo/real, a data-hora da medida e a calibracao em energia.
2.  Para cada linha gama tabelada do nuclideo de calibracao, localiza o pico,
    ajusta uma (ou varias, em caso de multipleto) gaussiana sobre fundo linear
    e calcula a area liquida com sua incerteza.
3.  Corrige a atividade da fonte da data de referencia do certificado para a
    data-hora de cada medida (decaimento exponencial).
4.  Calcula a eficiencia absoluta de pico por linha:

        eps(E) = N_liq / ( A * t_vivo * I_gama )

5.  Opcionalmente combina varios espectros (media ponderada ou aritmetica) e
    ajusta a curva classica log-log:

        ln(eps) = sum_{i=0..g} a_i * [ln(E/E0)]^i

    com pesos 1/u^2, devolvendo coeficientes, matriz de covariancia e banda
    de incerteza do ajuste.
6.  Grava CSV com os pontos, CSV com a curva amostrada, grafico e relatorio.

Uso tipico
----------
    python3 efic_hpge.py \
        --spe inferior.Spe --spe superior.Spe \
        --nuclideo Eu-152 \
        --atividade 4798 --data-ref 2018-06-20 \
        --combinar media \
        --grau 4 \
        --prefixo saida/eficiencia

Argumentos principais
---------------------
    --spe            caminho de um .Spe (repetir a opcao para varios)
    --rotulo         rotulo de exibicao para o .Spe anterior (opcional)
    --atividade      atividade da fonte na data de referencia, em Bq
    --u-atividade    incerteza relativa da atividade, em % (padrao 0)
    --data-ref       data de referencia do certificado (AAAA-MM-DD[ HH:MM])
    --meia-vida      meia-vida em anos (padrao: a do nuclideo na biblioteca)
    --combinar       'media' (aritmetica) | 'ponderada' | 'nenhuma'
    --grau           grau do polinomio em log-log (padrao 4)
    --i-minimo       intensidade minima da linha para entrar no ajuste, em %
    --max-incerteza  incerteza relativa maxima da area, em % (padrao 15)
    --prefixo        prefixo dos arquivos de saida

Convencoes e limitacoes (leia)
------------------------------
*  A eficiencia obtida e a eficiencia ABSOLUTA DE PICO (full-energy peak),
   valida APENAS para a geometria em que a fonte foi medida.
*  NAO ha correcao de coincidencia verdadeira (true coincidence summing).
   Para nuclideos com cascatas intensas, como o Eu-152, em geometria proxima
   isso deprime as areas dos picos em cascata e distorce a curva. Veja a nota
   impressa no relatorio.
*  NAO ha subtracao de espectro de fundo ambiental nem correcao de
   auto-atenuacao na fonte.
*  O tempo morto e tratado pelo uso do TEMPO VIVO, como de praxe.

Dependencias: numpy, scipy, matplotlib.
===============================================================================
"""

from __future__ import annotations

import argparse
import datetime as dt
import math
import os
import re
import sys

import numpy as np
from scipy.optimize import curve_fit

# ----------------------------------------------------------------------------
# 1. BIBLIOTECAS DE NUCLIDEOS
# ----------------------------------------------------------------------------
# Energia em keV, intensidade (probabilidade de emissao por decaimento) em %,
# e incerteza da intensidade em %. Valores de bibliotecas de referencia de uso
# corrente; CONFIRA-OS contra a biblioteca adotada no seu laboratorio
# (DDEP/LNHB, IAEA, NuDat) antes de publicar resultados.
#
# 'grupo' marca linhas que devem ser ajustadas juntas como multipleto: linhas
# com o mesmo rotulo de grupo entram em um unico ajuste multigaussiano.

BIBLIOTECA = {
    "Eu-152": {
        "meia_vida_anos": 13.517,
        "u_meia_vida_anos": 0.009,
        "linhas": [
            # (E_keV,    I_%,    u_I_%,  grupo)
            (121.7817,  28.41,  0.13,   None),
            (244.6974,   7.55,  0.04,   None),
            (295.9387,   0.440, 0.005,  None),
            (344.2785,  26.59,  0.12,   None),
            (367.7891,   0.859, 0.009,  None),
            (411.1165,   2.238, 0.010,  None),
            (443.9650,   2.800, 0.019,  None),
            (488.6792,   0.414, 0.004,  None),
            (564.0210,   0.492, 0.005,  None),
            (586.2648,   0.455, 0.005,  None),
            (656.4880,   0.146, 0.002,  None),
            (688.6700,   0.856, 0.006,  None),
            (778.9045,  12.97,  0.06,   None),
            (867.3800,   4.245, 0.023,  None),
            (919.3370,   0.419, 0.004,  "g919"),
            (926.3170,   0.271, 0.003,  "g919"),
            (1005.272,   0.646, 0.008,  None),
            (1085.837,  10.13,  0.06,   "g1086"),
            (1089.737,   1.734, 0.011,  "g1086"),
            (1112.076,  13.54,  0.06,   None),
            (1212.948,   1.415, 0.009,  None),
            (1299.142,   1.633, 0.010,  None),
            (1408.013,  20.85,  0.09,   None),
            (1457.643,   0.497, 0.004,  None),
            (1528.104,   0.279, 0.003,  None),
        ],
    },
    # Esqueleto para outras fontes de calibracao de uso comum.
    "Ba-133": {
        "meia_vida_anos": 10.551,
        "u_meia_vida_anos": 0.011,
        "linhas": [
            (53.1622,   2.14,  0.03, None),
            (79.6142,   2.65,  0.05, "g80"),
            (80.9979,  32.9,   0.3,  "g80"),
            (160.6120,  0.638, 0.005, None),
            (223.2368,  0.453, 0.004, None),
            (276.3989,  7.16,  0.05, None),
            (302.8508, 18.34,  0.13, None),
            (356.0129, 62.05,  0.19, None),
            (383.8485,  8.94,  0.06, None),
        ],
    },
    "Co-60": {
        "meia_vida_anos": 5.2711,
        "u_meia_vida_anos": 0.0008,
        "linhas": [
            (1173.228, 99.85, 0.03, None),
            (1332.492, 99.9826, 0.0006, None),
        ],
    },
    "Cs-137": {
        "meia_vida_anos": 30.05,
        "u_meia_vida_anos": 0.08,
        "linhas": [
            (661.657, 84.99, 0.20, None),
        ],
    },
}

DIAS_POR_ANO = 365.2422  # ano tropico medio; coerente com meia-vida em anos


# ----------------------------------------------------------------------------
# 2. LEITURA DO FORMATO .Spe
# ----------------------------------------------------------------------------

class Espectro:
    """Espectro lido de um arquivo .Spe da ORTEC."""

    def __init__(self, caminho, rotulo=None):
        self.caminho = caminho
        self.rotulo = rotulo or os.path.basename(caminho)
        self.spec_id = ""
        self.detector = ""
        self.aplicativo = ""
        self.data_medida = None      # datetime do INICIO da contagem
        self.tempo_vivo = None       # s
        self.tempo_real = None       # s
        self.canal_inicial = 0
        self.contagens = None        # np.ndarray de int
        self.cal_energia = None      # (a0, a1, a2, ...) -> E = sum a_i ch^i
        self.unidade_energia = "keV"
        self.cal_forma = None        # (b0, b1, b2) -> FWHM em canais
        self._ler()

    # -- leitura ------------------------------------------------------------
    def _ler(self):
        # Maestro grava em cp1252; alguns arquivos chegam ja transcodificados
        # para UTF-8 (com U+FFFD onde o acento se perdeu). Tenta UTF-8 e cai
        # para cp1252, substituindo o que nao decodificar.
        with open(self.caminho, "rb") as fh:
            cru = fh.read()
        for codec in ("utf-8", "cp1252", "latin-1"):
            try:
                texto = cru.decode(codec)
                break
            except UnicodeDecodeError:
                continue
        else:
            texto = cru.decode("latin-1", "replace")
        texto = texto.replace("�", "?")
        linhas = texto.replace("\r\n", "\n").split("\n")

        # Indexa as secoes: cada chave '$NOME:' aponta para as linhas seguintes
        secoes, chave_atual = {}, None
        for linha in linhas:
            if linha.startswith("$") and linha.rstrip().endswith(":"):
                chave_atual = linha.strip().rstrip(":")
                secoes[chave_atual] = []
            elif chave_atual is not None:
                secoes[chave_atual].append(linha)

        if "$DATA" not in secoes:
            raise ValueError("%s: secao $DATA ausente; nao e um .Spe valido"
                             % self.caminho)

        self.spec_id = (secoes.get("$SPEC_ID") or [""])[0].strip()

        for linha in secoes.get("$SPEC_REM", []):
            if linha.startswith("DETDESC#"):
                self.detector = linha.split("#", 1)[1].strip()
            elif linha.startswith("AP#"):
                self.aplicativo = linha.split("#", 1)[1].strip()

        # $DATE_MEA: MM/DD/AAAA HH:MM:SS
        if secoes.get("$DATE_MEA"):
            txt = secoes["$DATE_MEA"][0].strip()
            for fmt in ("%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M",
                        "%d/%m/%Y %H:%M:%S", "%m/%d/%y %H:%M:%S"):
                try:
                    self.data_medida = dt.datetime.strptime(txt, fmt)
                    break
                except ValueError:
                    continue
            if self.data_medida is None:
                raise ValueError("%s: nao entendi $DATE_MEA = %r"
                                 % (self.caminho, txt))

        # $MEAS_TIM: "<tempo_vivo> <tempo_real>"  (nesta ordem, na ORTEC)
        if secoes.get("$MEAS_TIM"):
            partes = secoes["$MEAS_TIM"][0].split()
            self.tempo_vivo = float(partes[0])
            self.tempo_real = float(partes[1]) if len(partes) > 1 else None

        # $DATA: primeira linha "<canal_inicial> <canal_final>", depois 1/linha
        dados = secoes["$DATA"]
        primeiro, ultimo = (int(x) for x in dados[0].split()[:2])
        self.canal_inicial = primeiro
        n = ultimo - primeiro + 1
        valores = []
        for linha in dados[1:]:
            linha = linha.strip()
            if not linha:
                continue
            valores.extend(int(float(x)) for x in linha.split())
            if len(valores) >= n:
                break
        if len(valores) < n:
            raise ValueError("%s: $DATA tem %d valores, esperados %d"
                             % (self.caminho, len(valores), n))
        self.contagens = np.asarray(valores[:n], dtype=float)

        # Calibracao em energia: $MCA_CAL e mais completa que $ENER_FIT
        if secoes.get("$MCA_CAL"):
            corpo = [l for l in secoes["$MCA_CAL"] if l.strip()]
            if len(corpo) >= 2:
                campos = corpo[1].split()
                coefs, unidade = [], self.unidade_energia
                for c in campos:
                    try:
                        coefs.append(float(c))
                    except ValueError:
                        unidade = c.strip()
                self.cal_energia = tuple(coefs)
                self.unidade_energia = unidade
        if self.cal_energia is None and secoes.get("$ENER_FIT"):
            self.cal_energia = tuple(float(x) for x in
                                     secoes["$ENER_FIT"][0].split())
        if self.cal_energia is None:
            raise ValueError("%s: sem calibracao em energia ($MCA_CAL/"
                             "$ENER_FIT); informe-a manualmente" % self.caminho)

        if secoes.get("$SHAPE_CAL"):
            corpo = [l for l in secoes["$SHAPE_CAL"] if l.strip()]
            if len(corpo) >= 2:
                self.cal_forma = tuple(float(x) for x in corpo[1].split())

    # -- utilidades ---------------------------------------------------------
    @property
    def n_canais(self):
        return len(self.contagens)

    @property
    def tempo_morto_pct(self):
        if not self.tempo_real:
            return float("nan")
        return 100.0 * (self.tempo_real - self.tempo_vivo) / self.tempo_real

    @property
    def data_meia_contagem(self):
        """Instante medio da contagem: o correto para corrigir decaimento."""
        return self.data_medida + dt.timedelta(seconds=0.5 * self.tempo_real)

    def energia(self, canal):
        """Energia a partir do numero de canal, pela calibracao do arquivo."""
        canal = np.asarray(canal, dtype=float)
        return sum(c * canal ** i for i, c in enumerate(self.cal_energia))

    def canal(self, energia_alvo):
        """Inverte a calibracao numericamente (monotonica na faixa util)."""
        chs = np.arange(self.n_canais, dtype=float)
        energias = self.energia(chs)
        return float(np.interp(energia_alvo, energias, chs))

    def fwhm_canais(self, canal):
        """FWHM em canais pela $SHAPE_CAL, se houver; senao estimativa."""
        if self.cal_forma:
            v = sum(c * canal ** i for i, c in enumerate(self.cal_forma))
            if v > 0.5:
                return v
        return 6.0


# ----------------------------------------------------------------------------
# 3. AJUSTE DE PICOS
# ----------------------------------------------------------------------------

def _multigaussiana(x, *p):
    """n gaussianas de sigma comum sobre fundo linear.

    p = [b0, b1, sigma, A1, mu1, A2, mu2, ...]
    A_i e a AREA de cada gaussiana (nao a altura).
    """
    b0, b1, sigma = p[0], p[1], p[2]
    y = b0 + b1 * x
    for k in range(3, len(p), 2):
        area, mu = p[k], p[k + 1]
        y = y + area / (sigma * math.sqrt(2.0 * math.pi)) * \
            np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    return y


class ResultadoPico:
    __slots__ = ("energia_nominal", "intensidade", "u_intensidade",
                 "canal_centroide", "energia_centroide", "fwhm_keV",
                 "area", "u_area", "area_soma", "u_area_soma",
                 "janela", "ok", "motivo")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))


def ajusta_multipleto(esp, linhas, n_fwhm_janela=5.0, max_iter=3):
    """Ajusta um grupo de linhas proximas em um unico multigaussiano.

    linhas: lista de tuplas (E_keV, I_%, u_I_%)
    Retorna lista de ResultadoPico, na mesma ordem de 'linhas'.
    """
    energias_nom = [l[0] for l in linhas]
    centros = [esp.canal(e) for e in energias_nom]
    fwhm0 = max(esp.fwhm_canais(c) for c in centros)
    sigma0 = fwhm0 / 2.3548

    resultados = [ResultadoPico(energia_nominal=l[0], intensidade=l[1],
                                u_intensidade=l[2], ok=False, motivo="")
                  for l in linhas]

    # Refina o centroide do pico mais intenso buscando o maximo local
    principal = int(np.argmax([l[1] for l in linhas]))
    busca = int(max(4, round(1.5 * fwhm0)))
    c0 = int(round(centros[principal]))
    lo, hi = max(0, c0 - busca), min(esp.n_canais, c0 + busca + 1)
    if hi - lo < 3:
        for r in resultados:
            r.motivo = "pico fora da faixa de canais"
        return resultados
    desloc = (lo + int(np.argmax(esp.contagens[lo:hi]))) - centros[principal]
    centros = [c + desloc for c in centros]

    sigma = sigma0
    popt = pcov = None
    for _ in range(max_iter):
        meia = max(4.0, n_fwhm_janela * sigma * 2.3548 / 2.0)
        ini = int(round(min(centros) - meia))
        fim = int(round(max(centros) + meia)) + 1
        ini, fim = max(0, ini), min(esp.n_canais, fim)
        if fim - ini < 2 * len(linhas) + 4:
            for r in resultados:
                r.motivo = "janela estreita demais"
            return resultados

        x = np.arange(ini, fim, dtype=float)
        y = esp.contagens[ini:fim]

        # Fundo pelas bordas da janela
        n_borda = max(3, int(0.15 * len(x)))
        xb = np.concatenate([x[:n_borda], x[-n_borda:]])
        yb = np.concatenate([y[:n_borda], y[-n_borda:]])
        b1, b0 = np.polyfit(xb, yb, 1)
        fundo = b0 + b1 * x
        liquido = np.clip(y - fundo, 0, None)

        p0 = [b0, b1, sigma]
        limites_inf = [-np.inf, -np.inf, 0.3 * sigma0]
        limites_sup = [np.inf, np.inf, 4.0 * sigma0]
        for c, (e, inten, _u) in zip(centros, linhas):
            janela_pico = (np.abs(x - c) <= 1.5 * sigma * 2.3548)
            area0 = max(liquido[janela_pico].sum(), 1.0)
            p0 += [area0, c]
            limites_inf += [0.0, c - 2.0 * sigma0 * 2.3548]
            limites_sup += [np.inf, c + 2.0 * sigma0 * 2.3548]

        # Pesos: desvio de Poisson do numero de contagens brutas
        sigma_y = np.sqrt(np.clip(y, 1.0, None))

        try:
            popt, pcov = curve_fit(_multigaussiana, x, y, p0=p0,
                                   sigma=sigma_y, absolute_sigma=True,
                                   bounds=(limites_inf, limites_sup),
                                   maxfev=40000)
        except Exception as exc:                      # noqa: BLE001
            for r in resultados:
                r.motivo = "ajuste nao convergiu (%s)" % type(exc).__name__
            return resultados

        novo_sigma = popt[2]
        centros = [popt[4 + 2 * i] for i in range(len(linhas))]
        if abs(novo_sigma - sigma) / sigma < 0.02:
            sigma = novo_sigma
            break
        sigma = novo_sigma

    erros = np.sqrt(np.clip(np.diag(pcov), 0, None))
    passo_keV = abs(esp.energia(centros[0] + 1) - esp.energia(centros[0]))
    fwhm_keV = popt[2] * 2.3548 * passo_keV

    # Soma de canais menos fundo ajustado, como verificacao independente
    meia = 2.5 * popt[2]
    for i, r in enumerate(resultados):
        area, mu = popt[3 + 2 * i], popt[4 + 2 * i]
        u_area = erros[3 + 2 * i]
        r.canal_centroide = mu
        r.energia_centroide = float(esp.energia(mu))
        r.fwhm_keV = fwhm_keV
        r.area = area
        r.u_area = u_area
        r.janela = (ini, fim)

        a, b = int(round(mu - meia)), int(round(mu + meia)) + 1
        a, b = max(0, a), min(esp.n_canais, b)
        bruto = esp.contagens[a:b].sum()
        xs = np.arange(a, b, dtype=float)
        fundo_int = (popt[0] + popt[1] * xs).sum()
        # Desconta a contribuicao das OUTRAS gaussianas do multipleto
        vizinhos = 0.0
        for j in range(len(resultados)):
            if j == i:
                continue
            aj, muj = popt[3 + 2 * j], popt[4 + 2 * j]
            vizinhos += (aj / (popt[2] * math.sqrt(2 * math.pi)) *
                         np.exp(-0.5 * ((xs - muj) / popt[2]) ** 2)).sum()
        r.area_soma = bruto - fundo_int - vizinhos
        r.u_area_soma = math.sqrt(max(bruto, 1.0) + max(fundo_int, 0.0))

        if area <= 0 or u_area <= 0 or not np.isfinite(u_area):
            r.motivo = "area ou incerteza invalida"
        else:
            r.ok = True

    return resultados


def ajusta_espectro(esp, nuclideo, n_fwhm_janela=5.0):
    """Ajusta todas as linhas tabeladas do nuclideo presentes no espectro."""
    linhas = BIBLIOTECA[nuclideo]["linhas"]
    e_min = float(esp.energia(0))
    e_max = float(esp.energia(esp.n_canais - 1))

    grupos, soltas = {}, []
    for (energia, inten, u_inten, grupo) in linhas:
        if not (e_min + 5 < energia < e_max - 5):
            continue
        if grupo:
            grupos.setdefault(grupo, []).append((energia, inten, u_inten))
        else:
            soltas.append([(energia, inten, u_inten)])

    resultados = []
    for conjunto in soltas + list(grupos.values()):
        resultados.extend(ajusta_multipleto(esp, conjunto,
                                            n_fwhm_janela=n_fwhm_janela))
    resultados.sort(key=lambda r: r.energia_nominal)
    return resultados


# ----------------------------------------------------------------------------
# 4. DECAIMENTO E EFICIENCIA
# ----------------------------------------------------------------------------

def atividade_na_data(atividade_ref, data_ref, data_medida,
                      meia_vida_anos, u_atividade_rel=0.0,
                      u_meia_vida_anos=0.0):
    """Atividade corrigida por decaimento, com incerteza propagada.

    Retorna (A, u_A, delta_t_dias, fator_decaimento).
    """
    dias = (data_medida - data_ref).total_seconds() / 86400.0
    t_meio_dias = meia_vida_anos * DIAS_POR_ANO
    lamb = math.log(2.0) / t_meio_dias
    fator = math.exp(-lamb * dias)
    atividade = atividade_ref * fator

    # u(A)/A = sqrt( (u(A0)/A0)^2 + (ln2 * t/T * u(T)/T)^2 )
    termo_t = 0.0
    if u_meia_vida_anos:
        termo_t = (math.log(2.0) * dias / t_meio_dias) * \
                  (u_meia_vida_anos / meia_vida_anos)
    u_rel = math.sqrt(u_atividade_rel ** 2 + termo_t ** 2)
    return atividade, atividade * u_rel, dias, fator


def eficiencia_por_linha(esp, resultados, atividade, u_atividade,
                         usar_soma=False):
    """eps = N / (A * t_vivo * I). Retorna lista de dicionarios."""
    saida = []
    for r in resultados:
        if not r.ok:
            saida.append({"energia": r.energia_nominal, "ok": False,
                          "motivo": r.motivo, "eficiencia": float("nan"),
                          "u_eficiencia": float("nan")})
            continue
        area = r.area_soma if usar_soma else r.area
        u_area = r.u_area_soma if usar_soma else r.u_area
        inten = r.intensidade / 100.0
        u_inten = r.u_intensidade / 100.0
        denom = atividade * esp.tempo_vivo * inten
        eps = area / denom
        u_rel = math.sqrt((u_area / area) ** 2 +
                          (u_atividade / atividade) ** 2 +
                          (u_inten / inten) ** 2)
        saida.append({
            "energia": r.energia_nominal,
            "energia_ajustada": r.energia_centroide,
            "desvio_keV": r.energia_centroide - r.energia_nominal,
            "fwhm_keV": r.fwhm_keV,
            "intensidade_pct": r.intensidade,
            "area": r.area, "u_area": r.u_area,
            "area_soma": r.area_soma, "u_area_soma": r.u_area_soma,
            "cps_liquido": area / esp.tempo_vivo,
            "eficiencia": eps,
            "u_eficiencia": eps * u_rel,
            "u_rel_area_pct": 100.0 * u_area / area,
            "ok": True, "motivo": "",
        })
    return saida


# ----------------------------------------------------------------------------
# 5. AJUSTE DA CURVA log-log
# ----------------------------------------------------------------------------

class CurvaEficiencia:
    """ln(eps) = sum a_i * [ln(E/E0)]^i, ajustado por minimos quadrados."""

    def __init__(self, energias, eficiencias, u_eficiencias, grau=3,
                 energia_ref=None, escalar_banda=True):
        energias = np.asarray(energias, float)
        eficiencias = np.asarray(eficiencias, float)
        u_eficiencias = np.asarray(u_eficiencias, float)

        if len(energias) < grau + 1:
            raise ValueError("pontos insuficientes (%d) para grau %d"
                             % (len(energias), grau))

        # Centrar ln(E) na media geometrica das energias deixa a matriz de
        # Vandermonde bem condicionada; sem isso, graus >= 5 degeneram.
        if energia_ref is None:
            energia_ref = float(np.exp(np.mean(np.log(energias))))

        self.grau = grau
        self.energia_ref = energia_ref
        self.escalar_banda = escalar_banda
        self.energias = energias
        self.eficiencias = eficiencias
        self.u_eficiencias = u_eficiencias

        x = np.log(energias / energia_ref)
        y = np.log(eficiencias)
        # u(ln eps) = u(eps)/eps
        u_y = u_eficiencias / eficiencias

        V = np.vander(x, grau + 1, increasing=True)
        W = np.diag(1.0 / u_y ** 2)
        VtW = V.T @ W
        self.cov = np.linalg.inv(VtW @ V)
        self.coefs = self.cov @ (VtW @ y)

        residuos = y - V @ self.coefs
        self.n_pontos = len(x)
        self.gl = self.n_pontos - (grau + 1)
        self.chi2 = float(residuos @ W @ residuos)
        self.chi2_red = self.chi2 / self.gl if self.gl > 0 else float("nan")
        self.residuos_rel_pct = 100.0 * (np.exp(residuos) - 1.0)
        self.desvio_rms_pct = float(np.sqrt(np.mean(self.residuos_rel_pct ** 2)))

        # Quando chi2/gl > 1 as incertezas declaradas nao explicam a dispersao
        # observada (coincidencia verdadeira, intensidades da biblioteca,
        # geometria). Inflar a covariancia por chi2/gl e a pratica usual para
        # que a banda do ajuste reflita a dispersao real dos pontos.
        self.fator_banda = 1.0
        if escalar_banda and self.gl > 0 and self.chi2_red > 1.0:
            self.fator_banda = math.sqrt(self.chi2_red)
            self.cov = self.cov * self.chi2_red

    def __call__(self, energia):
        e = np.asarray(energia, float)
        x = np.log(np.atleast_1d(e) / self.energia_ref)
        V = np.vander(x, self.grau + 1, increasing=True)
        y = np.exp(V @ self.coefs)
        return float(y[0]) if e.ndim == 0 else y

    def incerteza(self, energia):
        """Incerteza 1-sigma da curva ajustada (banda de confianca)."""
        e = np.asarray(energia, float)
        x = np.log(np.atleast_1d(e) / self.energia_ref)
        V = np.vander(x, self.grau + 1, increasing=True)
        var_ln = np.einsum("ij,jk,ik->i", V, self.cov, V)
        u = np.exp(V @ self.coefs) * np.sqrt(np.clip(var_ln, 0, None))
        return float(u[0]) if e.ndim == 0 else u

    def formula(self):
        termos = " + ".join("%+.6g*L^%d" % (c, i)
                            for i, c in enumerate(self.coefs))
        return "ln(eps) = %s,  L = ln(E[keV]/%.6g)" % (termos, self.energia_ref)


def combina_pontos(lista_de_pontos, modo="media"):
    """Combina conjuntos de pontos de varios espectros por energia nominal.

    modo = 'media'     -> media aritmetica simples (o que se pede ao tratar
                          duas posicoes como extremos de uma fonte extensa)
    modo = 'ponderada' -> media ponderada por 1/u^2
    """
    chaves = None
    for pontos in lista_de_pontos:
        validos = {p["energia"] for p in pontos if p["ok"]}
        chaves = validos if chaves is None else (chaves & validos)
    chaves = sorted(chaves or [])

    combinados = []
    for energia in chaves:
        vals, uncs = [], []
        for pontos in lista_de_pontos:
            p = next(q for q in pontos if q["energia"] == energia)
            vals.append(p["eficiencia"])
            uncs.append(p["u_eficiencia"])
        vals = np.array(vals)
        uncs = np.array(uncs)

        if modo == "ponderada" and np.all(uncs > 0):
            w = 1.0 / uncs ** 2
            eps = float(np.sum(w * vals) / np.sum(w))
            u_eps = float(1.0 / math.sqrt(np.sum(w)))
        else:
            eps = float(np.mean(vals))
            # Incerteza da media aritmetica de n valores correlacionados pela
            # mesma atividade: propagacao direta dos u_i.
            u_eps = float(np.sqrt(np.sum(uncs ** 2)) / len(vals))

        # Dispersao entre as posicoes: nao e incerteza estatistica, e a
        # variacao real de geometria. Guardada separadamente.
        dispersao = float(np.max(vals) - np.min(vals))
        combinados.append({
            "energia": energia, "eficiencia": eps, "u_eficiencia": u_eps,
            "valores": vals.tolist(), "dispersao_abs": dispersao,
            "dispersao_pct": 100.0 * dispersao / eps if eps else float("nan"),
            "ok": True, "motivo": "",
        })
    return combinados


# ----------------------------------------------------------------------------
# 6. SAIDAS
# ----------------------------------------------------------------------------

def grava_csv_pontos(caminho, blocos, combinados=None):
    cols = ["espectro", "energia_nominal_keV", "energia_ajustada_keV",
            "desvio_keV", "fwhm_keV", "intensidade_pct", "area_liquida",
            "u_area", "area_soma_canais", "u_area_soma", "cps_liquido",
            "eficiencia", "u_eficiencia", "u_rel_area_pct"]
    with open(caminho, "w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for rotulo, pontos in blocos:
            for p in pontos:
                if not p["ok"]:
                    continue
                fh.write("%s,%.4f,%.4f,%+.4f,%.3f,%.4f,%.1f,%.1f,%.1f,%.1f,"
                         "%.6g,%.6e,%.6e,%.3f\n" % (
                             rotulo, p["energia"], p["energia_ajustada"],
                             p["desvio_keV"], p["fwhm_keV"],
                             p["intensidade_pct"], p["area"], p["u_area"],
                             p["area_soma"], p["u_area_soma"],
                             p["cps_liquido"], p["eficiencia"],
                             p["u_eficiencia"], p["u_rel_area_pct"]))
        if combinados:
            for p in combinados:
                fh.write("MEDIA,%.4f,,,,,,,,,,%.6e,%.6e,\n"
                         % (p["energia"], p["eficiencia"], p["u_eficiencia"]))


def grava_csv_curva(caminho, curva, n=400):
    e = np.logspace(np.log10(curva.energias.min()),
                    np.log10(curva.energias.max()), n)
    eps = curva(e)
    u = curva.incerteza(e)
    with open(caminho, "w", encoding="utf-8") as fh:
        fh.write("energia_keV,eficiencia,u_eficiencia,u_rel_pct\n")
        for ei, epsi, ui in zip(e, eps, u):
            fh.write("%.4f,%.6e,%.6e,%.4f\n"
                     % (ei, epsi, ui, 100.0 * ui / epsi))


def grafico(caminho, blocos, combinados, curva, curvas_individuais=None,
            titulo=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(9.5, 7.6), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})

    cores = ["#1f4e79", "#8c2d04", "#2d6a4f", "#6a3d9a"]
    for i, (rotulo, pontos) in enumerate(blocos):
        e = [p["energia"] for p in pontos if p["ok"]]
        y = [p["eficiencia"] for p in pontos if p["ok"]]
        u = [p["u_eficiencia"] for p in pontos if p["ok"]]
        ax.errorbar(e, y, yerr=u, fmt="o", ms=4.5, lw=1,
                    color=cores[i % len(cores)], alpha=0.75,
                    label=rotulo, zorder=3)
        if curvas_individuais and rotulo in curvas_individuais:
            c = curvas_individuais[rotulo]
            eg = np.logspace(np.log10(min(e)), np.log10(max(e)), 300)
            ax.plot(eg, c(eg), "--", lw=1.0,
                    color=cores[i % len(cores)], alpha=0.55, zorder=2)

    ec = np.array([p["energia"] for p in combinados])
    yc = np.array([p["eficiencia"] for p in combinados])
    uc = np.array([p["u_eficiencia"] for p in combinados])
    ax.errorbar(ec, yc, yerr=uc, fmt="s", ms=6, lw=1.2, color="k",
                label="media das duas posicoes", zorder=5)

    eg = np.logspace(np.log10(ec.min()), np.log10(ec.max()), 400)
    yg, ug = curva(eg), curva.incerteza(eg)
    ax.plot(eg, yg, "-", lw=2, color="#b00020",
            label="ajuste log-log grau %d" % curva.grau, zorder=4)
    ax.fill_between(eg, yg - 2 * ug, yg + 2 * ug, color="#b00020",
                    alpha=0.15, lw=0, label="banda 2$\\sigma$ do ajuste",
                    zorder=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Eficiência absoluta de pico")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8.5, loc="lower left")
    if titulo:
        ax.set_title(titulo, fontsize=11)

    res = 100.0 * (yc / curva(ec) - 1.0)
    axr.axhline(0, color="k", lw=0.8)
    axr.errorbar(ec, res, yerr=100.0 * uc / curva(ec), fmt="s", ms=5,
                 color="k", lw=1)
    axr.axhspan(-curva.desvio_rms_pct, curva.desvio_rms_pct,
                color="#b00020", alpha=0.12, lw=0)
    axr.set_xscale("log")
    axr.set_xlabel("Energia (keV)")
    axr.set_ylabel("Resíduo (%)")
    axr.grid(True, which="both", alpha=0.3)

    fig.subplots_adjust(left=0.095, right=0.98, top=0.945, bottom=0.075)
    fig.savefig(caminho, dpi=145)
    plt.close(fig)


# ----------------------------------------------------------------------------
# 7. LINHA DE COMANDO
# ----------------------------------------------------------------------------

def _parse_data(txt):
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d",
                "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%d/%m/%Y"):
        try:
            return dt.datetime.strptime(txt, fmt)
        except ValueError:
            continue
    raise argparse.ArgumentTypeError("data invalida: %r" % txt)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Curva de eficiencia de HPGe a partir de espectros .Spe",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--spe", action="append", required=True, metavar="ARQUIVO",
                    help="arquivo .Spe (repita a opcao para varios)")
    ap.add_argument("--rotulo", action="append", default=[],
                    help="rotulo para cada --spe, na mesma ordem")
    ap.add_argument("--nuclideo", default="Eu-152",
                    choices=sorted(BIBLIOTECA), help="nuclideo da fonte")
    ap.add_argument("--atividade", type=float, required=True,
                    help="atividade na data de referencia, em Bq")
    ap.add_argument("--u-atividade", type=float, default=0.0,
                    help="incerteza relativa da atividade, em %%")
    ap.add_argument("--data-ref", type=_parse_data, required=True,
                    help="data de referencia do certificado")
    ap.add_argument("--meia-vida", type=float, default=None,
                    help="meia-vida em anos (padrao: biblioteca)")
    ap.add_argument("--combinar", default="media",
                    choices=["media", "ponderada", "nenhuma"])
    ap.add_argument("--grau", type=int, default=3,
                    help="grau do polinomio log-log")
    ap.add_argument("--varrer-grau", action="store_true",
                    help="testar graus 2..6 e imprimir o diagnostico")
    ap.add_argument("--i-minimo", type=float, default=1.0,
                    help="intensidade minima da linha para o ajuste, em %%")
    ap.add_argument("--i-forte", type=float, default=10.0,
                    help="intensidade que define 'linha forte' na verificacao "
                         "cruzada, em %%")
    ap.add_argument("--max-incerteza", type=float, default=15.0,
                    help="incerteza relativa maxima da area, em %%")
    ap.add_argument("--desvio-maximo", type=float, default=2.0,
                    help="desvio maximo centroide-nominal aceito, em keV")
    ap.add_argument("--usar-soma", action="store_true",
                    help="usar soma de canais em vez da area gaussiana")
    ap.add_argument("--prefixo", default="eficiencia",
                    help="prefixo dos arquivos de saida")
    args = ap.parse_args(argv)

    dados = BIBLIOTECA[args.nuclideo]
    t_meio = args.meia_vida or dados["meia_vida_anos"]
    u_t_meio = 0.0 if args.meia_vida else dados.get("u_meia_vida_anos", 0.0)

    rel = []
    rel.append("=" * 78)
    rel.append("CALIBRACAO DE EFICIENCIA - HPGe")
    rel.append("=" * 78)
    rel.append("Fonte .............. %s" % args.nuclideo)
    rel.append("Atividade de ref. .. %.4g Bq em %s  (u = %.2f %%)"
               % (args.atividade, args.data_ref.strftime("%d/%m/%Y"),
                  args.u_atividade))
    rel.append("Meia-vida .......... %.4f anos (%.1f d)"
               % (t_meio, t_meio * DIAS_POR_ANO))
    rel.append("")

    blocos, curvas_ind = [], {}
    for i, caminho in enumerate(args.spe):
        rotulo = args.rotulo[i] if i < len(args.rotulo) else None
        esp = Espectro(caminho, rotulo)

        atividade, u_atividade, dias, fator = atividade_na_data(
            args.atividade, args.data_ref, esp.data_meia_contagem,
            t_meio, args.u_atividade / 100.0, u_t_meio)

        rel.append("-" * 78)
        rel.append("ESPECTRO: %s" % esp.rotulo)
        rel.append("  ID ............... %s" % esp.spec_id)
        rel.append("  Detector ......... %s   [%s]"
                   % (esp.detector, esp.aplicativo))
        rel.append("  Inicio da medida . %s"
                   % esp.data_medida.strftime("%d/%m/%Y %H:%M:%S"))
        rel.append("  Tempo vivo/real .. %.0f s / %.0f s  (tempo morto %.2f %%)"
                   % (esp.tempo_vivo, esp.tempo_real, esp.tempo_morto_pct))
        rel.append("  Canais ........... %d" % esp.n_canais)
        rel.append("  Calibracao ....... E = %s  [%s]"
                   % (" + ".join("%.6e*ch^%d" % (c, k) for k, c
                                 in enumerate(esp.cal_energia)),
                      esp.unidade_energia))
        rel.append("  Decaimento ....... %.2f d desde a referencia; "
                   "fator %.5f" % (dias, fator))
        rel.append("  Atividade ........ %.2f Bq  (u = %.2f Bq)"
                   % (atividade, u_atividade))

        resultados = ajusta_espectro(esp, args.nuclideo)
        pontos = eficiencia_por_linha(esp, resultados, atividade,
                                      u_atividade, usar_soma=args.usar_soma)

        # Filtros de qualidade
        aceitos = []
        rel.append("")
        rel.append("  %9s %11s %8s %7s %11s %9s %8s  %s"
                   % ("E_nom", "E_ajust", "dE", "FWHM", "area", "u_a%",
                      "eps", "uso"))
        for p in pontos:
            if not p["ok"]:
                rel.append("  %9.3f %11s %8s %7s %11s %9s %8s  REJEITADA (%s)"
                           % (p["energia"], "-", "-", "-", "-", "-", "-",
                              p["motivo"]))
                continue
            motivos = []
            if p["intensidade_pct"] < args.i_minimo:
                motivos.append("I < %.2g%%" % args.i_minimo)
            if p["u_rel_area_pct"] > args.max_incerteza:
                motivos.append("u(area) > %.3g%%" % args.max_incerteza)
            if abs(p["desvio_keV"]) > args.desvio_maximo:
                motivos.append("dE > %.2g keV" % args.desvio_maximo)
            usada = not motivos
            if usada:
                aceitos.append(p)
            rel.append("  %9.3f %11.3f %+8.3f %7.3f %11.1f %9.2f %8.3e  %s"
                       % (p["energia"], p["energia_ajustada"], p["desvio_keV"],
                          p["fwhm_keV"], p["area"], p["u_rel_area_pct"],
                          p["eficiencia"],
                          "sim" if usada else "nao (%s)" % "; ".join(motivos)))
        rel.append("")
        rel.append("  linhas aceitas: %d de %d" % (len(aceitos), len(pontos)))
        blocos.append((esp.rotulo, aceitos))

        if len(aceitos) >= args.grau + 1:
            try:
                curvas_ind[esp.rotulo] = CurvaEficiencia(
                    [p["energia"] for p in aceitos],
                    [p["eficiencia"] for p in aceitos],
                    [p["u_eficiencia"] for p in aceitos], grau=args.grau)
            except Exception:                              # noqa: BLE001
                pass

    # -- combinacao ---------------------------------------------------------
    if args.combinar == "nenhuma" or len(blocos) == 1:
        combinados = [dict(p) for p in blocos[0][1]]
    else:
        combinados = combina_pontos([b[1] for b in blocos],
                                    modo=args.combinar)

    rel.append("=" * 78)
    rel.append("COMBINACAO: %s de %d espectro(s)" % (args.combinar, len(blocos)))
    rel.append("=" * 78)
    if len(blocos) == 2 and "valores" in (combinados[0] if combinados else {}):
        rel.append("%10s %13s %13s %13s %9s %9s"
                   % ("E (keV)", "eps[1]", "eps[2]", "media", "disp.%",
                      "razao"))
        razoes = []
        for p in combinados:
            v = p["valores"]
            razao = v[0] / v[1] if v[1] else float("nan")
            razoes.append(razao)
            rel.append("%10.3f %13.5e %13.5e %13.5e %9.2f %9.3f"
                       % (p["energia"], v[0], v[1], p["eficiencia"],
                          p["dispersao_pct"], razao))
        razoes = np.array(razoes)
        rel.append("")
        rel.append("  Razao eps[1]/eps[2] = %.3f +- %.3f (desvio-padrao entre "
                   "linhas)" % (razoes.mean(), razoes.std(ddof=1)))
        rel.append("  Uma razao praticamente constante em energia indica que os")
        rel.append("  dois espectros diferem so por angulo solido (distancia),")
        rel.append("  e nao por absorcao ou efeito dependente da energia.")

    energias_c = [p["energia"] for p in combinados]
    efic_c = [p["eficiencia"] for p in combinados]
    u_efic_c = [p["u_eficiencia"] for p in combinados]

    if args.varrer_grau:
        rel.append("")
        rel.append("ESCOLHA DO GRAU")
        rel.append("  %6s %10s %9s %11s" % ("grau", "chi2/gl", "rms %",
                                            "max|res| %"))
        for g in range(2, 7):
            if len(energias_c) < g + 2:
                continue
            try:
                c = CurvaEficiencia(energias_c, efic_c, u_efic_c, grau=g)
            except Exception:                                  # noqa: BLE001
                continue
            rel.append("  %6d %10.2f %9.2f %11.2f"
                       % (g, c.chi2_red, c.desvio_rms_pct,
                          float(np.max(np.abs(c.residuos_rel_pct)))))

    curva = CurvaEficiencia(energias_c, efic_c, u_efic_c, grau=args.grau)

    rel.append("")
    rel.append("AJUSTE DA CURVA MEDIA")
    rel.append("  %s" % curva.formula())
    for i, (c, u) in enumerate(zip(curva.coefs,
                                   np.sqrt(np.diag(curva.cov)))):
        rel.append("    a%d = %+.8e  +-  %.2e" % (i, c, u))
    rel.append("  pontos = %d   graus de liberdade = %d"
               % (curva.n_pontos, curva.gl))
    rel.append("  chi2 = %.2f   chi2/gl = %.2f" % (curva.chi2, curva.chi2_red))
    rel.append("  desvio rms dos residuos = %.2f %%" % curva.desvio_rms_pct)
    if curva.fator_banda > 1.0:
        rel.append("  covariancia inflada por chi2/gl (banda x %.2f), pois as"
                   % curva.fator_banda)
        rel.append("  incertezas declaradas nao explicam a dispersao observada")
    rel.append("")
    rel.append("  %10s %11s %9s  %s" % ("E (keV)", "residuo %", "u_curva %",
                                        "obs"))
    for e, r in zip(curva.energias, curva.residuos_rel_pct):
        uc = curva.incerteza(float(e))
        marca = ""
        if abs(r) > 3.0 * curva.desvio_rms_pct:
            marca = "<== discrepante"
        elif abs(r) > 2.0 * curva.desvio_rms_pct:
            marca = "<-- atencao"
        rel.append("  %10.3f %+11.2f %9.2f  %s"
                   % (e, r, 100.0 * uc / curva(float(e)), marca))

    # -- verificacao cruzada com as linhas de maior intensidade -------------
    # A dispersao dos pontos vem quase toda das linhas fracas (estatistica e
    # intensidade tabelada). Reajustar usando so as linhas intensas e comparar
    # as duas curvas e a melhor medida da exatidao real da interpolacao.
    limiar_forte = args.i_forte
    e_fortes = []
    for p in combinados:
        inten = next((l[1] for l in dados["linhas"]
                      if abs(l[0] - p["energia"]) < 1e-6), 0.0)
        if inten >= limiar_forte:
            e_fortes.append(p)
    if len(e_fortes) >= args.grau + 2:
        cf = CurvaEficiencia([p["energia"] for p in e_fortes],
                             [p["eficiencia"] for p in e_fortes],
                             [p["u_eficiencia"] for p in e_fortes],
                             grau=args.grau)
        rel.append("")
        rel.append("VERIFICACAO CRUZADA - ajuste so com linhas de I >= %.3g %%"
                   % limiar_forte)
        rel.append("  pontos = %d   chi2/gl = %.2f   rms dos residuos = %.2f %%"
                   % (cf.n_pontos, cf.chi2_red, cf.desvio_rms_pct))
        rel.append("  %10s %13s %13s %9s"
                   % ("E (keV)", "curva cheia", "so fortes", "dif %"))
        for e in (121.8, 300.0, 411.8, 661.7, 1173.2, 1332.5, 1408.0):
            if not (curva.energias.min() <= e <= curva.energias.max()):
                continue
            a, b = curva(float(e)), cf(float(e))
            rel.append("  %10.1f %13.5e %13.5e %+9.2f"
                       % (e, a, b, 100.0 * (b / a - 1.0)))
        rel.append("  A diferenca entre as duas curvas e a estimativa realista")
        rel.append("  da exatidao da interpolacao; o rms do ajuste completo")
        rel.append("  e inflado pelas linhas fracas e e conservador.")

    rel.append("")
    rel.append("  Eficiencia interpolada em energias de interesse:")
    for e in (59.5, 121.8, 300.0, 411.8, 661.7, 1173.2, 1332.5, 1408.0):
        if curva.energias.min() <= e <= curva.energias.max():
            rel.append("    %8.1f keV -> %.4e  +- %.2e"
                       % (e, curva(float(e)), curva.incerteza(float(e))))

    rel.append("")
    rel.append("AVISOS")
    rel.append("  * Sem correcao de coincidencia verdadeira (TCS). O %s tem"
               % args.nuclideo)
    rel.append("    cascatas intensas; em geometria proxima as areas das linhas")
    rel.append("    em cascata saem subestimadas e a curva fica deformada.")
    rel.append("  * Sem subtracao de fundo ambiental nem auto-atenuacao.")
    rel.append("  * Curva valida so para a geometria medida.")
    if args.u_atividade == 0.0:
        rel.append("  * u(atividade) = 0: informe --u-atividade com o valor do")
        rel.append("    certificado para que a incerteza fique completa.")

    texto = "\n".join(rel)
    print(texto)

    base = args.prefixo
    os.makedirs(os.path.dirname(base) or ".", exist_ok=True)
    with open(base + "_relatorio.txt", "w", encoding="utf-8") as fh:
        fh.write(texto + "\n")
    grava_csv_pontos(base + "_pontos.csv", blocos,
                     combinados if len(blocos) > 1 else None)
    grava_csv_curva(base + "_curva.csv", curva)
    grafico(base + ".png", blocos, combinados, curva, curvas_ind,
            titulo="Eficiência de pico — %s — fonte de %s"
                   % (Espectro(args.spe[0]).detector, args.nuclideo))

    print("\nArquivos gravados:")
    for sufixo in ("_relatorio.txt", "_pontos.csv", "_curva.csv", ".png"):
        print("  " + base + sufixo)
    return 0


if __name__ == "__main__":
    sys.exit(main())
