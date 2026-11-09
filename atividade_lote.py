#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 atividade_lote.py - Atividade de um lote de amostras irradiadas em conjunto,
                     a partir de espectros HPGe (.Chn ou .Spe) e de uma curva
                     de eficiencia levantada com efic_hpge.py
===============================================================================

O que faz
---------
1.  Reconstroi (ou le) a curva de eficiencia de pico do detector.
2.  Para cada espectro do lote, ajusta o pico do radionuclideo de interesse
    (gaussiana sobre fundo linear), calcula a area liquida e sua incerteza.
3.  Calcula a atividade no inicio da contagem, corrigindo o decaimento
    OCORRIDO DURANTE a aquisicao:

        A(t_c) = N * C_d / ( eps(E) * I_gama * t_vivo )
        C_d    = lambda*t_real / (1 - exp(-lambda*t_real))

4.  Transporta a atividade ao fim da irradiacao t0:

        A(t0) = A(t_c) * exp(+lambda * dt),   dt = t_c - t0

5.  Calcula a taxa de reacao por unidade de massa saturada, util quando as
    amostras sao monitores de fluxo:

        A(t0) / (1 - exp(-lambda * t_irr))  -> atividade de saturacao
        e, se houver massa, a atividade especifica em Bq/g.

6.  Calcula a Atividade Minima Detectavel (criterio de Currie, 95 %).
7.  Grava CSV com o orcamento de incerteza termo a termo e imprime um
    resumo com a dispersao do lote.

Uso
---
    python3 atividade_lote.py \
        --dir pasta_com_espectros \
        --nuclideo Au-198 \
        --fim-irradiacao "2026-10-27 11:00" \
        --tempo-irradiacao 120 \
        --spe-eu inferior.Spe --spe-eu superior.Spe \
        --atividade-eu 4798 --data-ref-eu 2018-06-20 --u-atividade-eu 3 \
        --massas massas.csv \
        --saida resultados_lote.csv

'massas.csv' e opcional: duas colunas, 'arquivo' (ou 'amostra') e 'massa_g'.

Alternativa a reconstruir a curva: passar os coeficientes prontos com
    --curva-coefs "a0,a1,a2,a3" --curva-eref 649.952 --curva-u-rel 3.0
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import glob
import math
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import efic_hpge as EF                                        # noqa: E402


# ---------------------------------------------------------------------------
# Biblioteca de radionuclideos de ativacao (energia keV, I, T1/2 em segundos)
# ---------------------------------------------------------------------------
NUCLIDEOS = {
    "Au-198": dict(E=411.8022, I=0.95620, u_I=0.00060,
                   T12=2.69410 * 86400, u_T12=0.00020 * 86400),
    "Na-24":  dict(E=1368.626, I=0.999936, u_I=0.000015,
                   T12=14.997 * 3600, u_T12=0.012 * 3600),
    "Mn-56":  dict(E=846.7638, I=0.98850, u_I=0.00030,
                   T12=2.5789 * 3600, u_T12=0.0010 * 3600),
    "K-42":   dict(E=1524.58, I=0.1808, u_I=0.0009,
                   T12=12.355 * 3600, u_T12=0.010 * 3600),
    "La-140": dict(E=1596.21, I=0.9540, u_I=0.0007,
                   T12=1.67855 * 86400, u_T12=0.00012 * 86400),
    "Sc-46":  dict(E=889.277, I=0.999840, u_I=0.000010,
                   T12=83.79 * 86400, u_T12=0.04 * 86400),
    "Co-60":  dict(E=1332.492, I=0.999826, u_I=0.000006,
                   T12=1925.28 * 86400, u_T12=0.14 * 86400),
    "Cr-51":  dict(E=320.0824, I=0.09910, u_I=0.00010,
                   T12=27.7010 * 86400, u_T12=0.0012 * 86400),
    "Fe-59":  dict(E=1099.245, I=0.5656, u_I=0.0018,
                   T12=44.495 * 86400, u_T12=0.009 * 86400),
    "Zn-65":  dict(E=1115.539, I=0.5004, u_I=0.0010,
                   T12=243.93 * 86400, u_T12=0.09 * 86400),
    "Br-82":  dict(E=776.517, I=0.8339, u_I=0.0030,
                   T12=35.282 * 3600, u_T12=0.007 * 3600),
}


# ---------------------------------------------------------------------------
# Leitura unificada: .Chn (binario) e .Spe (ASCII)
# ---------------------------------------------------------------------------
MESES = {m: i + 1 for i, m in enumerate(
    ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
     "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"])}


class EspectroChn:
    """Espectro .Chn da ORTEC, com a mesma interface usada por efic_hpge."""

    def __init__(self, caminho, rotulo=None):
        self.caminho = caminho
        self.rotulo = rotulo or os.path.splitext(os.path.basename(caminho))[0]
        with open(caminho, "rb") as fh:
            d = fh.read()

        if struct.unpack_from("<h", d, 0)[0] != -1:
            raise ValueError("%s: assinatura .Chn invalida" % caminho)

        segundos = d[6:8].decode("ascii", "replace")
        t_real, t_vivo = struct.unpack_from("<ii", d, 8)
        data_txt = d[16:24].decode("ascii", "replace")
        hora_txt = d[24:28].decode("ascii", "replace")
        self.canal_inicial, n = struct.unpack_from("<hh", d, 28)

        self.tempo_real = t_real * 0.02
        self.tempo_vivo = t_vivo * 0.02
        self.contagens = np.asarray(
            struct.unpack_from("<%di" % n, d, 32), dtype=float)
        self.n_canais = n

        ano = int(data_txt[5:7]) + (2000 if data_txt[7:8] == "1" else 1900)
        self.data_medida = dt.datetime(
            ano, MESES.get(data_txt[2:5].upper(), 1), int(data_txt[0:2]),
            int(hora_txt[0:2]), int(hora_txt[2:4]), int(segundos))

        self.cal_energia = (0.0, 1.0)
        self.cal_forma = None
        self.spec_id = self.detector = ""
        rodape = 32 + 4 * n
        if len(d) >= rodape + 32:
            f = d[rodape:]
            if struct.unpack_from("<h", f, 0)[0] in (-101, -102):
                self.cal_energia = struct.unpack_from("<3f", f, 4)
                self.cal_forma = struct.unpack_from("<3f", f, 16)
                if len(f) >= 384:
                    self.detector = self._pstring(f, 256)
                    self.spec_id = self._pstring(f, 320)

    @staticmethod
    def _pstring(buf, off):
        n = buf[off]
        return buf[off + 1:off + 1 + n].decode("ascii", "replace").strip()

    @property
    def data_meia_contagem(self):
        return self.data_medida + dt.timedelta(seconds=0.5 * self.tempo_real)

    def energia(self, canal):
        canal = np.asarray(canal, dtype=float)
        return sum(c * canal ** i for i, c in enumerate(self.cal_energia))

    def canal(self, energia_alvo):
        chs = np.arange(self.n_canais, dtype=float)
        return float(np.interp(energia_alvo, self.energia(chs), chs))

    def fwhm_canais(self, canal):
        if self.cal_forma:
            v = sum(c * canal ** i for i, c in enumerate(self.cal_forma))
            if v > 0.5:
                return v
        return 6.0


def abre_espectro(caminho, rotulo=None):
    if caminho.lower().endswith(".chn"):
        return EspectroChn(caminho, rotulo)
    return EF.Espectro(caminho, rotulo)


# ---------------------------------------------------------------------------
# Curva de eficiencia
# ---------------------------------------------------------------------------
class CurvaFixa:
    """Curva a partir de coeficientes conhecidos, com incerteza relativa."""

    def __init__(self, coefs, energia_ref, u_rel_pct):
        self.coefs = np.asarray(coefs, float)
        self.energia_ref = energia_ref
        self.u_rel = u_rel_pct / 100.0
        self.grau = len(self.coefs) - 1

    def __call__(self, e):
        L = math.log(float(e) / self.energia_ref)
        return math.exp(sum(c * L ** i for i, c in enumerate(self.coefs)))

    def incerteza(self, e):
        return self(e) * self.u_rel

    def formula(self):
        termos = " + ".join("%+.6g*L^%d" % (c, i)
                            for i, c in enumerate(self.coefs))
        return "ln(eps) = %s,  L = ln(E/%.6g)" % (termos, self.energia_ref)


def curva_de_spe(args):
    """Reconstroi a curva media a partir dos .Spe de Eu-152."""
    blocos = []
    for caminho in args.spe_eu:
        esp = EF.Espectro(caminho)
        A, uA, _, _ = EF.atividade_na_data(
            args.atividade_eu, args.data_ref_eu, esp.data_meia_contagem,
            args.meia_vida_eu, args.u_atividade_eu / 100.0, 0.009)
        res = EF.ajusta_espectro(esp, "Eu-152")
        pts = EF.eficiencia_por_linha(esp, res, A, uA)
        blocos.append([p for p in pts if p["ok"]
                       and p["intensidade_pct"] >= 1.0
                       and p["u_rel_area_pct"] <= 15
                       and abs(p["desvio_keV"]) <= 2])
    comb = (EF.combina_pontos(blocos, "media") if len(blocos) > 1
            else [dict(p) for p in blocos[0]])
    return EF.CurvaEficiencia([p["energia"] for p in comb],
                              [p["eficiencia"] for p in comb],
                              [p["u_eficiencia"] for p in comb],
                              grau=args.grau)


# ---------------------------------------------------------------------------
# Calculo por amostra
# ---------------------------------------------------------------------------
def processa(esp, nuc, curva, t0, t_irr_s, u_eps_extra_pct, massa_g=None):
    """Devolve um dicionario com tudo o que entra no relatorio."""
    res = EF.ajusta_multipleto(esp, [(nuc["E"], 100.0, 0.0)])[0]
    passo = abs(float(esp.energia(1)) - float(esp.energia(0)))

    lam = math.log(2.0) / nuc["T12"]
    x = lam * esp.tempo_real
    C_d = x / (1.0 - math.exp(-x))

    eps = curva(float(nuc["E"]))
    u_eps_ajuste = curva.incerteza(float(nuc["E"]))
    # Incerteza total da eficiencia: banda do ajuste + sistematicas (certificado
    # da fonte, exatidao da interpolacao, geometria). Somadas em quadratura.
    u_eps = math.hypot(u_eps_ajuste, eps * u_eps_extra_pct / 100.0)

    saida = {
        "amostra": esp.rotulo,
        "id": getattr(esp, "spec_id", ""),
        "inicio_contagem": esp.data_medida,
        "t_vivo_s": esp.tempo_vivo,
        "t_real_s": esp.tempo_real,
        "tempo_morto_pct": 100.0 * (esp.tempo_real - esp.tempo_vivo)
                           / esp.tempo_real,
        "eps": eps, "u_eps": u_eps,
        "C_decaimento": C_d,
        "ok": res.ok, "motivo": res.motivo or "",
    }

    # Fundo e AMD, calculados mesmo quando o pico nao e detectado
    canal = int(round(esp.canal(nuc["E"])))
    fwhm_ch = esp.fwhm_canais(canal)
    meia = max(2, int(round(1.25 * fwhm_ch)))
    a, b = max(0, canal - 3 * meia), min(esp.n_canais, canal + 3 * meia)
    esq = esp.contagens[a:canal - meia]
    dir_ = esp.contagens[canal + meia:b]
    n_lat = len(esq) + len(dir_)
    B = (esq.sum() + dir_.sum()) * (2.0 * meia) / n_lat if n_lat else 0.0
    L_D = 2.71 + 4.65 * math.sqrt(max(B, 0.0))
    saida["fundo_janela"] = B
    saida["AMD_Bq"] = L_D * C_d / (eps * nuc["I"] * esp.tempo_vivo)

    if not res.ok:
        return saida

    N, uN = res.area, res.u_area
    A_c = N * C_d / (eps * nuc["I"] * esp.tempo_vivo)

    dt_s = (esp.data_medida - t0).total_seconds()
    A_0 = A_c * math.exp(lam * dt_s)

    # Termos relativos de incerteza
    t_N = uN / N
    t_eps = u_eps / eps
    t_I = nuc["u_I"] / nuc["I"]
    t_T = abs(math.log(2.0) * dt_s / nuc["T12"]) * (nuc["u_T12"] / nuc["T12"])
    u_rel_c = math.sqrt(t_N ** 2 + t_eps ** 2 + t_I ** 2)
    u_rel_0 = math.sqrt(u_rel_c ** 2 + t_T ** 2)

    sat = 1.0 - math.exp(-lam * t_irr_s)
    A_sat = A_0 / sat

    saida.update({
        "centroide_keV": res.energia_centroide,
        "desvio_keV": res.energia_centroide - nuc["E"],
        "fwhm_keV": res.fwhm_keV,
        "area": N, "u_area": uN,
        "area_soma": res.area_soma, "u_area_soma": res.u_area_soma,
        "confere_soma_pct": 100.0 * (res.area_soma / N - 1.0),
        "cps": N / esp.tempo_vivo,
        "espera_h": dt_s / 3600.0,
        "fator_decaimento": math.exp(lam * dt_s),
        "A_contagem_Bq": A_c, "u_A_contagem_Bq": A_c * u_rel_c,
        "A_t0_Bq": A_0, "u_A_t0_Bq": A_0 * u_rel_0,
        "A_saturacao_Bq": A_sat, "u_A_saturacao_Bq": A_sat * u_rel_0,
        "u_rel_total_pct": 100.0 * u_rel_0,
        "contrib_area_pct": 100.0 * t_N,
        "contrib_eps_pct": 100.0 * t_eps,
        "contrib_I_pct": 100.0 * t_I,
        "contrib_T12_pct": 100.0 * t_T,
        "passo_keV": passo,
    })
    if massa_g:
        saida["massa_g"] = massa_g
        saida["A_t0_Bq_por_g"] = A_0 / massa_g
        saida["u_A_t0_Bq_por_g"] = A_0 * u_rel_0 / massa_g
    return saida


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _data(txt):
    for f in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d",
              "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%d/%m/%Y"):
        try:
            return dt.datetime.strptime(txt, f)
        except ValueError:
            continue
    raise argparse.ArgumentTypeError("data invalida: %r" % txt)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Atividade de um lote de amostras a partir de espectros "
                    "HPGe e de uma curva de eficiencia",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--dir", help="pasta com os espectros do lote")
    ap.add_argument("--espectro", action="append", default=[],
                    help="espectro avulso (repita a opcao)")
    ap.add_argument("--padrao", default="*.Chn",
                    help="mascara dos arquivos dentro de --dir")
    ap.add_argument("--nuclideo", default="Au-198", choices=sorted(NUCLIDEOS))
    ap.add_argument("--fim-irradiacao", type=_data, required=True,
                    help="instante do FIM da irradiacao (t0)")
    ap.add_argument("--tempo-irradiacao", type=float, required=True,
                    help="duracao da irradiacao, em minutos")
    # Curva: ou a partir dos .Spe, ou por coeficientes
    ap.add_argument("--spe-eu", action="append", default=[],
                    help="espectro .Spe da fonte de Eu-152 (repita)")
    ap.add_argument("--atividade-eu", type=float, default=4798.0)
    ap.add_argument("--u-atividade-eu", type=float, default=0.0,
                    help="incerteza da atividade do certificado, em %%")
    ap.add_argument("--data-ref-eu", type=_data,
                    default=dt.datetime(2018, 6, 20))
    ap.add_argument("--meia-vida-eu", type=float, default=13.517)
    ap.add_argument("--grau", type=int, default=3)
    ap.add_argument("--curva-coefs",
                    help="coeficientes a0,a1,... em vez de reconstruir")
    ap.add_argument("--curva-eref", type=float, default=649.952)
    ap.add_argument("--curva-u-rel", type=float, default=3.0,
                    help="incerteza relativa da curva fixa, em %%")
    ap.add_argument("--u-eps-extra", type=float, default=2.0,
                    help="sistematica adicional da eficiencia (exatidao da "
                         "interpolacao, geometria), em %%")
    ap.add_argument("--massas", help="CSV com colunas arquivo/amostra e massa_g")
    ap.add_argument("--saida", default="resultados_lote.csv")
    args = ap.parse_args(argv)

    # -- espectros do lote --------------------------------------------------
    caminhos = list(args.espectro)
    if args.dir:
        caminhos += sorted(glob.glob(os.path.join(args.dir, args.padrao)))
    if not caminhos:
        ap.error("nenhum espectro encontrado (use --dir ou --espectro)")

    # -- massas opcionais ---------------------------------------------------
    massas = {}
    if args.massas:
        with open(args.massas, newline="", encoding="utf-8-sig") as fh:
            for linha in csv.DictReader(fh):
                chave = (linha.get("arquivo") or linha.get("amostra") or
                         "").strip()
                try:
                    massas[chave] = float(
                        str(linha.get("massa_g", "")).replace(",", "."))
                except ValueError:
                    pass

    # -- curva --------------------------------------------------------------
    if args.curva_coefs:
        coefs = [float(c) for c in args.curva_coefs.split(",")]
        curva = CurvaFixa(coefs, args.curva_eref, args.curva_u_rel)
    elif args.spe_eu:
        curva = curva_de_spe(args)
    else:
        ap.error("informe --spe-eu (um ou dois) ou --curva-coefs")

    nuc = NUCLIDEOS[args.nuclideo]
    t_irr_s = args.tempo_irradiacao * 60.0

    print("=" * 78)
    print("ATIVIDADE DO LOTE - %s" % args.nuclideo)
    print("=" * 78)
    print("Linha analisada .... %.4f keV   I = %.4f %%"
          % (nuc["E"], 100 * nuc["I"]))
    print("Meia-vida .......... %.6g s (%.4f h)"
          % (nuc["T12"], nuc["T12"] / 3600))
    print("Fim da irradiacao .. %s" % args.fim_irradiacao.strftime(
        "%d/%m/%Y %H:%M:%S"))
    print("Duracao irradiacao . %.1f min  ->  fator de saturacao %.6f"
          % (args.tempo_irradiacao,
             1 - math.exp(-math.log(2) / nuc["T12"] * t_irr_s)))
    print("Curva .............. %s" % curva.formula())
    print("eps(%.1f keV) ...... %.5e  +- %.2e"
          % (nuc["E"], curva(float(nuc["E"])), curva.incerteza(float(nuc["E"]))))
    print("Sistematica extra .. %.2f %% somada em quadratura" % args.u_eps_extra)
    print()

    linhas = []
    for caminho in caminhos:
        try:
            esp = abre_espectro(caminho)
        except Exception as exc:                               # noqa: BLE001
            print("  !! %s: %s" % (os.path.basename(caminho), exc))
            continue
        massa = massas.get(os.path.basename(caminho)) or massas.get(esp.rotulo)
        linhas.append(processa(esp, nuc, curva, args.fim_irradiacao,
                               t_irr_s, args.u_eps_extra, massa))

    cab = ("%-26s %16s %7s %9s %8s %12s %12s %7s"
           % ("amostra", "inicio contagem", "espera", "area", "u_a%",
              "A(t_c) Bq", "A(t0) Bq", "u_tot%"))
    print(cab)
    print("-" * len(cab))
    validos = []
    for r in linhas:
        if not r["ok"]:
            print("%-26s %16s  pico nao ajustado (%s); AMD = %.2f Bq"
                  % (r["amostra"][:26],
                     r["inicio_contagem"].strftime("%d/%m %H:%M:%S"),
                     r["motivo"], r["AMD_Bq"]))
            continue
        validos.append(r)
        alerta = "  <== contagem ANTES de t0, confira a data" \
            if r["espera_h"] < 0 else ""
        print("%-26s %16s %6.2fh %9.1f %8.2f %12.2f %12.2f %7.2f%s"
              % (r["amostra"][:26],
                 r["inicio_contagem"].strftime("%d/%m %H:%M:%S"),
                 r["espera_h"], r["area"], r["contrib_area_pct"],
                 r["A_contagem_Bq"], r["A_t0_Bq"], r["u_rel_total_pct"],
                 alerta))

    if validos:
        a0 = np.array([r["A_t0_Bq"] for r in validos])
        print()
        print("Dispersao do lote em A(t0): media %.2f Bq, desvio-padrao "
              "%.2f Bq (%.2f %%)" % (a0.mean(), a0.std(ddof=1),
                                     100 * a0.std(ddof=1) / a0.mean()))
        print("                     minimo %.2f Bq, maximo %.2f Bq, "
              "razao max/min %.3f" % (a0.min(), a0.max(), a0.max() / a0.min()))
        if all("A_t0_Bq_por_g" in r for r in validos):
            esp_ = np.array([r["A_t0_Bq_por_g"] for r in validos])
            print("Atividade especifica: media %.4g Bq/g, desvio-padrao "
                  "%.4g Bq/g (%.2f %%)"
                  % (esp_.mean(), esp_.std(ddof=1),
                     100 * esp_.std(ddof=1) / esp_.mean()))
        print()
        print("Orcamento medio de incerteza: area %.2f %% | eficiencia %.2f %% "
              "| I_gama %.3f %% | T(1/2) %.3f %%"
              % (np.mean([r["contrib_area_pct"] for r in validos]),
                 np.mean([r["contrib_eps_pct"] for r in validos]),
                 np.mean([r["contrib_I_pct"] for r in validos]),
                 np.mean([r["contrib_T12_pct"] for r in validos])))

    # -- CSV ----------------------------------------------------------------
    campos = ["amostra", "id", "inicio_contagem", "t_vivo_s", "t_real_s",
              "tempo_morto_pct", "centroide_keV", "desvio_keV", "fwhm_keV",
              "area", "u_area", "area_soma", "confere_soma_pct", "cps",
              "eps", "u_eps", "C_decaimento", "espera_h", "fator_decaimento",
              "A_contagem_Bq", "u_A_contagem_Bq", "A_t0_Bq", "u_A_t0_Bq",
              "A_saturacao_Bq", "u_A_saturacao_Bq", "u_rel_total_pct",
              "contrib_area_pct", "contrib_eps_pct", "contrib_I_pct",
              "contrib_T12_pct", "massa_g", "A_t0_Bq_por_g",
              "u_A_t0_Bq_por_g", "fundo_janela", "AMD_Bq", "motivo"]
    with open(args.saida, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        for r in linhas:
            linha = dict(r)
            linha["inicio_contagem"] = r["inicio_contagem"].strftime(
                "%Y-%m-%d %H:%M:%S")
            w.writerow(linha)
    print("\nCSV gravado em: %s" % args.saida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
