#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
teste_sintetico.py - Teste de ponta a ponta das rotinas com dados sinteticos.

Gera, com valores conhecidos, espectros .Spe no formato ORTEC (duas posicoes
de uma fonte de Eu-152 e uma amostra com Au-198 e Na-24) e uma planilha UV-Vis,
roda efic_hpge.py, atividade_lote.py, varredura.py e uvvis.py e confere se os
valores recuperados concordam com os verdadeiros dentro das incertezas.

Nao substitui a verificacao com os espectros reais; serve para garantir que o
fluxo funciona na estrutura do repositorio e para detectar regressoes.

Uso (a partir da raiz do repositorio)
  python3 tests/teste_sintetico.py
"""
import csv
import datetime as dt
import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(RAIZ, "src")
sys.path.insert(0, SRC)
import efic_hpge as EF  # noqa: E402

RNG = np.random.default_rng(20261027)
N_CANAIS = 16384
CAL = (0.0688557, 0.1485534, -1.23033e-8)          # calibracao do sistema real
COEF = [-5.44871, -1.04138, -0.0341442, 0.0526425]  # curva do TCC
E0 = 649.952


def energia(ch):
    return CAL[0] + CAL[1] * ch + CAL[2] * ch ** 2


def eps_verdadeira(e):
    L = math.log(e / E0)
    return math.exp(sum(c * L ** i for i, c in enumerate(COEF)))


def fwhm_kev(e):
    return 0.85 + 0.00065 * e


def espectro(linhas, live, fundo_nivel=25.0):
    """linhas: lista de (E, contagens esperadas no fotopico)."""
    ch = np.arange(N_CANAIS, dtype=float)
    E = energia(ch)
    mu = fundo_nivel * np.exp(-E / 500.0) + 1.5
    for e, n in linhas:
        s = fwhm_kev(e) / 2.3548 / CAL[1]
        c0 = np.interp(e, E, ch)
        mu += n / (s * math.sqrt(2 * math.pi)) * np.exp(-0.5 * ((ch - c0) / s) ** 2)
    return RNG.poisson(mu)


def grava_spe(caminho, contagens, inicio, live, real, spec_id):
    with open(caminho, "w", encoding="utf-8") as fh:
        fh.write("$SPEC_ID:\n%s\n$SPEC_REM:\nDETDESC# SINTETICO\nAP# teste\n" % spec_id)
        fh.write("$DATE_MEA:\n%s\n" % inicio.strftime("%m/%d/%Y %H:%M:%S"))
        fh.write("$MEAS_TIM:\n%.0f %.0f\n" % (live, real))
        fh.write("$DATA:\n0 %d\n" % (N_CANAIS - 1))
        fh.write("\n".join("%8d" % c for c in contagens) + "\n")
        fh.write("$MCA_CAL:\n3\n%.6E %.6E %.6E keV\n" % CAL)
        fh.write("$SHAPE_CAL:\n3\n%.6E %.6E %.6E\n" % (0.85 / CAL[1], 0.00065, 0.0))


def roda(args, cwd):
    r = subprocess.run([sys.executable, "-I"] + args, cwd=cwd,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:])
        raise SystemExit("falhou: %s" % " ".join(args))
    return r.stdout


def main():
    tmp = tempfile.mkdtemp(prefix="teste_tcc_")
    ok = True
    try:
        # ---- 1. Eu-152 nas duas posicoes ----------------------------------
        A0, data_ref = 4798.0, dt.datetime(2018, 6, 20)
        inicio = dt.datetime(2026, 10, 27, 15, 36, 43)
        T12 = 13.517 * 365.2422 * 86400
        for nome, fator, live in (("inferior", 2 * 1.613 / 2.613, 1200.0),
                                  ("superior", 2 / 2.613, 1800.0)):
            meio = inicio + dt.timedelta(seconds=live / 2)
            A = A0 * math.exp(-math.log(2) * (meio - data_ref).total_seconds() / T12)
            linhas = [(e, A * live * I / 100 * eps_verdadeira(e) * fator)
                      for (e, I, _u, _g) in EF.BIBLIOTECA["Eu-152"]["linhas"]]
            grava_spe(os.path.join(tmp, "Eu152_%s.Spe" % nome),
                      espectro(linhas, live), inicio, live, live * 1.02, "Eu-152 " + nome)
        out = roda([os.path.join(SRC, "efic_hpge.py"),
                    "--spe", "Eu152_inferior.Spe", "--spe", "Eu152_superior.Spe",
                    "--nuclideo", "Eu-152", "--atividade", "4798", "--data-ref", "2018-06-20",
                    "--combinar", "media", "--grau", "3", "--prefixo", "efic/eficiencia"], tmp)
        curva_pts = list(csv.DictReader(open(os.path.join(tmp, "efic/eficiencia_pontos.csv"))))
        media = [p for p in curva_pts if p["espectro"] == "MEDIA"]
        z = [(float(p["eficiencia"]) - eps_verdadeira(float(p["energia_nominal_keV"])))
             / float(p["u_eficiencia"]) for p in media]
        fortes = [abs(float(p["eficiencia"]) / eps_verdadeira(float(p["energia_nominal_keV"])) - 1)
                  for p in media if float(p["energia_nominal_keV"]) in (121.7817, 344.2785, 1408.013)]
        print("[efic_hpge] %d linhas; |z| max. = %.1f; desvio max. nas linhas de 121, 344 e 1408 keV = %.1f%% (estatistica de contagem)"
              % (len(media), max(map(abs, z)), 100 * max(fortes)))
        ok &= len(media) >= 10 and max(map(abs, z)) < 3.5

        # ---- 2. amostra com Au-198 e Na-24 --------------------------------
        t0 = dt.datetime(2026, 10, 27, 11, 30)
        tc = dt.datetime(2026, 10, 28, 10, 12, 11)
        live, real = 900.0, 936.0
        verdade = {"Au-198": (445.3, 411.8022, 0.9562, 2.69410 * 86400),
                   "Na-24": (1085.1, 1368.626, 0.999936, 14.997 * 3600)}
        linhas = []
        for nuc, (A_t0, E, I, T) in verdade.items():
            lam = math.log(2) / T
            A_c = A_t0 * math.exp(-lam * (tc - t0).total_seconds())
            Cd = lam * real / (1 - math.exp(-lam * real))
            linhas.append((E, A_c * eps_verdadeira(E) * I * live / Cd))
        lam_na = math.log(2) / verdade["Na-24"][3]
        A_c_na = verdade["Na-24"][0] * math.exp(-lam_na * (tc - t0).total_seconds())
        linhas.append((1732.007, 0.07 * linhas[1][1]))   # escape duplo do 2754 keV
        os.makedirs(os.path.join(tmp, "lote"))
        grava_spe(os.path.join(tmp, "lote", "Amostra_10.Spe"), espectro(linhas, live, 40),
                  tc, live, real, "Amostra Controle")
        for nuc in ("Au-198", "Na-24"):
            roda([os.path.join(SRC, "atividade_lote.py"), "--dir", "lote", "--padrao", "*.Spe",
                  "--nuclideo", nuc, "--fim-irradiacao", "2026-10-27 11:30",
                  "--tempo-irradiacao", "120",
                  "--spe-eu", "Eu152_inferior.Spe", "--spe-eu", "Eu152_superior.Spe",
                  "--saida", "lote_%s.csv" % nuc], tmp)
            r = next(csv.DictReader(open(os.path.join(tmp, "lote_%s.csv" % nuc))))
            A, u = float(r["A_t0_Bq"]), float(r["u_A_t0_Bq"])
            z = (A - verdade[nuc][0]) / u
            print("[atividade_lote] %-6s A(t0) = %.1f +- %.1f Bq (verdadeiro %.1f; %.1f sigma)"
                  % (nuc, A, u, verdade[nuc][0], z))
            ok &= abs(z) < 3

        out = roda([os.path.join(SRC, "varredura.py"), "--dir", "lote", "--padrao", "*.Spe",
                    "--saida", "varredura.csv"], tmp)
        ach = list(csv.DictReader(open(os.path.join(tmp, "varredura.csv"))))
        ident = {a["identificacao"].split("  (")[0] for a in ach}
        print("[varredura] %d picos; identificados: %s" % (len(ach), sorted(ident)))
        ok &= {"Au-198", "Na-24", "Na-24, escape duplo de 2754 keV"} <= ident

        # ---- 3. UV-Vis ------------------------------------------------------
        from openpyxl import Workbook
        lamb = np.arange(200, 1001, 1.0)
        wb = Workbook(); ws = wb.active
        ws.append(["Leitura sintetica"]); ws.append([])
        ws.append(["Comprimento de onda (nm)", "Controle"])
        y = (0.036 + 0.030 * np.exp(-0.5 * ((lamb - 527) / 28) ** 2)
             + 0.6 * np.exp(-(lamb - 200) / 25) + RNG.normal(0, 0.0004, lamb.size))
        for l, v in zip(lamb, y):
            ws.append([float(l), float(v)])
        wb.save(os.path.join(tmp, "uv.xlsx"))
        roda([os.path.join(SRC, "uvvis.py"), "--planilha", "uv.xlsx", "--saida", "uv.csv"], tmp)
        u = next(csv.DictReader(open(os.path.join(tmp, "uv.csv"))))
        print("[uvvis] lambda_max = %s nm (verdadeiro 527), A_RPS = %s (verdadeiro 0,030)"
              % (u["lambda_max_nm"], u["A_RPS"]))
        ok &= abs(float(u["lambda_max_nm"]) - 527) <= 2 and abs(float(u["A_RPS"]) - 0.030) < 0.003
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\nRESULTADO: %s" % ("todos os testes passaram" if ok else "FALHA"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
