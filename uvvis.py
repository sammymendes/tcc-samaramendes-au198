#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analise de espectros UV-Vis de nanoparticulas de ouro em matrizes
biopolimericas, a partir da planilha exportada da leitora de placas.

Para cada amostra extrai:
  * lambda_max da banda de ressonancia plasmonica de superficie (RPS)
  * absorvancia no maximo e largura a meia altura da banda
  * razao A_RPS / A_450  -> diametro medio pelo metodo de Haiss (2007)
  * A_450                -> concentracao de particulas pelo mesmo metodo
  * indice de agregacao  A_RPS / A_700

Parametros adotados no TCC (secao 5.4), funcao parametros_tcc():
  * linha de base: media da absorvancia entre 900 e 960 nm
  * lambda_max: maximo do espectro suavizado (media movel de 7 pontos)
    entre 490 e 620 nm
  * razao pico/450 nm, razao pico/600 nm e indice de agregacao
    A(650 nm)/A(lambda_max), todos apos a subtracao da linha de base

A funcao analisa() e a versao original da analise (linha de base em 880-950
nm, indice A_700/A_RPS e diametro de Haiss). O diametro de Haiss NAO e
aplicavel a estas amostras, porque a absorcao dos biopolimeros se estende
ate 450 nm (TCC, secao 6.1.1); ele e mantido apenas como registro.

Uso
  python3 uvvis.py --planilha dados/brutos/Resultados_UV-VIS.xlsx \
                   [--saida resultados/uvvis_parametros.csv]

Referencia do metodo de dimensionamento:
  W. Haiss, N.T.K. Thanh, J. Aveyard, D.G. Fernig,
  "Determination of size and concentration of gold nanoparticles from
  UV-Vis spectra", Anal. Chem. 79 (2007) 4215-4221.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys

import numpy as np
from openpyxl import load_workbook

PLANILHA = "dados/brutos/Resultados_UV-VIS.xlsx"

# Coeficientes de Haiss et al. (2007) para AuNP esfericas
B1, B2 = 3.00, 2.20


def carrega(caminho=PLANILHA):
    wb = load_workbook(caminho, read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    linhas = list(ws.iter_rows(values_only=True))

    # Localiza a linha de cabecalho (a que comeca com 'Comprimento de onda')
    i_cab = next(i for i, r in enumerate(linhas)
                 if r and isinstance(r[0], str)
                 and r[0].strip().lower().startswith("comprimento de onda")
                 and r[1])
    nomes = [str(c).strip() for c in linhas[i_cab][1:] if c]

    lam, dados = [], [[] for _ in nomes]
    for r in linhas[i_cab + 1:]:
        if not r or not isinstance(r[0], (int, float)):
            continue
        lam.append(float(r[0]))
        for j in range(len(nomes)):
            v = r[1 + j]
            dados[j].append(float(v) if isinstance(v, (int, float)) else np.nan)

    return np.array(lam), nomes, {n: np.array(d) for n, d in zip(nomes, dados)}


def suaviza(y, n=7):
    """Media movel simples, para nao deixar o ruido deslocar o maximo."""
    k = np.ones(n) / n
    pad = n // 2
    yy = np.pad(y, pad, mode="edge")
    return np.convolve(yy, k, mode="valid")


def analisa(lam, a, faixa=(480.0, 700.0)):
    """Extrai os parametros da banda plasmonica."""
    s = suaviza(a)
    m = (lam >= faixa[0]) & (lam <= faixa[1])
    if not m.any():
        return None
    i = int(np.argmax(s[m]))
    lmax = float(lam[m][i])
    amax = float(s[m][i])

    def A(nm):
        return float(np.interp(nm, lam, s))

    # Linha de base no vermelho, onde nem o ouro nem o polimero absorvem
    base = float(np.mean(s[(lam >= 880) & (lam <= 950)]))

    # Largura a meia altura da banda, acima da linha de base
    meia = base + (amax - base) / 2.0
    esq = dire = float("nan")
    idx = int(np.argmin(np.abs(lam - lmax)))
    for k in range(idx, 0, -1):
        if s[k] <= meia:
            esq = float(np.interp(meia, [s[k], s[k + 1]], [lam[k], lam[k + 1]]))
            break
    for k in range(idx, len(lam) - 1):
        if s[k] <= meia:
            dire = float(np.interp(meia, [s[k], s[k - 1]], [lam[k], lam[k - 1]]))
            break
    fwhm = dire - esq if np.isfinite(esq) and np.isfinite(dire) else float("nan")

    a450 = A(450.0) - base
    aspr = amax - base
    razao = aspr / a450 if a450 > 0 else float("nan")

    # Diametro de Haiss: d = exp(B1 * (Aspr/A450) - B2)
    d = math.exp(B1 * razao - B2) if np.isfinite(razao) else float("nan")

    # Concentracao de particulas (cm^-1 de caminho optico):
    #   N = A450 * 1e14 / (d^2 * (-0.295 + 1.36*exp(-((d-96.8)/78.2)^2)))
    if np.isfinite(d) and d > 0:
        denom = d ** 2 * (-0.295 + 1.36 * math.exp(-(((d - 96.8) / 78.2) ** 2)))
        conc = a450 * 1e14 / denom if denom > 0 else float("nan")
    else:
        conc = float("nan")

    return {
        "lambda_max": lmax, "A_max": amax, "A_base": base,
        "A_spr": aspr, "A_450": a450, "razao": razao,
        "fwhm": fwhm, "borda_esq": esq, "borda_dir": dire,
        "d_nm": d, "N_part_mL": conc,
        "A_700": A(700.0) - base, "A_280": A(280.0),
        "agreg": (A(700.0) - base) / aspr if aspr > 0 else float("nan"),
        "satura_abaixo_nm": float(lam[np.argmax(a < 2.3)]) if (a >= 2.3).any() else None,
    }


def parametros_tcc(lam, a):
    """Parametros da banda plasmonica conforme a secao 5.4 do TCC."""
    s = suaviza(a)
    # linha de base: media do espectro SUAVIZADO entre 900 e 960 nm
    base = float(np.mean(s[(lam >= 900) & (lam <= 960)]))
    y = s - base
    m = (lam >= 490) & (lam <= 620)
    i = int(np.argmax(y[m]))
    lmax = float(lam[m][i])
    amax = float(y[m][i])
    A = lambda nm: float(np.interp(nm, lam, y))                # noqa: E731
    return {
        "lambda_max_nm": lmax,
        "A_RPS": amax,
        "A250": A(250.0),
        "razao_pico_450": amax / A(450.0) if A(450.0) != 0 else float("nan"),
        "razao_pico_600": amax / A(600.0) if A(600.0) != 0 else float("nan"),
        "indice_agregacao_650": A(650.0) / amax if amax != 0 else float("nan"),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--planilha", default=PLANILHA)
    ap.add_argument("--saida", default="resultados/uvvis_parametros.csv")
    ap.add_argument("--original", action="store_true",
                    help="imprime tambem a analise original (Haiss), so como registro")
    args = ap.parse_args(argv)

    lam, nomes, dados = carrega(args.planilha)
    print("Faixa: %.0f a %.0f nm, passo %.0f nm, %d amostras"
          % (lam.min(), lam.max(), lam[1] - lam[0], len(nomes)))
    print()
    cab = ("%-18s %8s %8s %9s %9s %9s"
           % ("amostra", "lmax", "A_RPS", "pico/450", "pico/600", "A650/Amax"))
    print(cab)
    print("-" * len(cab))
    res = {}
    for n in nomes:
        r = parametros_tcc(lam, dados[n])
        res[n] = r
        print("%-18s %8.1f %8.4f %9.3f %9.3f %9.3f"
              % (n, r["lambda_max_nm"], r["A_RPS"], r["razao_pico_450"],
                 r["razao_pico_600"], r["indice_agregacao_650"]))
    with open(args.saida, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["amostra", "lambda_max_nm", "A_RPS", "A250", "razao_pico_450",
                    "razao_pico_600", "indice_agregacao_650"])
        for n, r in res.items():
            w.writerow([n, "%.1f" % r["lambda_max_nm"], "%.4f" % r["A_RPS"],
                        "%.4f" % r["A250"], "%.3f" % r["razao_pico_450"],
                        "%.3f" % r["razao_pico_600"], "%.3f" % r["indice_agregacao_650"]])
    print("\nCSV: %s" % args.saida)

    if args.original:
        print("\nANALISE ORIGINAL (registro; diametro de Haiss nao aplicavel):")
        for n in nomes:
            r = analisa(lam, dados[n])
            print("%-18s %8.1f %8.4f %8.1f %8.3f %8.1f %9.3f"
                  % (n, r["lambda_max"], r["A_max"], r["fwhm"], r["razao"],
                     r["d_nm"], r["agreg"]))
    return lam, nomes, dados, res


if __name__ == "__main__":
    main()
