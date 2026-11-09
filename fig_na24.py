#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig_na24.py - Figura 28 do TCC: Na-24 nas dez amostras.

  (a) atividades de Au-198 e de Na-24 no fim da irradiacao, por amostra
  (b) fundo sob o fotopico de 411,8 keV do Au-198 em funcao da area do
      fotopico de 1368,6 keV do Na-24, com a regressao linear (secao 6.5.3)

Entradas (saidas de atividade_lote.py, geradas por executar_tudo.sh)
  resultados/lote_Au-198_atividades.csv
  resultados/lote_Na-24_atividades.csv
  dados/amostras.csv                    arquivo -> rotulo do TCC

Uso
  python3 src/fig_na24.py [--resultados resultados] [--saida resultados/fig28_na24.png]
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

COR_SIS = {"ALG": "#2a78d6", "CHI": "#eb6834", "GEL": "#1baf7a", "CTRL": "#52514e"}
COR_AU, COR_NA = "#eda100", "#6250d6"
TINTA, TINTA2, GRADE = "#2b2b29", "#52514e", "#e6e5e0"


def estilo():
    fam = "Liberation Serif" if any("Liberation Serif" in f.name
                                    for f in font_manager.fontManager.ttflist) else "DejaVu Serif"
    plt.rcParams.update({"font.family": fam, "font.size": 9,
                         "axes.edgecolor": TINTA2, "xtick.color": TINTA2,
                         "ytick.color": TINTA2, "axes.labelcolor": TINTA,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 9.5, "axes.titleweight": "bold",
                         "axes.titlelocation": "left"})


def br(v, _=None):
    return ("%g" % round(v, 6)).replace(".", ",").replace("-", "−")


def ler(caminho):
    with open(caminho, encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#") and l.strip()))


def carrega(resultados, dados):
    rot = {r["arquivo"]: r["rotulo_tcc"] for r in ler(os.path.join(dados, "amostras.csv"))}
    au = {r["amostra"]: r for r in ler(os.path.join(resultados, "lote_Au-198_atividades.csv"))}
    na = {r["amostra"]: r for r in ler(os.path.join(resultados, "lote_Na-24_atividades.csv"))}
    linhas = []
    for arq in sorted(rot):
        r = rot[arq]
        sist = "CTRL" if r.startswith("Controle") else r.split()[0].split("-")[1]
        linhas.append(dict(
            rotulo=r, sistema=sist,
            A_Au=float(au[arq]["A_t0_Bq"]), u_Au=float(au[arq]["u_A_t0_Bq"]),
            A_Na=float(na[arq]["A_t0_Bq"]), u_Na=float(na[arq]["u_A_t0_Bq"]),
            fundo=float(au[arq]["fundo_janela"]), area_Na=float(na[arq]["area"])))
    return linhas


def grafico(L, saida):
    estilo()
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.5), dpi=300,
                               gridspec_kw=dict(width_ratios=[1.15, 1]))

    # (a) barras agrupadas, escala logaritmica
    x = np.arange(len(L)); w = 0.38
    ek = dict(lw=0.8, capsize=2, ecolor=TINTA2)
    a.bar(x - w / 2 - 0.01, [r["A_Au"] for r in L], w, color=COR_AU,
          yerr=[r["u_Au"] for r in L], error_kw=ek, label="Au-198 (411,8 keV)")
    a.bar(x + w / 2 + 0.01, [r["A_Na"] for r in L], w, color=COR_NA,
          yerr=[r["u_Na"] for r in L], error_kw=ek, label="Na-24 (1368,6 keV)")
    a.set_yscale("log"); a.set_ylim(80, 2e4)
    a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%d" % v if v < 1e4 else "10 000"))
    a.set_xticks(x)
    a.set_xticklabels([r["rotulo"].replace("Au-", "").replace(" ", "\n", 1) for r in L],
                      fontsize=7)
    a.set_ylabel("Atividade no fim da irradiação (Bq)")
    a.set_title("(a) Au-198 e Na-24 por amostra")
    a.legend(loc="upper right", fontsize=7.5, frameon=False, ncol=2,
             handlelength=1.2, columnspacing=1.0)
    a.grid(axis="y", color=GRADE, lw=0.6); a.set_axisbelow(True)

    # (b) fundo sob o pico do ouro x area do Na-24
    xs = np.array([r["area_Na"] for r in L]); ys = np.array([r["fundo"] for r in L])
    inc, b0 = np.polyfit(xs, ys, 1)
    r2 = np.corrcoef(xs, ys)[0, 1] ** 2
    xx = np.linspace(0, 4500, 50)
    b.plot(xx, b0 + inc * xx, color="#9a9994", lw=1.2, zorder=1)
    desloc = {"ALG 0,5%": (-6, 6), "ALG 1,0%": (6, -6), "ALG 2,0%": (-6, 6),
              "CHI 0,6%": (-6, -7), "CHI 1,2%": (6, -6), "CHI 1,8%": (-6, 6),
              "GEL 0,05%": (6, -5), "GEL 0,11%": (6, 5), "GEL 0,21%": (-6, 5),
              "Controle": (6, 5)}
    for r in L:
        b.scatter(r["area_Na"], r["fundo"], s=30, color=COR_SIS[r["sistema"]], zorder=3,
                  edgecolor="white", linewidth=0.8)
    # rotulos: amostras isoladas direto; o agrupamento perto da origem em pilha
    perto = [r for r in L if r["area_Na"] < 800]
    for r in L:
        nome = r["rotulo"].replace("Au-", "")
        if r in perto:
            continue
        dx, dy = desloc[nome]
        b.annotate(nome, (r["area_Na"], r["fundo"]), xytext=(dx, dy), textcoords="offset points",
                   fontsize=7, color=TINTA, ha="left" if dx > 0 else "right", va="center")
    for k, r in enumerate(sorted(perto, key=lambda r: -r["fundo"])):
        b.annotate(r["rotulo"].replace("Au-", ""), (r["area_Na"], r["fundo"]),
                   xytext=(1700, 38 - 6.0 * k), fontsize=7, color=TINTA, va="center",
                   arrowprops=dict(arrowstyle="-", lw=0.5, color="#b5b4ae",
                                   shrinkA=1, shrinkB=3))
    b.text(0.03, 0.97, "fundo = %s + %s × área(Na-24)\nr² = %s"
           % (br(round(b0, 1)), br(round(inc, 4)), br(round(r2, 3))),
           transform=b.transAxes, ha="left", va="top", fontsize=7.5, color=TINTA)
    b.set_xlim(0, 4500); b.set_ylim(0, 120)
    b.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%d" % v))
    b.set_xlabel("Área do fotopico de 1368,6 keV do Na-24 (contagens)")
    b.set_ylabel("Fundo sob o fotopico de 411,8 keV (contagens)")
    b.set_title("(b) Fundo sob o pico do Au-198")
    b.grid(color=GRADE, lw=0.6); b.set_axisbelow(True)
    b.legend(handles=[Line2D([], [], marker="o", ls="", color=COR_SIS[k], label=n)
                      for k, n in (("ALG", "Au-ALG"), ("CHI", "Au-CHI"), ("GEL", "Au-GEL"),
                                   ("CTRL", "Controle"))],
             loc="lower right", fontsize=7, frameon=False, handletextpad=0.2)

    fig.tight_layout(w_pad=2.0)
    fig.savefig(saida, dpi=300)
    plt.close(fig)
    return b0, inc, r2


def main(argv=None):
    raiz = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resultados", default=os.path.join(raiz, "resultados"))
    ap.add_argument("--dados", default=os.path.join(raiz, "dados"))
    ap.add_argument("--saida", default=os.path.join(raiz, "resultados/fig28_na24.png"))
    args = ap.parse_args(argv)
    L = carrega(args.resultados, args.dados)
    b0, inc, r2 = grafico(L, args.saida)
    print("Figura 28 gravada em", args.saida)
    print("regressao: fundo(411,8 keV) = %.1f + %.4f x area(Na-24), r2 = %.3f" % (b0, inc, r2))


if __name__ == "__main__":
    main()
