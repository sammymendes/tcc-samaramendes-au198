#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figuras_reanalise.py - Figuras 27 e 29 do TCC, a partir da saida de
reanalise_diluicao.py (resultados/resultados_por_amostra.json).

  Figura 27  atividade de Au-198 por amostra: medida e corrigida pela diluicao
  Figura 29  razao em relacao ao controle por sistema: sem correcao, corrigida
             pela diluicao e normalizada pelo Na-24

Uso
  python3 src/figuras_reanalise.py [--resultados resultados]
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Patch

COR = {"ALG": "#2a78d6", "CHI": "#eb6834", "GEL": "#1baf7a", "CTRL": "#7a7974"}
virgula = lambda v, nd=0: ("%.*f" % (nd, v)).replace(".", ",")


def estilo():
    fam = "Liberation Serif" if any("Liberation Serif" in f.name
                                    for f in font_manager.fontManager.ttflist) else "DejaVu Serif"
    plt.rcParams.update({"font.family": fam, "font.size": 10,
                         "axes.edgecolor": "#52514e", "xtick.color": "#52514e",
                         "ytick.color": "#52514e", "axes.spines.top": False,
                         "axes.spines.right": False})


def figura27(R, ctrl, caminho):
    fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=300)
    w = 0.38
    for i, r in enumerate(R):
        c = COR[r["sistema"]]
        ax.bar(i - w / 2 - 0.01, r["Au_A"], w, facecolor="white", edgecolor=c, hatch="////",
               linewidth=1.0, yerr=r["Au_uA"], error_kw=dict(lw=0.8, capsize=2, ecolor="#52514e"))
        ax.bar(i + w / 2 + 0.01, r["A_corr"], w, color=c, yerr=r["uA_corr"],
               error_kw=dict(lw=0.8, capsize=2, ecolor="#52514e"))
        ax.text(i + w / 2 + 0.01, r["A_corr"] + r["uA_corr"] + 12, virgula(r["A_corr"]),
                ha="center", va="bottom", fontsize=7)
    ax.axhline(ctrl["Au_A"], ls="--", lw=0.9, color="#7a7974")
    ax.text(-0.45, ctrl["Au_A"] + 7, "controle (%s Bq)" % virgula(ctrl["Au_A"]),
            ha="left", va="bottom", fontsize=7.5, color="#52514e")
    ax.set_xticks(range(len(R)))
    ax.set_xticklabels([r["rotulo"].replace(" ", "\n", 1) for r in R], fontsize=7.5)
    ax.set_ylabel("Atividade de Au-198 no fim da irradiação (Bq)")
    ax.set_ylim(0, 1000)
    ax.grid(axis="y", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    ax.legend(handles=[Patch(facecolor="white", edgecolor="#52514e", hatch="////",
                             label="medida (por alíquota de 0,5 mL)"),
                       Patch(facecolor="#52514e", label="corrigida pela diluição ($\\times f_d$)")],
              loc="upper center", ncol=2, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 1.1))
    fig.tight_layout(); fig.savefig(caminho); plt.close(fig)


def figura29(R, caminho):
    sis = [("ALG", "Au-ALG"), ("CHI", "Au-CHI"), ("GEL", "Au-GEL")]
    met = [("razao_bruta", "medida, sem correção", "o", "#b5b3ab"),
           ("razao_corr", "corrigida pela diluição", "s", "#2a78d6"),
           ("R_rel", "normalizada pelo Na-24", "^", "#eb6834")]
    fig, ax = plt.subplots(figsize=(6.4, 2.9), dpi=300)
    for j, (s, nome) in enumerate(sis):
        g = [r for r in R if r["sistema"] == s]
        for k, (chave, rotulo, mk, c) in enumerate(met):
            y = j + (k - 1) * 0.22
            if chave == "R_rel" and s == "ALG":
                ax.text(0.42, y, "não aplicável (Na do alginato)", va="center",
                        fontsize=7.5, color="#52514e")
                continue
            v = [r[chave] for r in g]
            ax.plot([min(v), max(v)], [y, y], color=c, lw=1.6, solid_capstyle="round")
            ax.plot(np.mean(v), y, marker=mk, ms=7, color=c, mec="white", mew=1, ls="none",
                    label=rotulo if j == 1 else None)
            ax.text(max(v) + 0.04, y, virgula(np.mean(v), 2), va="center", fontsize=7.5)
    ax.axvline(1, color="#7a7974", ls="--", lw=0.9)
    ax.text(1.01, 2.5, "controle = 1", fontsize=7.5, color="#52514e", va="bottom")
    ax.set_yticks(range(3)); ax.set_yticklabels([n for _, n in sis]); ax.invert_yaxis()
    ax.set_ylim(2.6, -0.5); ax.set_xlim(0, 2.3)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: virgula(v, 1)))
    ax.set_xlabel("Razão em relação ao controle (média; barra = mín.–máx. das três concentrações)")
    ax.grid(axis="x", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    ax.legend(loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(); fig.savefig(caminho); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resultados", default="resultados")
    a = ap.parse_args()
    D = json.load(open(os.path.join(a.resultados, "resultados_por_amostra.json"), encoding="utf-8"))
    R = D["amostras"]; ctrl = next(r for r in R if r["sistema"] == "CTRL")
    estilo()
    figura27(R, ctrl, os.path.join(a.resultados, "fig27_atividade_corrigida.png"))
    figura29(R, os.path.join(a.resultados, "fig29_razoes_controle.png"))
    print("figuras gravadas em", a.resultados)


if __name__ == "__main__":
    main()
