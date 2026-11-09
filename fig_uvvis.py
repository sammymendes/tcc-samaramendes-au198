#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig_uvvis.py - Figura 25 do TCC: espectros UV-Vis das dez amostras.

  (a) espectros completos, com a linha de base de 900-960 nm subtraida
  (b) banda plasmonica normalizada no maximo
  (c) posicao (lambda_max) e assimetria (A(650)/A(lambda_max)) da banda
  (d) amplitude da banda confrontada com a atividade de Au-198 medida,
      com a coluna da placa indicada pelo marcador

Os espectros sao tratados por uvvis.parametros_tcc() (secao 5.4 do TCC):
media movel de 7 pontos e linha de base media entre 900 e 960 nm, ambas
sobre o espectro suavizado.

Entradas
  dados/brutos/Resultados_UV-VIS.xlsx   planilha da leitora de microplacas
  dados/amostras.csv                    poco -> rotulo do TCC
  dados/medidas_10_amostras.csv         atividade de Au-198 por amostra

Uso
  python3 src/fig_uvvis.py [--planilha ...] [--saida resultados/fig25_uvvis.png]
"""
import argparse
import csv
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uvvis  # noqa: E402

COR = {"ALG": "#2a78d6", "CHI": "#eb6834", "GEL": "#1baf7a", "CTRL": "#52514e"}
TRACO = {0: ":", 1: "--", 2: "-"}          # menor, intermediaria, maior concentracao
TINTA, TINTA2, GRADE = "#2b2b29", "#52514e", "#e6e5e0"
LIMIAR_AMPLITUDE = 0.007                   # abaixo disso a banda e proxima do ruido


def estilo():
    fam = "Liberation Serif" if any("Liberation Serif" in f.name
                                    for f in font_manager.fontManager.ttflist) else "DejaVu Serif"
    plt.rcParams.update({"font.family": fam, "font.size": 9,
                         "axes.edgecolor": TINTA2, "xtick.color": TINTA2,
                         "ytick.color": TINTA2, "axes.labelcolor": TINTA,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 9.5, "axes.titleweight": "bold",
                         "axes.titlelocation": "left",
                         "axes.formatter.use_locale": False})


def virgula(ax):
    """Separador decimal com virgula nos dois eixos."""
    from matplotlib.ticker import FuncFormatter

    def f(v, _):
        s = ("%g" % round(v, 6))
        return s.replace(".", ",").replace("-", "−")
    ax.xaxis.set_major_formatter(FuncFormatter(f))
    ax.yaxis.set_major_formatter(FuncFormatter(f))


def ler(caminho):
    with open(caminho, encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#") and l.strip()))


def amostras(dados):
    """Uma entrada por poco, com rotulo do TCC, sistema, nivel e atividade."""
    tab = ler(os.path.join(dados, "amostras.csv"))
    med = {r["poco_UV"]: r for r in ler(os.path.join(dados, "medidas_10_amostras.csv"))}
    out = {}
    for r in tab:
        rot = r["rotulo_tcc"]
        sist = "CTRL" if rot.startswith("Controle") else rot.split()[0].split("-")[1]
        nom = r["rotulo_nominal"]
        nivel = 2 if "6%" in nom else 1 if "3%" in nom else 0
        out[r["poco_UV"]] = dict(rotulo=rot, sistema=sist, nivel=nivel,
                                 coluna=r["poco_UV"][1:],
                                 A_Au=float(med[r["poco_UV"]]["Au198_A_t0_Bq"]))
    return out


def curto(rot):
    return rot.replace("Au-", "").replace("Controle", "Controle")


def grafico(lam, esp, info, saida):
    estilo()
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 6.6), dpi=300)
    (a, b), (c, d) = axs
    ordem = sorted(info, key=lambda p: (["ALG", "CHI", "GEL", "CTRL"].index(info[p]["sistema"]),
                                        info[p]["nivel"]))

    # (a) espectros completos
    for p in ordem:
        s, i = esp[p], info[p]
        ctrl = i["sistema"] == "CTRL"
        a.plot(lam, s["y"], color=COR[i["sistema"]], ls="-" if ctrl else TRACO[i["nivel"]],
               lw=1.6 if ctrl else 1.1)
    a.axvspan(490, 620, color="#f3ead2", lw=0, zorder=0)
    a.text(555, 0.295, "banda plasmônica", ha="center", va="top", fontsize=7.5,
           color=TINTA2)
    a.annotate("absorção dos biopolímeros\n(sai da escala abaixo de 250 nm)", xy=(238, 0.22),
               xytext=(320, 0.19), ha="left", va="center", fontsize=7.5, color=TINTA2,
               arrowprops=dict(arrowstyle="-", lw=0.6, color=TINTA2))
    a.set_xlim(200, 900); a.set_ylim(0, 0.30)
    a.set_xlabel("Comprimento de onda (nm)")
    a.set_ylabel("Absorvância (linha de base subtraída)")
    a.set_title("(a) Espectros completos")

    # (b) banda normalizada
    for p in ordem:
        s, i = esp[p], info[p]
        if s["amax"] < LIMIAR_AMPLITUDE:
            continue
        ctrl = i["sistema"] == "CTRL"
        b.plot(lam, s["y"] / s["amax"], color=COR[i["sistema"]],
               ls="-" if ctrl else TRACO[i["nivel"]], lw=1.6 if ctrl else 1.1)
    b.axvline(527, color="#9a9994", lw=0.7, zorder=0)
    b.text(525.5, -0.06, "527 nm", fontsize=7, color=TINTA2, va="bottom", ha="right")
    for p, x, (tx, ty) in (("A03", 600, (560, 1.17)), ("E02", 655, (700, 0.80))):
        b.annotate(f"{curto(info[p]['rotulo'])}: banda alargada,\nmáximo em {esp[p]['lmax']:.0f} nm",
                   xy=(x, float(np.interp(x, lam, esp[p]['y'] / esp[p]['amax']))),
                   xytext=(tx, ty), fontsize=7, color=TINTA, va="center",
                   arrowprops=dict(arrowstyle="-", lw=0.6, color=TINTA2))
    b.set_xlim(420, 780); b.set_ylim(-0.1, 1.38)
    b.set_xlabel("Comprimento de onda (nm)")
    b.set_ylabel("Absorvância normalizada no máximo")
    b.set_title("(b) Forma da banda plasmônica")
    b.text(0.99, 0.995, "amostras com A(λ$_{máx}$) < 0,007 omitidas", transform=b.transAxes,
           ha="right", va="top", fontsize=6.5, color=TINTA2)

    # (c) posicao e assimetria
    desloc_c = {"B03": (6, -2), "D03": (6, -2), "E02": (-6, -10), "A03": (-6, 4),
                "C03": (6, -2)}
    # agrupamento em 526-528 nm: rotulos empilhados a direita, com linhas de chamada
    pilha = {"C02": 0.46, "B02": 0.38, "D02": 0.30, "E03": 0.22, "A02": 0.14}
    for p in ordem:
        s, i = esp[p], info[p]
        fraca = s["amax"] < LIMIAR_AMPLITUDE
        c.scatter(s["lmax"], s["ind650"], s=34, marker="o", zorder=3,
                  facecolor="white" if fraca else COR[i["sistema"]],
                  edgecolor=COR[i["sistema"]], linewidth=1.3)
        if p in pilha:
            c.annotate(curto(i["rotulo"]), (s["lmax"], s["ind650"]), xytext=(530.2, pilha[p]),
                       fontsize=6.8, color=TINTA, ha="left", va="center",
                       arrowprops=dict(arrowstyle="-", lw=0.5, color="#9a9994",
                                       shrinkA=1, shrinkB=4))
            continue
        dx, dy = desloc_c[p]
        c.annotate(curto(i["rotulo"]), (s["lmax"], s["ind650"]), xytext=(dx, dy),
                   textcoords="offset points", fontsize=6.8, color=TINTA,
                   ha="left" if dx > 0 else "right", va="center")
    c.axhline(0, color="#9a9994", lw=0.6, zorder=0)
    c.set_xlim(523, 540); c.set_ylim(-0.35, 1.1)
    c.set_xlabel("λ$_{máx}$ da banda plasmônica (nm)")
    c.set_ylabel("A(650) / A(λ$_{máx}$)")
    c.set_title("(c) Posição e assimetria da banda")

    # (d) amplitude x atividade de Au-198
    desloc_d = {"A02": (-6, 7), "B03": (6, 0), "C02": (-6, -7), "D03": (6, 0),
                "E02": (6, 0), "A03": (6, -3), "B02": (6, 0), "C03": (6, 0),
                "D02": (6, 0), "E03": (6, 4)}
    for p in ordem:
        s, i = esp[p], info[p]
        fraca = s["amax"] < LIMIAR_AMPLITUDE
        d.scatter(i["A_Au"], s["amax"], s=34, zorder=3,
                  marker="s" if i["coluna"] == "02" else "o",
                  facecolor="white" if fraca else COR[i["sistema"]],
                  edgecolor=COR[i["sistema"]], linewidth=1.3)
        dx, dy = desloc_d[p]
        d.annotate(curto(i["rotulo"]), (i["A_Au"], s["amax"]), xytext=(dx, dy),
                   textcoords="offset points", fontsize=6.8, color=TINTA,
                   ha="left" if dx > 0 else "right", va="center")
    d.axhline(LIMIAR_AMPLITUDE, color="#9a9994", lw=0.7, ls="--", zorder=0)
    d.text(1000, LIMIAR_AMPLITUDE + 0.0006, "0,007: próximo do ruído", ha="right",
           fontsize=6.5, color=TINTA2)
    d.set_xlim(0, 1000); d.set_ylim(0, 0.042)
    d.set_xlabel("Atividade de Au-198 medida, A(t$_0$) (Bq)")
    d.set_ylabel("Amplitude da banda, A(λ$_{máx}$)")
    d.set_title("(d) Amplitude da banda e Au-198")
    d.legend(handles=[Line2D([], [], marker="s", ls="", color=TINTA2, label="coluna 02 da placa"),
                      Line2D([], [], marker="o", ls="", color=TINTA2, label="coluna 03 da placa")],
             loc="center right", bbox_to_anchor=(1.0, 0.64), fontsize=7, frameon=False,
             handletextpad=0.3)

    for ax in axs.flat:
        ax.grid(color=GRADE, lw=0.6); ax.set_axisbelow(True)
        virgula(ax)

    # legenda comum: cor = sistema, traco = concentracao
    h = [Line2D([], [], color=COR[k], lw=1.6, label=n) for k, n in
         (("ALG", "Au-ALG"), ("CHI", "Au-CHI"), ("GEL", "Au-GEL"), ("CTRL", "Controle"))]
    h += [Line2D([], [], color=TINTA2, lw=1.1, ls=TRACO[k], label=n) for k, n in
          ((0, "menor concentração"), (1, "concentração intermediária"),
           (2, "maior concentração"))]
    h += [Line2D([], [], marker="o", ls="", mfc="white", mec=TINTA2, label="A(λ$_{máx}$) < 0,007")]
    fig.legend(handles=h, loc="lower center", ncol=4, fontsize=7.5, frameon=False,
               bbox_to_anchor=(0.5, 0.0), columnspacing=1.4, handlelength=2.4)
    fig.tight_layout(rect=(0, 0.075, 1, 1), h_pad=1.6, w_pad=1.6)
    fig.savefig(saida, dpi=300)
    plt.close(fig)


def main(argv=None):
    raiz = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--planilha", default=os.path.join(raiz, "dados/brutos/Resultados_UV-VIS.xlsx"))
    ap.add_argument("--dados", default=os.path.join(raiz, "dados"))
    ap.add_argument("--saida", default=os.path.join(raiz, "resultados/fig25_uvvis.png"))
    args = ap.parse_args(argv)

    lam, nomes, A = uvvis.carrega(args.planilha)
    info = amostras(args.dados)
    esp = {}
    for nome in nomes:
        poco = re.search(r"\((\w+)\)", nome).group(1)
        s = uvvis.suaviza(A[nome])
        y = s - float(np.mean(s[(lam >= 900) & (lam <= 960)]))
        par = uvvis.parametros_tcc(lam, A[nome])
        esp[poco] = dict(y=y, lmax=par["lambda_max_nm"], amax=par["A_RPS"],
                         ind650=par["indice_agregacao_650"])
    os.makedirs(os.path.dirname(args.saida) or ".", exist_ok=True)
    grafico(lam, esp, info, args.saida)
    a_ = np.array([info[p]["A_Au"] for p in esp]); b_ = np.array([esp[p]["amax"] for p in esp])
    print("Figura 25 gravada em", args.saida)
    print("correlacao amplitude x A(Au-198): r = %.2f (n = %d)" % (np.corrcoef(a_, b_)[0, 1], len(a_)))


if __name__ == "__main__":
    main()
