#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reanalise_diluicao.py - Correcao das atividades de Au-198 pela diluicao das
formulacoes, normalizacao pelo Na-24 e comparacao com a atividade teorica
(secoes 5.3.4, 5.7, 6.4 a 6.7 do TCC).

Entradas
  dados/medidas_10_amostras.csv  atividades de Au-198 e Na-24 no fim da
                                 irradiacao (saida de atividade_lote.py e
                                 varredura.py), uma linha por amostra
  dados/lote_Na24_atividades.csv saida de atividade_lote.py para o Na-24
                                 (usada para a incerteza das areas do Na-24)
  dados/formulacoes.csv          preparo: rotulo corrigido, sistema, massa de
                                 biopolimero, volume final e fator de diluicao

Saidas (pasta resultados/)
  tabela7_atividade_corrigida.csv   Tabela 7 do TCC
  tabela8_na24_tracador.csv         Tabela 8
  tabela9_atividade_teorica.csv     Tabela 9
  tabela10_medias_por_sistema.csv   Tabela 10
  resultados_por_amostra.json       todos os valores, para figuras_reanalise.py
  relatorio_reanalise.txt           resumo impresso na tela

Modelo
  Correcao pela diluicao      A_corr = f_d * A(t0),  f_d = V_final / V_AuNP
  Normalizacao pelo Na-24     R = A(Au-198)/A(Na-24);  A_Na = R * A_ctrl(Na-24)
  Equacao de ativacao         A = N (sigma0*phi_th + I0*phi_epi)(1 - exp(-lambda*t_irr))

As incertezas das razoes entre amostras incluem apenas a contagem (areas
liquidas): eficiencia, intensidade de emissao e meia-vida sao comuns a todas
as amostras e se cancelam nas razoes.

Uso
  python3 src/reanalise_diluicao.py [--dados dados] [--saida resultados]
        [--massa-haucl4-mg 5.0] [--fluxo 1e9]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os

import numpy as np

# --------------------------------------------------------------------------
# Constantes fisicas e dados nucleares
# --------------------------------------------------------------------------
NA = 6.02214076e23            # mol^-1
M_AU = 196.96657              # g/mol
M_HAUCL4_3H2O = 393.83        # g/mol (HAuCl4.3H2O)
M_NA = 22.98977               # g/mol
M_CITRATO = 294.10            # g/mol (citrato trissodico di-hidratado)
FRACAO_NA_ALGINATO = M_NA / 198.11   # alginato de sodio puro (unidade monomerica)

T12_AU = 2.69410 * 86400.0    # s   (NNDC)
T12_NA = 14.997 * 3600.0      # s   (NNDC)
SIGMA0_AU = 98.65e-24         # cm^2 (Mughabghab, 2006)
I0_AU = 1550e-24              # cm^2
SIGMA0_NA = 0.530e-24         # cm^2
I0_NA = 0.311e-24             # cm^2

# Preparo (secao 5.2 e 5.3)
V_SINTESE_ML = 60.0           # 50 mL de HAuCl4 + 10 mL de citrato
M_CITRATO_SINTESE_G = 0.100   # 10 mL de citrato a 1% (m/v)
V_ALIQUOTA_ML = 0.5
V_AUNP_ML = 10.0
T_IRR_S = 7200.0


def ler_csv(caminho):
    with open(caminho, encoding="utf-8") as fh:
        linhas = [l for l in fh if not l.startswith("#") and l.strip()]
    return list(csv.DictReader(linhas))


def fator_saturacao(t12, t_irr=T_IRR_S):
    return 1.0 - math.exp(-math.log(2.0) / t12 * t_irr)


def fmt(x, nd=1):
    return ("%.*f" % (nd, x)).replace(".", ",")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dados", default="dados")
    ap.add_argument("--saida", default="resultados")
    ap.add_argument("--massa-haucl4-mg", type=float, default=5.0,
                    help="massa de HAuCl4.3H2O por sintese (padrao: nominal, 5,0 mg)")
    ap.add_argument("--fluxo", type=float, default=1e9,
                    help="fluxo termico de referencia, cm^-2 s^-1 (padrao 1e9)")
    args = ap.parse_args()
    os.makedirs(args.saida, exist_ok=True)

    med = {r["amostra"]: r for r in ler_csv(os.path.join(args.dados, "medidas_10_amostras.csv"))}
    form = ler_csv(os.path.join(args.dados, "formulacoes.csv"))
    lote_na = {r["id"].replace("Amostra ", "", 1).strip(): r
               for r in ler_csv(os.path.join(args.dados, "lote_Na24_atividades.csv"))}

    R = []
    for f in form:
        m = med[f["amostra"]]
        r = {
            "rotulo": f["rotulo"], "sistema": f["sistema"],
            "massa_mg": float(f["massa_biopolimero_mg"]), "f_d": float(f["f_d"]),
            "Au_area": float(m["Au198_area"]), "Au_u_area": float(m["Au198_u_area"]),
            "Au_A": float(m["Au198_A_t0_Bq"]), "Au_uA": float(m["Au198_u_A_t0_Bq"]),
            "Au_AMD": float(m["Au198_AMD_Bq"]),
            "Na_area": float(m["Na24_area"]), "Na_u_area": float(lote_na[f["amostra"]]["u_area"]),
            "Na_A": float(m["Na24_A_t0_Bq"]), "Na_uA": float(m["Na24_u_A_t0_Bq"]),
            "inicio": m["inicio_contagem"], "t_real": float(m["t_real_s"]),
        }
        R.append(r)
    ctrl = next(r for r in R if r["sistema"] == "CTRL")
    urc = lambda r, k: r[k + "_u_area"] / r[k + "_area"]   # incerteza relativa de contagem

    # ---- correcao pela diluicao e normalizacao pelo Na-24 ------------------
    for r in R:
        r["A_corr"] = r["f_d"] * r["Au_A"]
        r["uA_corr"] = r["f_d"] * r["Au_uA"]
        r["razao_bruta"] = r["Au_A"] / ctrl["Au_A"]
        r["razao_corr"] = r["A_corr"] / ctrl["Au_A"]
        r["u_razao_corr"] = r["razao_corr"] * math.hypot(urc(r, "Au"), urc(ctrl, "Au"))
        r["Na_corr"] = r["f_d"] * r["Na_A"]
        r["R"] = r["Au_A"] / r["Na_A"]
        r["uR"] = r["R"] * math.hypot(urc(r, "Au"), urc(r, "Na"))
    for r in R:
        r["R_rel"] = r["R"] / ctrl["R"]
        r["u_R_rel"] = r["R_rel"] * math.sqrt(urc(r, "Au") ** 2 + urc(r, "Na") ** 2
                                              + urc(ctrl, "Au") ** 2 + urc(ctrl, "Na") ** 2)
        r["A_Na"] = r["R"] * ctrl["Na_A"]
        r["tracador_aplicavel"] = r["sistema"] != "ALG"
        r["fd_pelo_Na"] = ctrl["Na_A"] / r["Na_A"]

    # ---- atividade teorica (equacao de ativacao) ---------------------------
    m_au_sintese = args.massa_haucl4_mg * 1e-3 * M_AU / M_HAUCL4_3H2O   # g
    conc_au = m_au_sintese / V_SINTESE_ML                                # g/mL
    S_au = fator_saturacao(T12_AU)
    S_na = fator_saturacao(T12_NA)
    for r in R:
        m_al = conc_au * V_ALIQUOTA_ML / r["f_d"]
        N = m_al / M_AU * NA
        r["m_Au_ug"] = m_al * 1e6
        r["A_teorica"] = N * SIGMA0_AU * args.fluxo * S_au
        r["med_sobre_teo"] = r["Au_A"] / r["A_teorica"]

    # sodio do citrato no controle
    n_na = M_CITRATO_SINTESE_G / M_CITRATO * 3 / V_SINTESE_ML * V_ALIQUOTA_ML   # mol
    A_na_teo = n_na * NA * SIGMA0_NA * args.fluxo * S_na
    fator_na_citrato = ctrl["Na_A"] / A_na_teo
    # sodio do alginato: inclinacao de A(Na-24) contra a massa de polimero
    alg = [r for r in R if r["sistema"] == "ALG"]
    x = np.array([r["massa_mg"] for r in alg]); y = np.array([r["Na_A"] for r in alg])
    w = 1.0 / np.array([r["Na_uA"] for r in alg]) ** 2
    X = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(X.T @ np.diag(w) @ X)
    p = cov @ X.T @ np.diag(w) @ y
    m_na_por_mg = 1e-3 * FRACAO_NA_ALGINATO * V_ALIQUOTA_ML / V_AUNP_ML    # g de Na por mg de alginato
    A_na_alg_teo = m_na_por_mg / M_NA * NA * SIGMA0_NA * args.fluxo * S_na
    fator_na_alg = p[1] / A_na_alg_teo
    phi_eff_au = ctrl["Au_A"] / ctrl["A_teorica"] * args.fluxo
    phi_eff_na = 0.5 * (fator_na_citrato + fator_na_alg) * args.fluxo

    # ---- medias por sistema -------------------------------------------------
    medias = {}
    for s in ("ALG", "CHI", "GEL", "CTRL"):
        g = [r for r in R if r["sistema"] == s]
        mm = lambda k: float(np.mean([r[k] for r in g]))
        medias[s] = {"A_medida": mm("Au_A"), "A_corr": mm("A_corr"),
                     "A_Na": mm("A_Na") if s != "ALG" else None,
                     "Na_corr": mm("Na_corr")}
        for k in ("A_medida", "A_corr", "A_Na"):
            v = medias[s][k]
            medias[s]["razao_" + k] = v / ctrl["Au_A"] if v is not None else None

    # ---- gravacao ------------------------------------------------------------
    def grava(nome, cab, linhas):
        with open(os.path.join(args.saida, nome), "w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh, delimiter=";")
            wr.writerow(cab); wr.writerows(linhas)

    grava("tabela7_atividade_corrigida.csv",
          ["amostra", "A_R", "u_A_R", "A_t0_Bq", "u_Bq", "AMD_Bq", "f_d", "A_corr_Bq", "u_corr_Bq",
           "A_corr_sobre_ctrl", "u_razao"],
          [[r["rotulo"], fmt(r["Au_area"]), fmt(r["Au_u_area"]), fmt(r["Au_A"]), fmt(r["Au_uA"]),
            fmt(r["Au_AMD"]), fmt(r["f_d"]), fmt(r["A_corr"]), fmt(r["uA_corr"]),
            fmt(r["razao_corr"], 3), fmt(r["u_razao_corr"], 3)] for r in R])
    grava("tabela8_na24_tracador.csv",
          ["amostra", "A_Na24_t0_Bq", "u_Bq", "fd_x_A_Bq", "fd_x_A_sobre_ctrl", "R_Au_Na", "u_R",
           "R_sobre_Rctrl", "u", "A_Au_normalizada_Bq"],
          [[r["rotulo"], fmt(r["Na_A"]), fmt(r["Na_uA"]), fmt(r["Na_corr"]),
            fmt(r["Na_corr"] / ctrl["Na_A"], 3), fmt(r["R"], 4), fmt(r["uR"], 4),
            fmt(r["R_rel"], 3) if r["tracador_aplicavel"] else "n.a.",
            fmt(r["u_R_rel"], 3) if r["tracador_aplicavel"] else "n.a.",
            fmt(r["A_Na"]) if r["tracador_aplicavel"] else "n.a."] for r in R])
    grava("tabela9_atividade_teorica.csv",
          ["amostra", "m_Au_aliquota_ug", "A_teorica_Bq", "A_medida_Bq", "medida_sobre_teorica"],
          [[r["rotulo"], fmt(r["m_Au_ug"]), fmt(r["A_teorica"]), fmt(r["Au_A"]),
            fmt(r["med_sobre_teo"], 2)] for r in R])
    grava("tabela10_medias_por_sistema.csv",
          ["sistema", "A_medida_Bq", "razao", "A_corr_Bq", "razao", "A_Na_Bq", "razao"],
          [[s, fmt(v["A_medida"]), fmt(v["razao_A_medida"], 2), fmt(v["A_corr"]),
            fmt(v["razao_A_corr"], 2), fmt(v["A_Na"]) if v["A_Na"] else "n.a.",
            fmt(v["razao_A_Na"], 2) if v["A_Na"] else "n.a."] for s, v in medias.items()])
    with open(os.path.join(args.saida, "resultados_por_amostra.json"), "w", encoding="utf-8") as fh:
        json.dump({"amostras": R, "medias": medias}, fh, ensure_ascii=False, indent=1)

    A = [r["Au_A"] for r in R]; C = [r["A_corr"] for r in R]
    rel = [
        "REANALISE - CORRECAO PELA DILUICAO, TRACADOR Na-24 E ATIVIDADE TEORICA",
        "",
        "Medias por sistema (razao em relacao ao controle):",
        "  %-6s %10s %8s %10s %8s %10s %8s" % ("", "medida", "razao", "corrig.", "razao", "Na-24", "razao"),
    ]
    for s, v in medias.items():
        rel.append("  %-6s %10.1f %8.2f %10.1f %8.2f %10s %8s" % (
            s, v["A_medida"], v["razao_A_medida"], v["A_corr"], v["razao_A_corr"],
            "%.1f" % v["A_Na"] if v["A_Na"] else "n.a.",
            "%.2f" % v["razao_A_Na"] if v["A_Na"] else "n.a."))
    rel += [
        "",
        "Razao entre extremos: medida %.2f -> corrigida %.2f (entre formulacoes: %.2f)" % (
            max(A) / min(A), max(C) / min(C),
            max(r["A_corr"] for r in R if r["sistema"] != "CTRL")
            / min(r["A_corr"] for r in R if r["sistema"] != "CTRL")),
        "Coeficiente de variacao: medida %.0f%% -> corrigida %.0f%%" % (
            100 * np.std(A, ddof=1) / np.mean(A), 100 * np.std(C, ddof=1) / np.mean(C)),
        "Fator de diluicao da gelatina pelo Na-24: %s (nominal 3,5)" % ", ".join(
            "%.2f" % r["fd_pelo_Na"] for r in R if r["sistema"] == "GEL"),
        "",
        "Atividade teorica (fluxo termico %.1e, sem epitermico): controle %.1f Bq; medido/teorico = %.2f" % (
            args.fluxo, ctrl["A_teorica"], ctrl["med_sobre_teo"]),
        "Fluxo efetivo pelo Au-198 do controle: %.2e cm^-2 s^-1" % phi_eff_au,
        "Na-24 do citrato (controle): teorico %.1f Bq, medido %.1f Bq, fator %.2f" % (
            A_na_teo, ctrl["Na_A"], fator_na_citrato),
        "Na-24 do alginato: inclinacao %.2f +- %.2f Bq/mg (teorico %.2f Bq/mg, fator %.2f); intercepto %.0f +- %.0f Bq" % (
            p[1], math.sqrt(cov[1, 1]), A_na_alg_teo, fator_na_alg, p[0], math.sqrt(cov[0, 0])),
        "Fluxo termico efetivo pelo Na-24: %.2e cm^-2 s^-1" % phi_eff_na,
        "",
        "Observacao: sem replicas, as diferencas entre formulacoes nao sao submetidas a",
        "teste de significancia; as incertezas das razoes sao apenas de contagem.",
    ]
    texto = "\n".join(rel)
    print(texto)
    with open(os.path.join(args.saida, "relatorio_reanalise.txt"), "w", encoding="utf-8") as fh:
        fh.write(texto + "\n")


if __name__ == "__main__":
    main()
