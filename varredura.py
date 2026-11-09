#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Varredura de picos em espectros HPGe: localiza TODOS os picos
estatisticamente significativos e tenta identifica-los contra uma
biblioteca de linhas gama de ativacao, de fundo natural e de raios X.

Uso:
    python3 varredura.py --dir pasta --padrao "*.Spe" [--sigma 4]
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import efic_hpge as EF                                        # noqa: E402
from atividade_lote import EspectroChn                        # noqa: E402


# ---------------------------------------------------------------------------
# Biblioteca de identificacao: (energia keV, rotulo, origem)
# ---------------------------------------------------------------------------
LINHAS = [
    # --- raios X e aniquilacao ---
    (66.99, "Hg Ka2 (CI do Au-198)", "raio X"),
    (68.89, "Hg Ka1 (CI do Au-198)", "raio X"),
    (70.82, "Hg Ka1'", "raio X"),
    (77.35, "Au Ka1 (fluorescencia)", "raio X"),
    (80.13, "Hg Kb1", "raio X"),
    (82.5,  "Hg Kb2", "raio X"),
    (72.8,  "Pb Ka2 (blindagem)", "raio X"),
    (75.0,  "Pb Ka1 (blindagem)", "raio X"),
    (84.9,  "Pb Kb1 (blindagem)", "raio X"),
    (511.0, "aniquilacao", "511"),

    # --- produto principal ---
    (411.802, "Au-198", "ativacao"),
    (675.884, "Au-198 (fraca, 0,80%)", "ativacao"),
    (1087.68, "Au-198 (fraca, 0,16%)", "ativacao"),
    (158.379, "Au-199", "ativacao"),
    (208.205, "Au-199", "ativacao"),

    # --- ativacao comum em matriz organica / sais ---
    (1368.626, "Na-24", "ativacao"),
    (2754.007, "Na-24", "ativacao"),
    # Picos de escape do gama de 2754 keV do Na-24. Acima de ~1022 keV a
    # criacao de pares passa a contribuir: se um ou os dois fotons de
    # aniquilacao escapam do cristal, a energia depositada fica 511 ou
    # 1022 keV abaixo da energia do gama.
    (2243.007, "Na-24, escape simples de 2754 keV", "escape"),
    (1732.007, "Na-24, escape duplo de 2754 keV", "escape"),
    (857.626, "Na-24, escape simples de 1368 keV", "escape"),
    (346.626, "Na-24, escape duplo de 1368 keV", "escape"),
    (1642.43, "Cl-38", "ativacao"),
    (2167.40, "Cl-38", "ativacao"),
    (1524.58, "K-42", "ativacao"),
    (846.764, "Mn-56", "ativacao"),
    (1810.72, "Mn-56", "ativacao"),
    (2113.09, "Mn-56", "ativacao"),
    (559.10, "As-76", "ativacao"),
    (657.05, "As-76", "ativacao"),
    (554.35, "Br-82", "ativacao"),
    (619.11, "Br-82", "ativacao"),
    (698.37, "Br-82", "ativacao"),
    (776.52, "Br-82", "ativacao"),
    (827.83, "Br-82", "ativacao"),
    (1044.08, "Br-82", "ativacao"),
    (1317.47, "Br-82", "ativacao"),
    (1474.88, "Br-82", "ativacao"),
    (320.08, "Cr-51", "ativacao"),
    (1115.54, "Zn-65", "ativacao"),
    (1099.25, "Fe-59", "ativacao"),
    (1291.60, "Fe-59", "ativacao"),
    (1173.23, "Co-60", "ativacao"),
    (1332.49, "Co-60", "ativacao"),
    (889.28, "Sc-46", "ativacao"),
    (1120.55, "Sc-46", "ativacao"),
    (1596.21, "La-140", "ativacao"),
    (487.02, "La-140", "ativacao"),
    (564.24, "Sb-122", "ativacao"),
    (657.76, "Ag-110m", "ativacao"),
    (884.68, "Ag-110m", "ativacao"),
    (442.90, "I-128", "ativacao"),
    (1297.09, "Ar-41", "ativacao"),
    (159.38, "Ca-47 / Ti-47", "ativacao"),
    (1297.09, "Ca-47", "ativacao"),
    (765.79, "Nb-95 / Mo-99", "ativacao"),
    (140.51, "Tc-99m (de Mo-99)", "ativacao"),
    (228.19, "Mo-99", "ativacao"),
    (739.50, "Mo-99", "ativacao"),

    # --- fundo natural ---
    (1460.82, "K-40 (fundo)", "fundo"),
    (2614.51, "Tl-208 (fundo)", "fundo"),
    (583.19, "Tl-208 (fundo)", "fundo"),
    # 510.77 keV do Tl-208 omitida: nao se distingue de 511 keV (aniquilacao)
    # nesta resolucao, e sem 583 e 2614 keV a atribuicao ao Tl-208 nao se sustenta
    (860.56, "Tl-208 (fundo)", "fundo"),
    (238.63, "Pb-212 (fundo)", "fundo"),
    (300.09, "Pb-212 (fundo)", "fundo"),
    (727.33, "Bi-212 (fundo)", "fundo"),
    (295.22, "Pb-214 (fundo)", "fundo"),
    (351.93, "Pb-214 (fundo)", "fundo"),
    (609.31, "Bi-214 (fundo)", "fundo"),
    (768.36, "Bi-214 (fundo)", "fundo"),
    (934.06, "Bi-214 (fundo)", "fundo"),
    (1120.29, "Bi-214 (fundo)", "fundo"),
    (1238.11, "Bi-214 (fundo)", "fundo"),
    (1377.67, "Bi-214 (fundo)", "fundo"),
    (1764.49, "Bi-214 (fundo)", "fundo"),
    (2204.21, "Bi-214 (fundo)", "fundo"),
    (911.20, "Ac-228 (fundo)", "fundo"),
    (968.97, "Ac-228 (fundo)", "fundo"),
    (338.32, "Ac-228 (fundo)", "fundo"),
    (209.25, "Ac-228 (fundo)", "fundo"),
    (186.21, "Ra-226 / U-235 (fundo)", "fundo"),
    (63.29, "Th-234 (fundo)", "fundo"),
    (92.6, "Th-234 (fundo)", "fundo"),
    (477.60, "Be-7 (cosmogenico)", "fundo"),
    (661.657, "Cs-137 (contaminacao)", "fundo"),
    (1173.0, "Ge(n,g) / captura", "fundo"),
    (139.68, "Ge-75m (captura no proprio Ge)", "fundo"),
    (198.39, "Ge-71m (captura no proprio Ge)", "fundo"),
    (596.0, "Ge-74(n,n') (fundo rapido)", "fundo"),
    (691.0, "Ge-72(n,n') (fundo rapido)", "fundo"),
]


def identifica(energia, tol=1.2):
    """Devolve os candidatos da biblioteca dentro de 'tol' keV."""
    cands = [(abs(energia - e), rot, org, e) for (e, rot, org) in LINHAS
             if abs(energia - e) <= tol]
    cands.sort()
    return cands


# ---------------------------------------------------------------------------
# Busca de picos
# ---------------------------------------------------------------------------
def busca_picos(esp, n_sigma=4.0, ch_min=60, ch_max=None, offset_keV=0.0):
    """Varre o espectro procurando excessos locais significativos.

    Para cada canal, compara a soma numa janela central de +-0.6 FWHM com o
    fundo estimado por duas janelas laterais de mesma largura. Marca onde o
    excesso supera n_sigma vezes o desvio de Poisson. Depois agrupa canais
    vizinhos e ajusta uma gaussiana em cada grupo.
    """
    c = esp.contagens
    n = len(c)
    ch_max = min(ch_max or n - 60, n - 60)

    soma = np.concatenate([[0.0], np.cumsum(c)])

    def janela(a, b):
        a = np.clip(a, 0, n)
        b = np.clip(b, 0, n)
        return soma[b] - soma[a]

    chs = np.arange(ch_min, ch_max)
    larg = np.array([max(2, int(round(0.6 * esp.fwhm_canais(float(k)))))
                     for k in chs])

    centro = janela(chs - larg, chs + larg + 1)
    esq = janela(chs - 3 * larg, chs - larg)
    dire = janela(chs + larg + 1, chs + 3 * larg + 1)
    n_lat = 4 * larg
    fundo = (esq + dire) * (2 * larg + 1) / np.maximum(n_lat, 1)

    liquido = centro - fundo
    var = centro + fundo * ((2 * larg + 1) / np.maximum(n_lat, 1))
    sig = liquido / np.sqrt(np.maximum(var, 1.0))

    marcados = chs[sig >= n_sigma]
    if len(marcados) == 0:
        return []

    # Agrupa canais marcados que estejam proximos
    grupos, atual = [], [marcados[0]]
    for k in marcados[1:]:
        if k - atual[-1] <= max(4, int(esp.fwhm_canais(float(k)))):
            atual.append(k)
        else:
            grupos.append(atual)
            atual = [k]
    grupos.append(atual)

    achados = []
    for g in grupos:
        k0 = g[int(np.argmax([c[k] for k in g]))]
        e_nom = float(esp.energia(k0))
        r = EF.ajusta_multipleto(esp, [(e_nom, 100.0, 0.0)])[0]
        if not r.ok or r.area < 1:
            continue
        i = list(chs).index(g[int(len(g) // 2)]) if g[len(g) // 2] in chs else 0
        achados.append({
            "canal": r.canal_centroide,
            "energia": r.energia_centroide + offset_keV,
            "energia_bruta": r.energia_centroide,
            "fwhm_keV": r.fwhm_keV,
            "area": r.area,
            "u_area": r.u_area,
            "sig": float(sig[i]),
            "cps": r.area / esp.tempo_vivo,
        })

    # remove duplicatas proximas
    achados.sort(key=lambda a: a["energia"])
    limpos = []
    for a in achados:
        if limpos and abs(a["energia"] - limpos[-1]["energia"]) < 0.8:
            if a["area"] > limpos[-1]["area"]:
                limpos[-1] = a
            continue
        limpos.append(a)
    return limpos


def main(argv=None):
    ap = argparse.ArgumentParser(description="Varredura de picos em espectros")
    ap.add_argument("--dir", required=True)
    ap.add_argument("--padrao", default="*.Spe")
    ap.add_argument("--sigma", type=float, default=4.0)
    ap.add_argument("--offset", type=float, default=None,
                    help="correcao de energia em keV; por padrao, calculada "
                         "pelo desvio do pico do Au-198")
    ap.add_argument("--saida", default="varredura.csv")
    args = ap.parse_args(argv)

    caminhos = sorted(glob.glob(os.path.join(args.dir, args.padrao)))
    todos, por_energia = [], {}

    for caminho in caminhos:
        esp = (EspectroChn(caminho) if caminho.lower().endswith(".chn")
               else EF.Espectro(caminho))
        rotulo = os.path.splitext(os.path.basename(caminho))[0]

        # Recalibracao em energia usando linhas conhecidas presentes no
        # proprio espectro. Com DOIS ancoras bem separados (Au-198 em 411,8 e
        # Na-24 em 1368,6) a correcao e linear e vale em toda a faixa; com um
        # so, vira um deslocamento constante que subcorrige o extremo alto.
        ancoras = []
        for e_ref in (411.8022, 1368.626):
            r = EF.ajusta_multipleto(esp, [(e_ref, 100.0, 0.0)])[0]
            if r.ok and r.u_area / max(r.area, 1) < 0.15:
                ancoras.append((r.canal_centroide, e_ref))

        if args.offset is not None:
            recal = lambda e, ch: e + args.offset          # noqa: E731
            nota = "offset fixo %+.3f keV" % args.offset
        elif len(ancoras) >= 2:
            (ch1, e1), (ch2, e2) = ancoras[0], ancoras[-1]
            incl = (e2 - e1) / (ch2 - ch1)
            const = e1 - incl * ch1
            recal = lambda e, ch: const + incl * ch        # noqa: E731
            nota = ("2 ancoras (%.1f e %.1f keV)" % (e1, e2))
        elif ancoras:
            ch1, e1 = ancoras[0]
            d = e1 - float(esp.energia(ch1))
            recal = lambda e, ch: e + d                    # noqa: E731
            nota = "1 ancora (%.1f keV), offset %+.3f keV" % (e1, d)
        else:
            recal = lambda e, ch: e                        # noqa: E731
            nota = "sem recalibracao"

        picos = busca_picos(esp, n_sigma=args.sigma)
        for p in picos:
            p["energia"] = recal(p["energia_bruta"], p["canal"])
        print("=" * 78)
        print("%s   [%s]   t_vivo %.0f s   recalibracao: %s"
              % (rotulo, getattr(esp, "spec_id", ""), esp.tempo_vivo, nota))
        print("  %9s %7s %9s %9s %6s  %s"
              % ("E (keV)", "FWHM", "area", "cps", "sig", "identificacao"))
        for p in picos:
            cands = identifica(p["energia"])
            if cands:
                d, rot, org, e_ref = cands[0]
                ident = "%s  (%.2f keV, %+.2f)" % (rot, e_ref, p["energia"] - e_ref)
            else:
                ident = "NAO IDENTIFICADA"
            print("  %9.2f %7.3f %9.1f %9.4f %6.1f  %s"
                  % (p["energia"], p["fwhm_keV"], p["area"], p["cps"],
                     p["sig"], ident))
            p["amostra"] = rotulo
            p["identificacao"] = ident
            p["origem"] = cands[0][2] if cands else "?"
            todos.append(p)
            chave = round(p["energia"] / 1.5) * 1.5
            por_energia.setdefault(chave, []).append(rotulo)
        print("  total: %d picos" % len(picos))

    print("\n" + "=" * 78)
    print("LINHAS PRESENTES EM MAIS DE UMA AMOSTRA")
    print("=" * 78)
    for chave in sorted(por_energia):
        amostras = sorted(set(por_energia[chave]))
        if len(amostras) >= 2:
            cands = identifica(chave, tol=2.0)
            rot = cands[0][1] if cands else "NAO IDENTIFICADA"
            print("  %8.1f keV  em %2d/%d amostras  %s"
                  % (chave, len(amostras), len(caminhos), rot))

    with open(args.saida, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "amostra", "energia", "energia_bruta", "canal", "fwhm_keV",
            "area", "u_area", "cps", "sig", "identificacao", "origem"],
            extrasaction="ignore")
        w.writeheader()
        w.writerows(todos)
    print("\nCSV: %s" % args.saida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
