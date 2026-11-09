#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conferir.py - Confere se o processamento dos dados brutos (pasta resultados/)
reproduz as tabelas arquivadas em dados/, que sao as usadas no TCC.

Verifica:
  1. coeficientes da curva de eficiencia (Tabela 5 do TCC);
  2. areas e atividades de Au-198 (dados/medidas_10_amostras.csv);
  3. areas e atividades de Na-24 (dados/lote_Na24_atividades.csv);
  4. parametros UV-Vis (dados/uvvis_parametros.csv).

Uso:  python3 src/conferir.py      (codigo de saida 1 se algo nao conferir)
"""
import csv
import re
import sys

COEF_TCC = [-5.44871, -1.04138, -0.0341442, 0.0526425]
falhas = []


def ler(caminho):
    with open(caminho, encoding="utf-8") as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#") and l.strip()))


def confere(nome, ref, novo, tol_abs=0.0, tol_rel=0.0):
    ok = abs(ref - novo) <= max(tol_abs, tol_rel * abs(ref))
    if not ok:
        falhas.append(f"{nome}: arquivado {ref:.6g}, recalculado {novo:.6g}")
    return ok


# 1. Curva de eficiencia
txt = open("resultados/eficiencia_Eu152_relatorio.txt", encoding="utf-8").read()
bloco = txt[txt.index("AJUSTE DA CURVA MEDIA"):]
coef = [float(x) for x in re.findall(r"a\d = ([-+0-9.e]+)", bloco)[:4]]
n1 = sum(confere(f"eficiencia a{i}", c, x, tol_rel=1e-5)
         for i, (c, x) in enumerate(zip(COEF_TCC, coef)))
print(f"[{'OK' if n1 == 4 else 'FALHA'}] curva de eficiencia: {n1}/4 coeficientes")

# 2. Au-198 (tabela consolidada arredondada a 0,1)
med = {r["amostra"]: r for r in ler("dados/medidas_10_amostras.csv")}
au = {r["id"].replace("Amostra ", "", 1).strip(): r
      for r in ler("resultados/lote_Au-198_atividades.csv")}
n2 = 0
for k, m in med.items():
    a = au[k]
    n2 += all([
        confere(f"Au-198 {k} area", float(m["Au198_area"]), float(a["area"]), tol_abs=0.051),
        confere(f"Au-198 {k} A(t0)", float(m["Au198_A_t0_Bq"]), float(a["A_t0_Bq"]), tol_abs=0.051),
        confere(f"Au-198 {k} u(A)", float(m["Au198_u_A_t0_Bq"]), float(a["u_A_t0_Bq"]), tol_abs=0.051),
    ])
print(f"[{'OK' if n2 == len(med) else 'FALHA'}] Au-198: {n2}/{len(med)} amostras")

# 3. Na-24 (arquivo completo, sem arredondamento)
na_ref = {r["amostra"]: r for r in ler("dados/lote_Na24_atividades.csv")}
na_new = {r["amostra"]: r for r in ler("resultados/lote_Na-24_atividades.csv")}
n3 = 0
for k, r in na_ref.items():
    n3 += all(confere(f"Na-24 {k} {c}", float(r[c]), float(na_new[k][c]), tol_rel=1e-6)
              for c in ("area", "u_area", "A_t0_Bq", "u_A_t0_Bq"))
print(f"[{'OK' if n3 == len(na_ref) else 'FALHA'}] Na-24: {n3}/{len(na_ref)} amostras")

# 4. UV-Vis (arquivado com 3 a 4 casas decimais)
uv_ref = {r["poco"]: r for r in ler("dados/uvvis_parametros.csv")}
uv_new = {re.search(r"\((\w+)\)", r["amostra"]).group(1): r
          for r in ler("resultados/uvvis_parametros.csv")}
cols = ("lambda_max_nm", "A_RPS", "A250", "razao_pico_450", "razao_pico_600",
        "indice_agregacao_650")
n4 = 0
for p, r in uv_ref.items():
    n4 += all(confere(f"UV-Vis {p} {c}", float(r[c]), float(uv_new[p][c]), tol_abs=6e-4)
              for c in cols)
print(f"[{'OK' if n4 == len(uv_ref) else 'FALHA'}] UV-Vis: {n4}/{len(uv_ref)} pocos")

if falhas:
    print("\nDiferencas encontradas:")
    for f in falhas:
        print("  " + f)
    sys.exit(1)
print("\nOs dados brutos reproduzem as tabelas arquivadas em dados/.")
