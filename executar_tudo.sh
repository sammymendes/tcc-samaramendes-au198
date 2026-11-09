#!/usr/bin/env bash
# Reproduz o tratamento de dados do TCC, na ordem do Capitulo 6.
# Requer os dados brutos em dados/brutos/ (ver dados/brutos/LEIA-ME.txt).
set -euo pipefail
cd "$(dirname "$0")"
B=dados/brutos
mkdir -p resultados

echo "== 1. Curva de eficiencia (Eu-152, duas posicoes no microtubo)"
python3 src/efic_hpge.py \
  --spe "$B/Eu152_inferior.Spe" --rotulo "posicao inferior (menor distancia)" \
  --spe "$B/Eu152_superior.Spe" --rotulo "posicao superior (maior distancia)" \
  --nuclideo Eu-152 --atividade 4798 --data-ref 2018-06-20 \
  --combinar media --grau 3 --varrer-grau \
  --prefixo resultados/eficiencia_Eu152

for NUC in Au-198 Na-24; do
  echo "== 2. Atividade do lote: $NUC"
  python3 src/atividade_lote.py --dir "$B/amostras" --padrao "*.Spe" \
    --nuclideo "$NUC" --fim-irradiacao "2026-10-27 11:30" --tempo-irradiacao 120 \
    --spe-eu "$B/Eu152_inferior.Spe" --spe-eu "$B/Eu152_superior.Spe" \
    --u-eps-extra 2 --saida "resultados/lote_${NUC}_atividades.csv"
done

echo "== 3. Varredura de picos e identificacao (4 sigma)"
python3 src/varredura.py --dir "$B/amostras" --padrao "*.Spe" --sigma 4 \
  --saida resultados/varredura_picos.csv

echo "== 4. UV-Vis"
python3 src/uvvis.py --planilha "$B/Resultados_UV-VIS.xlsx" \
  --saida resultados/uvvis_parametros.csv

echo "== 5. Conferencia com as tabelas arquivadas em dados/"
python3 src/conferir.py

echo "== 6. Reanalise (diluicao, Na-24, atividade teorica) e figuras"
python3 src/reanalise_diluicao.py
python3 src/figuras_reanalise.py

echo "== 7. Figuras 25 (UV-Vis) e 28 (Na-24)"
python3 src/fig_uvvis.py
python3 src/fig_na24.py
