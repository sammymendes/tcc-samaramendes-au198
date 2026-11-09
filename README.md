# Ativação neutrônica de nanopartículas de ouro dispersas em biopolímeros

Rotinas computacionais em Python e dados do Trabalho de Conclusão de Curso
**"Ativação neutrônica de nanopartículas de ouro dispersas em biopolímeros:
produção de Au-198 e medida de atividade por espectrometria gama"**, de
Samara Mendes Matos (Bacharelado em Física Médica, Instituto de Física, UFRJ, 2026).

O trabalho mede, por espectrometria gama com detector HPGe, a atividade de Au-198
produzida no canal J9 do reator Argonauta (IEN/CNEN-RJ) em suspensões de
nanopartículas de ouro (AuNPs) preparadas com alginato, quitosana e gelatina, e
mostra que a maior parte das diferenças de atividade entre as formulações vem da
diluição da suspensão no preparo.

## Estrutura

```
tcc-au198-ativacao/
├── src/
│   ├── ler_chn.py               leitor de espectros ORTEC .Chn (binário), exporta CSV
│   ├── efic_hpge.py             curva de eficiência com fonte de Eu-152 (.Spe) → Figura 26, Tabelas 4 e 5
│   ├── atividade_lote.py        atividade de Au-198 ou Na-24 por amostra, em lote → Tabela 7 (medida), Tabela 8
│   ├── varredura.py             varredura de picos a 4σ e identificação de radionuclídeos → Tabela 6
│   ├── uvvis.py                 parâmetros da banda plasmônica a partir da planilha da leitora → Tabela 2
│   ├── conferir.py              confere se os dados brutos reproduzem as tabelas de dados/
│   ├── reanalise_diluicao.py    correção pela diluição, traçador Na-24 e atividade teórica → Tabelas 7 a 10
│   ├── figuras_reanalise.py     Figuras 27 e 29
│   ├── fig_uvvis.py             Figura 25 (espectros UV-Vis)
│   └── fig_na24.py              Figura 28 (Na-24 e fundo sob o pico do Au-198)
├── dados/
│   ├── medidas_10_amostras.csv  tabela consolidada por amostra (Au-198, Na-24, UV-Vis)
│   ├── lote_Na24_atividades.csv saída de atividade_lote.py para o Na-24
│   ├── uvvis_parametros.csv     parâmetros UV-Vis usados no TCC
│   ├── formulacoes.csv          preparo: rótulo, sistema, massa de biopolímero, volume final, f_d
│   ├── amostras.csv             correspondência arquivo de espectro ↔ amostra
│   └── brutos/                  dados brutos originais (ver LEIA-ME.txt)
│       ├── Eu152_inferior.Spe   fonte de Eu-152 na posição inferior (1200 s)
│       ├── Eu152_superior.Spe   fonte de Eu-152 na posição superior (1800 s)
│       ├── amostras/            Amostra_01.Spe … Amostra_10.Spe (900 s cada)
│       └── Resultados_UV-VIS.xlsx  espectros UV-Vis (Multiskan SkyHigh, placa µDrop Duo, 10 µL, caminho óptico de 0,5 mm)
├── tests/teste_sintetico.py     teste de ponta a ponta com dados sintéticos
├── resultados/                  saídas geradas
├── executar_tudo.sh             reproduz todo o tratamento a partir dos dados brutos
├── requirements.txt
└── LICENSE
```

## Fluxo de tratamento

| Etapa | Script | Entrada | Saída | Seção do TCC |
|---|---|---|---|---|
| 1. Leitura dos espectros | `ler_chn.py` (e leitores internos de `.Spe`/`.Chn`) | `.Chn`/`.Spe` | contagens, tempos, calibração em energia | 5.6 |
| 2. Calibração em eficiência | `efic_hpge.py` | 2 × `.Spe` de Eu-152 | pontos ε(E), curva ln ε = Σ aᵢ [ln(E/E₀)]ⁱ, relatório, gráfico | 5.6.1, 6.2 |
| 3. Atividade por amostra | `atividade_lote.py` | 10 espectros + curva | A(t_c), A(t₀), AMD, orçamento de incerteza | 5.6.2, 6.4, 6.5 |
| 4. Identificação | `varredura.py` | 10 espectros | linhas a 4σ, picos de escape, recalibração por dois pontos | 6.3 |
| 5. UV-Vis | `uvvis.py` | planilha da leitora | λ_máx, A_RPS, pico/450, pico/600, A(650)/A(λ_máx) | 5.4, 6.1 |
| 5b. Conferência | `conferir.py` | `resultados/` e `dados/` | OK/FALHA por tabela | — |
| 6. Reanálise | `reanalise_diluicao.py` | `dados/*.csv` | Tabelas 7 a 10, fluxo efetivo | 5.7, 6.4 a 6.7 |
| 7. Figuras da reanálise | `figuras_reanalise.py` | saída da etapa 6 | Figuras 27 e 29 | 6.4, 6.7 |
| 8. Figuras 25 e 28 | `fig_uvvis.py`, `fig_na24.py` | planilha UV-Vis; saídas da etapa 3 | Figura 25; Figura 28 e regressão do fundo | 6.1, 6.5 |

As etapas 2 a 5 partem dos dados brutos em `dados/brutos/`. A conferência compara
o resultado com as tabelas arquivadas em `dados/` (as usadas no TCC): a curva de
eficiência, as atividades de Au-198 e Na-24 e os parâmetros UV-Vis saem idênticos.
As etapas 6 e 7 dependem apenas dos arquivos em `dados/` e reproduzem os números
do Capítulo 6.

## Como executar

Requer Python 3.9 ou superior.

```bash
pip install -r requirements.txt

# todo o tratamento, a partir dos dados brutos em dados/brutos/
./executar_tudo.sh

# só a reanálise e as figuras (a partir dos dados consolidados já incluídos)
python3 src/reanalise_diluicao.py
python3 src/figuras_reanalise.py

# teste das rotinas com dados sintéticos de valores conhecidos
python3 tests/teste_sintetico.py
```

Parâmetros usados no TCC: fonte de Eu-152 de 4798 Bq em 20/06/2018; fim da
irradiação em 27/10/2026 às 11h30; irradiação de 120 min; curva de grau 3 pela
média das duas posições; componente sistemática de 2% na eficiência
(`--u-eps-extra 2`). A incerteza do certificado do Eu-152 não foi incluída
(`--u-atividade-eu 0`), porque o certificado não foi localizado.

`reanalise_diluicao.py` aceita `--massa-haucl4-mg` (padrão: 5,0 mg, valor nominal)
e `--fluxo` (padrão: 10⁹ cm⁻²·s⁻¹). Ao informar a massa efetivamente pesada, a
atividade teórica é recalculada.

### Observações sobre as rotinas

- `uvvis.py` contém duas análises. A função `parametros_tcc()` é a usada no TCC
  (seção 5.4). A função `analisa()` é a versão inicial, mantida como registro e
  impressa com `--original`; o diâmetro de Haiss que ela calcula não é aplicável
  a estas amostras (seção 6.1.1). No método do TCC, a linha de base (média entre
  900 e 960 nm) é tomada sobre o espectro já suavizado, como o máximo da banda.
- Em `varredura.py`, a linha de 510,77 keV do Tl-208 foi retirada da biblioteca:
  ela se confundia com a aniquilação (511 keV), e nenhuma outra linha do Tl-208
  aparece nos espectros.
- `fig_uvvis.py` e `fig_na24.py` foram reescritos para o repositório, com os
  rótulos corrigidos das amostras; os originais não foram preservados.

## Modelo de medida

- Atividade no início da contagem: A(t_c) = A_R·C_d / (ε·I_γ·T_v), com
  C_d = λT_r / (1 − e^(−λT_r)).
- Atividade no fim da irradiação: A(t₀) = A(t_c)·e^(λΔt).
- Correção pela diluição: A_corr = f_d·A(t₀), com f_d = V_final / 10 mL.
- Traçador: R = A(Au-198)/A(Na-24); A_Na = R·A_controle(Na-24). Não se aplica ao
  alginato, que contém sódio próprio.
- Atividade teórica: A = N(σ₀φ_th + I₀φ_epi)(1 − e^(−λt_irr)), com σ₀ = 98,65 b e
  I₀ ≈ 1550 b para o Au-197.

Dados nucleares: NuDat 3 (NNDC). Seções de choque: Mughabghab (2006).

## Limitações

Cada formulação foi preparada e medida uma única vez; as incertezas das razões
entre amostras são apenas de contagem e não incluem a variabilidade de preparo,
pipetagem, posicionamento no reator e geometria de contagem. Ver a seção 6.8 do TCC.

## Como citar

MATOS, Samara Mendes. **Ativação neutrônica de nanopartículas de ouro dispersas em
biopolímeros**: produção de Au-198 e medida de atividade por espectrometria gama.
2026. Trabalho de Conclusão de Curso (Bacharelado em Física Médica) – Instituto de
Física, Universidade Federal do Rio de Janeiro, Rio de Janeiro, 2026.

Repositório: <https://github.com/sammymendes/tcc-samaramendes-au198>
