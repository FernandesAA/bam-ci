# Especificacao — pre-processamento do BAM 1.2 subsazonal

## Objetivo

Reescrever em Python o pre-processamento que hoje roda em GrADS/Fortran e produz
as **analises** usadas como condicao inicial / forçante do BAM 1.2 subsazonal.
O resultado vive no pacote `bam-ci` (modulo `bamci`).

## Onde roda

| item | valor |
|---|---|
| host de producao | `jaci.cptec.inpe.br` |
| usuario | `ioper` |
| diretorio | `/p/projetos/ioper/models/BAM-1.2-subsaz` |
| espelho local | `~/dsh/Bam-1.2-subsaz` (vazio em 2026-11-07) |

## O que ja existe no pacote

1. `bamci.ctl` — descritor GrADS completo (dset, undef, options, xdef/ydef/zdef/tdef, vars).
2. `bamci.grads` — leitura do binario float32 (sequencial e um arquivo por tempo),
   `yrev`, conversao de `undef` para NaN, leitura por (variavel, tempo, nivel).
3. `bamci.vertical` — `p(k) = a(k)*p0 + b(k)*ps`, interpolacao linear em log(p)
   por coluna com extrapolacao constante no topo e na superficie.
4. CLI minima (`bamci ctl`, `bamci dump`) e 14 testes.

## O que falta definir (com o codigo original em maos)

1. **Entradas**: quais analises (ECMWF/NCEP/GSI?), formatos, grades e datas;
   quais variaveis entram e quais sao derivadas (umidade, geopotencial, etc.).
2. **Niveis do modelo**: coeficientes `a`/`b` do hibrido ou sigma puro, `p0`.
3. **Transformacoes**: interpolacao horizontal (se houver mudanca de grade),
   correcoes de superficie, mascaramento, unidades.
4. **Saida**: formato binario do estado inicial aceito pelo BAM (layout, ordem,
   cabecalho, precisao) e nomenclatura dos arquivos.
5. **Especificidades subsazonais**: qual e a parte do pre-processamento que muda
   em relacao ao BAM operacional deterministico.
6. **Criterio de equivalencia**: comparacao bit-a-bit, por tolerancia ou por
   estatisticas dos campos contra a saida Fortran/GrADS.

## Principio de port

Cada etapa do codigo original vira uma funcao pura (arrays entram, arrays saem),
com teste de equivalencia contra a saida do programa antigo. O CLI orquestra as
etapas; nada de estado global implicito.
