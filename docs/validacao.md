# Validacao contra a producao do CPTEC/INPE

## Condicao inicial — 2026-11-05 00Z (GFS)

Comparacao com `GANL01R2026110500S.unf.TQ0126L042`:

| grandeza | resultado |
|---|---|
| tamanho / registros | identicos (11.055.828 B, 3 + 4x42) |
| cabecalho sigma | identico (a0 = 1,0; a-1 = 0,0) |
| correlacao da topografia | **0,994** |
| pressao de superficie | rms **9,4 hPa** (~1%) |
| temperatura virtual k=1 / k=21 | rms **1,06 K** / **0,76 K** |
| divergencia / vorticidade | rms 2,7e-5 / 3,7e-5 s-1 |

Perfil por nivel (media global, K):

| k | p (hPa) | bam-ci | producao |
|---|---|---|---|
| 1 | 996,0 | 280,09 | 280,18 |
| 21 | — | 247,08 | 246,95 |
| 41 | 7,65 | 233,72 | 235,35 |
| 42 | 2,04 | **254,16** | 236,20 |

O nivel 42 e o unico discrepante: o GFS tem **254,28 K** em 2 hPa, valor que o
bam-ci reproduz; o arquivo antigo fica 18 K abaixo, num valor compativel com ~5-6 hPa
(indicacao de extrapolacao/indice errado no ultimo nivel do codigo legado).

## SST — 2026-11-05 (RTG/NCEP, dado pontual)

| grandeza | resultado |
|---|---|
| tamanho | identico (589.824 B) |
| mascara oceano/terra | **identica** bit a bit (24.469 pontos de terra) |
| SST | rms **0,58 K**, corr **0,9990**, vies +0,13 K, max 7,6 K |

O residuo e compativel com interpolacao bilinear x interpolacao por area
(`AreaInterpolation.f90` no legado). Teste de orientacao de grade (4 combinacoes)
mostrou que a SST usa longitudes **crescentes** (rms 0,58 K) e nao decrescentes
(rms 2,49 K).

## Como reproduzir

```bash
bamci check   arquivo.unf --mend 126
bamci compare nosso.unf producao.unf --mend 126
```
