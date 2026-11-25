# bam-ci

Pre-processamento do **BAM 1.2 subsazonal** (modelo global do CPTEC/INPE) em Python.

O pacote cobre as etapas que hoje rodam em GrADS + Fortran:

1. **Condicao inicial espectral** a partir de analises GRIB (ERA5, FNL, GFS);
2. **SST semanal** no formato que o modelo le (OISST v2.1 ou RTG/NCEP);
3. **Ensaio MB09** de perturbacao de condicoes iniciais por EOF (15 membros).

```
GRIB (ERA5 / FNL / GFS)
   -> grade gaussiana da truncacao pedida (T126 -> 384x192, T213 -> 640x320, ...)
   -> NewPs + interpolacao vertical sigma L<kmax>
   -> derivadas espectrais (divergencia / vorticidade)
   -> FFT + Legendre (filtro de Hoskins)
   -> GANL<CASE><AAAAMMDDHH>S.unf.T<Q|L><mend>L<kmax>
```

A resolucao e um argumento (`-r TQ0213L64`, `--mend 254 --kmax 42`); TQ126L042
continua sendo o default.

## Indice

- [Instalacao](instalacao.md)
- [Uso](uso.md) — CLI e submissao PBS
- [Formatos de arquivo](formatos.md) — IC, SST e DeltaSigma
- [Validacao](validacao.md) — numeros medidos contra a producao
- [MB09: especificacao do ensaio](mb09.md)
- [Artefatos do ensaio 2026-11-19](artefatos_ensaio_20261119.md)
- [Literatura: perturbacao de condicoes iniciais](literatura-perturbacoes.md)
- [API](api.md) — modulos e funcoes principais
- [Problemas conhecidos](troubleshooting.md)
- [Especificacoes do legado](SPEC.md) — extraidas do codigo Fortran/GrADS

## Estado

| etapa | situacao |
|---|---|
| condicao inicial (TQ126L42, sigma) | validada contra a producao (topografia corr 0.994, ps rms 9,4 hPa, T rms ~1 K) |
| outras resolucoes (T213L64, T254L42, TL95L42) | pipeline ponta a ponta com tamanho de arquivo exato |
| SST semanal (RTG/OISST) | validada (mascara identica, SST rms 0.58 K, corr 0.9990) |
| DeltaSigma (espessuras sigma) | 7 arquivos em `DeltaSigma/` (L009..L160) acompanham o pacote, soma 1 |
| ensaio MB09 (EOF) | nucleo + driver BAM + metricas + shell 0-9; presets oper/ext1/trop/ext2/etsa, 15 membros |
| testes | 387 (+1 skip), sem dependencia de rede |

## Licenca

MIT. Ver [LICENSE](../LICENSE) e [CITATION.cff](../CITATION.cff).