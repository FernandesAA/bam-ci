# bam-ci

Pre-processamento do **BAM 1.2 subsazonal** (modelo global do CPTEC/INPE) em Python.

O pacote cobre as duas etapas que hoje rodam em GrADS + Fortran:

1. **Condicao inicial espectral** a partir de analises GRIB (ERA5, FNL, GFS);
2. **SST semanal** no formato que o modelo le (OISST v2.1 ou RTG/NCEP).

```
GRIB (ERA5 / FNL / GFS)
   -> grade gaussiana TQ126 (384 x 192)
   -> NewPs + interpolacao vertical sigma L42
   -> derivadas espectrais (divergencia / vorticidade)
   -> FFT + Legendre (filtro de Hoskins)
   -> GANL<CASE><AAAAMMDDHH>S.unf.TQ0126L042   (11.055.828 B)
```

## Indice

- [Instalacao](instalacao.md)
- [Uso](uso.md) — CLI e submissao PBS
- [Formatos de arquivo](formatos.md) — IC e SST
- [Validacao](validacao.md) — numeros medidos contra a producao
- [API](api.md) — modulos e funcoes principais
- [Problemas conhecidos](troubleshooting.md)
- [Especificacoes do legado](SPEC.md) — extraidas do codigo Fortran/GrADS

## Estado

| etapa | situacao |
|---|---|
| condicao inicial (TQ126L42, sigma) | validada contra a producao (topografia corr 0.994, ps rms 9,4 hPa, T rms ~1 K) |
| SST semanal (RTG/OISST) | validada (mascara identica, SST rms 0.58 K, corr 0.9990) |
| testes | 89, sem dependencia de rede |

## Licenca

MIT. Ver [LICENSE](../LICENSE) e [CITATION.cff](../CITATION.cff).
