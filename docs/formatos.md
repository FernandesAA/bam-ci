# Formatos de arquivo

## Condicao inicial — GANL<CASE><AAAAMMDDHH>S.unf.TQ0126L042

Fortran **UNFORMATTED sequencial**, big endian, marcadores de 4 bytes por registro,
`imax = 384`, `jmax = 192`, `kmax = 42`, `ncoeff = (126+1)(126+2) = 16.256`.

| registro | conteudo |
|---|---|
| 1 | `ForecastDay` (i4), `TimeOfDay` (r4), `DateInp(4)`, `DateCur(4)`, `SigI(43)`, `SigL(42)` |
| 2 | topografia (coeficientes) |
| 3 | ln(ps/10) — ps em hPa, ou seja ln(cbar) |
| 4 .. K+3 | temperatura virtual (K) |
| K+4 .. 3K+3 | divergencia e vorticidade (s-1) alternadas por nivel |
| 3K+4 .. 4K+3 | umidade especifica (kg/kg) |

Ordenacao das datas: **hora, mes, dia, ano**. ForecastDay=0 e TimeOfDay=0 na CI.

### Ordem dos coeficientes

Pares (cosseno, seno) com **m externo** e n interno:

```
p(m,n) = m(M+1) - m(m-1)/2 + (n-m)      real em 2p, imaginario em 2p+1
```

Normalizacao 4pi: `pln(0,0) = 1/sqrt(2)`, `pln = Pbar/sqrt(2)` (m=0) e `Pbar/2` (m>0),
com `Cte = 1/sqrt(2)` ligando coeficiente e media (`c(0,0) = sqrt(2) * media`).
O seno de m=0 e sempre nulo.

### Grade

Gaussiana com **192 latitudes** (89,28 -> -89,28) e **384 longitudes**, na convencao
do modelo: **longitudes decrescentes** (`[..., ::-1]`).

## SST semanal — SSTWeekly<AAAAMMDD>.G00192

Fortran **UNFORMATTED DIRECT**, recl = 192*384*4, sem marcadores, 2 registros
(`589.824` bytes):

| registro | conteudo |
|---|---|
| 1 | mascara oceano/terra: **+1 oceano, -1 terra** |
| 2 | SST em K, **-999** sobre terra |

A grade e a mesma 192x384, porem com **longitudes crescentes** (0 -> 360) — diferente
da CI. A mascara de terra tem 24.469 pontos.

## DeltaSigma.LNNN — espessuras sigma

Texto puro, formato Fortran `(5F9.6)` (5 valores por linha, campo de 9 colunas):
NNN **espessuras de camada** em coordenada sigma, uma por camada, ordenadas da
superficie para o topo. Nao ha cabecalho nem contagem de valores — o `NNN` do
nome e a unica indicacao do numero de camadas, e o pipeline exige que ele seja
igual ao `kmax` pedido.

Propriedades (verificadas em `tests/test_deltasigma_files.py`):

- a soma das NNN espessuras e 1 (tolerancia 1e-4 nos arquivos do pacote);
- todos os valores sao positivos;
- `bamci.levels.read_delta_sigma` normaliza a soma quando ela se afasta de 1
  por mais de 1e-5, e `sigma_interfaces` devolve as interfaces
  `0, d1, d1+d2, ..., 1`.

```bash
bamci run --case GFS --date 2026110500 -r TQ0213L64 \
          --grib gfs.grib --delta-sigma DeltaSigma/DeltaSigma.L064 --out-dir model/datain
```

### Procedencia (1.1.4)

Copia **byte-identica** de
`/p/projetos/ioper/fix/sazonal/bam_SeasonalClima-1.2_oper/pre/datain` (jaci),
mantendo os nomes originais: `DeltaSigma.L009`, `.L018`, `.L028`, `.L042`,
`.L064`, `.L096` e `.L160` (83 a 1472 bytes). Nada foi reescrito ou arredondado.

**Nota sobre o L064**: na origem `DeltaSigma.L064` soma **0.999976** — residuo de
2.4e-5 do arquivo operacional, que e identico ao `DeltaSigma.L064_NEW-42to64` da
mesma pasta. E o unico dos sete acima de 1e-6; como o leitor normaliza, a CI
gerada nao muda. O `DeltaSigma.L064_ORIG` daquela pasta (soma exata 1) **nao**
foi copiado; os demais `DeltaSigma.L064_*` (`_NEW`, `_NEW-28to64`) tambem nao.

## Filtro espectral

A producao aplica o filtro de Hoskins (`SmthPerCut=0.12`, `Iter=10`):

```
w(n) = exp(cx * [n(n+1)]^2) ^ Iter,   cx = -ln(red)/rm^2,
rm = 2M(2M+1), rn = (M-1)M, red = per_cut^(-rm^2/(rn^2*Iter))
```
