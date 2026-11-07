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

## Filtro espectral

A producao aplica o filtro de Hoskins (`SmthPerCut=0.12`, `Iter=10`):

```
w(n) = exp(cx * [n(n+1)]^2) ^ Iter,   cx = -ln(red)/rm^2,
rm = 2M(2M+1), rn = (M-1)M, red = per_cut^(-rm^2/(rn^2*Iter))
```
