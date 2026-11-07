# Especificacao das ENTRADAS — dados meteorologicos do BAM 1.2 subsazonal

Documento derivado por leitura do codigo legado (GrADS + shell + Fortran de consumo).
Espelho local read-only: `reference/GetICnAGCM_subsaeasonal/` (copia de
`/oper/scripts/subsaz/Bruno/GetICnAGCM_subsaeasonal` na jaci).

**Convencao de citacao:** todos os caminhos sao relativos a
`reference/GetICnAGCM_subsaeasonal/`, no formato `arquivo:linha`.
Citacoes de `Sources/` sao do codigo Fortran que **consome** os binarios e por isso
sao a fonte de verdade para o layout (item 4).

Este documento cobre **somente as entradas** (dados meteorologicos). A geracao do
estado inicial espectral e o contrato do modelo BAM ficam em `SPEC.md`.

---

## 0. Resumo executivo

O pre-processamento legado nao le dados meteorologicos diretamente: ele converte
GRIB em **um unico binario GrADS float32 big-endian** por data, com 6 variaveis
(topografia, pressao ao nivel do mar, T, u, v, q), e o Fortran `ICnAGCM` le esse
binario por acesso direto e o interpola para a grade espectral do modelo.

Fluxo comum a todas as fontes:

```
GRIB (por mes ou por data)
   -> gribmap -> indice .gmp                        [] GrADSScripts/GRIBMap.ksh:27-28
   -> .ctl    (descreve o GRIB: xdef/ydef/zdef/tdef) [] GrADSScripts/Set*ctl.gs
   -> Get*ICn.gs: open, display, set gxout fwrite,
      set fwrite -be -> <PREFIXO>ATM<AAAAMMDDHH>.dat [] GrADSScripts/Get*ICn*.gs
   -> Fortran ICnAGCM/SpectralICn.f90: leitura direta
                                     [] Sources/ICnAGCM/SpectralICn.f90:94-148
```

Fontes existentes: **ERA-Interim** (1.5 e 0.75), **ERA-40** (web e interno),
**FNL grib1**, **FNL grib2**, **GFS 0.25**. **ERA5 nao existe no legado** (ver §8).

---

## 1. Tabela geral por fonte

| fonte | Case | grade (Im x Jm) | dxy | Km (T,u,v) | Kq (q) | dir de dados |
|---|---|---|---|---|---|---|
| ERA-Interim 1.5 | `EIT` | 240 x 121 | 1.5 | 37 | 37 | `DataAnl/EIT` |
| ERA-Interim 0.75 | `EIH` | 480 x 241 | 0.75 | 37 | 37 | `DataAnl/EIH` |
| ERA-40 web | `E4W` | 144 x 73 | 2.5 | 23 | 23 | `DataAnl/E4W` |
| ERA-40 interno CPTEC | `E40` | 320 x 161 | 1.125 | 23 | 23 | `DataAnl/E40` |
| FNL grib1 | `FNL` | 360 x 181 | 1.0 | 26 | **21** | `DataAnl/FNL` |
| FNL grib2 | `FNL` | 360 x 181 | 1.0 | 26 | **21** | `DataAnl/FNL` |
| GFS 0.25 | `GFS` | 1440 x 721 | 0.25 | 26 | **21** | `DataAnl/GFS` |

Dimensoes em `Sources/ICnAGCM/InputParameters.f90:161-204` (bloco `SELECT CASE (Case)`).
O diretorio de entrada e montado como `DirMainInp//'DataAnl/'//Case//'/'`
(`Sources/ICnAGCM/InputParameters.f90:212`), coerente com a criacao dos subdiretorios
em `Configure.ksh:26,92-110`.

Nomes de arquivo de saida (`DataInp`), todos lidos pelo Fortran:

| fonte | binario lido | prefixo escrito pelo .gs |
|---|---|---|
| E40 | `ERA40CATM<AAAAMMDDHH>.dat` | `ERA40CATM` (`GrADSScripts/GetERA40CICn.gs:90`) |
| E4W | `ERA40WATM<AAAAMMDDHH>.dat` | `ERA40WATM` (`GrADSScripts/GetERAWebICn.gs:123`) |
| EIT | `ERAINTATM<AAAAMMDDHH>.dat` | `ERAINTATM` (`GrADSScripts/GetERAWebICn.gs:102`) |
| EIH | `ERAINTATM<AAAAMMDDHH>.dat` | `ERAINTATM` (`GrADSScripts/GetERAHWebICn.gs:102`) |
| FNL | `NMCFNLATM<AAAAMMDDHH>.dat` | `NMCFNLATM` (`GrADSScripts/GetFnlICgrib2.gs:74`) |
| GFS | `NMCGFSATM<AAAAMMDDHH>.dat` | `NMCGFSATM` (`GrADSScripts/GetNMCGFSFnlICn.gs:74`) |

---

## 2. ERA-Interim (EIT 1.5 e EIH 0.75), via web/ECMWF

### 2.1 Arquivos GRIB esperados

- **Um arquivo GRIB por variavel e por mes**, nao por data:
  `ERAINTATM<VAR><AAAA><MM>.grb`, com `VAR` em `{TMP, UME, UZN, VMD, OZN}`
  (`GrADSScripts/SetERAIntctl.gs:83`; idem `SetERAInHctl.gs:83`).
- Mais um arquivo de superficie por mes, `ERAINTSFCVAR<AAAA><MM>.grb`
  (`GrADSScripts/SetERAIntctl.gs:111`), e um **invariante** sem data,
  `ERAINTSFCINVARIANT.grb` (`GrADSScripts/SetERAIntctl.gs:137`).
- Caminho fisico: `DirMain/DataAnl/EIT` (EIH vai para `DataAnl/EIH`:
  `GrADSScripts/SetERAInHctl.gs:86,114,140`).
- Horas: 00/06/12/18 UTC — validado em `GrADSScripts/GetERAWebICn.gs:39-41`;
  o driver sempre usa `00z` (`GrADSScripts/roda_EIT.sh:40`).
- Anos aceitos: `1979 <= yy <= 2017` (`GrADSScripts/SetERAIntctl.gs:55`);
  `GetERAHWebICn.gs:23-24` idem.
- Contagem de tempos: `tdef <lm> linear <time> 06hr` com `lm = dias_do_mes*4`
  (`GrADSScripts/SetERAIntctl.gs:102`; `gera_ctl.gs:177-179`).
  O argumento `dtm` e **sempre forcado a 6 h** por um bug de condicao
  (`GrADSScripts/SetERAIntctl.gs:11-13`: `if (dtm != 06 | dtm != 12 | dtm != 24) dtm=06`
  e sempre verdadeiro) — para o port, assumir passo de 6 h.

### 2.2 Como os .ctl sao gerados

Sete descritores por execucao — **cinco de atmosfera (um por variavel)**, um de
superficie e um invariante — escritos por concatenacao de `echo` dentro do .gs, e
um indice `.gmp` por `gribmap`. O nome do .ctl atmosferico carrega o nome da
variavel (`SetERAIntctl.gs:84`) e cada um descreve **uma** variavel (`vars 1`):

| .ctl | dset | vars | referencia |
|---|---|---|---|
| `ERAINTATM<VAR>_<ini>_<fim>_06.ctl` | `^ERAINTATM<VAR><AAAA><MM>.grb` | 1 | `SetERAIntctl.gs:84-105` |
| `ERAINTSFCVAR<ini>_<fim>_06.ctl` | `^ERAINTSFCVAR<AAAA><MM>.grb` | 4 | `SetERAIntctl.gs:112-133` |
| `ERAINTSFCINVARIANT.ctl` | `^ERAINTSFCINVARIANT.grb` | 1 | `SetERAIntctl.gs:138-156` |

Conteudo tipico (EIT), literal de `GrADSScripts/SetERAIntctl.gs:90-105`:

```
dset ^ERAINTATMTMP201001.grb
dtype grib
index ^ERAINTATMTMP_2010010100_2010013118_06.gmp
title ERA Interim at ECMWF
options yrev
undef 9.99e30
xdef 240 linear   0.0 1.5
ydef 121 linear -90.0 1.5
zdef  37 levels 1000 975 ... 3 2 1
tdef 124 linear 00z01jan2010 06hr
vars 1
temp 37 130,100 Temperature [K]
endvars
```

- `options yrev` (linha 94): o GRIB vem ordenado norte->sul e o GrADS inverte para
  `ydef -90 +1.5` (sul->norte). **Isso e o que fixa a orientacao do .dat final** (§7).
- `undef 9.99e30` (linha 95) — valor de ausencia do ECMWF.
- `xdef 240 linear 0.0 1.5` / `ydef 121 linear -90.0 1.5` (linhas 96-97).
  EIH: `480 linear 0.0 0.75` / `241 linear -90.0 0.75` (`SetERAInHctl.gs:96-97`).
- `zdef 37 levels` com a lista **1000 975 950 925 900 875 850 825 800 775 750 700 650
  600 550 500 450 400 350 300 250 225 200 175 150 125 100 70 50 30 20 10 7 5 3 2 1**
  (`SetERAIntctl.gs:98-101`). Ordem: superficie -> topo (z=1 e 1000 hPa).
- Descritor de variavel no formato `nome <n_niveis> <codigo_grib>,<tipo_nivel>`:
  `temp 37 130,100 Temperature [K]` (`SetERAIntctl.gs:104`).
- O indice e gerado por `GrADSScripts/GRIBMap.ksh:21-28`: carrega
  `../Configure.ksh`, faz `cd ${DirData}/${subdir}` e roda
  `${DirGrADS}/gribmap -v -i <ctl>` filtrando linhas `MATCH`/`NOOOO`.

### 2.3 Variaveis lidas do GRIB (lista completa)

Codigos ECMWF (tabela 128), tipo de nivel 100 = isobarico em hPa
(`GrADSScripts/SetERAIntctl.gs:65-69`):

| # | nome no .ctl | codigo GRIB | niveis | unidade origem | unidade destino |
|---|---|---|---|---|---|
| 1 | `temp` | 130 | 37 | K | K |
| 2 | `umes` | 133 (q especifica) | 37 | kg/kg | kg/kg |
| 3 | `uvel` | 131 | 37 | m/s | m/s |
| 4 | `vvel` | 132 | 37 | m/s | m/s |
| 5 | `ozon` | 203 (ozone mass mixing ratio) | 37 | kg/kg | kg/kg |
| 6 | `psnm` | 151 | 0 (superficie) | Pa | **hPa** (/100) |
| 7 | `pslc` | 134 | 0 | Pa | hPa (/100) — so no binario de ozonio |
| 8 | `sstk` | 34 | 0 | K | K |
| 9 | `sice` | 31 | 0 | 0..1 | 0..1 |
| 10 | `gpsf` | 129 (geopotencial) | 0 | m2/s2 | **m** (/9.80665) |

Superficie: `GrADSScripts/SetERAIntctl.gs:129-132`; invariante:
`SetERAIntctl.gs:155`.

### 2.4 Derivadas e correcoes

O ERA-Interim **nao tem nenhuma derivacao de umidade nem de temperatura**: `q` e
lida pronta (codigo 133) e `T` e a temperatura absoluta (codigo 130). As unicas
transformacoes sao conversao de unidade:

- `lpnm='psnm.5/100'` — Pa -> hPa (`GrADSScripts/GetERAWebICn.gs:147`).
- `ltop='gpsf.6(t=1)/9.80665'` — geopotencial -> altura geopotencial em m
  (`GrADSScripts/GetERAWebICn.gs:148`). `(t=1)` e necessario porque o invariante
  tem `tdef 1 linear 12z01jan1989 1mo` (`SetERAIntctl.gs:153`).
- **Nenhuma correcao de vies/ajuste**: `GrADSScripts/corre.gs` e
  `GrADSScripts/corre_fnl.gs` sao apenas buffers de uma linha reescritos pelo driver
  (`roda_EIT.sh:43,56`; `roda_FNL_grib1.sh:43,47`). Nao ha conteudo cientifico neles.
- A temperatura virtual **nao** e calculada no GrADS: o .dat guarda T seca e o
  Fortran aplica `TempOut = TempOut*(1 + cTv*q)` com `cTv = Rv/Rd - 1 = 0.608`
  (`Sources/ICnAGCM/SpectralICn.f90:160`; `Sources/ICnAGCM/InputParameters.f90:341`).

---

## 3. FNL grib1 (NCEP Final Analysis)

### 3.1 Arquivos GRIB

- O driver varre `/scratchin/.../GetICnAGCM/DataAnl/FNL/*.grib1`
  (`GrADSScripts/roda_FNL_grib1.sh:3`) — **um arquivo por data de analise**.
- O `.ctl` aponta para um **template** por data e hora:
  `dset ^fnl_<AAAA><MM>%d2_%h2_00` (`GrADSScripts/SetNMCFnlctl.gs:71,78`) com
  `options yrev template` (`:82`). Ano e mes ficam fixos na string; dia e hora
  variam por tempo (`%d2`, `%h2`).
  > Observacao: o template nao inclui a extensao `.grib1` que o driver usa no
  > `ls`. Divergencia a confirmar com o script de download (ausente do espelho).
- Hora: o driver sempre usa `00z` (`GrADSScripts/roda_FNL_grib1.sh:41`);
  `hh` validado em 00/06/12/18 (`SetNMCFnlctl.gs:44-46`).
- Anos aceitos: `1999 <= yy <= 2015` (`SetNMCFnlctl.gs:56`);
  `GetNMCFnlICn.gs:59` idem.
- `tdef <lm> linear <time> 06hr` com `lm = dias_do_mes*4` (o mesmo bug de `dtm`
  forca 6 h: `SetNMCFnlctl.gs:11-13,122-124`). Na pratica so o tempo 1 casa com o
  arquivo existente.

### 3.2 .ctl

Escrito em `DirMain/DataAnl/FNL` (`GrADSScripts/SetNMCFnlctl.gs:74`), conteudo
literal de `:78-101`:

```
dset ^fnl_200712%d2_%h2_00
dtype grib
index ^fnl_2007120100_2007123118_06.gmp
title NCEP Final Analysis at NCAR
options yrev template
undef -2.56e30
xdef 360 linear   0.0 1.0
ydef 181 linear -90.0 1.0
zdef  26 levels 1000 975 950 925 900 850 800 750 700 650
                 600 550 500 450 400 350 300 250 200 150
                 100  70  50  30  20  10
tdef 124 linear 00z01dec2007 06hr
vars 10
temp 26 11,100   Temperature [K]
umrl 21 52,100   Relative Humidity [%]
uvel 26 33,100   Zonal Wind [m/s]
vvel 26 34,100   Meridional Wind [m/s]
psnm  0  2,102,0 Mean Sea Level Pressure [Pa]
pslc  0  1,  1,0 Surface Pressure Pressure [Pa]
tsfc  0 11,  1,0 Surface Temperature [K]
sice  0 91,  1,0 Sea Ice Cover [ice=1;no ice=0]
topo  0  7,  1,0 Topography [m]
lsmk  0 81,  1,0 Land Sea Mask [land=1;sea=0]
endvars
```

- `undef -2.56e30` e o valor de ausencia do NCEP (linha 83) — **diferente do ECMWF**.
- `zdef` em **hPa** (linhas 86-88), 26 niveis.
- `umrl` declarada em apenas **21 niveis** (linha 92) porque a FNL nao publica UR
  acima de 100 hPa.
- Descritor `nome <n_niveis> <codigo_grib>,<tipo_nivel>,<nivel>` (formato grib1).
- `psnm` usa `2,102,0` = codigo 2, tipo de nivel 102 (MSL), nivel 0.
  `pslc` usa `1,1,0` = codigo 1, tipo 1 (superficie).

### 3.3 Variaveis e unidades

Codigos NCEP grib1 (tabela 2), `GrADSScripts/SetNMCFnlctl.gs:91-100`:

| nome .ctl | codigo GRIB | tipo/nivel | niveis | origem | destino |
|---|---|---|---|---|---|
| `temp` | 11 | 100 | 26 | K | K |
| `umrl` | 52 | 100 | **21** | % (UR) | **q kg/kg** (Tetens, §3.4) |
| `uvel` | 33 | 100 | 26 | m/s | m/s |
| `vvel` | 34 | 100 | 26 | m/s | m/s |
| `psnm` | 2 | 102,0 | 0 | Pa | **hPa** (/100) |
| `pslc` | 1 | 1,0 | 0 | Pa | hPa (so no binario de ozonio) |
| `tsfc` | 11 | 1,0 | 0 | K | K (usada so como TSM) |
| `sice` | 91 | 1,0 | 0 | 0..1 | 0..1 |
| `topo` | 7 | 1,0 | 0 | m | m (**sem** conversao) |
| `lsmk` | 81 | 1,0 | 0 | 0/1 | 0/1 (mascara de oceano para a TSM) |

### 3.4 Umidade especifica a partir de UR (formula exata)

`GrADSScripts/GetNMCFnlICn.gs:202-221` (funcao `SpecifHumidity`), com o racional
documentado no proprio arquivo:

```
Tetens:  Es(T) = Eo * exp[ a*(T-To)/(T-b) ]
         Eo = 6.1078 hPa ; To = 273.16 K ; a = 17.2693882 ; b = 35.86 K
UR:      umrl = 100*E/Es   =>  E = (umrl/100)*Es        [hPa]
q:       umes = Eps*E / [ lev - (1-Eps)*E ]   Eps = 0.622 ; lev em hPa
```

Linhas literais do GrADS (`GetNMCFnlICn.gs:217-220`):

```
'set z 1 21'
'define E=(umrl/100)*6.1078*exp(17.2693882*(temp-273.16)/(temp-35.86))'
'define umes=0.622*E/(lev-(1-0.622)*E)'
'set z 1'
```

- `lev` e o valor do nivel de pressao do `zdef` (em hPa no grib1).
- A expressao e sempre aplicada nos 21 niveis 1000..100 hPa, na ordem do `zdef`.
- Nao ha limite superior de q nem correcao de saturacao; nao ha tratamento de
  UR > 100 %.

### 3.5 Escrita do binario

`GrADSScripts/GetNMCFnlICn.gs:120-155`: `set x 1 360` / `set y 1 181`
(`:122-123`), `set time <data>`, `set gxout fwrite`,
`set fwrite -be <arq>.dat` (`:126`), e a sequencia:

1. `display topo` (`:129`) — conversao nenhuma;
2. `display psnm/100` (`:132`, `lpnm='psnm/100'` em `:78`);
3. para `n=1..3` (`temp, uvel, vvel`) e `k=1..26`: `set z k` + `display`
   (`:134-145`);
4. para `k=1..21`: `set z k` + `display umes` (`:147-153`).

---

## 4. FNL grib2

### 4.1 Arquivos e .ctl

Mesma logica do grib1, com tres diferencas estruturais:

| item | grib1 | grib2 | referencia |
|---|---|---|---|
| `dtype` | `grib` | `grib2` | `SetFnl_grib2.gs:79` |
| `zdef` | hPa (1000..10) | **Pa** (100000..1000) | `SetFnl_grib2.gs:86-88` |
| descritor var | `cod,tipo,nivel` | `tipo,cod` + `disc,cat,num` | `SetFnl_grib2.gs:91-100` |

`GrADSScripts/SetFnl_grib2.gs:86-100`:

```
zdef  26 levels 100000 97500 95000 ... 2000 1000
vars 10
temp 26,100 0,0,0   Temperature [K]
umrl 21,100 0,1,1   Relative Humidity [%]
uvel 26,100 0,2,2   Zonal Wind [m/s]
vvel 26,100 0,2,3   Meridional Wind [m/s]
psnm  0,101,0 0,3,1  Mean Sea Level Pressure [Pa]
pslc  0,1,0   0,3,0  Surface Pressure Pressure [Pa]
tsfc  0,1,0  0,0,0  Surface Temperature [K]
sice  0,1,0  10,2,0  Sea Ice Cover [ice=1;no ice=0]
topo  0,1,0 0,3,5   Topography [m]
lsmk  0,1,0 2,0,0  Land Sea Mask [land=1;sea=0]
```

O primeiro campo do descritor (`26,100`) e `n_niveis,nivel_type` do GrADS 2.x para
grib2; o segundo e `disciplina,categoria,numero` do GRIB2. `topo` e
`0,3,5` = Geopotential height em **m**, portanto **sem** `/9.80665` (diferente do
ERA, que traz geopotencial em m2/s2).

O driver passa `dtm=6` e `ntimes=1` (`GrADSScripts/roda_FNL_grib2.sh:43`), logo
`tdef` tem **1 tempo** (nao 124): `if (ntimes > 0) lm = ntimes`
(`GrADSScripts/SetFnl_grib2.gs:29-33`). Arquivos esperados:
`.../DataAnl/FNL/*.grib2` (`roda_FNL_grib2.sh:3`); anos `2000..2015`
(`SetFnl_grib2.gs:56`).

### 4.2 A CORRECAO DA PRESSAO EM Pa (ponto critico)

Duas versoes do escritor existem:

- `GrADSScripts/GetNMCFnlICn.gs:219` -> `umes=0.622*E/(lev-(1-0.622)*E)`
  (**grib1**, `lev` em hPa) — usada por `roda_FNL_grib1.sh:47`.
- `GrADSScripts/GetFnlICgrib2.gs:219` -> `umes=0.622*E/((lev/100)-(1-0.622)*E)`
  (**grib2**, `lev` em Pa) — usada por `roda_FNL_grib2.sh:47`.

O `README` do diretorio registra exatamente isso:
`GrADSScripts/README:11` — *"nele (GetFnlICgrib2.gs) a pressao tem que ser dividida
por 100 na equacao da umidade especifica porque os niveis de pressao estao em Pa."*

Todas as demais constantes de Tetens permanecem em hPa (Eo=6.1078, a=17.2693882,
b=35.86), porque `E` e `Eps*E` tambem estao em hPa.

> Existe uma terceira copia `GrADSScripts/GetNMCFnlICn_grib2.gs` (identica a
> `GetFnlICgrib2.gs`, tambem com `lev/100`) que **nao e chamada por nenhum driver**:
> `roda_FNL_grib2.sh` chama `GetFnlICgrib2.gs`. Manter apenas uma no port.

### 4.3 Escrita do binario

Identica ao grib1: mesmo Im/Jm (360x181), mesmo KmaxInp/KmaxInq (26/21), mesmo
`set fwrite -be`, mesma ordem de `display` (`GrADSScripts/GetFnlICgrib2.gs:100-155`).
O unico byte que muda em relacao ao grib1 e o **conteudo de q** (por causa de
`lev/100`).

---

## 5. GFS 0.25 (grib2)

Implementado nos .gs mas **sem driver shell** no espelho (nao ha `roda_GFS*.sh`;
`GrADSScripts/README:3-7` so cita EIT, FNL grib1 e FNL grib2).

- `.ctl`: `GrADSScripts/SetNMCFnlGFSctl_grib2.gs:78-101`, mesma estrutura do FNL grib2,
  `dtype grib2`, `options yrev template`, `undef -2.56e30`, grade
  `xdef 1440 linear 0.000000 0.250000` / `ydef 721 linear -90.000000 0.25`
  (`:85-86`), `zdef 26 levels` em Pa (`:87-89`), `vars 10` com os mesmos
  `(disc,cat,num)` do FNL grib2 (`:91-101`).
- Nome do GRIB: `gfs.t<HH>z.pgrb2.0p25.f000.<AAAA><MM><DD><HH>.grib2` (`:71`), em
  `/scratchout/oper/tempo/externos/Download/FORECAST/GFS_025gr/<AAAA><MM><DD>/`
  (`:75`) e `dset` e **caminho absoluto** (sem `^`, `:79`) — unica fonte que nao
  usa caminho relativo ao `DirMain`.
- Escritor: `GrADSScripts/GetNMCGFSFnlICn.gs` — identico ao FNL grib2 exceto
  `set x 1 1440` / `set y 1 721` (`:122-123`) e `xdef/ydef` 0.25 na saida
  (`:105-106`). Umidade com `lev/100` (`:219`).
- Anos: `2000..2025` (`SetNMCFnlGFSctl_grib2.gs:56`; `GetNMCGFSFnlICn.gs:59`).
- **Natureza do dado**: `f000` de previsao (nao analise), 1 tempo por arquivo.

---

## 6. ERA-40 (duas variantes; so para contexto historico)

- **E4W** (web/ECMWF, 144x73, 23 niveis): `GrADSScripts/SetERA40Wctl.gs:132-136`;
  `zdef 23 levels 1000 925 850 775 700 600 500 400 300 250 200 150 100 70 50 30
  20 10 7 5 3 2 1`; anos `1957..2002` (`:91`); atm/sfc/inv como no ERA-Interim, mas
  a sfc so tem `psnm` (`:163-165`) e o invariante usa `gpsf 0 129,109,1` (`:187`).
  Escritor `GetERAWebICn.gs:122-142` (`case='E4W'`).
- **E40** (interno CPTEC, `/rede/pnud/ERA40/`, 320x161, 1.125): nao gera .ctl —
  ja existe `e4_oper_an_pl.ctl`, `e4_oper_an_sfc.ctl`, `e4_invariant_an_sfc.ctl`
  (`GrADSScripts/GetERA40CICn.gs:87-89`), e usa nomes de variavel nativos
  (`t u v q`, `:91`). Le `direra` de `Configure.ksh` (`:83-85`).
- `GrADSScripts/GetE40PsTvZs.gs` produz um binario **diferente** (`ERA40?PST`), com
  **temperatura virtual** ja calculada (`tqvar='temp.1*(1+0.608*umes.2)'`, `:126`) em
  apenas **13 niveis** 1000..100 (`:157-158`). E consumido por
  `Sources/PsLcERA/GetERAPsLc.f90:29-32` (topo, psnm, 13x Tvir), nao pelo `ICnAGCM`.

---

## 7. LAYOUT EXATO DO BINARIO .dat (item 4)

### 7.1 Contrato lido pelo Fortran

`Sources/ICnAGCM/SpectralICn.f90:94-148`:

```fortran
INQUIRE (IOLENGTH=LRec) WorkInp          ! WorkInp e REAL(KIND=r4) (ImaxInp,JmaxInp)
OPEN (UNIT=nfgrd, FILE=...//DataInp, STATUS='OLD', &
      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=LRec)
nRec=1
READ (UNIT=nfgrd, REC=nRec) WorkInp ; TopoInp(:,JmaxInp:1:-1)=...   ! rec 1
nRec=nRec+1
READ (UNIT=nfgrd, REC=nRec) WorkInp ; PsnmInp(:,JmaxInp:1:-1)=...   ! rec 2
DO k=1,KmaxInp ; nRec=nRec+1 ; READ(...) ; TempInp(:,:,k)   ! KmaxInp registros
DO k=1,KmaxInp ; nRec=nRec+1 ; READ(...) ; UVelInp(:,:,k)   ! KmaxInp registros
DO k=1,KmaxInp ; nRec=nRec+1 ; READ(...) ; VVelInp(:,:,k)   ! KmaxInp registros
DO k=1,KmaxInq ; nRec=nRec+1 ; READ(...) ; SpHuInp(:,:,k,1) ! KmaxInq registros
```

Consequencias (todas obrigatorias no port):

1. **Acesso direto, sem cabecalho e sem registro de controle.** Nao ha bytes de
   contagem Fortran: `RECL = 4*ImaxInp*JmaxInp` bytes (IOLENGTH de um array `r4`),
   e cada registro e um campo 2D completo.
2. **Precisao:** `REAL(KIND=r4) = SELECTED_REAL_KIND(6)` = IEEE float32 de 4 bytes
   (`Sources/ICnAGCM/InputParameters.f90:28`). O `-be` do GrADS
   (`GetERAWebICn.gs:206`, `GetNMCFnlICn.gs:126`, `GetFnlICgrib2.gs:126`,
   `GetNMCGFSFnlICn.gs:126`, `GetIntSST.gs:149`, `GetTrcOzn.gs:207`) e
   **big-endian**.
3. **Ordem dos registros (fixa):**
   `topo, psnm, T(k=1..KmaxInp), u(k=1..KmaxInp), v(k=1..KmaxInp), q(k=1..KmaxInq)`.
4. **Ordem de k dentro de cada variavel = ordem do zdef**, isto e,
   superficie -> topo (k=1 e o primeiro nivel listado = 1000 hPa).
5. **Ordem das linhas: SUL -> NORTE.** Evidencia dupla:
   - o `.ctl` da saida declara `ydef <Jm> linear -90.0 <dxy>`
     (`GetERAWebICn.gs:179`, `GetNMCFnlICn.gs:106`, `GetFnlICgrib2.gs:106`) e o .gs
     escreve `set y 1 <Jm>` nessa ordem (`GetERAWebICn.gs:203`);
   - o Fortran **inverte** ao ler: `TopoInp(:,JmaxInp:1:-1)`
     (`Sources/ICnAGCM/SpectralICn.f90:99,107,116,126,136,146`).
   Ou seja: registro 1 contem a linha do polo sul; o ultimo, o polo norte.
6. **x mais rapido (row-major)**: `WorkInp(ImaxInp,JmaxInp)` e o arquivo e
   x-contiguo -> em Python `np.frombuffer(...).reshape(Jmax, Imax)` sem transpor.
7. **A longitude comeca em 0** (`xdef 240 linear 0.0 1.5`), grade regular
   (nao gaussiana) e **sem coluna duplicada em 360**.
8. **Endianness x plataforma:** o `OPEN` do Fortran **nao** tem `CONVERT=`. O codigo
   so funciona porque o `-be` do GrADS casava com uma maquina big-endian (SX6,
   `Configure.ksh:44,53`). Em x86 (little-endian) o legado produziria lixo — logo o
   port **deve** declarar explicitamente `>f4` e nao depender da plataforma.

### 7.2 Como o .ctl da saida declara o binario

Exemplo ERA-Interim (`GrADSScripts/GetERAWebICn.gs:174-198`):

```
dset ^ERAINTATM2010010100.dat
title ERA Interim at ECMWF
undef 1e+20
options big_endian
xdef 240 linear 0.000 1.5
ydef 121 linear -90.0 1.5
zdef  37 levels 1000 975 ... 3 2 1
tdef   1 linear 00z01jan2010 6hr
vars   6
topo   0 99 Topography [m]
psnm   0 99 Mean Sea Level Pressure [hPa]
temp  37 99 Absolute Temperature [K]
uvel  37 99 Zonal Wind [m/s]
vvel  37 99 Meridional Wind [m/s]
umes  37 99 Specific Humidity [kg/kg]
endvars
```

Pontos de atencao para o parser (`bam-ci/src/bamci/ctl.py`):

- **`options big_endian`** (linha 177): unica opcao declarada, **sem `yrev`** — o
  `.dat` ja esta em ordem sul->norte e o eixo `ydef` cresce de -90 para +90.
  Nao confundir com o `yrev` dos `.ctl` de GRIB de entrada.
- **`undef 1e+20`** (linha 176) no ICn — diferente do `9.99e30` (ERA) e `-2.56e30`
  (FNL) dos GRIB. Ja a TSM usa `-9.99e+08` (`GetIntFnlSST.gs:86`, ecoado em `:147`).
- **`tdef 1`** sempre: um arquivo por data (`GetERAWebICn.gs:190`).
- **Descritor de variavel na forma de contagem implicita**: `temp 37 99 ...` significa
  *37 niveis = os 37 primeiros do zdef*, **sem listar os valores**; e
  `topo 0 99 ...` significa *0 niveis = campo de superficie*; `99` e o campo de
  unidades do GrADS (codigo de variavel), nao um nivel.
  O parser atual (`ctl.py:245-257`) interpreta `tokens[1]` como contagem e exige os
  valores explicitos. **Verificado**: rodando `parse_ctl` sobre um .ctl ICn do EIT,
  o resultado e `CtlError: linha 12: niveis insuficientes em 'temp 37 99 ...'`.
  Pior: `topo 0 99 ...` passa, mas com `levels=()` o `Var` significa "todos os
  niveis do zdef" (`ctl.py:105`), ou seja a topografia seria lida como 37 niveis em
  vez de 1 — erro **silencioso**. Precisa aceitar a forma implicita:
  `niveis = zdef.values()[:count]` e `count == 0` = campo de superficie.
- No FNL/GFS a umidade tem **menos niveis que o zdef**: `umes 21 99 ...`
  (`GetNMCFnlICn.gs:117`) => `z=1..21` = 1000..100 hPa. O layout por variavel,
  portanto, nao e uniforme (26/26/26/21).

### 7.3 Tamanhos (verificacao de sanidade)

| fonte | Im x Jm | bytes/registro | registros | total |
|---|---|---|---|---|
| EIT | 240 x 121 | 116 160 | 2 + 37*3 + 37 = **150** | 17 424 000 |
| EIH | 480 x 241 | 462 720 | 150 | 69 408 000 |
| E4W | 144 x 73 | 42 048 | 2 + 23*3 + 23 = **94** | 3 952 512 |
| FNL (grib1/grib2) | 360 x 181 | 260 640 | 2 + 26*3 + 21 = **101** | 26 324 640 |
| GFS | 1440 x 721 | 4 152 960 | 101 | 419 448 960 |

Formula: `nbytes = 4*Im*Jm*(2 + 3*KmaxInp + KmaxInq)`.

### 7.4 Pos-processamento dentro do Fortran (nao replicar no escritor)

- Extrapolacao vertical de q quando `KmaxInq < KmaxInp` (FNL/GFS): `MoistVertExtrap`
  (`Sources/ICnAGCM/SpectralICn.f90:149-152`). O `.dat` guarda **apenas** os 21
  niveis; os 5 niveis do topo sao preenchidos no Fortran.
- `Tv = T*(1 + 0.608*q)` (`:160`).
- `pres(k) = a(k) + b(k)*ps` com `a_hybr` em Pa e `ps` em hPa
  (`:198-199`; nota `:313` "a_hybr(k)*100 convert back from mb to Pa").

---

## 8. Binarios auxiliares (SST/TSM e ozonio)

### 8.1 TSM/SST — `<PREFIXO>SST<AAAAMMDDHH>.dat`

- Escritores: `GrADSScripts/gera_bin_tsm.gs` (usado por `roda_EIT.sh:56`),
  `GrADSScripts/GetIntFnlSST.gs` (FNL, 3 casos EIT/EIH/FNL),
  `GrADSScripts/GetIntSST.gs` (variante EIT/FNL).
- Conteudo: **2 registros** — rec 1 `sstk` (K), rec 2 `sice` (0..1);
  ctl `vars 2`, `zdef 1 levels 1000`, `tdef 1 linear <time> 6hr`, `options big_endian`
  (`gera_bin_tsm.gs:136-141`).
- Leitura Fortran: `Sources/SST/SSTFill.f90:29-30`, mesma inversao
  `JmaxInp:1:-1` (linhas 33, 39, 40).
- **Para o FNL a TSM e derivada da temperatura de pele com mascaramento de
  continente**: `lvar='const(maskout(tsfc,1-2*lsmk),<Und>,-u) const(maskout(sice,1-2*lsmk),<Und>,-u)'`
  (`GetIntFnlSST.gs:125`) — `maskout` anula o continente e `const(...,-u)` repoe o
  valor de ausencia. Para ERA usa `sstk sice` diretos (`:97,111`).
- Valores de ausencia: `GetIntFnlSST.gs:86` `Und=-9.99e+08` (todos os casos);
  `GetIntSST.gs:86` `Und=9.99e30`; `gera_bin_tsm.gs:86` `Und=-9.99e+08`.
  O Fortran espera `Undef=-9.99e8` para EIT e `9.99e30` para EIH/FNL
  (`Sources/SST/InputParameters.f90:93-107`) — **inconsistente entre si**; na
  pratica `SSTFill.f90:41,45,51,55` (`WHERE (SSTIn > 400) SSTIn=Undef`) resgata o
  caso. No port, gravar `undef` de forma unica e coerente.
- Dimensoes por caso: `GetIntFnlSST.gs:91-127` (EIT 240x121/1.5, EIH 480x241/0.75,
  FNL 360x181/1.0); `Sources/SST/InputParameters.f90:92-108`.
- **O GFS nao tem fonte de TSM** — herda a do ERA (ver §10 e `README:9`).

### 8.2 Ozonio — `<PREFIXO>OZN<AAAAMMDDHH>.dat`

- Escritores: `GrADSScripts/GetTrcOzn.gs` (casos E40/EIT/E4W) e
  `GrADSScripts/gera_bin_oz.gs` (E40/EIT/EIH/E4W, versao mais nova).
- Conteudo: **1 + Km registros** — rec 1 `pslc` em **hPa**, rec 2.. `ozon` (razao de
  mistura em massa, kg/kg) nivel a nivel (`GetTrcOzn.gs:197-198,209-218`).
- ctl: `vars 2` (`pslc 0 99`, `ozon <Km> 99`), `undef 1e+20`, `options big_endian`
  (`GetTrcOzn.gs:179-199`).
- Leitura Fortran: `Sources/Ozone/Ozone.f90:29-36` (pslc primeiro, depois Kmax
  niveis), `Undef=9.99E30` (`Sources/Ozone/InputParameters.f90:153`), casos apenas
  `E40/E4W/EIT` (`:113-135`).
- **Nao existe ozonio para FNL nem para GFS** — `README:9` diz que nos hindcasts
  FNL a TSM e o ozonio foram mantidos do ERA-Interim.
- **Bug legado em `gera_bin_oz.gs:207`**: `if (case = 'EIT' & case = 'EIH')` e
  sempre falso, entao o `.ctl` sai com `zdef 23` mesmo quando 37 registros foram
  escritos (EIT). `GetTrcOzn.gs:185` tem a condicao correta (`if (case = 'EIT')`).
  **Nao reproduzir o bug**.
- `roda_EIT.sh` **nao** regenera ozonio: as chamadas a `gera_bin_oz.gs` estao
  comentadas (`roda_EIT.sh:53-54,59`).

### 8.3 PsLc (so ERA-40, para o executavel PsLcERA)

- `GrADSScripts/GetE40PsTvZs.gs`: 3 variaveis, `topo`, `psnm` (hPa) e `tvir`
  (temperatura virtual) em 13 niveis 1000..100 (`:157-158,171-173`).
- Leitura: `Sources/PsLcERA/GetERAPsLc.f90:29-32`; `Kmax=13`, `PresRef`
  `1000..100` (`Sources/PsLcERA/InputParameters.f90:13,60`).

---

## 9. Ordem de execucao dos drivers

### 9.1 roda_EIT.sh (ERA-Interim)

1. Le `lista_pci` (uma data `AAAAMMDD` por linha, ex. `19990106`,
   `GrADSScripts/lista_pci:1`) e converte o mes para 3 letras
   (`roda_EIT.sh:4-38`).
2. Escreve um `corre.gs` de uma linha chamando `gera_ctl.gs` e roda
   `grads -lbc corre.gs` (`:43-45`) -> gera os 3 `.ctl` + `.gmp` do ERA-Interim
   (atmosfera/superficie/invariante). `gera_ctl.gs` define a **funcao
   `SetERAIntctl`**, ou seja e uma copia renomeada de `SetERAIntctl.gs`.
3. Reescreve `corre.gs` chamando `gera_bin_tsm.gs <data> EIT` e roda
   `grads -lbc` (`:56-57`) -> gera `ERAINTSST<data>.dat/.ctl`.
4. Etapas **comentadas** (nao executam): `junta_atm.gs` (`:47-51`, geraria
   `ERAINTATM<data>.dat`, mas nao e chamada em producao) e `gera_bin_oz.gs`
   (`:53-54,59`).
   > **O binario atmosferico do ERA-Interim nao e gerado por roda_EIT.sh.** Ele
   > deve ser produzido por chamada manual a `GetERAWebICn.gs <data> EIT`
   > (documentado em `GetERAWebICn.gs:256`) ou por `junta_atm.gs` (identico,
   > `junta_atm.gs:146-227`).
5. Diretorios: dados e saida ambos em `DirMain/DataAnl/EIT`
   (`GetERAWebICn.gs:107-108`); `DirMain` vem de `Configure.ksh:24` e e lido via
   arquivo `dirmain` (`Configure.ksh:60`; `SetERAIntctl.gs:72-76`).
6. Listas de datas: `lista_pci` (1043 linhas, 1999-2018) e `lista`/`lista2`
   (esquemas de datas por mes para outros experimentos,
   `GrADSScripts/lista:1-194`, `lista2:1-124`).

### 9.2 roda_FNL_grib1.sh

1. `lista=$(ls -1 .../DataAnl/FNL/*.grib1)` (`:3`).
2. Extrai ano/mes/dia do **caminho absoluto por posicao fixa de caractere**
   (`ano=${file:73:4}`, `mes=${file:77:2}`, `dia=${file:79:2}`, `:8,77,79`) e
   converte o mes para 3 letras (`:10-36`). Fragil: depende do comprimento exato
   do prefixo do diretorio.
3. `SetNMCFnlctl.gs <data>` via `grads -lbc` (`:43-45`) -> `.ctl` + `.gmp`.
4. `GetNMCFnlICn.gs <data>` via `grads -lbc` (`:47-49`) ->
   `NMCFNLATM<data>.dat/.ctl`.
5. Sem TSM nem ozonio: `README:9` (mantidos do ERA-Interim).

### 9.3 roda_FNL_grib2.sh

Identico, com `*.grib2` (`:3`), `SetFnl_grib2.gs <data> 6 1` (`:43`) e
`GetFnlICgrib2.gs <data>` (`:47`). O terceiro argumento `1` fixa `lm=1`
(um tempo no `.ctl`).

### 9.4 Chamada ao executavel

Os scripts que chamam o binario Fortran (`RunICnAGCM.ksh`, `RunSSTAGCM.ksh`,
`RunOznAGCM.ksh`, `RunERAPsLc.ksh`, `RunCheckICn.ksh`) **nao estao no espelho**:
`RunScripts/bash`, `RunScripts/ksh` e `RunScripts/scripts2` estao vazios (sao
gerados por `configure:134-177` a partir de `backup/RunScripts.tar.gz`, tambem
ausente). O que se sabe com certeza:

- o executavel recebe um namelist `ICnAGCMNML`
  (`Sources/ICnAGCM/InputParameters.f90:144`), com `Case`, `DateICn` (AAAAMMDDHH),
  `DirMainInp`/`DirMainOut` (`:141-142`);
- `DirICnInp = DirMainInp//'DataAnl/'//Case//'/'` (`:212`) e
  `DirICnOut = DirMainOut//'model/datain/'` (`:213`);
- o executavel roda a partir do diretorio de producao, coerente com
  `DirMainInp='./ '` (`:141`).

Essa e a **lacuna de rastreabilidade** mais importante do item 5: para reproduzir
o pipeline completo e preciso obter `RunScripts/*/RunICnAGCM.ksh` da jaci.

---

## 10. Diferencas entre fontes que afetam o resultado

| # | aspecto | ERA-Interim | ERA-40 | FNL grib1 | FNL grib2 | GFS |
|---|---|---|---|---|---|---|
| 1 | umidade de origem | q pronta (133) | q pronta | **UR (52) -> Tetens** | UR (0,1,1) -> Tetens **com lev/100** | idem grib2 |
| 2 | niveis de q | 37/37 | 23/23 | 26/**21** | 26/**21** | 26/**21** |
| 3 | zdef, unidade | hPa | hPa | hPa | **Pa** | **Pa** |
| 4 | topografia | geopotencial m2/s2 **/9.80665** | geopotencial **/9.80665** | `topo` em m | `0,3,5` em m | `0,3,5` em m |
| 5 | pressao MSL | 151 Pa **/100** | 151 Pa **/100** | 2 Pa **/100** | 0,3,1 Pa **/100** | idem |
| 6 | undef do GRIB | 9.99e30 | 9.99e30 | -2.56e30 | -2.56e30 | -2.56e30 |
| 7 | n. de tempos no tdef | ~124 (mes) | ~124 | ~124 | **1** | 1 |
| 8 | TSM | `sstk` (34) | `sstk` | **`tsfc` mascarada por `lsmk`** | idem | **inexistente** |
| 9 | ozonio | 203 | 203 (E4W) / `o3` (E40) | **inexistente** | **inexistente** | **inexistente** |
| 10 | grade | 240x121 / 480x241 | 144x73 / 320x161 | 360x181 | 360x181 | 1440x721 |
| 11 | interpolacao horizontal | Fortran | Fortran | Fortran | Fortran | Fortran |

Impactos qualitativos a registrar no port:

- **(1)/(2)** FNL e GFS tem **q so ate 100 hPa**; os 5 niveis superiores sao
  extrapolados pelo Fortran (`MoistVertExtrap`), nao pelo pre-processamento.
  Qualquer reimplementacao que gere q em 26 niveis muda o resultado.
- **(3)/(2)** o caso grib2 usa `lev/100`; se o port reusar a formula do grib1 com
  `lev` em Pa, q sai ~100x errado (o denominador `lev - 0.378*E` domina).
- **(4)** topografia: no ERA a conversao e obrigatoria; no NCEP nao. Aplicar
  `/9.80665` no NCEP dividiria a topografia por 9.8.
- **(7)** com `tdef 1` (FNL grib2/GFS) o `set time` do GrADS cai sempre em t=1;
  com `tdef 124` (ERA/FNL grib1) e o indice do tempo que seleciona a data.
- **(8)** TSM do FNL e uma **aproximacao** (temperatura de pele mascarada), nao uma
  analise de TSM; e o `maskout(...,1-2*lsmk)` inverte a mascara (oceano=1).
- **(11)** a grade de entrada nao e a grade do modelo em nenhum caso: a
  interpolacao horizontal e feita pelo Fortran (`SpectralICn.f90:100-104`), nao
  pelo GrADS. O `.dat` mantem a grade nativa da fonte.

---

## 11. O que e COMUM e o que e ESPECIFICO (para ERA5 / FNL / GFS)

### 11.1 Comum a todas as fontes (nao parametrizar)

1. **Contrato binario** do §7.1: float32 big-endian, acesso direto, registros
   `topo, psnm, T^KmaxInp, u^KmaxInp, v^KmaxInp, q^KmaxInq`, sem cabecalho,
   x-contiguo, linha 1 = polo sul, k=1 = 1000 hPa.
2. **Semantica das 6 variaveis**: topografia em m; psnm em **hPa**; T seca em K;
   u, v em m/s; q em kg/kg.
3. **Cadeia** `GRIB -> .ctl + gribmap -> .dat` e a nomenclatura
   `<PREFIXO>ATM<AAAAMMDDHH>.dat` + `.ctl` companion.
4. **Convencao de data** `AAAAMMDDHH` e hora em {00,06,12,18} UTC; mes em 3 letras
   minusculas na sintaxe GrADS `hhzddmclyyyy`.
5. **`options yrev`** nos `.ctl` de GRIB (fonte norte->sul) e **ausencia de
   `yrev`** no `.ctl` do `.dat` (que ja esta sul->norte).
6. **Temperatura virtual e interpolacao vertical/horizontal ficam no Fortran** —
   nunca no pre-processamento.
7. **`tdef 1`** no `.ctl` do `.dat` (um arquivo por data).

### 11.2 Especifico de cada fonte (parametrizar)

| item | onde parametrizar |
|---|---|
| grade (Im, Jm, dxy) e ordem dos niveis | por fonte (§1) |
| codigos GRIB e tipo de nivel (hPa vs Pa) | tabela por fonte (§2.3, §3.3, §4.1) |
| obtencao de q (direta vs Tetens) e divisor de `lev` | flag por fonte (§3.4, §4.2) |
| conversao de topografia (m2/s2 vs m) | flag por fonte |
| `undef` de entrada | por fonte |
| existencia de TSM e de ozonio, e a regra de fallback | por fonte (§8, §9.4) |
| numero de tempos no `.ctl` de entrada | por fonte |
| nome/padrao do GRIB (mensal vs por data; template vs caminho absoluto) | por fonte |

### 11.3 ERA5 — nao existe no legado; o que falta

Nao ha **nenhuma** referencia a ERA5 no espelho (busca por `era5|ERA5|ERA-5` em toda
a arvore nao retorna nada; os limites de ano param em 2017/2018:
`SetERAIntctl.gs:55`, `gera_bin_tsm.gs:24`, `gera_bin_oz.gs:31`,
`junta_atm.gs:24`). Para suportar ERA5 e preciso decidir/verificar:

1. **Formato de aquisicao**: GRIB1 (o CDS ainda entrega, com os mesmos codigos
   ECMWF 130/131/132/133/203/151/134/34/31/129) ou GRIB2/NetCDF. Se GRIB2, o
   `dtype` e os descritores mudam como em `SetFnl_grib2.gs`.
2. **Niveis**: ERA5 tem os mesmos 37 niveis 1000..1 hPa — a `zdef` do EIT deve
   servir sem mudanca; falta confirmar que a lista de 37 e identica.
3. **Ozonio (codigo 203)**: confirmar se o ERA5 de producao traz
   `ozone_mass_mixing_ratio` em todos os niveis/anos; caso nao traga, o ozonio
   tera de vir de outra fonte (o legado ja faz esse fallback para FNL, `README:9`).
4. **TSM/gelo**: ERA5 fornece `sst` (34) e `siconc`/`ci` (31) — mapear o nome
   correto (`ci` no ERA5 e concentracao de gelo, equivalente ao 31 do ERA-Interim).
5. **Grade**: 0.25 grau nativo (1440x721) exigiria nova entrada na tabela do §1 e
   no `SELECT CASE` do Fortran (`InputParameters.f90:161-204`) — ou seja **exige
   mudanca no Fortran**, nao so no pre-processamento. Alternativa: regridar para
   0.75/1.5 antes de escrever o `.dat` e reusar `EIH`/`EIT`.
6. **Janela temporal**: ERA5 vai ate o presente; remover o guarda `yy > 2017`.

### 11.4 GFS — o que falta

1. **Driver shell**: `SetNMCFnlGFSctl_grib2.gs` e `GetNMCGFSFnlICn.gs` existem, mas
   nao ha `roda_GFS*.sh`. E preciso escrever o orquestrador (no port isso e trivial,
   mas o legado nao tem referencia de ordem de passos para GFS).
2. **`PressureLevels.NMCGFS`**: o Fortran pede esse arquivo em
   `Sources/ICnAGCM/InputParameters.f90:202`; ele vive em `DataLevels/`, que **nao
   esta no espelho** (`Configure.ksh:71-75` cria um symlink para
   `$DirHome/DataLevels`; `configure:108-110` o extrai de
   `backup/DataLevels.tar.gz`, ausente). Sem ele nao se sabe a ordem/valores
   exatos de `PresInp` para GFS.
3. **TSM e ozonio**: GFS `f000` nao traz SST nem ozonio. O legado resolve herdando
   do ERA-Interim (`README:9`); o port precisa explicitar essa regra de fallback
   por fonte, nao deixa-la implicita.
4. **`f000` e previsao, nao analise**: o nome `gfs.tHHz.pgrb2.0p25.f000`
   (`SetNMCFnlGFSctl_grib2.gs:71`) e o campo de passo 0 do ciclo; para um
   hindcast subsazonal isso e aceitavel como analise, mas deve ser documentado
   como diferenca em relacao a FNL (que e analise).
5. **Grade 0.25 (1440x721)** ja esta prevista no Fortran
   (`InputParameters.f90:198-203`), logo o GFS e a fonte com **menos** mudancas
   estruturais necessarias — falta apenas o orquestrador e a TSM/ozonio.
6. **Caminho absoluto** em `dset` (`SetNMCFnlGFSctl_grib2.gs:75,79`): o port deve
   receber esse diretorio como configuracao, nao como constante.

---

## 12. Defeitos e fragilidades do legado (nao reproduzir)

1. `dtm` sempre forcado a 6 h por `if (dtm != 06 | dtm != 12 | dtm != 24)`
   (`SetERAIntctl.gs:11-13` e copias) — o passo de 12 h/24 h nunca e usado.
2. `gera_bin_oz.gs:207`: `if (case = 'EIT' & case = 'EIH')` sempre falso ->
   `.ctl` de ozonio com `zdef 23` para dados de 37 niveis.
3. `gera_bin_tsm.gs:86` usa `Und=-9.99e+08` tambem para EIT, enquanto
   `GetIntSST.gs:86` usa `9.99e30` e o Fortran espera valores diferentes por caso
   (`Sources/SST/InputParameters.f90:97,102,107`).
4. Template do FNL sem extensao (`SetNMCFnlctl.gs:71`) vs glob `*.grib1`
   (`roda_FNL_grib1.sh:3`).
5. Data extraida por offset fixo de caractere no caminho absoluto
   (`roda_FNL_grib1.sh:8`, `roda_FNL_grib2.sh:8`).
6. `corre.gs`/`corre_fnl.gs` escritos com aspas desbalanceadas
   (`roda_EIT.sh:43`, `roda_FNL_grib1.sh:43`): `echo "gera_ctl.gs' '00z... ; 'quit' "`.
   O separador script/argumento e a sequencia literal `' '`.
7. Cursor de arquivo implicito: os `Get*` leem a lista de `.ctl` com
   `read(<arquivo>)` repetido, dependendo do ponteiro interno do GrADS avancar
   (`GetERAWebICn.gs:109-120`). No port, usar indice explicito.
8. Copias divergentes do mesmo escritor: `GetNMCFnlICn.gs` /
   `GetNMCFnlICn_grib2.gs` / `GetFnlICgrib2.gs`; `GetIntSST.gs` /
   `GetIntFnlSST.gs` / `gera_bin_tsm.gs`; `GetTrcOzn.gs` / `gera_bin_oz.gs`;
   `junta_atm.gs` / `GetERAWebICn.gs`; `gera_ctl.gs` / `SetERAIntctl.gs`.
9. `GrADSScripts/rmctl` apaga `./CTL/*` com `rm -f` (`rmctl:3`) — sem impacto
   cientifico, mas e estado global compartilhado entre execucoes.

---

## 13. Checklist de aceitacao para o port (entradas)

Para cada fonte suportada, o novo codigo deve:

- [ ] listar as variaveis GRIB com **codigo, tipo de nivel, niveis e unidades**,
      conforme as tabelas dos §2.3/§3.3/§4.1;
- [ ] derivar q por Tetens com `(Eo,To,a,b)=(6.1078, 273.16, 17.2693882, 35.86)` e
      `Eps=0.622`, usando `lev` em hPa — **dividindo o zdef por 100 quando a
      fonte declara pressao em Pa** (FNL grib2, GFS, e ERA5 se vier em Pa);
- [ ] converter psnm Pa->hPa e topografia m2/s2->m **apenas** quando a fonte
      entrega geopotencial;
- [ ] escrever o `.dat` com `nbytes = 4*Im*Jm*(2 + 3*KmaxInp + KmaxInq)`,
      `>f4`, ordem `topo, psnm, T, u, v, q`, linha 1 = polo sul, k=1 = 1000 hPa;
- [ ] emitir o `.ctl` companion com `undef`, `options big_endian`, `xdef/ydef`
      da grade nativa, `zdef` na ordem superficie->topo, `tdef 1`, e descritores na
      **forma implicita** (`var 0 ...` para superficie, `var <n> ...` para volume);
- [ ] validar hora em {00,06,12,18}, mes valido, data dentro do periodo da fonte e
      existencia do GRIB esperado, com erro claro em vez de silencio;
- [ ] declarar explicitamente, por fonte, se TSM e ozonio vem da propria fonte ou
      sao herdados (regra do `README:9`).

---

## 14. Mapa de arquivos citados

| arquivo | papel |
|---|---|
| `GrADSScripts/README` | regras de uso e a nota sobre Pa no grib2 (linha 11) |
| `GrADSScripts/roda_EIT.sh` | orquestrador ERA-Interim |
| `GrADSScripts/roda_FNL_grib1.sh` | orquestrador FNL grib1 |
| `GrADSScripts/roda_FNL_grib2.sh` | orquestrador FNL grib2 |
| `GrADSScripts/gera_ctl.gs` | copia de `SetERAIntctl.gs` usada pelo driver EIT |
| `GrADSScripts/SetERAIntctl.gs` | .ctl do ERA-Interim 1.5 |
| `GrADSScripts/SetERAInHctl.gs` | .ctl do ERA-Interim 0.75 |
| `GrADSScripts/SetERA40Wctl.gs` | .ctl do ERA-40 web |
| `GrADSScripts/SetNMCFnlctl.gs` | .ctl do FNL grib1 |
| `GrADSScripts/SetFnl_grib2.gs` | .ctl do FNL grib2 |
| `GrADSScripts/SetNMCFnlGFSctl_grib2.gs` | .ctl do GFS 0.25 |
| `GrADSScripts/GetERAWebICn.gs` | .dat atmosferico ERA-Interim 1.5 / ERA-40 web |
| `GrADSScripts/GetERAHWebICn.gs` | .dat atmosferico ERA-Interim 0.75 |
| `GrADSScripts/GetERA40CICn.gs` | .dat atmosferico ERA-40 interno |
| `GrADSScripts/junta_atm.gs` | duplicata de `GetERAWebICn.gs` |
| `GrADSScripts/GetNMCFnlICn.gs` | .dat FNL grib1 (Tetens em hPa) |
| `GrADSScripts/GetFnlICgrib2.gs` | .dat FNL grib2 (Tetens com `lev/100`) |
| `GrADSScripts/GetNMCFnlICn_grib2.gs` | duplicata nao usada de `GetFnlICgrib2.gs` |
| `GrADSScripts/GetNMCGFSFnlICn.gs` | .dat GFS 0.25 |
| `GrADSScripts/GetIntSST.gs`, `GetIntFnlSST.gs`, `gera_bin_tsm.gs` | TSM/gelo |
| `GrADSScripts/GetTrcOzn.gs`, `gera_bin_oz.gs` | ozonio |
| `GrADSScripts/GetE40PsTvZs.gs` | PsLc/Tvir do ERA-40 (13 niveis) |
| `GrADSScripts/GRIBMap.ksh` | indice `.gmp` por `gribmap` |
| `GrADSScripts/lista`, `lista2`, `lista_pci` | listas de datas de experimentos |
| `GrADSScripts/corre.gs`, `corre_fnl.gs` | buffers de uma linha (sem ciencia) |
| `Configure.ksh`, `configure` | `DirMain`, `DirGrADS`, criacao de `DataAnl/*` |
| `RunScripts/lista` | datas `AAAAMMDDHH` (1043 linhas, 1999-2018) |
| `Sources/ICnAGCM/SpectralICn.f90` | **leitor do .dat** (linhas 94-148) |
| `Sources/ICnAGCM/InputParameters.f90` | `ImaxInp/JmaxInp/KmaxInp/KmaxInq` por caso |
| `Sources/SST/*`, `Sources/Ozone/*`, `Sources/PsLcERA/*` | leitores dos binarios auxiliares |
