# Núcleo numérico Fortran do gerador de IC do BAM 1.2 subsazonal

Especificação derivada **somente da leitura do código** em
`reference/GetICnAGCM_subsaeasonal/` (cópia read-only de
`/oper/scripts/subsaz/Bruno/GetICnAGCM_subsaeasonal` na jaci).
Citações no formato `arquivo:linha`, caminhos relativos a essa raiz; índices
Fortran 1-based. Nada foi obtido executando o legado (não há binário de
referência no espelho): tudo vem de leitura de código + **verificação numérica
independente em Python/numpy** nos pontos marcados ✅.
Documentos irmãos: `spec_inputs.md` (GRIB→binário GrADS), `spec_output.md`
(envelope do arquivo de IC). Aqui: o **núcleo numérico**.

---

## 0. Fluxo de `SpectralICn` (ordem significativa)

```
InitParameters          ICnAGCM/InputParameters.f90:100-377
GetArrays               ICnAGCM/InputArrays.f90:49-94
Init{Linear|Area}Interp ICnAGCM/SpectralICn.f90:48-52
InitFFT / InitLegendre  ICnAGCM/Fourier.f90:24-38 / Legendre.f90:34-115
topografia (espectral opcional + em grade)      SpectralICn.f90:58-104
leitura+interp. horiz. Topo,Psnm,T(k),u(k),v(k),q(k<=KmaxInq)  :94-159
MoistVertExtrap (se KmaxInq<KmaxInp)            :149-152
T -> T_vir  (TempOut *= 1+cTv*q)                :160
NewSigma                                        :173-174
pressoes dos niveis de saida (hibrido/sigma)    :182-225
NewPs -> PsOut, LnPsOut                         :233-236
recalcula pressoes com PsOut                    :238-270
GrADS entrada                                   :272-279
T_vir -> T  (TempOut /= 1+cTv*q)                :281
VertSigmaInter (u,v,T,q)                        :282-284
T -> T_vir  (TempSig *= 1+cTv*q)                :285
GrADS saida                                     :287-294
DectoSpherHarm+TransSpherHarm (Topo,LnPs,T,v,q) :296-304
UVtoDivgVort+TransSpherHarm (Div,Vort)          :306-308
escrita da IC (_hibrid ou _sigma)               :310-354
```

`TempOut` muda de significado 3× (T real → virtual → real → virtual) e `PresSig`
é recalculado depois de `NewPs`: a ordem não é cosmética.

---

## 1. Namelist e parâmetros

### 1.1 Namelist (`InputParameters.f90:111-112`, lido de stdin, `:144`)

| nome | default | significado | citação |
|---|---|---|---|
| `MendOut` | 62 | truncamento triangular de saída | `:131` |
| `KmaxOut` | 28 | nº de camadas do modelo | `:132` |
| `DateICn` | `'yyyymmddhh'` | data/hora da análise | `:135,146-149` |
| `Iter` | 10 | iterações do *smoothing* de topografia (**morto**, §9) | `:128` |
| `SmthPerCut` | 0.12 | corte do *smoothing* (**morto**, §9) | `:129` |
| `Case` | `'E40'` | base de dados (§1.4) | `:136-139` |
| `Grid` | `'Q'` | entra **só** no nome dos arquivos | `:140,288` |
| `SigmaLevel` | `.false.` | `.false.`=híbrido sigma-p; `.true.`=sigma puro | `:95,183` |
| `GrADS`/`GrADSOnly` | T/F | diagnósticos / parar após eles | `:133-134,272-278` |
| `DirMainInp`/`Out` | `'./ '` | raízes de entrada/saída | `:141-142` |

### 1.2 Dimensões derivadas (`:306-324`)

`Mend1Out=Mend+1`, `Mend2Out=Mend+2`, `Mend3Out=Mend+3`;
`Mnwv2Out=Mend1*Mend2` = **nº de coeficientes reais**; `Mnwv0Out=Mnwv2/2` = nº de
complexos; `Mnwv3Out=Mnwv2+2*Mend1`, `Mnwv1Out=Mnwv3/2` (dimensões do vetor
`pln`); `ImxOut=ImaxOut+2` (folga do FFT), `JmaxhfOut=JmaxOut/2`,
`KmaxOutp=KmaxOut+1`; `NTracers=1`, `Kdim=1`, `ICaseDec=1`, `ICaseRec=-1`
(constantes); `MFactorFourier=64`, `MTrigsFourier=3*ImaxOut/2`;
`ForecastDay=0`, `TimeOfDay=0`.
**`LinearInterp = (ImaxInp <= ImaxOut)`** (`:226-230`) é o único critério
bilinear × área.
Para `MendOut=62`, `KmaxOut=28`: `ImaxOut,JmaxOut = 192,96`; `Imx=194`;
`Jmaxhf=48`; `Mend1..3 = 63,64,65`; `Mnwv2=4032`, `Mnwv0=2016`, `Mnwv1=2079`,
`Mnwv3=4158`; `KmaxOutp=29`; `dx=1.875°`.

### 1.3 Constantes (`:329-350`, r8)

`EMRad=6.37E6` m; `Grav=9.80665`; `Rd=287.05`; `Rv=461.50`; `Cp=1004.6`;
`Lc=2.5E6`; `Gamma=-6.5E-3` K/m; `GEps=1e-9`; `cTv=Rv/Rd-1=0.607752`;
`RoCp=Rd/Cp=0.285736`; `RoCp1=RoCp+1`; `RoCpR=1/RoCp=3.50`;
`Cte=1/√2=0.7071067812`; `Undef=-999999.0`; `Lon0=0`, `Lat0=90`;
`PolarMean=.false.`; `MaskInput=1`; `FlagInput=(T,T,T,T,F)`;
`FlagOutput=(T,T,F,T,T)`. Semântica das flags: `AreaInterpolation.f90:19-32`
(1) começa no polo norte; (2) em Greenwich; (3) latitudes no centro da caixa;
(4) longitudes no centro; (5) grade gaussiana.
⇒ entrada = regular lat/lon centrada em 90N/0E; saída = gaussiana com latitudes
na **borda norte** da caixa (`:372`).

### 1.4 `Case` (`:161-204`) e arquivos de níveis (`:258-273`)

| Case | ImaxInp | JmaxInp | KmaxInp | KmaxInq | DataPrs | interp. |
|---|---|---|---|---|---|---|
| E40 | 320 | 161 | 23 | 23 | PressureLevels.ERA40C | área |
| E4W | 144 | 73 | 23 | 23 | PressureLevels.ERA40W | **bilinear** |
| EIT | 240 | 121 | 37 | 37 | PressureLevels.ERAINT | área |
| EIH | 480 | 241 | 37 | 37 | PressureLevels.ERAINH | área |
| FNL | 360 | 181 | 26 | 21 | PressureLevels.NMCFNL | área |
| GFS | 1440 | 721 | 26 | 21 | PressureLevels.NMCGFS | área |

Níveis: `DataPrs` com `(8F7.1)` → `KmaxInp` pressões (hPa) superfície→topo;
`DataSig` com `(5F9.6)` → `KmaxOut` espessuras sigma; `DataHyb` (lista livre
`a b`) → `KmaxOut+1` pares de **interfaces** em Pa, e `a_hybr` é convertido a
hPa no próprio laço (`a/100`, `:266-273`); `b_hybr` não é convertido.
Na escrita `a` volta a Pa (`HybIOut=a*100`, `SpectralICn.f90:320`).
A soma de `DelSigmaOut` é apenas impressa com o resíduo (`:284-285`); o resíduo
vai para o último nível (§7.1).

### 1.5 Caminhos e nomes (quirks obrigatórias)

`InputParameters.f90:212-219, 287-300`:

```
DirICnInp = TRIM(DirMainInp)//'DataAnl/'//Case//'/ '   ! branco no fim
DirICnOut = TRIM(DirMainOut)//'model/datain/ '
DirVert   = TRIM(DirMainInp)//'/DataLevels/ '          ! barra dupla
DirTop    = TRIM(DirMainInp)//'/DataTop/'
DataSig   = 'DeltaSigma.L'  //I3.3(KmaxOut)     -> DeltaSigma.L028
DataHyb   = 'HybridLevels.L'//I3.3(KmaxOut)     -> HybridLevels.L028
DataTop   = 'Topography.TQ' //I4.4(MendOut)      -> Topography.TQ0062
DataTopGrid='Topography.G'  //I5.5(JmaxOut)      -> Topography.G00096
TruncOut  = 'T'//Grid//I4.4(Mend)//'L'//I3.3(Kmax)  -> TQ0062L028 (10 chars)
DataOut   = 'GANL'//Case//DateICn//'S.unf.'//TruncOut
```

`SpectralICn` anexa `'_hibrid'` ou `'_sigma'` (`:311,315`). `TRIM` só é
aplicado a `DirMainInp/Out`; no port use `pathlib` mas **documente** a diferença.

---

## 2. Precisão e unidades

Cálculo em **r8**; r4 só em I/O/buffers (`InputParameters.f90:27-29`,
`InputArrays.f90:17-23`). Arquivos GrADS são r4: leitura
`REAL(WorkInp,r8)` (`SpectralICn.f90:99..146`), escrita `REAL(...,r4)`
(`:337-353`). **Não há** `-r8`/`-i8` nas flags
(`Sources/makefiles/Makefile.linux_gnu:12`, `Makefile.linux_intel:12`).
Unidades: pressão em **hPa** (`a_hybr` inclusive); `LnPsOut = ln(ps_hPa/10)` =
ln(cbar) (`Utils.f90:156`; confere com `10*EXP(Cte*LnPsSC(1))` em
`SpectralICn.f90:372` e `CheckICn.f90:36`); T em K; q em kg/kg; vento m/s;
topografia m.

---

## 3. Geometria

### 3.1 `GetImaxJmax` (`InputParameters.f90:380-446`)

`Lfft=40000` (`:390`); gera os `2^a3^b5^c` com `a>=1` (laço de `n2` começa em 1,
`:403-416`) e `<=Lfft`, ordena por *bubble sort* (`:419-432`);
`Imax =` menor candidato `>= 3*Mend+1` (`:434-440`); `Jmax=Imax/2`, arredondado
para par (`:441-442`). Para `Mend=62`: 187 → `192=2^6·3` ⇒ 192, 96 ✅.

### 3.2 Grades regulares (`LinearInterpolation.f90:68-111`)

`Lon(i)=Lon0+(i-1)·360/Imax`; `Lat(j)=Lat0-(j-1)·dy`, `dy=2·Lat0/(Jmax-1)`.
⇒ `Lat(1)=+90`, `Lat(Jmax)=-90`: **orientação norte→sul** (coerente com
`Lat0=90` e `FlagInput(1)=.true.`).

### 3.3 Gaussianas (`Legendre.f90:216-263`)

```
prec=100·EPSILON(r8); scale=2/Jmax²; dradz=atan(1)/Jmax; rad=0
para j=1..Jmax/2:  drad=dradz
   repita: p2=P_Jmax(rad); repita {p1=p2; rad+=drad; p2=P_Jmax(rad)} até sign(p1)/=sign(p2)
           se drad<=prec pare; rad-=drad; drad*=0.25
   colrad(j)=rad
   wgt(j)=scale·(1-cos²rad)/P_{Jmax-1}(rad)²      ! :258-259
   rcs2(j)=1/sin²rad = 1/cos²(lat)                 ! :260
glat(j)=90-colrad(j)·180/π ; glat(Jmax-j+1)=-glat(j)   ! :81-83
```

✅ `Σwgt=1` (hemisfério). Com `Jmax=96`: `cos(colrad)` e `wgt` coincidem com
`numpy.polynomial.legendre.leggauss(96)` a `7e-15`/`4e-14`;
`glat(1)=88.57216851400727°`, `glat(48)=0.9326299678376557°`.
Há **três** cópias do cálculo (`Legendre.f90:216`, `LinearInterpolation.f90:114-149`,
`AreaInterpolation.f90:465-498`) com tolerâncias diferentes → ver A2 (§14).

### 3.4 Norte→sul no I/O

O binário GrADS tem a **primeira linha no sul** (`ydef Jm linear -90.0 Dxy` +
`fwrite -be`, `GrADSScripts/junta_atm.gs`). Todo o I/O inverte as linhas:
`TempInp(:,JmaxInp:1:-1,k)=REAL(WorkInp,r8)` (`SpectralICn.f90:116`; idem
`:99,107,126,136,146`), `SSTFill.f90:33,38-40`, `Ozone.f90:33,36`,
`TopoNewOut(i,j)=WorkOut(i,JmaxOut+1-j)` (`:85-89`), eixos `.ctl` com
`glat(JmaxOut:1:-1)` (`:549,672`). **Regra:** internamente `j=1↔+90°` sempre; a
inversão mora no I/O. Longitude não inverte (`i=1↔0°E`); exceção: saída do SST
(`Lon0Out=0.5, Lat0Out=89.5`).

---

## 4. Interpolação horizontal bilinear (`LinearInterpolation.f90`)

Usada se `ImaxInp<=ImaxOut`. Entrada regular, saída **gaussiana**
(`:31-51`), sem consultar as flags.

Pesos (`:180-294`): para cada `i` de saída, `iLon=MINLOC(|LonIn-LonOut(i)|)`
(empate → primeiro); `dLon<0` → Lower=`iLon-1` (ou `Idim` com
`LowerLonIn=LonIn(Idim)-360`), `=0` → iguais, `>0` → Upper=`iLon+1` (ou 1 com
`+360`) (`:196-224`). Periodicidade por **±360°** nos vizinhos, sem coluna extra.
Para `j`: `jLat=MINLOC(|LatIn-LatOut(j)|)`; **clamp** `LatOutA=min(max(LatOut,
LatIn(Jdim)),LatIn(1))` (`:233-234`); mesmo esquema de vizinhos (`:236-260`).

| caso | LL | LU | RL | RU |
|---|---|---|---|---|
| `ddx=0,ddy=0` | 1 | 0 | 0 | 0 |
| `ddx=0` | `1-dby` | `dby` | 0 | 0 |
| `ddy=0` | `1-dbx` | 0 | `dbx` | 0 |
| geral | `1-dbx-dby+dbx·dby` | `dby-dbx·dby` | `dbx-dbx·dby` | `dbx·dby` |

(`dbx=dx/ddx, dby=dy/ddy`; `:266-292`.)
Aplicação (`:297-397`): 4 cantos válidos → bilinear (`:333-337`); 0 →
`Undef`; 1 → cópia; 2 ou 3 → **média renormalizada** pela soma dos pesos válidos
(`:345-391`). A soma implementada é `LL·f_LL+LU·f_LU+RL·f_RL+RU·f_RU`
(ordem de arredondamento própria).
⚠️ `Undef=-999999` (`InputParameters.f90:350`) mas o GrADS grava `1e20`
(`junta_atm.gs`) ⇒ **na prática nada é indefinido** e todos os testes `/=Undef`
passam. Decisão A3 (§14).

---

## 5. Interpolação horizontal por área (`AreaInterpolation.f90`)

Usada se `ImaxInp>ImaxOut`. Remapeamento de 1ª ordem **renormalizado** (média
ponderada por sobreposição, não fluxo), sem máscara (`MaskInput=1`) e sem média
polar (`PolarMean=.false.`). Coordenadas em `sin(lat)` (lat) e rad (lon);
`eps=EPSILON(r4)` (`:164-166`) — tolerância em **precisão simples**.

Bordas de latitude (`:170-209`), `joi=JdimIn+JdimOut+2`, `FlagInput(5)=.false.`:

```
dlat=π/(JdimIn-1); dof=-(π+dlat)/2                     ! :198-199
W(joi+j)=sin(dof+dlat·(j-1)), j=2..JdimIn              ! :204-206
W(1+joi)=-1 ; W(JdimIn+1+joi)=+1                       ! :207-208
```

(lista crescente de `sin(lat)` com as `JdimIn+1` **bordas**).
Saída gaussiana (`:211-242`): `GaussianLatitudes(JdimOut/2)` (`:465-498`) dá as
colatitudes; `gLats(j)=90-rad·Gcol` e `gLats(JdimOut-j+1)=-gLats(j)`
(`:223-226`); bordas `sin((drltm+drltp)/2)` conforme `j<=lath`, `j>lath+1` ou
equador, com `W(1+joo)=-1`, `W(lath+1+joo)=0`, `W(JdimOut+1+joo)=+1`,
`joo=2·JdimIn+JdimOut+3`.
Fusão (`:269-310`): união ordenada; a cada elemento novo

```
wtlat(j3-1)   = W(j3)-W(j3-1)                       ! espessura      :276/288/299
mplat(j3-1,1) = j1-1 ; se FlagInput(1) -> JdimIn+2-j1               ! :277-278
mplat(j3-1,2) = j2-1 ; se FlagOutput(1)-> JdimOut+2-j2              ! :279-280
```

(`+2-j` = conversão "crescente em sin(lat)" → índice norte→sul); `latd=j3-2`.
Longitudes (`:314-460`): `W(i+ioi)=(dof+(i-1))·2π/IdimIn`, `dof=0.5,ici=0` se
centro (`FlagInput(4)`) senão `dof=0,ici=1`; fusão análoga com
`mplon(i3-1,1)=i1-ici` (deslocado de `IdimIn/2` se `FlagInput(2)=.false.`);
fecho do círculo em `:406-459` com `wtlon(i3-1)=2π+W(1)-W(i3-1)` (`:447`) ⇒
`SUM(wtlon)=2π`.

Aplicação (`DoAreaInterpolation`, `:82-143`):

```
FieldOutput=0; work=0
para j=1..latd, i=1..lond:
   se MaskInput(lni,lti)/=0 e FieldInput(lni,lti)/=Undef:
      FieldOutput(lno,lto) += FieldInput(lni,lti)·wlt·wtlon(i)
      work(lno,lto)        += wlt·wtlon(i)
FieldOutput /= work   (Undef se work==0)                 ! :128-136
se PolarMean: média zonal nas duas colunas de polo          ! :138-141
```

`j` externo e `i` interno (`:111-126`) fixam o arredondamento. Constantes são
preservadas exatamente, mas o operador **não** é o adjunto do conservativo
clássico (divide por `work`).

---

## 6. Interpolação vertical

### 6.1 `VertSigmaInter` (`VerticalInterpolation.f90:18-181`)

Arrays `(Im·Jm,K)` com coluna plana `c=i+(j-1)·Im` (`i` mais rápido; chamada em
`SpectralICn.f90:282-284`). `p1,p2` em hPa na ordem **superfície→topo**.

1. `z=-ln(p)` (`:93-113`) — crescente porque `p` decresce (exigência do
   `rsearch`, que aborta se `z1(i,1)>z1(i,Km1)`, `:358-368`).
2. `c1(:,:,1..4)=(u,v,t,q)`; tracers em `c1(:,:,3+n)` (`:96-108`).
3. `terp3(Im·Jm,Km1,Km2,3+Nt,z1,c1,z2,c2)` (`:120`).
4. Pós-processamento (`dz=z2(i,k)-z1(i,1)`, `:126-176`):
   - `u2=c2(:,1)`, `v2=c2(:,2)` sem correção;
   - `dz>=0`: `t2=c2(:,3)`, `q2=c2(:,4)`;
   - `dz<0` (abaixo do nível mais baixo, p>p1(1)):
     * `1<k<Km2` → **1-2-1**: `t2=0.25c2(k-1,3)+0.5c2(k,3)+0.25c2(k+1,3)`
       (`:144`, idem `q2` com coluna 4);
     * `k=1` → `t2=t1(1)`, `q2=q1(1,1)` (`:146-148`);
     * `k=Km2` e `Km2==Km1` → `t2=t1(k)`, `q2=q1(k,1)` (assimetria: usa `t1(k)`);
     * `k=Km2` e `Km2/=Km1` → gradiente constante (`:162-165`):

```
t2 = t1(1)·EXP(dltdz·dz)                        dltdz  = Gamma·Rd/Grav = -0.19029
q2 = q1(1,1)·EXP(dlpvdrt·(1/t2-1/t1(1)) - dz)   dlpvdrt= -Lc/Rv = -5417.1
```

       com `dz<0`: aquece para baixo a -6.5 K/km e mantém **UR constante**;
   - tracers `n>=2`: `q2(:,:,n)=c2(:,:,3+n)`.

Em produção `KmaxOut=28` vs `KmaxInp∈{23,26,37}` ⇒ só o ramo `Km2/=Km1` ocorre,
mas o ramo `Km2==Km1` existe (blocos `:140-152` quase idênticos ao `ELSE`).

### 6.2 `terp3`/`rsearch`/`bsrch` (`:184-417`)

`rsearch` devolve `l2(i,k2)` (`:307-328`): `0` = abaixo de `z1(i,1)`; `Km1` =
acima de `z1(i,Km1)`; senão `z1(l)<=z<z1(l+1)`, com `rc=0` quando `z` coincide
com um nó. `bsrch` (`:375-417`): `y<x(1)`→`indx=1,rc=1`; `y>x(n)`→`indx=n+1,rc=1`;
varredura linear para igualdade e para intervalo — não é busca binária, é O(n)
por ponto (no port use `np.searchsorted` respeitando as 4 regras).
Interpolação (`:235-281`): `k1=0` → `q1(:,1,:)`; `k1=Km1` → `q1(:,Km1,:)`;
`k1=1` ou `Km1-1` → **linear** (`ffa=(z2s-z1b)/(z1a-z1b)`,
`ffb=(z2s-z1a)/(z1b-z1a)`); `1<k1<Km1-1` → **Lagrange cúbico** com os 4 nós
`k1-1..k1+2`:

```
ffa=(z2s-z1b)/(z1a-z1b)·(z2s-z1c)/(z1a-z1c)·(z2s-z1d)/(z1a-z1d)   (idem b,c,d)
q2s = min(max(Σ ffi·q_i, min(q_b,q_c)), max(q_b,q_c))            ! monotonicidade :275-276
```

Pesos calculados por coluna `i` e reusados para todos os `Nm` campos.

### 6.3 Pressões de saída (`SpectralICn.f90:182-270`)

`gPresaux(k)=a_hybr(k)+b_hybr(k)·Ps` para `k=1..KmaxOut+1` (`:196-201,239-245`).
`IDSLInp=2` é **constante local** (`:40`) ⇒ sempre
`PresSig(k)=0.5·(gPresaux(k)+gPresaux(k+1))` (média **aritmética** das
interfaces, `:205-211,246-251`). O ramo "Phillips"
`((p_k^{RoCp1}-p_{k+1}^{RoCp1})/(RoCp1(p_k-p_{k+1})))^{RoCpR}` é código morto
(`:219-220,256-257`). Se `SigmaLevel`: `PresSig=Ps·SigLayerOut(k)`
(`:185-191,263-269`).
⚠️ No caso híbrido (default) a média de Exner de `NewSigma` **não** é usada; ela só
importa com `SigmaLevel=.true.` (A4).

---

## 7. Utils

### 7.1 `NewSigma` (`Utils.f90:18-53`)

```
SigInter(1)=1 ; SigInter(k+1)=SigInter(k)-DelSigmaOut(k) ; SigInter(KmaxOutp)=0  ! :38-42
SigLayer(k) = ((σ_k^{RoCp1}-σ_{k+1}^{RoCp1})/(RoCp1(σ_k-σ_{k+1})))^{1/RoCp}     ! :48-51
```

`SigLayer^{RoCp}` = média em σ de `σ^{RoCp}=(p/ps)^{RoCp}` (média de Exner da
camada). `SigInter(KmaxOutp)=0` forçado absorve o resíduo de `ΣDelSigma`.

### 7.2 `NewPs` (`Utils.f90:56-158`)

Chamada: `NewPs(ImaxOut,JmaxOut,KmaxInp,TopoPrs,TopoOut,PsInp,TempOut,PresOut,
PsOut,LnPsOut)` (`SpectralICn.f90:233-236`), com **`TopoPrs≡0`** (`:162`) e
`TempOut` = temperatura **virtual** (`:160`); `PresOut(i,j,k)=PresInp(k)`
(uniforme, `:163-169`). `goRd=Grav/Rd=0.0341635`.

```
(a) k=1, GammaLR=Gamma:                                              ! :84-107
    Zu = gTopoInp - Tvir(1)/GammaLR·((ps/p1)^(-GammaLR/goRd)-1)
    se TopoOut<=Zu: ps_out = p1·(1+GammaLR/Tvir(1)·(TopoOut-Zu))^(-goRd/GammaLR)   |G|>GEps
                    ps_out = p1·EXP(-goRd/Tvir(1)·(TopoOut-Zu))                    senão
    senão: ps_out=0 (SENTINELA) ; ls++
(b) k=2..KmaxInp, só onde ps_out==0:                                 ! :110-139
    GammaLR = -goRd·ln(Tvir(k-1)/Tvir(k))/ln(p(k-1)/p(k))
    |G|>GEps: Zu -= Tvir(k)/GammaLR·((p(k-1)/p(k))^(-GammaLR/goRd)-1)
    senão:    Zu += Tvir(k)/goRd·ln(p(k-1)/p(k))
    se TopoOut<=Zu: recalcula ps_out com p(k),Tvir(k) ; ls--
(c) ls>0 (topografia acima do topo): k=KmaxInp, G=0                     ! :142-153
    ps_out = p_top·EXP(-goRd/T_top·(TopoOut-Zu))
(d) LnPsOut = ln(ps_out/10)   → ln(cbar)                                ! :156
```

Como `GammaLR<0` e `p(k-1)>p(k)`, `Zu` **cresce** com `k` ⇒ acha o primeiro nível
acima da topografia ✅ (`Zu(1)≈+108 m` com ps=1013.25, p1=1000, T=280 K).
Consequências de `gTopoInp=0`: o "solo" é `z=0` (referencial fictício, não a
topografia interpolada); `PsInp` é **pressão ao nível do mar** (`psnm`/hPa em
`junta_atm.gs`) reprojetada hidrostaticamente para `TopoOut`; `Zu` persiste entre
as iterações de (b) (é a caminhada vertical por coluna); o sentinela 0 é ambíguo
com um valor físico nulo (impossível).
`GetERAPsLc` (`PsLcERA/GetERAPsLc.f90:67-139`) é a mesma rotina em **r4**, com
`PresRef` fixo de 13 níveis e **sem** o `/10` (saída em hPa).

### 7.3 `MoistVertExtrap` (`Utils.f90:161-308`)

Só se `KmaxInq<KmaxInp` (`SpectralICn.f90:149-152`), **antes** da interp. horizontal.

1. Tabela LIMS `SpHuStrat(8,19)` por `RESHAPE(...,(/8,19/))` ⇒ preenchimento
   coluna-maior: coluna 1 = 90N, linhas 1..8 = 80,70,…,10 hPa; valores ×1e-6
   (tabela em ppmv de massa → kg/kg, `:191-222`).
2. `Lat=Lat0·(1-2(j-1)/(JmaxInp-1))` (`:227`) — assume grade **regular**.
3. Saturação (`:233-237`): `Fac=6.1078·EXP(17.2693882·(T-273.16)/(T-35.86))`
   [hPa]; `q_sat=0.622·Fac/(p-0.378·Fac)`.
4. `q_sat·0.15 ≤ q ≤ q_sat` nos `KmaxInq` níveis baixos (`:243-248`).
5. `Dy=(90-Lat)/10+1`, `jDy=MIN(INT(Dy),8)`, `Dy-=jDy`; interpola em latitude:
   `q80 = S(1,jDy)(1-Dy)+S(1,jDy+1)Dy`, `q10 = S(8,·)` (`:252-256`).
   ⚠️ Para latitudes ao sul de -20° o índice satura e `Dy>1` (extrapolação, não
   interpolação) — bug latente, afeta só `FNL`/`GFS`.
6. `xLamb=(ln q80-ln q(KmaxInq))/(ln80-ln p(KmaxInq))`; para `p(k)>=80`:
   `q=q(KmaxInq)·(p/p(KmaxInq))^xLamb` (`:258-272`).
7. `80>p>=10`: bilinear na tabela em `p/10` e latitude, `Dp=8-p/10`,
   `kDp=MIN(INT(Dp),8)`, ordem de soma `(kDp,jDy),(kDp+1,jDy+1),(kDp,jDy+1),
   (kDp+1,jDy)` (`:274-284`). Acima de 10 hPa: `q=q10` constante (`:288-290`).
8. Reimpõe `q≤q_sat` e grava de volta em todos os níveis (inclusive o clamp de RH
   dos `KmaxInq` primeiros) — `:294-303`.

---

## 8. Transformada espectral

### 8.1 FFT (`Fourier.f90`, FFT99 do NCAR)

Convenção declarada em `:190-197`:

```
f(x)=f(0)+f(1)e^{-ix}+f(-1)e^{+ix}+…+f(max)e^{-max·ix} ;  a(2m+1)=Re f(m), a(2m+2)=Im f(m)
```

`isign=-1` grade→onda, `isign=+1` onda→grade (`:229-233`); o `1/n` só existe em
`fft99b` (`:411`). Portanto, para `n` pontos reais:

```
A(m) = (1/n)Σ_j f_j e^{+2πijm/n} = np.fft.ifft(f)[m]        (grade→onda)
f_j  = A(0) + 2Σ_{m=1}^{n/2-1} Re[A(m)e^{-2πijm/n}] + A(n/2)(-1)^j
```

No port: `np.fft.ifft(f,axis=0)[:Mend+1]` na ida e `np.fft.fft(A,axis=0)` na
volta. Um cosseno unitário dá `a(2m+1)=0.5` (não é a convenção do `rfft`×2).
Chamada: `fft991(fna,gna,inc=1,jump=Imx,n=Imax,lot=Ldim,isign=-1)`
(`Decomposition.f90:62-63`); as 2 posições de folga (`Imax+1,Imax+2`) **devem ser
zero** (`:52-55`) e são zeradas também na volta (`Fourier.f90:306-312`).
`trigs(3·Imax/2)` e `ifax(64)` (`Fourier.f90:45-46`) bastam para `Imax=192`.

### 8.2 `DectoSpherHarm` (`Decomposition.f90:27-71`)

```
para j=1..Jmaxhf: fn=f(:,j,:) ; fs=f(:,Jmax+1-j,:) ; zerar i>Imax
   fft991(fn) ; fft991(fs) ; SymAsy(fn,fs)          ! fn=sim, fs=anti :581-632
   FourierToSpecCoef(fn,fs,qf,Ldim,j)               ! :511-578
```

`SymAsy`: `a=a_n+a_s`, `b=a_n-a_s` (só `i=1..Imax`). `FourierToSpecCoef` por
nível `k`: `l=0`; para `nn=1..Mend1` (`nn-1=n-m`), `Mmax=2(Mend2-nn)`
(coeficientes **reais** da coluna); `s(mm)=fp` se `n-m` par e `fm` se ímpar;
`fln(:,k)+=s(:)·qlnw(:,lat)` com `qlnw=pln·wgt` (`:96,574`).
A paridade `n-m` é a correta (`P_n^m(-μ)=(-1)^{n-m}P_n^m`) e `fln` acumula a
quadratura sobre as latitudes (laço de latitudes externo).

### 8.3 Legendre: numeração, recursão e normalização (`Legendre.f90`)

Numeração (`:49-77`):

```
la0(mm,nn): nn=1..Mend1, mm=1..Mend2-nn   -> coeficientes complexos (Mnwv0)
la1(mm,nn): nn=1..Mend1,mm=1..Mend2-nn ; depois nn=2..Mend2, mm=1..Mend3-nn -> vetor pln (Mnwv1)
(m,n) <-> (mm=m+1, nn=n-m+1)
```

Ordinal do coeficiente complexo na saída de `DectoSpherHarm`:

```
p(m,n) = (n-m)(Mend+1) - (n-m)(n-m-1)/2 + m + 1            OK (1..2016)
```

`epslon` (`:176-213`): zeros na 1ª coluna e, para `nn>=2`,

```
eps[la1(mm,nn)] = sqrt(((n-1)²-m²)/(4(n-1)²-1))   ( = eps_{n-1}^m )
```

`pln2` (`:322-396`): `sinlat=cos(colrad)=μ`, `coslat=sin(colrad)`,
`x(mm)=sqrt(2mm+1)`, `y(mm)=sqrt(1+0.5/mm)`, `rthf=sqrt(0.5)=Cte`;

```
pln[mm] = rthf·prod_{j<mm}(coslat·y(j))       ! diagonal n=m, mm=m+1  :376-381
pln[mm+Mend1] = x(mm)·sinlat·pln[mm]          ! n=m+1                :383-385
pln[la1(mm,nn)] = (sinlat·pln[la1(mm,nn-1)] - eps[la1(mm,nn-1)]·pln[la1(mm,nn-2)])
                  / eps[la1(mm,nn)]                                   :386-394
```

✅ **Normalização verificada** contra `scipy.special.lpmv` (sem Condon-Shortley):
- `m=0`: `pln = sqrt(2n+1)·P_n /sqrt(2)`;
- `m>0`: `pln = sqrt(2(2n+1)(n-m)!/(n+m)!)·P_n^m /2`.
Isto é o harmônico real 4π-normalizado dividido por sqrt(2) (m=0) ou 2 (m>0), o
que o torna **ortonormal na quadratura gaussiana**:

```
sum_{j=1}^{Jmax/2} 2·wgt(j)·pln(m,n,j)² = 1     para todo (m,n)      OK (1±2e-12)
```

Consequências: `Cte=1/sqrt(2)` converte `a(0,0)` em média global (a(0,0)=sqrt(2)·média;
impressões usam `Cte·a`, `SpectralICn.f90:367,372,378,395`, `CheckICn.f90:32,36,45`);
a inversa (`SpecCoefToFourierScal/Vect`, `:833-978`, **não usada** aqui) soma
`Σ qln·fln` sem pesos porque `wgt` já foi embutido na ida.

### 8.4 Derivadas pseudo-espectrais (`plnder`, `:399-508`)

```
raa=wgt/a ; wcsa=rcs2·wgt/a ; rcs2=1/cos²(lat)                 ! :455-456
x[la1(mm,nn)] = m (nn=1) ; = n-1 (nn>=2)                        ! :458-468
der(la0(mm,nn)) = x[la1(mm,nn+1)]·eps[la1(mm,nn)]·pln[la1(mm,nn-1)]
                - x[la1(mm,nn)]·eps[la1(mm,nn+1)]·pln[la1(mm,nn+1)]   ! :469-479
der(la0(mm,1))  = -x[la1(mm,1)]·eps[la1(mm,2)]·pln[la1(mm,2)]        ! :480-482
dpln=der ; der=wcsa·der ; plnwcs=wcsa·m·pln                          ! :483-498
pln[la1(mm,nn)] *= x[la1(mm,nn)]·x[la1(mm,nn+1)]·raa   (laplaciano)  ! :499-506
```

⚠️ `der`/`dpln`/`plnwcs` têm comprimento `Mnwv0` e são indexados em ordem
**`la0`**, enquanto `pln` usa `la1` (coincidem só na 1ª coluna) — erro fácil.
✅ `der` cru = `cos(lat)·dP/dlat` (razão constante em todos os `(m,n)` testados,
3 latitudes) ⇒ `der = wgt/(a·cos²φ)·cosφ·dP/dφ` é o operador meridional do
divergente/vorticidade na forma de **pseudo-vento** `(U,V)=(u,v)cosφ`;
`plnwcs = wgt·m·P/(a·cos²φ)` é o zonal.

### 8.5 `UVtoDivgVort` (`Decomposition.f90:74-133`) e o sinal de D/zeta

```
coslat(j)=SIN(colrad)=cos(lat)                                          ! :93
para j=1..Jmaxhf:
  un=u(:,j,:)·coslat ; us=u(:,Jmax+1-j,:)·coslat ; idem vn,vs             ! :97-100
  fft991 em un,us,vn,vs ; SymAsy(un,us) ; SymAsy(vn,vs)                  ! :114-123
  GetConv(us,un,vs,vn,qdiv,j) ; GetVort(us,un,vs,vn,qrot,j)              ! :124-125
qdiv = -qdiv                                                             ! :131
```

`GetConv` (`:635-731`), por nível: (i) termo meridional com `qder`:
`s=bm` se `n-m` par, senão `bp` (`:691-705`), `fln+=s·qder`; (ii) termo zonal com
`qlnwcs`: para `mm=1,3,5,…`, `s(2l-1)=a(mm+1)`, `s(2l)=-a(mm)` (`a` = `ap` se
`n-m` par, senão `am`) ⇒ multiplicação por `-i·m` (`:710-725`).
`GetVort` (`:734-830`) espelha: meridional com `am/ap` (`:789-806`), zonal com
`s(2l-1)=-a(mm+1)`, `s(2l)=+a(mm)` sobre `bp/bm` ⇒ `+i·m` (`:809-827`).
⚠️ **A1 (§14).** As fórmulas fixam `qrot`/`qdiv` a menos do sinal **relativo**
entre os dois termos (depende da integração por partes com `1/cos²φ` e da
convenção da FFT). `qdiv=-qdiv` sugere que `GetConv` produz convergência (o
comentário em `:128-131` confirma a intenção), mas `GetVort` não é negado.
Portar literalmente e validar com o teste de rotação sólida (§12).

### 8.6 `TransSpherHarm` e a ordenação do arquivo (`Decomposition.f90:136-183`)

Permutação de pares (real em `2p-1`, imaginário em `2p`, `:157-158,173-174`).
Com `ICase=1` (= `ICaseDec`): destino `l` (mm=1..Mend1, nn=1..Mend2-mm — m
externo) ← origem `la0(mm,nn)`. Logo a ordenação **final do arquivo** é

```
l(m,n) = m(Mend+1) - m(m-1)/2 + (n-m+1)        [1..Mnwv0]        OK
```

e a **interna** (saída de `DectoSpherHarm`, (n-m) externo) é `p(m,n)`.
Mend=62: `(0,0)→p=1,l=1`; `(0,1)→64,2`; `(1,1)→2,64`; `(0,62)→2016,63`;
`(62,62)→63,2016`; `(30,40)→616,1466`. ✅ bijetora (2016).
⇒ `Topography.TQ0062` é lido **sem** reordenação (`SpectralICn.f90:69-71`) e
escrito **sem** reordenação (`:331-335`): já está na ordem `l(m,n)`, a mesma da IC.

### 8.7 Fatores para auditoria

| fator | origem | citação |
|---|---|---|
| `1/n` | `fft99b` (só ida) | `Fourier.f90:411` |
| `1/sqrt2` (m=0), `1/2` (m>0) | `pln2` | `Legendre.f90:370-394` |
| `wgt` | `glats` | `Legendre.f90:259` |
| `Cte=1/sqrt2` | coeficiente↔média (impressão) | `SpectralICn.f90:367` |
| `cosφ` em u,v | `UVtoDivgVort` | `Decomposition.f90:93-100` |
| `1/a`, `1/cos²φ` | `plnder` | `Legendre.f90:455-456` |

---

## 9. *Smoothing* espectral (`SmoothCoef`, `SpectralICn.f90:412-461`)

```
rm=(2MendC)(2MendC+1) ; rn=(MendC-1)MendC ; red=SmthPerCut^(-(rm²/rn²)/Iter)
cx=-ln(red)/rm² ; ck=(EXP(cx·((n)(n+1))²))^Iter ; Weight(2mn-1)=Weight(2mn)=ck
qOut(mn,k) *= Weight(mn)
```

Filtro de Hoskins (`rmn=mns(mns-1)=n(n+1)`).
⚠️ **Nunca é chamada** (busca textual: nome só em `:412,432,438,461`;
`Iter`/`SmthPerCut` só no namelist e aqui) ⇒ parâmetros mortos e *smoothing* de
topografia **não aplicado** nesta versão. Implementar e não chamar, salvo A5.

---

## 10. SST e ozônio (resumo)

**SST** (`Sources/SST`): namelist `SSTAGCMNML` (`InputParameters.f90:64`);
entrada = 2 registros r4 (TSM em K, fração de gelo 0–1) de `ImaxInp×JmaxInp`
(`SSTFill.f90:24-31`) com linhas invertidas; saída **fixa** 360×180, 1 registro
r4 em K, grade descentrada `Lon0Out=0.5, Lat0Out=89.5, DLat=1`
(`InputParameters.f90:110-140`), `.ctl` com `yrev big_endian`.
Algoritmo (`SSTFill.f90:33-108`): `INTERP_LAND` 2× (preenchimento de buracos na
grade de entrada pelos vizinhos válidos, pesos = razão de **índices**;
`:257-359`, com código morto após `RETURN` em `:365-404` a não portar);
filtros físicos `SST>400→Undef, SST<271.40→271.40, 0<=SIce<=1, 1<SIce<2→1,
SIce>2→Undef`; bilinear para 360×180 (`SST/HorizontalInterpolation.f90`, igual a
§4, com o bug latente de `:48-52`: ramo `GridIn=='Gauss'` escreve em `LatOut`);
repetição dos filtros; **preenchimento zonal** por latitude com máscara `CMask`
em string e interpolação linear entre vizinhos válidos (`:57-105`; latitude toda
inválida → `SSTSeaIce`); gelo final `0.5<=SIce<=1.01 → SST=271.40` (`:107-108`).

**Ozônio** (`Sources/Ozone`): namelist `ERAOzonNML` inclui a grade de saída
(default 192×96×28, `:78-93`); entrada = R1 pressão de superfície (hPa) +
R2…R(KmaxInp+1) ozônio r4 (`Ozone.f90:28-38`); saída `KmaxOut` registros r4
(`:58-65`), nome `OZON<Case><DateICn>S.grd.G<J5><K3>`, unidade **kg/kg** de
razão de mistura (`:81-82`). Bilinear (`GridOut='Gauss'`; aqui o ramo gaussiano
está correto) e depois `DoVerticalInterpolation`
(`Ozone/VerticalInterpolation.f90:18-129`) = mesmo `terp3` **sem** tratamento de
temperatura: abaixo do domínio usa 1-2-1 e `c1(i,1,n)`; acima, quando
`KmaxOut/=KmaxInp`, usa `(c2+c1(KmaxInp))/2` (`:116-123`). Pressões de saída
`Sigma(k)·Psfc` (sigma puro, `Ozone.f90:50-52`), ordem superfície→topo.

**`ps_lc`** (`PsLcERA/GetERAPsLc.f90`): mesma `NewPs` em **r4**, `PresRef` fixo
de 13 níveis (`PsLcERA/InputParameters.f90:60-62`), `TopoRef=0`, saída em hPa
(sem `/10`).

---

## 11. Plataforma (endianness, RECL, precisão)

| item | comportamento | citação |
|---|---|---|
| precisão | r4/r8 explícitos; sem `-r8`/`-i8` | `Sources/makefiles/Makefile.linux_gnu:12` |
| endianness | **nativa**; nenhum `CONVERT=`/`-fconvert` no código | grep em `Sources/**/*.f90` |
| entrada | GrADS `fwrite -be` → big-endian r4 (`options big_endian`) | `GrADSScripts/junta_atm.gs` |
| saída IC | UNFORMATTED sequencial, payload r4 ⇒ marcadores de 4 bytes + endianness nativa | `SpectralICn.f90:310-318,337-353` |
| GrADS | DIRECT com `RECL=INQUIRE(IOLENGTH=…)` ⇒ RECL em **bytes** (gfortran/ifort) | `:471-473,594-596` |
| topo espectral | SEQUENTIAL unformatted, `Mnwv2` r4 | `:61-71` |
| topo em grade | DIRECT, 1 registro r4 `ImaxOut×JmaxOut` | `:73-84` |
| entrada 3D | DIRECT, `RECL=4·ImaxInp·JmaxInp` bytes | `:94-98` |

⚠️ (1) Se o compilador da jaci usar RECL em *palavras* (Cray antigo — há
`Makefile.cray_cray`/`Makefile.gnu_cray` na árvore), os arquivos existentes têm
RECL diferente: confirmar a plataforma de origem antes de fixar
`recl_bytes=4*Imax*Jmax`. (2) Marcadores de registro exigem `dtype='>f4'` explícito
se a referência for big-endian e o port rodar em x86. (3) `IOLENGTH` de um array
r4 dá 4·N: trocar o buffer para r8 dobra o RECL silenciosamente (é por isso que
`SSTFill.f90:110` usa `REAL(SST,r4)` **dentro** do `IOLENGTH`).

---

## 12. Verificações numéricas executadas (Python/numpy) e testes

1. `GetImaxJmax(62)=(192,96)`; `Imx=194`; `Jmaxhf=48` ✅
2. `Mnwv2=4032`, `Mnwv0=2016`, `Mnwv1=2079`, `Mnwv3=4158` ✅
3. `glats` vs `leggauss(96)`: `|dmu|<7.2e-15`, `|dw|<4.2e-14`, `sum w=1.0000000000007` ✅
4. `pln2` vs `lpmv`: razão exata `1/sqrt2` (m=0) e `1/2` (m>0) ✅
5. Ortonormalidade `sum_half 2·wgt·pln²=1` nos 2016 pares (`|erro|<2.2e-12`) ✅
6. `der` cru = `cosφ·dP/dφ` (3 latitudes × 8 pares) ✅
7. `p(m,n)` e `l(m,n)` conferidos contra os laços de `la0`/`TransSpherHarm` ✅
8. `Zu(1)≈+108 m` em `NewPs` (ps=1013.25, p1=1000, T=280) — sinal correto ✅

### Testes de aceitação recomendados

1. **Rotação sólida** (decisivo): `u=u0·cosφ, v=0` ⇒ `zeta=2u0·sinφ/a`, `D=0`;
   checar `qdiv≈0` e `qrot` no harmônico `(0,1)` com o sinal de `+2u0/a`.
2. **Divergência pura**: `v=v0·cosφ, u=0` ⇒ `D=-2v0·sinφ/a`, `zeta=0`.
3. **Ida-e-volta**: campo de um único harmônico → só o coeficiente `(m,n)`
   (com os fatores `sqrt2`/`2` previstos) e resíduos <1e-12.
4. **Média global**: `Cte·a(0,0)` = média de área a 1e-12.
5. **Vertical**: `T=T0-6.5e-3·z` com `p` hidrostática → cúbica a 1e-10 no
   interior, exatamente linear nos dois intervalos de borda.
6. **`NewPs` identidade**: `gTopoOut==gTopoInp` e `gPsfcInp` = pressão do nível 1
   ⇒ `gPsfcOut==gPsfcInp`.
7. **Bilinear vs área**: campo linear em lat/lon é reproduzido exatamente pelos
   dois no interior.

---

## 13. Receita de port

1. Grade de saída por `GetImaxJmax` (§3.1) — **não** `3*Mend+1`.
2. Gaussianas: `leggauss(JmaxOut)`, `mu=nós[::-1][:Jmax/2]`,
   `wgt=pesos[::-1][:Jmax/2]`, `colrad=arccos(mu)`, `rcs2=1/(1-mu²)`
   (diferença esperada ~1e-14 vs. legado).
3. Legendre: recursão exata de §8.3 (única forma de casar com `pln`) ou
   `lpmv(m,n,mu)·(-1)^m·sqrt((2n+1)(n-m)!/(n+m)!)·(1/sqrt2 se m=0 senão 1/2)`;
   cuidado com Condon-Shortley e com underflow (`pln(62,62)≈8e-100` a 88.57°).
4. FFT: `ifft`/`fft` com `axis=0` e sem as folgas de 2 pontos.
5. Preservar a ordem dos laços nos acúmulos (latitudes externas em
   `FourierToSpecCoef`; `j` externo em `DoAreaInterpolation`; `k2/n/i` em `terp3`).
6. Payload r4, cálculo r8 (o legado só calcula em r4 em `GetERAPsLc`).
7. Funções puras; índices `p(m,n)`/`l(m,n)` pré-computados uma vez (vetores 2016).
8. Constantes de grade (`Flag*`, `MaskInput`, `PolarMean`) como argumentos default.

---

## 14. Ambiguidades que exigem decisão do autor

**A1. Sinal relativo de D e zeta (§8.5).** `qdiv=-qdiv` mas `qrot` não é negado; a
identidade de integração por partes decide o sinal relativo. Portar literalmente
(bug-compatível) e validar com rotação sólida, ou corrigir? Recomendo portar e
registrar o teste.

**A2. Qual rotina de latitudes gaussianas é canônica (§3.3).** Três cópias com
tolerâncias distintas (`Legendre.f90:216`, `LinearInterpolation.f90:114`,
`AreaInterpolation.f90:465`), diferenças ~1e-14 rad; a de `LinearInterpolation`
alimenta a interpolação e a de `Legendre` a transformada. Recomendo `leggauss`
com teste de 1e-12 em `glat`.

**A3. Indefinidos na entrada (§4).** Fortran compara com `-999999`, GrADS grava
`1e20` ⇒ hoje nada é indefinido. (a) preservar; (b) tratar `1e20`/`-999999` como
indefinido e usar a renormalização de §4; (c) NaN. Maior impacto numérico do
documento; pode melhorar a física.

**A4. Pressão de camada no híbrido (§6.3).** `IDSLInp=2` = média **aritmética**
das interfaces; a média de Exner de `NewSigma` só vale com `SigmaLevel=.true.`.
Manter (fidelidade) ou unificar (consistência)? Diferença típica <0.1 hPa nas
camadas baixas, mas muda `PresSig` e a interpolação vertical.

**A5. `SmoothCoef` não é chamado (§9).** Implementar e não chamar, ou ativar o
*smoothing* de topografia (mudaria `TopoSC` e a consistência com
`Topography.TQ0062`)?

**A6. Origem da topografia em `NewPs` (§7.2).** `TopoPrs=0` é uma referência
fictícia; se `Topography.TQ0062` existe, `TopoOut` é substituído pela topografia
espectral (`SpectralICn.f90:227-230`) **antes** de `NewPs`. Manter `gTopoInp=0`
ou usar a topografia interpolada?

**A7. RECL/endianness dos arquivos de referência (§11).** Sem saber a plataforma
(Cray/SX6 big-endian com RECL em palavras × x86 em bytes) o teste numérico pode
passar e o arquivo ainda ser ilegível pelo modelo.

**A8. `Undef` do SST em `EIH`/`FNL` (§4, §10).** Fortran `9.99E30` × script
`-9.99e+08`. Provavelmente o Fortran está errado e o valor é `-9.99e8` em todos
os casos.

---

## 15. Índice de citações

| regra | citação |
|---|---|
| namelist/defaults | `ICnAGCM/InputParameters.f90:111-144` |
| `Case`, níveis | `:161-204`, `:258-273` |
| dimensões, constantes | `:306-324`, `:329-350` |
| flags de grade | `:362-374`; `AreaInterpolation.f90:19-32` |
| `GetImaxJmax` | `:380-446` |
| orientação N→S | `SpectralICn.f90:85-89,99,107,116,126,136,146` |
| bilinear | `LinearInterpolation.f90:68-111,180-397` |
| área | `AreaInterpolation.f90:82-143,146-498` |
| vertical | `VerticalInterpolation.f90:18-417` |
| pressões de saída | `SpectralICn.f90:182-270` |
| `NewSigma`/`NewPs`/`MoistVertExtrap` | `Utils.f90:18-53,56-158,161-308` |
| FFT | `Fourier.f90:184-315` |
| decomposição | `Decomposition.f90:27-183` |
| Legendre (`pln2`,`epslon`,`glats`,`plnder`,`GetConv`,`GetVort`) | `Legendre.f90:176-830` |
| inversa (não usada) | `Legendre.f90:833-978` |
| `SmoothCoef` (não chamada) | `SpectralICn.f90:412-461` |
| SST / ozônio / `ps_lc` | `SST/*`, `Ozone/*`, `PsLcERA/*` |
