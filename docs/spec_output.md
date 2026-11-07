# Especificacao da SAIDA: arquivo de condicao inicial (nfcft) e diagnosticos GrADS

Documento de referencia para o port em Python (pacote **bam-ci**) do gerador de
condicoes iniciais do BAM 1.2 subsazonal (CPTEC/INPE). Cobre o arquivo binario de
condicao inicial escrito pelo programa Fortran **SpectralICn** (unidade **nfcft**),
os diagnosticos **.GrADS/.GrADS.ctl**, o verificador **CheckICn** e as datas.

## 0. Fontes e convencoes

Todas as citacoes **arquivo:linha** usam estes caminhos:

| abreviatura | caminho real |
|---|---|
| **ICnAGCM/SpectralICn.f90** | reference/GetICnAGCM_subsaeasonal/Sources/ICnAGCM/SpectralICn.f90 |
| **ICnAGCM/InputParameters.f90** | reference/GetICnAGCM_subsaeasonal/Sources/ICnAGCM/InputParameters.f90 |
| **ICnAGCM/InputArrays.f90** | reference/GetICnAGCM_subsaeasonal/Sources/ICnAGCM/InputArrays.f90 |
| **ICnAGCM/Legendre.f90** | reference/GetICnAGCM_subsaeasonal/Sources/ICnAGCM/Legendre.f90 |
| **ICnAGCM/Fourier.f90** | reference/GetICnAGCM_subsaeasonal/Sources/ICnAGCM/Fourier.f90 |
| **ICnAGCM/Decomposition.f90** | reference/GetICnAGCM_subsaeasonal/Sources/ICnAGCM/Decomposition.f90 |
| **ICnAGCM/Utils.f90** | reference/GetICnAGCM_subsaeasonal/Sources/ICnAGCM/Utils.f90 |
| **CheckICn/CheckICn.f90** | reference/GetICnAGCM_subsaeasonal/Sources/CheckICn/CheckICn.f90 |
| **CheckICn/InputParameters.f90** | reference/GetICnAGCM_subsaeasonal/Sources/CheckICn/InputParameters.f90 |
| **CheckICn/Makefile** | reference/GetICnAGCM_subsaeasonal/Sources/CheckICn/Makefile |
| **PDF** | reference/GetICnAGCM_subsaeasonal/doc/Reanalysis2CPTECAGCM.pdf (paginas do texto extraido com pdftotext) |
| **bamci/grads.py**, **bamci/ctl.py** | bam-ci/src/bamci/grads.py, bam-ci/src/bamci/ctl.py |

Convencoes: posicoes de array Fortran sao **1-based**; deslocamentos em bytes sao
**0-based**; "registro" = registro de arquivo Fortran UNFORMATTED sequencial.
Nada neste documento foi obtido executando o legado (nao ha arquivo de IC de
exemplo no espelho): tudo vem da leitura do codigo e do PDF.

## 1. Envelope binario do arquivo de IC

### 1.1 Formato

O arquivo e aberto como **UNFORMATTED sequencial** e cada **WRITE** gera **um
registro** (ICnAGCM/SpectralICn.f90:310-318, 322, 327, 337-353):

    IF(.not.SigmaLevel)THEN
    OPEN (UNIT=nfcft, FILE=TRIM(DirICnOut)//TRIM(DataOut)//'_hibrid', STATUS='REPLACE', &
                      FORM='UNFORMATTED')
    ELSE
    OPEN (UNIT=nfcft, FILE=TRIM(DirICnOut)//TRIM(DataOut)//'_sigma', STATUS='REPLACE', &
                      FORM='UNFORMATTED')
    END IF

Consequencias obrigatorias para o port:

1. **IEEE 32 bits, big endian** ("FORTRAN sequential IEEE 32 bits big endian file",
   PDF p.1). O legado rodava em maquina big endian (Configure.ksh:44,53-56 e
   Sources/makefiles/Makefile.SUPER-UX_SX6); os makefiles Linux/Intel **nao** tem
   opcao de conversao de byte order (Sources/makefiles/Makefile.linux_gnu:13,
   Makefile.tupay_intel:13), ou seja, num build little-endian o proprio legado
   geraria arquivo little-endian — o Python deve **forcar big endian sempre**.
2. **Marcadores de registro de 4 bytes** (comprimento do registro) antes e depois
   de cada registro. E o que a nota 2 do PDF (p.2) chama de "2 control words: one is
   at the beginning and the other at the end of each record"; a formula de tamanho
   total do PDF usa exatamente "+2" palavras por registro. Todos os registros sao
   muito menores que 2 GB, logo nao ha marcador de 8 bytes.
3. **Um WRITE = um registro**; nao existe registro extra, cabecalho magico,
   checksum, indice ou metadado de dimensao. O leitor (o BAM) precisa saber
   **mend** e **kmax** de antemao — na pratica eles vem do **nome do arquivo**
   (secao 6, item 7) e do namelist do modelo.
4. Tipos: **i4** = SELECTED_INT_KIND(9) = inteiro de 4 bytes;
   **r4** = SELECTED_REAL_KIND(6) = real IEEE de 4 bytes
   (ICnAGCM/InputParameters.f90:27-29). Nao ha r8 no arquivo: todas as escritas
   passam por REAL(...,r4) (ICnAGCM/SpectralICn.f90:320-321, 332-352).
5. Ordem dentro de cada registro: **ordem de memoria Fortran** (primeiro indice
   varia mais rapido). Para vetores 1-D isso e trivial; para campos 2-D (apenas no
   .GrADS) significa longitude mais rapida que latitude (secao 3).
6. Arquivo aberto com STATUS='REPLACE' — sobrescreve sem avisar.

### 1.2 Numero de registros e blocos

K = KmaxOut (numero de camadas), M = MendOut (truncamento). Sao **3 + 4*K**
registros, na ordem:

| bloco | registros (1-based) | qtd | conteudo | tamanho do payload |
|---|---|---|---|---|
| cabecalho | 1 | 1 | ForecastDay..(a,b) ou (sigma) | 48+8K (hibrido) / 44+8K (sigma) |
| topografia | 2 | 1 | TopoSC (m) | 4*Mnwv2 |
| ln(ps) | 3 | 1 | LnPsSC (cbar) | 4*Mnwv2 |
| temperatura virtual | 4 .. 3+K | K | TvirSC(:,k) (K) | 4*Mnwv2 cada |
| divergencia/vorticidade | 4+K .. 3+3K | 2K | DivgSC(:,k), VortSC(:,k) **alternados** por nivel | 4*Mnwv2 cada |
| umidade especifica | 4+3K .. 3+4K | K | SpHuSC(:,k) (kg/kg) | 4*Mnwv2 cada |

com **Mnwv2 = (M+1)*(M+2)** (ICnAGCM/InputParameters.f90:306-311). Para
M=62/K=28: registros 4..31 (T), 32..87 (div/vort, 32=div k=1, 33=vort k=1, ...,
86=div k=28, 87=vort k=28), 88..115 (q). Total 115 registros.

Ordem exata dos blocos (ICnAGCM/SpectralICn.f90:337-353):

    WRITE (UNIT=nfcft) WorkSC                     ! topografia
    WRITE (UNIT=nfcft) WorkSC                     ! ln(ps)
    DO k=1,KmaxOut                                ! temperatura virtual
       WRITE (UNIT=nfcft) WorkSC
    END DO
    DO k=1,KmaxOut                                ! divergencia e vorticidade
       WRITE (UNIT=nfcft) WorkSC                  ! divergencia
       WRITE (UNIT=nfcft) WorkSC                  ! vorticidade
    END DO
    DO k=1,KmaxOut                                ! umidade especifica
       WRITE (UNIT=nfcft) WorkSC
    END DO

A divergencia e escrita **antes** da vorticidade, dentro do mesmo laco de nivel —
nao ha dois blocos separados. O PDF descreve o mesmo (PDF p.1-2, "Records from
kmax+4 to 3*kmax+3").

### 1.3 Registro 1 — cabecalho (completo)

Escrito em ICnAGCM/SpectralICn.f90:322 (hibrido) e 327 (sigma):

    IF ( .not.SigmaLevel ) THEN
       HybIOut=REAL(a_hybr * 100._r4, r4)
       HybLOut=REAL(b_hybr, r4)
       WRITE (UNIT=nfcft) ForecastDay, TimeOfDay, DateInp, DateCur, HybIOut, HybLOut
    ELSE
       SigIOut=REAL(SigInterOut, r4)
       SigLOut=REAL(SigLayerOut, r4)
       WRITE (UNIT=nfcft) ForecastDay, TimeOfDay, DateInp, DateCur,SigIOut,SigLOut
    END IF

Layout (payload, antes dos marcadores):

| offset | tamanho | campo | tipo | conteudo |
|---|---|---|---|---|
| 0 | 4 | ForecastDay | i4 | dia de previsao; **0** na condicao inicial (ICnAGCM/InputParameters.f90:326) |
| 4 | 4 | TimeOfDay | r4 | hora do dia de previsao; **0.0** na CI (ICnAGCM/InputParameters.f90:327) |
| 8 | 16 | DateInp(1..4) | 4 x i4 | **hora, mes, dia, ano** (nesta ordem!) |
| 24 | 16 | DateCur(1..4) | 4 x i4 | mesma ordem; = DateInp na CI |
| 40 | 4*(K+1) | HybIOut ou SigIOut | r4 | hibrido: **a(k) em Pa**; sigma: interfaces sigma (K+1) |
| 40+4*(K+1) | 4*(K+1) hib. / 4*K sigma | HybLOut ou SigLOut | r4 | hibrido: **b(k)** (K+1); sigma: camadas (K) |

Tamanhos do payload: **hibrido = 48 + 8K** bytes; **sigma = 44 + 8K** bytes.
Para K=28: 272 e 268 bytes. Note a **assimetria**: no caso hibrido **b tambem tem
K+1 valores** (niveis de interface), no caso sigma as camadas tem apenas K
(ICnAGCM/InputArrays.f90:54-56: SigIOut(KmaxOutp), SigLOut(KmaxOut),
HybIOut(KmaxOutp), HybLOut(KmaxOutp); ICnAGCM/InputParameters.f90:209-210
aloca a_hybr e b_hybr com KmaxOut+1).

O PDF (p.1) documenta apenas o caso sigma: "Record 1: ForecastDay, TimeOfDay,
DateInp(4), DateCur(4), SigIOut(kmax+1), SigLOut(kmax)" com tamanho
((11+2*kmax)+2)*nbyte. Isso **nao** vale para o arquivo _hibrid (o PDF esta
desatualizado nesse ponto).

#### 1.3.1 Caso hibrido (_hibrid)

- a e b vem de **DataLevels/HybridLevels.L<kkk>**, K+1 linhas, leitura
  list-directed de dois reais por linha, com **a convertido de Pa para mbar
  (a/100)** e depois **reconvertido para Pa (a*100)** na escrita
  (ICnAGCM/InputParameters.f90:216-217, 266-273; ICnAGCM/SpectralICn.f90:319-321).
  Logo: o arquivo contem **a em Pa** e **b adimensional** — o valor numerico e o
  mesmo do arquivo de niveis, com erro de ida-e-volta de 1e-2 (float64)
  irrelevante em float32.
- O modelo reconstroi p(nivel) = a(k) + b(k)*ps (ICnAGCM/SpectralICn.f90:198-199,
  242-243).
- Uso interno das pressoes de camada: **IDSLInp = 2 fixo** (ICnAGCM/SpectralICn.f90:40)
  = media aritmetica das duas interfaces (linhas 203-211, 246-251). O arquivo nao
  leva essa informacao: apenas a e b. O modelo repete a sua propria escolha.
- Nome do arquivo: ver secao 6.1 (sufixo _hibrid).

#### 1.3.2 Caso sigma (_sigma)

- Sigma de **interface** SigIOut = SigInterOut(K+1): SigInterOut(1)=1.0,
  SigInterOut(k+1)=SigInterOut(k)-DelSigmaOut(k), SigInterOut(K+1)=0.0
  (ICnAGCM/Utils.f90:37-42). Ou seja: **K+1 valores decrescentes, de 1 a 0**.
- Sigma de **camada** SigLOut = SigLayerOut(K) **nao** e o ponto medio
  aritmetico: e o valor "tipo Phillips"
  ((s_k^(Rd/Cp+1) - s_(k+1)^(Rd/Cp+1)) / ((Rd/Cp+1)*(s_k - s_(k+1))))^(Cp/Rd)
  (ICnAGCM/Utils.f90:44-51). Reproduzir a formula exatamente, com Rd e Cp do
  legado (Rd=287.05, Cp=1004.6; ICnAGCM/InputParameters.f90:334-337).
- DelSigmaOut vem de **DataLevels/DeltaSigma.L<kkk>** lido com formato (5F9.6)
  (ICnAGCM/InputParameters.f90:216, 262-264). O legado **sempre** le os dois
  arquivos (DeltaSigma e HybridLevels), independente de SigmaLevel
  (ICnAGCM/InputParameters.f90:262-273).
- ATENCAO: no caso sigma o registro 1 tem 4 bytes **menos** que no hibrido. Um
  leitor que assume o layout errado desalinha todos os registros seguintes.

### 1.4 Tamanho total

    bytes = 8 + L1 + (2+4K) * (8 + 4*Mnwv2)

com L1 = 48+8K (hibrido) ou 44+8K (sigma); o 8 externo e o par de marcadores do
registro 1 e (8 + 4*Mnwv2) e um registro de dados completo (marcador + payload +
marcador). Formula equivalente do PDF (p.2, caso sigma, com o typo "nbyte2"):

    tsICn = ((11+2*kmax)+2 + (2+4*kmax)*((mend+1)*(mend+2)+2)) * 4

Exemplo numerico **M=62, K=28** (defaults do namelist,
ICnAGCM/InputParameters.f90:131-132): Mnwv2 = 63*64 = **4032** coeficientes =
16128 bytes por registro; 114 registros de dados = 1838592 bytes.

| | _hibrid | _sigma |
|---|---|---|
| payload do registro 1 | 272 | 268 |
| inicio do registro 2 (offset) | 280 | 276 |
| passo entre registros de dados | 16136 | 16136 |
| tamanho total do arquivo | **1839784** | **1839780** |

Deslocamento do registro r (1-based, r>=2), caso hibrido:
offset(r) = 280 + (r-2)*16136; o payload comeca em offset(r)+4.

### 1.5 Escritor minimo (pseudocodigo Python)

```python
import struct
import numpy as np

def _rec(payload: bytes) -> bytes:                 # marcadores de 4 bytes BE
    n = struct.pack(">i", len(payload))
    return n + payload + n

def _f4(a) -> bytes:                               # IEEE binary32 big endian
    return np.ascontiguousarray(np.asarray(a, dtype=">f4")).tobytes()

def _i4(a) -> bytes:
    return np.ascontiguousarray(np.asarray(a, dtype=">i4")).tobytes()

def write_icn(path, *, mend, kmax, forecast_day, time_of_day,
              date_inp, date_cur, vertical, topo, lnps, tvir, divg, vort, sphu):
    mnwv2 = (mend + 1) * (mend + 2)
    for arr in (topo, lnps):
        assert np.asarray(arr).size == mnwv2
    for arr in (tvir, divg, vort, sphu):
        assert np.asarray(arr).size == kmax * mnwv2
    head = (_i4([forecast_day]) + _f4([time_of_day]) +
            _i4(date_inp) + _i4(date_cur) + _f4(vertical["a_or_sig_i"]) +
            _f4(vertical["b_or_sig_l"]))
    with open(path, "wb") as fh:                   # um write por registro
        fh.write(_rec(head))
        fh.write(_rec(_f4(topo)))
        fh.write(_rec(_f4(lnps)))
        for k in range(kmax):
            fh.write(_rec(_f4(tvir[k])))
        for k in range(kmax):
            fh.write(_rec(_f4(divg[k])))
            fh.write(_rec(_f4(vort[k])))
        for k in range(kmax):
            fh.write(_rec(_f4(sphu[k])))
```

Teste de aceitacao do envelope, sem o BAM: tamanho exato (secao 1.4), cada
marcador igual ao payload, e leitura de volta com
np.frombuffer(..., dtype=">f4") reproduzindo os mesmos arrays.

## 2. Convencao dos coeficientes espectrais

### 2.1 Tamanho e triangularidade

- Cada registro de dados tem exatamente **Mnwv2 = (M+1)*(M+2)** reais
  (ICnAGCM/InputParameters.f90:306-311), isto e **(M+1)*(M+2)/2 = Mnwv0 pares**
  (parte real = cosseno, parte imaginaria = seno).
- Truncamento **triangular**: pares (m,n) com 0 <= m <= n <= M. Para M=62:
  4032 reais = 2016 pares.
- **Nao ha coeficiente extra, padding, mascara ou cabecalho** dentro do registro:
  nem o par (m=0, n=0) e omitido, nem o seno de m=0 (que e identicamente zero) e
  suprimido. Arrays internos maiores existem apenas para a transformada de campos
  vetoriais (Mnwv1/Mnwv3, ICnAGCM/InputParameters.f90:310-312) e **nao** vao para
  o arquivo.

### 2.2 Ordenacao (m,n): a numeracao la0

A ordem do arquivo e a numeracao **la0** de ICnAGCM/Legendre.f90:49-65, usada em
ICnAGCM/Decomposition.f90:144-181 para converter para a numeracao la1 (por grau)
usada pela base de Legendre. Em termos de (m,n):

- o indice de par p (1-based) e
  **p(m,n) = d*(M+1) - d*(d-1)/2 + m + 1**, com **d = n - m**;
- dentro de uma diagonal d, m cresce;
- a diagonal d tem M+1-d elementos;
- real em **2p-1**, imaginario em **2p** (Fortran); em Python (0-based):
  i_real = 2*(p-1), i_imag = 2*(p-1)+1.

Para M=62 (ilustrativo, mostra que a diagonal e o criterio, nao o grau):

| (m,n) | d | p |
|---|---|---|
| (0,0) | 0 | 1 |
| (1,1) | 0 | 2 |
| (62,62) | 0 | 63 |
| (0,1) | 1 | 64 |
| (1,2) | 1 | 65 |
| (0,62) | 62 | 2016 |

Ou seja: **os primeiros 63 pares sao os setoriais (m=n=0..62)**; o segundo
coeficiente do modo zonal (m=0, n=1) esta na posicao **64**, nao na 2. Qualquer
port que assuma "m=0 contiguo no inicio" produz arquivo invalido.

Equivalente vetorizado:

```python
p = np.empty((M + 1, M + 1), dtype=int)   # indice de par (1-based)
for mm in range(M + 1):                     # ordem = m
    for nn in range(mm, M + 1):             # grau n >= m
        d = nn - mm
        p[mm, nn] = d * (M + 1) - d * (d - 1) // 2 + mm + 1
# par p -> reais[2*(p-1)] (cosseno), reais[2*(p-1)+1] (seno)
```

### 2.3 Normalizacao e fator Cte

O arquivo guarda **o array interno do modelo**, saida de
DectoSpherHarm/TransSpherHarm (ICnAGCM/SpectralICn.f90:296-308), sem qualquer
rescalonamento na escrita. As invariantes verificaveis:

- Funcoes de Legendre normalizadas com **P(0,0) = 1/sqrt(2)**
  (ICnAGCM/Legendre.f90:370, rthf=SQRT(0.5_r8); primeira linha de pln2,
  ICnAGCM/Legendre.f90:376-381).
- Pesos gaussianos com soma 2 (ICnAGCM/Legendre.f90:216-263: latitudes e pesos de
  Gauss), de modo que a media global de um campo f e **0.5 * soma_j w_j * f_j**.
- A FFT grade->onda inclui o fator 1/N (ICnAGCM/Fourier.f90:411, 418-419), logo
  para um campo constante de valor c o primeiro coeficiente e **c*sqrt(2)**.
- Consequencia pratica, usada pelo proprio verificador: **media global =
  c1 * Cte**, com **Cte = 1/sqrt(2) = 0.70710678**
  (CheckICn/CheckICn.f90:31-32, 36, 44-45, 51-52, 61-63;
  ICnAGCM/InputParameters.f90:348). Inversamente **c1 = sqrt(2) * media**.
- O seno de m=0 (posicoes 2p com m=0) e **exatamente 0.0** e serve de teste de
  sanidade do empacotamento.

Nao existe, no arquivo nem no codigo, fator de conversao adicional: se o port
errar a normalizacao de Legendre, o BAM aceita o arquivo (o envelope continua
valido) mas os campos reconstruidos ficam errados por um fator global. O
criterio de equivalencia tem de ser por coeficiente (secao 6.1).

### 2.4 Unidades e grandezas por campo

| registro | campo | unidade no arquivo | evidencia |
|---|---|---|---|
| 2 | TopoSC | **m** (coeficientes espectrais da topografia) | ICnAGCM/SpectralICn.f90:365-368; PDF p.1 |
| 3 | LnPsSC | **ln(ps/10 hPa)** = ln(ps em cbar); ps=10*exp(c) em hPa | ICnAGCM/Utils.f90:155-156; CheckICn/CheckICn.f90:36 |
| 4..3+K | TvirSC | **K**, temperatura **virtual** T*(1+cTv*q), cTv=Rv/Rd-1 | ICnAGCM/SpectralICn.f90:285, 160; ICnAGCM/InputParameters.f90:341 |
| 4+K..3+3K | DivgSC | **s^-1**, divergencia (sinal = -convergencia) | ICnAGCM/Decomposition.f90:128-131; PDF p.1 |
| 4+K..3+3K | VortSC | **s^-1**, vorticidade, sem troca de sinal | ICnAGCM/Legendre.f90:734-830 |
| 4+3K..3+4K | SpHuSC | **kg/kg** (o relatorio multiplica por 1000 para g/kg) | ICnAGCM/SpectralICn.f90:392-397; PDF p.1 |
| 1 | a / b | **Pa** / adimensional | ICnAGCM/SpectralICn.f90:313, 320-321 |
| 1 | sigma | adimensional, 1 -> 0 (interface), camadas por Phillips | ICnAGCM/Utils.f90:37-51 |

O PDF (p.1) diz "logarithm of surface pressure in **kPa**". O codigo mostra que e
**cbar**: gLnPsOut = LOG(gPsfcOut/10.0) com ps em hPa (ICnAGCM/Utils.f90:155-156),
e o verificador imprime 10*EXP(c1*Cte) como hPa (CheckICn/CheckICn.f90:36).
1 cbar = 10 hPa. Errar esse fator desloca a pressao em 10x.

### 2.5 Como o BAM le esse arquivo

O fonte do modelo nao esta neste espelho; o espelho do leitor e o **CheckICn**,
que reproduz a sequencia de leitura do BAM (CheckICn/CheckICn.f90:14-18, 30, 34,
42-64) e o PDF (p.1-2). O leitor:

1. abre por nome exato (nao procura, nao valida sufixo);
2. le o registro 1 em um unico READ com a lista
   ForecastDay, TimeOfDay, DateInitial(4), DateCurrent(4), SigInter(K+1), SigLayer(K);
3. le 1 registro de topografia, 1 de ln(ps), K de Tvir, 2K alternados
   (div, vort), K de SpHu, cada um em um array de Mnwv2 reais;
4. **nao** tem cabecalho de dimensoes: se M ou K do arquivo nao coincidirem com o
   que o modelo espera, o resultado e erro de I/O (registro curto/EOF) ou lixo
   silencioso (byte order errado);
5. usa os arrays como estado inicial: os diagnosticos do gerador
   (ICnAGCM/SpectralICn.f90:356-398) sao calculados exatamente sobre esses arrays,
   e sao os valores de referencia para o port.

## 3. Arquivo .GrADS (diagnostico)

### 3.1 Quando e escrito

Se GrADS=.TRUE. (default, ICnAGCM/InputParameters.f90:133):
GetGrADSInp e chamado **antes** da interpolacao vertical (sobre o dado de
entrada ja interpolado horizontalmente para a grade do modelo, nos niveis de
pressao de entrada; ICnAGCM/SpectralICn.f90:272-279) e GetGrADSOut **depois**
(estado final nos niveis sigma do modelo; ICnAGCM/SpectralICn.f90:287-294). Se
GrADSOnly=.TRUE., o programa **para** depois dos diagnosticos, sem escrever o
arquivo de IC (ICnAGCM/SpectralICn.f90:278, 293) — util para testar a cadeia.

### 3.2 Binario

- Arquivos: DirICnOut + DataOut + '.GrADS' (ICnAGCM/SpectralICn.f90:472-473) e
  DirICnOut + DataOut + '.GrADSOut' (ICnAGCM/SpectralICn.f90:595-596).
- **ACCESS='DIRECT'** com RECL = INQUIRE(IOLENGTH) gWorkOut
  (ICnAGCM/SpectralICn.f90:471-473, 594-596), onde gWorkOut e REAL*4
  (ImaxOut,JmaxOut) (ICnAGCM/InputArrays.f90:22-23, 83). Portanto:
  **float32 cru, sem marcadores de registro, sem cabecalho**, um campo 2-D por
  registro, 4*ImaxOut*JmaxOut bytes por registro (73728 bytes para
  ImaxOut=192, JmaxOut=96).
- **Big endian** explicito no descritor: options yrev big_endian
  (ICnAGCM/SpectralICn.f90:544, 667).
- Ordem dentro do registro: Fortran column-major, **longitude mais rapida**,
  depois latitude, com **j=1 = Polo Norte** (glat(1)=+90;
  ICnAGCM/Legendre.f90:80-83). O yrev do .ctl diz ao GrADS que a primeira linha
  do arquivo e o **ultimo** nivel do ydef — e o ydef e escrito de sul para norte
  (glat(JmaxOut:1:-1), ICnAGCM/SpectralICn.f90:548-549), fechando a consistencia:
  **no disco as linhas vao de norte para sul**.
- Registros, nesta ordem (ICnAGCM/SpectralICn.f90:480-515 para a entrada,
  603-638 para a saida):

| # | .GrADS (entrada) | .GrADSOut (saida) |
|---|---|---|
| 1 | topografia (m) | topografia (m) |
| 2 | pressao de superficie PsOut (**hPa**, nao log) | idem |
| 3 .. 2+KmaxInp | T virtual da analise (K), KmaxInp niveis | 3 .. 2+KmaxOut: T virtual nos niveis sigma |
| 3+KmaxInp .. 2+2KmaxInp | umidade especifica (kg/kg) | idem |
| 3+2KmaxInp .. 2+3KmaxInp | vento zonal (m/s) | idem |
| 3+3KmaxInp .. 2+4KmaxInp | vento meridional (m/s) | idem |

Total de registros = 2 + 4*Kmax (KmaxInp no arquivo de entrada, KmaxOut no de
saida); bytes = 4*ImaxOut*JmaxOut*(2+4*Kmax). Divergencia e vorticidade **nao**
sao gravadas (blocos comentados, ICnAGCM/SpectralICn.f90:491-500, 614-623);
ozonio e tracers tambem estao comentados (516-532, 639-655).

### 3.3 Descritor .GrADS.ctl

Escrito como texto formatado (ICnAGCM/SpectralICn.f90:535-583 e 658-706):

    dset ^<DataOut>.GrADS
    options yrev big_endian
    undef 1e20
    xdef <ImaxOut> linear 0.0  <360/ImaxOut>          ! formato 2F12.6
    ydef <JmaxOut> levels                             ! 6F10.3, glat(Jmax:1:-1)
    zdef <KmaxInp> levels                             ! 6F10.3, PresInp (hPa)
    tdef 1 linear <Tdef> 1dy                          ! Tdef = 'hhzddmmmyyyy'
    vars 6                                            ! ou 6+NTracers (morta)
    topo   0 99 Topography        TQ0062L028 (m)
    pslc   0 99 Surface Pressure  TQ0062L028 (hPa)
    tvir  37 99 Virt Temperature  TQ0062L028 (K)
    umes  37 99 Specific Humidity TQ0062L028 (kg/kg)
    uvel  37 99 Zonal Wind        TQ0062L028 (m/s)
    vvel  37 99 Meridional Wind   TQ0062L028 (m/s)
    endvars

Detalhes que o port precisa copiar:

- dset usa apenas o **nome base** (o .ctl deve morar no mesmo diretorio do
  binario); ^ = relativo ao .ctl.
- a ordem das variaveis no vars **e** a ordem dos registros (topo, pslc, tvir,
  umes, uvel, vvel) — inclusive o tvir declarado depois de pslc.
- zdef de entrada = **PresInp** (hPa, lido de PressureLevels.*,
  ICnAGCM/InputParameters.f90:258-260); zdef de saida = **1000*SigLayerOut**
  (aproximacao p = sigma*1000 hPa, ICnAGCM/SpectralICn.f90:674) — nao e a pressao
  real de cada coluna.
- TruncOut no texto e a string de truncamento (TQ0062L028,
  ICnAGCM/InputParameters.f90:287-290).
- undef 1e20 e decorativo: nao ha valor indefinido nos campos gravados.
- bug cosmetico: a mensagem de erro do .GrADSOut imprime .GrADS.ctl
  (ICnAGCM/SpectralICn.f90:658-665) e vars suporta 6+NTracers sem escrever os
  tracers (553-581) — com NTracers=1 (fixo, ICnAGCM/InputParameters.f90:318) nao
  ha efeito.

### 3.4 Para que serve

- Inspecao visual rapida, no GrADS, do que entrou (.GrADS: analise interpolada
  horizontalmente, nos niveis de pressao de entrada) e do que saiu (.GrADSOut:
  estado final nos niveis sigma, **o mesmo** que gerou os coeficientes). E o unico
  "dump" legivel da cadeia legada.
- **Nao** e lido pelo BAM. Portanto, para o port, o .GrADS e opcional e de menor
  prioridade; se implementado, deve ser byte-compativel para permitir comparacao
  visual lado a lado com o legado.
- Serve tambem como oraculo de teste do port: .GrADSOut e a representacao em
  grade dos mesmos campos que foram transformados em coeficientes.

## 4. CheckICn: o que valida

Programa separado, compilado de CheckICn/CheckICn.f90 + CheckICn/InputParameters.f90
(CheckICn/Makefile:9-13; o alvo se chama SST por copy-paste, mas e o default).

### 4.1 O que ele faz

1. Monta o nome do arquivo a partir do namelist ChkICnNML (Mend, Kmax, DateICn,
   Case, Grid, DirICn, DirOut; CheckICn/InputParameters.f90:41) no padrao
   GANL<Case><yyyymmddhh>S.unf.T<Grid><mmmm>L<kkk>
   (CheckICn/InputParameters.f90:61-64).
2. Abre TRIM(DirICn)//NameICn como UNFORMATTED SEQUENTIAL OLD
   (CheckICn/CheckICn.f90:14-15).
3. Le o registro 1 em ForecastDay, TimeOfDay, DateInitial, DateCurrent, SigInter,
   SigLayer e **imprime** (CheckICn/CheckICn.f90:17-27).
4. Le topografia e imprime c1*Cte (media global, m) e MIN/MAX de qWork(2:)
   (CheckICn/CheckICn.f90:29-32).
5. Le ln(ps) e imprime a pressao media global 10*EXP(c1*Cte) em hPa e
   MIN/MAX dos demais coeficientes, rotulando-os "kPa"
   (CheckICn/CheckICn.f90:34-39).
6. Le K temperaturas virtuais, 2K divergencia/vorticidade e K umidades,
   imprimindo por nivel c1*Cte e MIN/MAX (a umidade e multiplicada por 1000,
   ou seja, o "g/kg" do rotulo; CheckICn/CheckICn.f90:41-64).

### 4.2 Tolerancias

**Nao existe nenhuma tolerancia no codigo.** CheckICn nao compara com nada, nao
calcula residuo, nao tem IF ... STOP, nao devolve codigo de erro: e um
**impressor** de estatisticas. Toda a validacao numerica e humana (comparar as 3
colunas com o que se espera).

### 4.3 O que ele considera erro

Somente erro de **I/O/Formato** (sem IOSTAT em nenhum ponto, logo o runtime
aborta):

- arquivo ausente ou nome montado errado (OPEN sem IOSTAT, linhas 14-15);
- registro mais curto que a lista de leitura, p.ex. Kmax do namelist maior que o
  do arquivo (o Fortran aborta ao ler alem do fim do registro);
- arquivo truncado / EOF antes do fim (Mend do namelist maior que o do arquivo,
  ou arquivo de outro truncamento);
- byte order/real size errados **nao** geram erro: produzem numeros absurdos
  impressos como se fossem validos (risco real ao comparar IC gerada em x86).

### 4.4 Limitacoes que o port deve contornar

1. **Nome truncado**: NameICn e CHARACTER(LEN=33)
   (CheckICn/InputParameters.f90:16) e o nome montado tem exatamente 33 caracteres
   (4+3+10+6+10). Como o gerador **sempre** acrescenta _hibrid ou _sigma
   (ICnAGCM/SpectralICn.f90:311, 315), o nome real tem ~40 caracteres e o
   CheckICn **nao consegue abrir o arquivo de verdade** — so abre copia/atual
   ligacao simbolica com o nome truncado. O port deve ter seu proprio leitor por
   nome completo.
2. **So o caso sigma**: as declaracoes sao SigInter(Kmax+1), SigLayer(Kmax)
   (CheckICn/InputParameters.f90:22-23, 80-83). Num arquivo _hibrid a leitura do
   registro 1 e "curta" (o Fortran permite ler menos itens que o registro e
   descarta o resto), de modo que ele **nao da erro** e imprime, sob os rotulos
   SigInter/SigLayer, os valores de **a** e dos K primeiros **b** — diagnostico
   silenciosamente errado.
3. **MIN/MAX poluidos**: MINVAL(qWork(2:)) inclui os slots de seno de m=0, que
   sao exatamente 0. Para campos de coeficientes positivos (tipico de Tvir) o
   minimo impresso sai **0**, nao o menor coeficiente fisico. Nao usar isso como
   criterio de aceitacao.
4. Se Mend/Kmax do namelist nao baterem com o arquivo, ver 4.3.

## 5. Datas, ForecastDay, TimeOfDay, DateInp/DateCur

### 5.1 Entrada de configuracao

Namelist ICnAGCMNML (ICnAGCM/InputParameters.f90:111-112, 144): MendOut,
KmaxOut, DateICn, Iter, SmthPerCut, Case, Grid, SigmaLevel, GrADS, GrADSOnly,
DirMainInp, DirMainOut. Defaults: MendOut=62, KmaxOut=28, DateICn='yyyymmddhh',
Case='E40', Grid='Q', SigmaLevel=.FALSE., GrADS=.TRUE., GrADSOnly=.FALSE.
(ICnAGCM/InputParameters.f90:131-142).

- DateICn e uma string de **10 digitos yyyymmddhh** lida por posicao fixa
  (I4)/(I2) (ICnAGCM/InputParameters.f90:146-149), **sem validacao** de faixa:
  mes 13 ou hora 99 passam e viram lixo nos arrays.
- DateInp(1)=Hour, DateInp(2)=Month, DateInp(3)=Day, DateInp(4)=Year
  (ICnAGCM/InputParameters.f90:155-158) — **hora primeiro, ano por ultimo**. Nao
  e a ordem ISO nem a do ECMWF.
- DateCur=DateInp (ICnAGCM/InputParameters.f90:159): na condicao inicial os dois
  sao identicos; o PDF explica que DateCur e a data corrente da previsao e so
  diverge depois, dentro do modelo (PDF p.1, nota 1).
- TimeGr = 'hh z dd mmm yyyy' montado com MonChar
  (ICnAGCM/InputParameters.f90:114-115, 150-154) e usado **apenas** na linha
  tdef do .GrADS.ctl (ICnAGCM/SpectralICn.f90:273-278, 288-293, 552, 675;
  formato de 12 caracteres '  z         ' com hour(1:2), day(4:5), mon(6:8),
  year(9:12)).
- ForecastDay e TimeOfDay **nao** estao no namelist: sao fixos
  ForecastDay=0_i4, TimeOfDay=0.0_r4
  (ICnAGCM/InputParameters.f90:326-327). Ou seja, todo arquivo de IC produzido
  por este gerador tem 0 e 0.0 no registro 1.
- Year/Month/Day/Hour sao extraidos de DateICn, mas **nada** os usa para escolher
  o dado de entrada: o nome do dado de entrada e
  DataInp='ERAINTATM'//DateICn//'.dat' etc. (ICnAGCM/InputParameters.f90:161-204).
- Bug de relatorio: a linha "Current Date" imprime DateInp
  (ICnAGCM/SpectralICn.f90:358-359).

### 5.2 O que o BAM consome

O modelo usa ForecastDay/TimeOfDay como **ponto de partida do relogio da
previsao** e DateInp/DateCur como data do estado. Para o gerador de CI, os
valores corretos sao sempre 0 / 0.0 / data da analise em ambos. O port deve
escrever os 4 inteiros na ordem hora-mes-dia-ano e **nao** inventar valores
diferentes para DateCur.

### 5.3 Efeito nos diagnosticos

tdef 1 linear <hh>z<dd><mmm><yyyy> 1dy — a data aparece no eixo T dos .ctl;
um valor errado aqui nao invalida nada para o BAM, mas os .GrADS so serao
comparaveis ao legado se a data for a mesma.

## 6. Checklist do que o port Python DEVE reproduzir

Envelope e tipos

1. Abrir em 'wb' e escrever **exatamente 3+4K registros**, cada um com marcador
   int32 BE antes e depois; nenhum byte a mais (nem metadado, nem padding).
2. Todos os reais em **float32 big endian**; ForecastDay e as 8 datas em
   **int32 big endian**; TimeOfDay como float32.
3. Tamanho final tem de bater com a secao 1.4 (1839784 bytes para
   M=62/K=28/_hibrid; 1839780 para _sigma).

Registro 1

4. Ordem exata: ForecastDay, TimeOfDay, DateInp[4]=(hour,month,day,year),
   DateCur[4], vertical.
5. Hibrido: **2*(K+1)** valores = a em **Pa** (nao hPa) seguido de b
   adimensional. Sigma: **K+1** interfaces (1 -> 0) seguido de **K** camadas
   calculadas pela formula de Phillips com kappa=Rd/Cp do legado.
6. Na duvida entre hPa e Pa: o legado le a em Pa, divide por 100, e multiplica por
   100 de volta na escrita (round-trip) — o port nao deve reintroduzir um /100.

Nomes e sufixos

7. Nome = DirICnOut + DataOut + '_hibrid' ou '_sigma', com
   DataOut = 'GANL' + Case + yyyymmddhh + 'S.unf.' + 'T' + Grid + M(4) + 'L' + K(3)
   (ICnAGCM/InputParameters.f90:287-294). Exemplo (Case='EIT', Grid='Q', M=62,
   K=28, 1999-01-06 00Z): GANLEIT1999010600S.unf.TQ0062L028_hibrid.
   O sufixo **nao** e configuravel: o port deve mante-lo (e confirmar contra o BAM
   qual e o nome efetivamente aberto — pendencia P1).
8. Nao truncar, nao renomear, nao comprimir; diretorio de destino
   DirICnOut = DirMainOut + 'model/datain/' (ICnAGCM/InputParameters.f90:213).

Fisica/unidades

9. LnPs = ln(ps/10) com ps em **hPa** (cbar). Nunca ln(Pa) nem ln(hPa).
10. Tvir = T*(1+cTv*q), cTv = Rv/Rd - 1 com Rd=287.05, Rv=461.50
    (ICnAGCM/InputParameters.f90:335-336, 341), aplicado **antes** da
    interpolacao vertical e da transformada (ICnAGCM/SpectralICn.f90:160, 281,
    285): o arquivo leva temperatura **virtual**.
11. Divergencia = **-convergencia** (ICnAGCM/Decomposition.f90:128-131).
12. Umidade em **kg/kg**; topografia em **m**; ventos nunca vao ao arquivo (so
    div/vort).
13. Nao suavizar a topografia: SmoothCoef existe mas **nunca e chamada**
    (ICnAGCM/SpectralICn.f90:412-461), embora Iter=10 e SmthPerCut=0.12 sejam
    lidos do namelist (ICnAGCM/InputParameters.f90:128-129). O port deve manter o
    comportamento (sem suavizacao) para produzir o mesmo arquivo.

Espectral

14. Um array de Mnwv2 = (M+1)(M+2) reais por registro, na numeracao la0
    (secao 2.2); c1 = sqrt(2)*media; seno de m=0 = 0.
15. Calculos em float64 e conversao para float32 **apenas na escrita** (o legado
    faz exatamente isso: arrays r8, REAL(...,r4) no WRITE).
16. Grade de saida gaussiana: ImaxOut = menor inteiro 5-suave >= 3*M+1,
    JmaxOut = ImaxOut/2 arredondado para par (ICnAGCM/InputParameters.f90:380-446;
    M=62 -> 192x96). As latitudes/pesos gaussianos vem de glats
    (ICnAGCM/Legendre.f90:216-263) e definem tanto a transformada quanto o ydef
    do .GrADS.

Diagnosticos

17. .GrADS/.GrADS.ctl (e .GrADSOut): float32 BE cru, sem marcadores,
    linhas N->S, 2+4*Kmax registros (secao 3). Opcional, mas e o oraculo visual.
18. Reproduzir o relatorio de ICnAGCM/SpectralICn.f90:356-398 (media de topo,
    pressao media, por nivel: Tvir, div, vort, q) — e o que o operador compara
    hoje para aceitar a IC. Proposta de comando no pacote: bamci report <arq>
    (equivalente ao CheckICn, mas sem o limite de 33 caracteres).

### 6.1 Tolerancias propostas para o port (o legado nao tem nenhuma)

| grandeza | criterio | justificativa |
|---|---|---|
| envelope (tamanho, marcadores, ordem dos blocos) | **exato** | nao ha margem: o leitor e sequencial |
| ForecastDay, datas | **exato** | inteiros |
| TimeOfDay | exato (0.0) | fixo no legado |
| sigma/b | abs <= 1e-6 | float32, valores O(1) |
| a (Pa) | abs <= 1e-2 Pa | ida-e-volta Pa->mbar->Pa |
| coeficientes (topo, lnps, T, div, vort, q) | rel <= 1e-6 **ou** abs <= 1e-10, o maior | eps float32 ~1.2e-7; caminho numerico diferente (FFT/quadratura) |
| c1*Cte (medias do relatorio) | rel <= 1e-6 | mesma base acima |
| pressao media global (hPa) | abs <= 1e-3 hPa | 10*EXP(c1*Cte) |
| Tvir (K) | abs <= 1e-3 K | idem |
| q (kg/kg) | abs <= 1e-9 | idem |

Se o port replicar a mesma cadeia numerica (mesma FFT, mesma quadratura
gaussiana, mesma interpolacao vertical), a comparacao **bit-a-bit** dos registros
depois do byte-swap e viavel e deve ser o teste de ouro; as tolerancias acima sao
o fallback quando a cadeia for reimplementada com outras bibliotecas.

## 7. Riscos e pendencias

- **P1 (bloqueante para "o BAM aceita")**: o nome exato que o modelo abre
  (..._hibrid vs ..._sigma, e se ele espera o sufixo) nao pode ser verificado
  neste espelho — o fonte do BAM nao esta aqui, so o pre-processamento. Precisa
  ser confirmado no espelho/jaci de Bam-1.2-subsaz (somente leitura).
- **P2**: o PDF esta desatualizado em dois pontos: registro 1 do caso hibrido
  (escreve 2*(K+1) coeficientes, nao sigma) e unidade de ln(ps) (cbar, nao kPa).
  Seguir o codigo, nao o PDF.
- **P3**: CheckICn nao valida nada e nem consegue abrir o arquivo real
  (LEN=33); nao trate-lo como criterio de aceitacao. O criterio deve ser o
  relatorio + comparacao de coeficientes.
- **P4**: ordem de linhas dos .GrADS. No legado o binario vai de norte para sul e
  o .ctl declara ydef de sul para norte com options yrev. Em bam-ci,
  bamci/grads.py:142-144 **inverte** as linhas quando yrev esta presente
  enquanto bamci/ctl.py:158-162 **inverte** o eixo de latitudes — as duas
  inversoes juntas podem duplicar a correcao. Para escrever .GrADS, gravar as
  linhas N->S cruas e nao aplicar inversao; resolver a inconsistencia antes de
  usar round-trip como teste.
- **P5**: no .GrADS (entrada) o eixo z tem **KmaxInp** niveis de pressao de
  entrada (37 para EIT), nao KmaxOut; no .GrADSOut tem KmaxOut. Nao confundir ao
  comparar com o estado final.
- **P6**: IDSLInp=2 e constante (media das interfaces) — se o BAM usar outra
  regra para as pressoes de camada, a comparacao de Tvir/q por nivel tera um vies
  sistematico que **nao** e bug do envelope.
- **P7**: DataTop/DataTopGrid (ICnAGCM/InputParameters.f90:218-219,
  ICnAGCM/SpectralICn.f90:58-92): se DirTop/Topography.TQ<mend> existir, os
  coeficientes de topografia do arquivo de IC sao **lidos desse arquivo**, e so a
  grade nova (Topography.G<jmax>) e re-transformada para o .GrADS/NewPs. O port
  precisa decidir se reproduz esse caminho ou sempre deriva a topografia do dado
  de entrada; os dois produzem arquivos diferentes.
- **P8**: codigo morto no legado que nao deve ser portado como funcionalidade:
  SmoothCoef (nunca chamada), NTracers>1 (variaveis no .ctl sem dados),
  ozonio/tracers no .GrADS (comentados), Undef=-999999 (nunca usado,
  ICnAGCM/InputParameters.f90:350).
