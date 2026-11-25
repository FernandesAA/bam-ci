# Artefatos do ensaio MB09 — 2026-11-19

Inventario da rodada de 2026-11-19 do ensaio MB09 (T126L42) e do bloqueio
T213L64. Coleta **somente leitura** na JACI em 2026-11-25; caminhos absolutos
a partir de `BASE` (o mesmo diretorio em que os comandos foram executados).

```
BASE = /p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2
```

Tamanhos em bytes (`ls -l`; diretorios com `du -sb`). `md5` de 32 hex;
`nao recuperado` quando o arquivo nao estava acessivel na coleta. Onde o
`md5` esta no `work/IC15_126/md5.txt`, a tabela repete o valor.

## 1. Entradas

| artefato | caminho absoluto | tamanho (B) | md5 |
|---|---|---:|---|
| ERA5 combinado | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/DataAnl/era5_2026111900.grib` | 394442116 | `a8a050df22f7aa6b62d8d324e56a56bf` |
| ERA5 niveis de pressao | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/DataAnl/era5_pl_2026111900.grib` | 384168780 | `ddfa169231eb6dd07c6c1d753508a311` |
| ERA5 superficie | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/DataAnl/era5_sl_2026111900.grib` | 10273336 | `8539ba2238acc6e5e1e37b10485c33db` |
| SST semanal T213 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/SSTWeekly20261119.G00320` | 1638400 | `51603178e9e18d2c1f99d74671a1eb1e` |

## 2. Condicoes iniciais de controle e aleatorias (T126L42)

| artefato | caminho absoluto | tamanho (B) | md5 |
|---|---|---:|---|
| IC controle (GANLNMC) | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/bam126/GANLNMC2026111900S.unf.TQ0126L042` | 11055828 | `657b4627077ea4d086b6c6bf266fe1bb` |
| aleatoria 01 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/rand126/mb09rand_01.unf.TQ0126L042` | 11055828 | `bb9f686452028cc546c4191827a7e6bb` |
| aleatoria 02 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/rand126/mb09rand_02.unf.TQ0126L042` | 11055828 | `0162ebf210a2fcbb07be6eb768c250e9` |
| aleatoria 03 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/rand126/mb09rand_03.unf.TQ0126L042` | 11055828 | `f3f666c9175fae27ddd4312390a50174` |
| aleatoria 04 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/rand126/mb09rand_04.unf.TQ0126L042` | 11055828 | `3977adc88ba83b718393da60580e42dd` |
| aleatoria 05 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/rand126/mb09rand_05.unf.TQ0126L042` | 11055828 | `334941a8b2f3ceb173ea13e2e835f1a3` |
| aleatoria 06 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/rand126/mb09rand_06.unf.TQ0126L042` | 11055828 | `4ccb6c9f5f73f6e533f3fa0e4662f50f` |
| aleatoria 07 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/rand126/mb09rand_07.unf.TQ0126L042` | 11055828 | `a3c43d60a4f85f15d2f1173bd57c6bb3` |
| manifesto das aleatorias | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/rand126/random.json` | 994 | `nao recuperado` |

## 3. Rodadas de 36 h (4 ranks, OMP_NUM_THREADS=2)

Binario usado: `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/model/exec/ParModel_MPI_O2_fix3`
(`md5 a5e2e569f553a7b99579f69c219d9bb5`, flags de compilacao `-O2 -fp-model=precise`).
Saidas em `<run>/out` (`.dir` + `.fct`; 3 h por segmento).

| rodada | caminho absoluto (`.../out`) | tamanho do diretorio (B) | job PBS | arquivos | md5 |
|---|---|---:|---|---:|---|
| controle | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/ens126/ctrl/out` | 1389481949 | `848830` | 37 | `nao recuperado` (diretorio) |
| membro 01 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/ens126/m01/out` | 1389481214 | `848831` | 37 | `nao recuperado` (diretorio) |
| membro 02 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/ens126/m02/out` | 1389481928 | `848832` | 37 | `nao recuperado` (diretorio) |
| membro 03 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/ens126/m03/out` | 1312053284 | `848833` | 35 | `nao recuperado` (diretorio) |
| membro 04 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/ens126/m04/out` | 1389480072 | `848834` | 37 | `nao recuperado` (diretorio) |
| membro 05 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/ens126/m05/out` | 1312053456 | `848835` | 35 | `nao recuperado` (diretorio) |
| membro 06 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/ens126/m06/out` | 1389480144 | `848836` | 37 | `nao recuperado` (diretorio) |
| membro 07 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/ens126/m07/out` | 1079762445 | `848854` | 29 | `nao recuperado` (diretorio) |

As rodadas com menos arquivos sao as que pararam antes de 36 h; o relatorio
`RELATORIO_TESTE.md` registra `qdel`/paradas. A rodada de controle passou de
36 h e foi encerrada por `qdel`.

## 4. EOF e as 15 condicoes iniciais

| artefato | caminho absoluto | tamanho (B) | md5 |
|---|---|---:|---|
| EOFs (M=11, 6..36 h, ETSA, `q_profile=artigo_l28`) | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/mb09/eof_TQ0126L042.npz` | 3488953292 | `74a04212fc613c06cdc8969c3334cbed` |
| controle (copia) | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_000.unf` | 11055828 | `657b4627077ea4d086b6c6bf266fe1bb` |
| +p1 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_001.unf` | 11055828 | `9fac6055f3583fa50bebe82ff38e4cce` |
| -p1 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_002.unf` | 11055828 | `46cfb1516f67fc7aa8befb67873f5770` |
| +p2 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_003.unf` | 11055828 | `3a8ccfed2238a9c0385a033fc8badcab` |
| -p2 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_004.unf` | 11055828 | `631d411f24f764fb307826854afadc01` |
| +p3 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_005.unf` | 11055828 | `d7ba0e0e1efdee1cf63b8c007ac1b928` |
| -p3 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_006.unf` | 11055828 | `299225eaf527febdad91adc53f548805` |
| +p4 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_007.unf` | 11055828 | `0122d07fd5b9036afd49593805a9a5df` |
| -p4 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_008.unf` | 11055828 | `07c329baa1a46ba01ff69aa8e31dbff9` |
| +p5 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_009.unf` | 11055828 | `74170c08cb0ffb8fba6761739a8eb9cb` |
| -p5 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_010.unf` | 11055828 | `cb13aa9a654046ed68ac190df8c074fa` |
| +p6 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_011.unf` | 11055828 | `40fda52796cce18e38b2a12dbd5092fb` |
| -p6 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_012.unf` | 11055828 | `cce30f4cd6db3d5f36361ddea936bf1a` |
| +p7 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_013.unf` | 11055828 | `fca9988c965d112393246a3bd8d3b5f3` |
| -p7 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/mb09_014.unf` | 11055828 | `abd5e53ea29f6711cf00730c13a4316d` |
| lista de md5 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/md5.txt` | 915 | `nao recuperado` |
| manifesto do apply | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/members.json` | 2901 | `nao recuperado` |
| log do apply | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/work/IC15_126/apply.log` | 405 | `nao recuperado` |

O apply usou `--q-positivity 3`; o `q_min` final de cada membro fica em
`members.json` e e ~`1e-8` kg/kg. O controle `mb09_000.unf` e copia exata do
IC de controle.

## 5. Relatorios e binarios

| artefato | caminho absoluto | tamanho (B) | md5 |
|---|---|---:|---|
| relatorio do teste | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/RELATORIO_TESTE.md` | 9585 | `8c8d2a734f48977b8abe7e92750014d2` |
| log do build O2 | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/build/BUILD_O2.log` | 49660 | `e41f76f147750a0cd461b6b9d61c1ef0` |
| copia/procedencia | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/COPIA.md` | 991 | `eafa8aa753217344317bf4d48ddf6a06` |
| origem da suite | `/p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2/pre/ORIGEM.md` | 665 | `6413610f7676c234ef5b5c6b352a8a53` |

Trilha dos binarios O2 (mesmo diretorio `model/exec`):

| binario | md5 | tamanho (B) | observacao |
|---|---|---:|---|
| `ParModel_MPI_O2` | `a5a49feafc3437a4429fc4e3dd38e12b` | 34215720 | sem patches |
| `ParModel_MPI_O2_fix` | `50d8009c584dc07142bea1110ad27580` | 34213960 | patches 1-4 |
| `ParModel_MPI_O2_fix2` | `22092f864c65e10fb23afa667779d9d2` | 34218632 | + patch 5 (guard) |
| `ParModel_MPI_O2_fix3` | `a5e2e569f553a7b99579f69c219d9bb5` | 34218568 | + patch 6; **usado no ensaio** |
| `ParModel_MPI_O2_fix4` | `8df71c9290aaa2709af1c4fd87b20342` | 34218880 | diagnostico NaN |
| `ParModel_MPI_O2_fix5` | `253515c8e2a5eaaabb31909597bd8d2e` | 34220376 | diagnostico NaN |
| `ParModel_MPI_O2_fix6` | `3afebd1a5a32df1006dc98c76758e445` | 34225344 | comparacoes NaN-safe |
| `ParModel_MPI_O2_fix7` | `d25433243e84e138d20ad7a0b6a810eb` | 34234576 | + patch 7 (`nCols`) |
| `ParModel_MPI_O2_i3d` | `ffc5282eae585c4c9921c706a61aa4f1` | 34216040 | variante 3D |

Flags de compilacao do ensaio: `-O2 -fp-model=precise` (alem de
`-convert big_endian -g -traceback -align array64byte -mcmodel=large -qopenmp`).

## 6. Bloqueio T213L64 (nao resolvido)

Sintoma: `NaN/Outside` em `tonpsadiabat`; o diagnostico mostra `thte=NaN` com
`prs/PSFC` validos e `tmk=NaN` em **k=2 e k=4, i=1..3 de blocos nCols=6**;
`qvp ~1e-22` e `prs` valido.

7 patches aplicados (arquivos `.orig*` e diffs em `build/BUILD_O2.log`):

| # | patch | local |
|---|---|---|
| 1 | `klev = 1` (nunca inicializado no ramo 2D) | `Utils/PhysicalFunctions.f90:3348` |
| 2 | `IF (klev > kMax) klev = kMax` | `Utils/PhysicalFunctions.f90:3383` |
| 3 | `klcl = 1` (evita `zrel(klcl)` OOB) | `Utils/PhysicalFunctions.f90:3474` |
| 4 | `kLevel = 150` em vez de `STOP` | `Utils/PhysicalFunctions.f90:3551` |
| 5 | guard `IF (.NOT.(ethmax >= 0.0_r8))` + fallback `ethpari` | `Utils/PhysicalFunctions.f90:3365` (fallback em `:3464`) |
| 6 | guard `IF (.NOT.(totprs > 0.0_r8))` (NaN-safe) | `Utils/PhysicalFunctions.f90:3422` |
| 7 | `nCols = MIN(ibMaxPerJB(latco), iMax)` em `Convection.f90` | `Physics/Convection/Convection.f90:824` e `:846` |

Descartado como causa/solucao: rebuild do modelo, NMI, SST com `-999`,
registros dos estaticos e PostGrib. O bloqueio persiste.

Proximo passo tecnico: `PRINT iMax x ibMaxPerJB(latco)` para confirmar a
largura real do bloco e dump cru de `t3/q3/prsi` **antes da conveccao** no
primeiro `k/i` afetado.

## 7. Como reproduzir

```bash
cd /p/projetos/ioper/models/BAM-1.2-subsaz/bam-1.2

# EOFs direto dos .fct do BAM (sem PostGrib)
python -m bamci.mb09 build --source bamfct \
  --control work/ens126/ctrl/out \
  --members work/ens126/m01/out work/ens126/m02/out work/ens126/m03/out \
            work/ens126/m04/out work/ens126/m05/out work/ens126/m06/out \
            work/ens126/m07/out \
  --window 6:36:3 --analysis-time 2026111900 --n-modes 7 \
  --source-delta-sigma DeltaSigma/DeltaSigma.L028 \
  --delta-sigma DeltaSigma/DeltaSigma.L042 --preset etsa \
  --out work/mb09/eof_TQ0126L042.npz

# 15 ICs finais com positividade de q
python -m bamci.mb09 apply --eof work/mb09/eof_TQ0126L042.npz \
  --ic work/bam126/GANLNMC2026111900S.unf.TQ0126L042 \
  --out-dir work/IC15_126 --prefix mb09 --resolution TQ0126L042 \
  --n-modes 7 --q-positivity 3
```

Para conferir os md5: `md5sum work/IC15_126/*.unf` e `cat work/IC15_126/md5.txt`.
