# Perturbacao de condicoes iniciais — referencias

Documento de referencia bibliografica sobre metodos de perturbacao de condicoes
iniciais (ICs) para modelos numericos de previsao de tempo e clima. Escrito para
apoiar a construcao de ensaios perturbados do BAM 1.2 subsazonal, em particular a
linhagem **Mendonca & Bonatti (2009)** ja implementada em
`/p/projetos/ioper/models/oensMB09_bam_xc50` (jaci).

Todos os registros bibliograficos deste documento foram conferidos em Crossref
(titulo, revista, volume, paginas, ano, autores). Onde o texto integral nao estava
acessivel, o abstract foi obtido de registro oficial indexado e a fonte e indicada.

> **Aviso de escopo.** Ate o fechamento desta versao nao foi possivel confirmar se
> existe publicacao especifica do CPTEC/INPE sobre perturbacao de ICs em regime
> **subsazonal** (30–60 dias). A literatura abaixo e de previsao de tempo de medio
> prazo (1–15 dias) e de ensemble em geral. Ver secao 6.

---

## 1. O metodo EOF de perturbacao

### 1.1 Referencia central (a que o `oensMB09_bam_xc50` segue)

**Mendonca, A. M. B.; Bonatti, J. P. (2009).** Experiments with EOF-Based
Perturbation Methods and Their Impact on the CPTEC/INPE Ensemble Prediction System.
*Monthly Weather Review*, **137**(4), 1438–1459.
DOI: [10.1175/2008MWR2581.1](https://doi.org/10.1175/2008MWR2581.1)

Abstract (registro Crossref, verbatim):

> The impact of modifications of the perturbation method based on empirical
> orthogonal functions (EOFs) used operationally upon the ensemble prediction system
> (EPS) at the Center for Weather Prediction and Climate Studies/National Institute
> for Space Research (CPTEC/INPE) is evaluated. The main changes proposed in this
> study are to apply the EOF method to perturb the midlatitudes, apply additional
> perturbations to the surface pressure (P) and specific humidity (Q) fields, and
> compute regional perturbations over South America. The impact of these
> modifications in the characteristics of the initial perturbations and in the
> quality of the EPS forecasts is investigated. The EPS forecasts are evaluated
> through average statistical scores over the period 15 December 2004–15 February
> 2005. The statistical scores used in the evaluation are pattern anomaly
> correlation, root-mean-square error, ensemble spread, Brier skill score, and
> perturbation versus error correlation analysis (PECA). Results indicate that with
> the inclusion of perturbations on P and Q, EOF-based perturbations acquire a more
> baroclinic structure. It is also observed that the simultaneous application of
> additional perturbations both in the extratropics and to the P and Q fields
> improves the performance of CPTEC EPS and enhances the quality of forecast
> perturbations. Moreover, regional EOF-based perturbations computed over South
> America have positive impact on the ensemble forecasts over the target region.

Resumo complementar (registro INPE/OAI-PMH do ePrint
`sid.inpe.br/mtc-m17@80/2007/10.11.13.28`), que descreve melhor o **algoritmo**:

> [...] selection of most unstable perturbations based on the linear growth rate of
> EOF coefficients; to use a simplified version of the breeding of growing modes
> method (Toth and Kalnay, 1993) for generation perturbed initial conditions. The
> impacts of these modifications are evaluated under the hypothesis of perfect model
> for both an operational suite during January 2005 and a case study for the first
> documented hurricane over South Atlantic, occurred in March 2004.

### 1.2 Metodo-base de Zhang & Krishnamurti

**Zhang, Z.; Krishnamurti, T. N. (1999).** A Perturbation Method for Hurricane
Ensemble Predictions. *Monthly Weather Review*, **127**(4), 447–469.
DOI: [10.1175/1520-0493(1999)127%3C0447:APMFHE%3E2.0.CO;2](https://doi.org/10.1175/1520-0493(1999)127%3C0447:APMFHE%3E2.0.CO;2)

Origem do "EOF method" do CPTEC: constroi os modos a partir de um conjunto de
diferencas de previsao, nao de uma serie temporal longa de analises.

Trabalhos anteriores dos mesmos autores, que formam a base do ensemble de furacao:

- **Zhang, Z.; Krishnamurti, T. N. (1997).** Ensemble Forecasting of Hurricane
  Tracks. *BAMS*, **78**(12), 2785–2795.
  DOI: [10.1175/1520-0477(1997)078%3C2785:EFOHT%3E2.0.CO;2](https://doi.org/10.1175/1520-0477(1997)078%3C2785:EFOHT%3E2.0.CO;2)
- **Krishnamurti, T. N. et al. (1997).** Physical Initialization and Hurricane
  Ensemble Forecasts. *Weather and Forecasting*, **12**(3), 503–514.
  DOI: [10.1175/1520-0434(1997)012%3C0503:PIAHEF%3E2.0.CO;2](https://doi.org/10.1175/1520-0434(1997)012%3C0503:PIAHEF%3E2.0.CO;2)

### 1.3 EOF aplicado sobre bred vectors (melhora substancial)

**Magnusson, L.; Kallen, E.; Nycander, J. (2008).** Initial state perturbations in
ensemble forecasting. *Nonlinear Processes in Geophysics*, **15**(5), 751–759.
DOI: [10.5194/npg-15-751-2008](https://doi.org/10.5194/npg-15-751-2008) — acesso aberto.

Abstract (texto integral, Copernicus):

> Due to the chaotic nature of atmospheric dynamics, numerical weather prediction
> systems are sensitive to errors in the initial conditions. To estimate the
> forecast uncertainty, forecast centres produce ensemble forecasts based on
> perturbed initial conditions. How to optimally perturb the initial conditions
> remains an open question and different methods are in use. One is the singular
> vector (SV) method, adapted by ECMWF, and another is the breeding vector (BV)
> method (previously used by NCEP). In this study we compare the two methods with a
> modified version of breeding vectors in a low-order dynamical system (Lorenz-63).
> We calculate the Empirical Orthogonal Functions (EOF) of the subspace spanned by
> the breeding vectors to obtain an orthogonal set of initial perturbations for the
> model. We will also use Normal Mode perturbations. Evaluating the results, we
> focus on the fastest growth of a perturbation. The results show a **large
> improvement for the BV-EOF perturbations compared to the non-orthogonalised BV**.

### 1.4 Evolucao espaco-temporal: bred vs Lyapunov vs singular vectors

**Pazo, D.; Rodriguez, M. A.; Lopez, J. M. (2010).** Spatio-temporal evolution of
perturbations in ensembles initialized by bred, Lyapunov and singular vectors.
*Tellus A*, **62**(1), 10–23.
DOI: [10.1111/j.1600-0870.2009.00419.x](https://doi.org/10.1111/j.1600-0870.2009.00419.x)

Compara os tres tipos de perturbacao nao apenas no instante inicial, mas ao longo
da integracao — util para justificar o horizonte de crescimento em previsao
subsazonal.

---

## 2. Metodos classicos (a linha de base contra a qual a EOF e medida)

| Metodo | Referencia | Ideia | Custo |
|---|---|---|---|
| **Breeding (BV)** | Toth & Kalnay (1993); Toth & Kalnay (1997) | integra a perturbacao e rescala periodicamente; amostra os modos de crescimento rapido sem adjunto | barato |
| **Singular vectors (SV)** | Molteni & Palmer (1993); Buizza & Palmer (1995) | maior amplificacao linear em norma de energia | exige adjunto |
| **Aleatorio / Monte Carlo** | Houtekamer & Derome (1995); Houtekamer et al. (1996) | ruido com estrutura de covariancia | barato, pouco eficiente |
| **Normal modes** | ver Magnusson et al. (2008) | modos normais lineares | barato; crescimento fraco |
| **Lag-average / multi-model** | ECMWF; Krishnamurti et al. (1999) | diversidade por defasagem temporal ou por modelo | trivial |

### 2.1 Breeding of growing modes

- **Toth, Z.; Kalnay, E. (1993).** Ensemble Forecasting at NMC: The Generation of
  Perturbations. *BAMS*, **74**(12), 2317–2330.
  DOI: [10.1175/1520-0477(1993)074%3C2317:EFANTG%3E2.0.CO;2](https://doi.org/10.1175/1520-0477(1993)074%3C2317:EFANTG%3E2.0.CO;2)
- **Toth, Z.; Kalnay, E. (1997).** Ensemble Forecasting at NCEP and the Breeding
  Method. *Monthly Weather Review*, **125**(12), 3297–3319.
  DOI: [10.1175/1520-0493(1997)125%3C3297:EFANAT%3E2.0.CO;2](https://doi.org/10.1175/1520-0493(1997)125%3C3297:EFANAT%3E2.0.CO;2)
- **Tracton, M. S.; Kalnay, E. (1993).** Operational Ensemble Prediction at the
  National Meteorological Center: Practical Aspects. *Weather and Forecasting*,
  **8**(3), 379–398.
  DOI: [10.1175/1520-0434(1993)008%3C0379:OEPATN%3E2.0.CO;2](https://doi.org/10.1175/1520-0434(1993)008%3C0379:OEPATN%3E2.0.CO;2)

### 2.2 Singular vectors

- **Molteni, F.; Palmer, T. N. (1993).** Predictability and finite-time instability
  of the northern winter circulation. *QJRMS*, **119**(510), 269–298.
  DOI: [10.1002/qj.49711951004](https://doi.org/10.1002/qj.49711951004)
- **Buizza, R.; Palmer, T. N. (1995).** The Singular-Vector Structure of the
  Atmospheric Global Circulation. *JAS*, **52**(9), 1434–1456.
  DOI: [10.1175/1520-0469(1995)052%3C1434:TSVSOT%3E2.0.CO;2](https://doi.org/10.1175/1520-0469(1995)052%3C1434:TSVSOT%3E2.0.CO;2)
- **Hartmann, D. L.; Buizza, R.; Palmer, T. N. (1995).** Singular Vectors: The
  Effect of Spatial Scale on Linear Growth of Disturbances. *JAS*, **52**(22),
  3885–3894.
  DOI: [10.1175/1520-0469(1995)052%3C3885:SVTEOS%3E2.0.CO;2](https://doi.org/10.1175/1520-0469(1995)052%3C3885:SVTEOS%3E2.0.CO;2)
- Estrategia operacional do ECMWF (documento de referencia, nao artigo):
  [The ECMWF initial ensemble perturbation strategy](https://www.ecmwf.int/en/about/media-centre/media-resources/tim-palmer-festschrift/The-ECMWF-initial-ensemble-perturbation-strategy)

### 2.3 Perturbacao aleatoria com estrutura de covariancia

- **Houtekamer, P. L.; Derome, J. (1995).** Methods for Ensemble Prediction.
  *Monthly Weather Review*, **123**(7), 2181–2196.
  DOI: [10.1175/1520-0493(1995)123%3C2181:MFEP%3E2.0.CO;2](https://doi.org/10.1175/1520-0493(1995)123%3C2181:MFEP%3E2.0.CO;2)
- **Houtekamer, P. L.; Lefaivre, L.; Derome, J.; Ritchie, H.; Mitchell, H. L.
  (1996).** A System Simulation Approach to Ensemble Prediction. *Monthly Weather
  Review*, **124**(6), 1225–1242.
  DOI: [10.1175/1520-0493(1996)124%3C1225:ASSATE%3E2.0.CO;2](https://doi.org/10.1175/1520-0493(1996)124%3C1225:ASSATE%3E2.0.CO;2)

### 2.4 Superensemble multivariado

- **Krishnamurti, T. N. et al. (1999).** Improved Weather and Seasonal Climate
  Forecasts from Multimodel Superensemble. *Science*, **285**(5433), 1548–1550.
  DOI: [10.1126/science.285.5433.1548](https://doi.org/10.1126/science.285.5433.1548)

---

## 3. Fontes INPE/CPTEC em portugues

Material de divulgacao e ePrints, uteis para citacao em relatorio tecnico em
portugues (nao passaram por revisao por pares no mesmo nivel dos artigos acima):

- **Mendonca, A. M. B.; Bonatti, J. P. (2007).** *Experiments with EOF-Based
  perturbation and breeding of growing modes methods and their impacts on
  CPTEC/INPE ensemble prediction system.* ePrint INPE,
  `sid.inpe.br/mtc-m17@80/2007/10.11.13.28`. Versao em PDF no acervo
  ([link](http://mtc-m16b.sid.inpe.br/col/sid.inpe.br/mtc-m17@80/2007/10.11.13.28/doc/v3.pdf)).
- **Mendonca, A. M. B. (2006).** *Impacto dos modos de perturbacao extratropicais na
  previsao de tempo.* CBMET ([PDF](http://mtc-m16b.sid.inpe.br/col/sid.inpe.br/mtc-m15@80/2006/11.07.13.52/doc/Mendon%e7a.Impacto.pdf)).
- **Coutinho, M. et al. (2004).** *Implementacao da tecnica de "Ensemble Prediction"
  no CPTEC — resultados preliminares* ([PDF](http://mtc-m16b.sid.inpe.br/col/cptec.inpe.br/walmeida/2004/07.05.15.14/doc/Coutinho_Implementacao%20da%20tecnica%20de%20.pdf)).
- **Mendonca, A. M. B.** *Ensemble global weather forecasting at CPTEC.* WMO
  CBS/ERA ([PDF](https://old.wmo.int/extranet/pages/prog/www/DPS/TC-DPFS-2002/Papers-Posters/Topic1-Mendonca.pdf)).
- **Bastarz, C. et al. (2016).** *Sistema de previsao por conjunto global no CPTEC:
  estado atual e planos futuros.*
  ([PDF](http://mtc-m21b.sid.inpe.br/col/sid.inpe.br/mtc-m21b/2016/11.29.18.08.53/doc/bastarz_sistema_poster.pdf)).
- **Experiments with EOF-Based Perturbation Method to Ensemble Weather Forecasting
  in Middle Latitudes** (2005), pre-print
  ([PDF](http://mtc-m16b.sid.inpe.br/col/cptec.inpe.br/adm_conf/2005/10.31.12.25/doc/1829-1832.pdf)).

---

## 4. Leitura critica aplicada ao BAM subsazonal

Tres pontos que a literatura deixa claros e que afetam diretamente o desenho do
ensaio perturbado:

### 4.1 O metodo EOF e estacionario e linear

Zhang & Krishnamurti constroem os modos a partir de **um** conjunto de diferencas de
previsao (no MB09: 11 diferencas de 48 h de uma unica data). Em escala subsazonal
(30–60 dias) a instabilidade dominante nao e a baroclinica de 48 h: os modos que
crescem mais em 48 h podem nao ser os que crescem na janela de 30 dias. O resultado
de Magnusson et al. (2008) — **ortogonalizar o substrato melhora muito** — sugere
que aplicar a mesma ortogonalizacao (Gram-Schmidt nos modos, ou EOF sobre o conjunto
de BV) e uma melhoria barata antes de partir para adjunto/singular vectors.

### 4.2 Amplitude fisica por variavel e necessaria

O MB09 escala cada variavel com uma STD fisica tabelada (T, vento, ps, q) e aplica
o padrao com sinal oposto nos pares (+/-). Mendonca & Bonatti (2009) mostram que
**acrescentar P e Q** faz a perturbacao adquirir "estrutura mais baroclinica" — ou
seja, perturbar apenas T e vento deixa a estrutura vertical incompleta. Para o BAM
em niveis sigma, vale verificar se a perturbacao de q e coerente com a de T
(balanco hidrostatico / umidade relativa), e nao apenas aplicada em paralelo.

### 4.3 Corte rigido de regiao e agressivo

O MB09 aplica 3 faixas e 2 caixas com corte **abrupto**. Isso produz descontinuidades
nos gradientes que o modelo pode amplificar artificialmente. A pratica usual e
compor padroes regionais **somando** contribuicoes suavemente ponderadas, em vez de
mascarar um padrao global com uma funcao degrau.

---

## 5. Mapa metodo -> parametro de implementacao

Referencia cruzada entre a literatura e as chaves/flags do `seasst` e do
`oensMB09_bam_xc50`:

| Concept da literatura | Onde aparece no codigo | Observacao |
|---|---|---|
| EOF de diferencas de previsao (Zhang & Krishnamurti) | base EOFT do `oensMB09_bam_xc50` | 11 diferencas de 48 h de uma data |
| Aplicar apenas a EOF 1 | idem | limita a diversidade; ver 4.1 |
| STD fisica tabelada por variavel | `--amplitude-model physical` | alternativa: `variance` |
| Pares (+/-) | `--pairing symmetric` | alternativa: `independent` |
| Selecao pelo crescimento linear dos coeficientes | nao implementado como filtro | candidato a melhoria |
| Ortogonalizacao tipo BV-EOF (Magnusson 2008) | **nao implementado** | candidato a melhoria |
| Regioes com corte (3 faixas + 2 caixas) | fixo no MB09 | candidato a suavizacao |
| Numero de modos | `--apply-modes N` | hoje N = 1 |
| Vento como EOF complexa u + iv | idem | consistente com SV do ECMWF |

---

## 6. Lacunas desta pesquisa

- **Nao foi confirmada** publicacao do CPTEC/INPE sobre perturbacao de ICs para
  previsao **subsazonal** especificamente. A literatura CPC/CPTEC sobre subsazonal
  trata mais de SST, TSM e ensembles multi-modelo do que do metodo de perturbacao.
  Se existir (por exemplo tese do PG-MET do INPE), precisa ser buscada no repositorio
  de teses do INPE, que nao esta indexado pelas ferramentas usadas aqui.
- **Textos integrais nao lidos**: Mendonca & Bonatti (2009) e Zhang & Krishnamurti
  (1999) estao sob paywall ou bloqueados por CDN; os abstracts foram obtidos de
  registros oficiais (Crossref e OAI-PMH do INPE) e sao consistentes entre si. As
  afirmacoes da secao 4 sao inferencias a partir desses abstracts e do codigo
  `oensMB09_bam_xc50`, **nao** citacoes literais do corpo dos artigos.
- **Nao cobertos aqui**: perturbacoes estocasticas de esquema fisico (SPPT, SKEB) e
  perturbacao de tendencia/superficie, que sao ortogonais ao tema "condicao inicial"
  e merecem documento proprio.

---

## 7. BibTeX

```bibtex
@article{MendoncaBonatti2009,
  author  = {Mendon{\c c}a, Antonio Marcos B. and Bonatti, Jos{\'e} Paulo},
  title   = {Experiments with {EOF}-Based Perturbation Methods and Their Impact on
             the {CPTEC}/{INPE} Ensemble Prediction System},
  journal = {Monthly Weather Review},
  volume  = {137}, number = {4}, pages = {1438--1459}, year = {2009},
  doi     = {10.1175/2008MWR2581.1}
}

@article{ZhangKrishnamurti1999,
  author  = {Zhang, Z. and Krishnamurti, T. N.},
  title   = {A Perturbation Method for Hurricane Ensemble Predictions},
  journal = {Monthly Weather Review},
  volume  = {127}, number = {4}, pages = {447--469}, year = {1999},
  doi     = {10.1175/1520-0493(1999)127<0447:APMFHE>2.0.CO;2}
}

@article{Magnusson2008,
  author  = {Magnusson, L. and K{\"a}ll{\'e}n, E. and Nycander, J.},
  title   = {Initial state perturbations in ensemble forecasting},
  journal = {Nonlinear Processes in Geophysics},
  volume  = {15}, number = {5}, pages = {751--759}, year = {2008},
  doi     = {10.5194/npg-15-751-2008}
}

@article{TothKalnay1993,
  author  = {Toth, Z. and Kalnay, E.},
  title   = {Ensemble Forecasting at {NMC}: The Generation of Perturbations},
  journal = {Bulletin of the American Meteorological Society},
  volume  = {74}, number = {12}, pages = {2317--2330}, year = {1993},
  doi     = {10.1175/1520-0477(1993)074<2317:EFANTG>2.0.CO;2}
}

@article{TothKalnay1997,
  author  = {Toth, Z. and Kalnay, E.},
  title   = {Ensemble Forecasting at {NCEP} and the Breeding Method},
  journal = {Monthly Weather Review},
  volume  = {125}, number = {12}, pages = {3297--3319}, year = {1997},
  doi     = {10.1175/1520-0493(1997)125<3297:EFANAT>2.0.CO;2}
}

@article{MolteniPalmer1993,
  author  = {Molteni, F. and Palmer, T. N.},
  title   = {Predictability and finite-time instability of the northern winter
             circulation},
  journal = {Quarterly Journal of the Royal Meteorological Society},
  volume  = {119}, number = {510}, pages = {269--298}, year = {1993},
  doi     = {10.1002/qj.49711951004}
}

@article{BuizzaPalmer1995,
  author  = {Buizza, R. and Palmer, T. N.},
  title   = {The Singular-Vector Structure of the Atmospheric Global Circulation},
  journal = {Journal of the Atmospheric Sciences},
  volume  = {52}, number = {9}, pages = {1434--1456}, year = {1995},
  doi     = {10.1175/1520-0469(1995)052<1434:TSVSOT>2.0.CO;2}
}

@article{Houtekamer1996,
  author  = {Houtekamer, P. L. and Lefaivre, L. and Derome, J. and Ritchie, H.
             and Mitchell, H. L.},
  title   = {A System Simulation Approach to Ensemble Prediction},
  journal = {Monthly Weather Review},
  volume  = {124}, number = {6}, pages = {1225--1242}, year = {1996},
  doi     = {10.1175/1520-0493(1996)124<1225:ASSATE>2.0.CO;2}
}

@article{Pazo2010,
  author  = {Paz{\'o}, D. and Rodr{\'i}guez, M. A. and L{\'o}pez, J. M.},
  title   = {Spatio-temporal evolution of perturbations in ensembles initialized by
             bred, {Lyapunov} and singular vectors},
  journal = {Tellus A},
  volume  = {62}, number = {1}, pages = {10--23}, year = {2010},
  doi     = {10.1111/j.1600-0870.2009.00419.x}
}
```
