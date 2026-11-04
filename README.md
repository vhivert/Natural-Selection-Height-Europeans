# Converging evidence of positive selection at height-associated loci in Europe

Scripts for the research project "Converging evidence of positive selection at height-associated loci in Europe" (Hivert et al., bioRxiv 2026, [doi:10.64898/2026.01.27.702172](https://doi.org/10.64898/2026.01.27.702172)).

The Qx script extends the Qx test of polygenic adaptation of Berg & Coop (2014) to adjust for population-level covariates such as ancient admixture. It is an edited version of [Jeremy Berg's code](https://github.com/jjberg2/PolygenicAdaptationCode).

## Contents

| File | Description |
| --- | --- |
| `functions_Qx_covar.R` | Extended from Jeremy Berg's code. Functions to compute Qx, covariate-adjusted Qx and per-population contributions, with asymptotic and empirical p-values |
| `AM_simulations.R` | New code for this study. Simulations of the effect of assortative mating on Qx |

## Requirements

- R (the analyses in the paper used R 4.4.3)
- R packages `plyr`, `data.table` and `MASS`
- Perl, to run `sampleSNPs.pl`. This script is not distributed here. Download it from the [`Scripts` folder of Jeremy Berg's repository](https://github.com/jjberg2/PolygenicAdaptationCode/tree/master/Scripts), make it executable (`chmod +x sampleSNPs.pl`), and give its folder as `scripts_path`.

## Quick start

Example data for all input files except the covariate file are available in the [`Example` folder of Jeremy Berg's repository](https://github.com/jjberg2/PolygenicAdaptationCode/tree/master/Example). See the note on comment lines under `match.pop.file` below.

```r
source("functions_Qx_covar.R")

result <- PolygenicAdaptationFunction(
  gwas.data.file = "height.gwas.txt",
  freqs.file = "height.freqs.txt",
  env.var.data.files = list("populations.env.txt"),
  match.pop.file = "match.file.txt",
  full.dataset.file = "genome.file.txt",
  covar.file = "ancestry.PCs.txt",     # NULL for the standard Qx test
  path = "Height",
  scripts_path = "Scripts",
  match.categories = c("MAF", "ldscore"),
  match.bins = list(seq(0, 0.5, 0.02),
                    c(0.4363, 13.4999, 20.9151, 32.0346, 305.1991002)),
  cov.SNPs.per.cycle = 20000,
  cov.cycles = 5,
  null.phenos.per.cycle = 1000,
  null.cycles = 10,
  load.cov.mat = FALSE,
  sim.null = TRUE,
  check.allele.orientation = TRUE,
  PCProj = TRUE
)
```

A similar example is at the end of `functions_Qx_covar.R`. All file names should be given either with a full path or with a path relative to the directory from which the function is called.

## Method

Let $Z$ be the vector of genetic values of $M$ populations, where $Z_i = \sum_l \beta_l p_{li}$ for $L$ GWAS SNPs with effect $\beta_l$ and frequency $p_{li}$. Let $\epsilon_l$ be the mean frequency of SNP $l$ across populations, and $V = \sum_l \beta_l^2 \epsilon_l (1 - \epsilon_l)$.

**Covariance matrix.** Genome-wide SNPs matched to the GWAS SNPs are used to estimate the $M \times M$ matrix

$$U = \frac{1}{n-1} \sum_l s_l s_l^\top, \qquad s_{li} = \frac{p_{li} - \epsilon_l}{\sqrt{\epsilon_l (1 - \epsilon_l)}}.$$

The rows of $U$ sum to zero, so the code uses a $(M-1) \times M$ matrix $T$ that mean-centers the populations and drops the last one in sorted order. The result does not depend on which population is dropped.

**Qx.** With $F = T U T^\top$,

$$Q_X = \frac{(TZ)^\top F^{-1} (TZ)}{V},$$

which follows a $\chi^2_{M-1}$ distribution under neutrality.

**Covariate adjustment.** Let $W$ be the matrix of $J$ mean-centered covariates, without the dropped population. Genetic values are regressed on the covariates by generalized least squares,

$$\hat\gamma = (W^\top F^{-1} W)^{-1} W^\top F^{-1} TZ, \qquad r = TZ - W\hat\gamma,$$

and $Q_{X,res} = r^\top F^{-1} r / V$ follows a $\chi^2_{M-J-1}$ distribution under neutrality.

**Population contributions.** With $Z_c = 2(Z - \bar{Z})$ and $F^+$ the Moore-Penrose inverse of the $M \times M$ centered covariance matrix, the contribution of population $i$ is

$$\frac{Z_{c,i} (F^+ Z_c)_i}{4V}.$$

Contributions sum to $Q_X$ and can be negative. The same formula is applied to the residuals and to the fitted values of the full $M$-population regression.

**Null distribution.** For each GWAS SNP, matched SNPs are drawn from the same joint bin of the matching variables and given the effect of that GWAS SNP, signed to apply to the minor allele. All statistics are then computed again on each null replicate.

## Input files

Files have a header and are whitespace-delimited. Frequencies in all files must refer to the same allele. Because `sampleSNPs.pl` always keeps lines whose first field contains an uppercase S, N or P, as a way of keeping the header, SNP identifiers must not contain these letters (rsIDs are fine).

**gwas.data.file** contains the effect size information for the SNPs identified in the GWAS. Required columns are `SNP`, `EFF` and `FRQ`, where `EFF` and `FRQ` refer to the same allele. Columns `A1` and `A2` are required if `check.allele.orientation = TRUE`.

**freqs.file** contains the frequency information from the test populations for the GWAS SNPs, with one line per SNP and population. Required columns are `SNP`, `CLST` (population) and `FRQ`, plus `A1` and `A2` if `check.allele.orientation = TRUE`. Any match category other than `MAF` must also be given in this file, with the same value on every line of a SNP. In Jeremy Berg's code these were read from `match.pop.file`. The file must contain exactly the same SNPs as `gwas.data.file`, and every SNP in every population.

**env.var.data.files** is a list of files with columns `CLST`, `ENV` and `REG`, giving the environmental variable and the integer grouping for the regional Z scores. At least one file is required, even if you are only interested in the Qx test. In this case you can fill in arbitrary integers for `ENV` and `REG` and ignore those results. These statistics were not used in our study.

**match.pop.file** contains the matching variables for all SNPs that may be used for the covariance matrix or for the null genetic values, from a single population. Required columns are `SNP`, `FRQ` and any match category other than `MAF`. Ideally, this file comes from the GWAS population. The file is read with `data.table::fread`, so remove the comment lines found at the top of Jeremy Berg's example file.

**full.dataset.file** contains the frequencies of all SNPs in all populations. Required columns are `SNP`, `CLST`, `A1`, `A2` and `FRQ`. Each SNP in `match.pop.file` must be present in every test population. Other populations may be present and are ignored.

**covar.file** (optional) contains a column `CLST` followed by $J$ numeric columns, with one line per test population and exactly the same populations as the other files. **Each covariate must be mean-centered across the test populations before loading**, as the function does not center them. Otherwise `Qx.res` is wrong and no longer equals the sum of `Qx.res.contrib`. The paper used six principal components of ancient admixture proportions. Use `NULL` or `"NULL"` for the standard Qx test.

```
CLST    PC1      PC2
POP1   0.021   -0.004
POP2  -0.013    0.017
POP3  -0.008   -0.013
```

## Arguments

**path** is the output folder, relative to the directory from which you run the analysis.

**scripts_path** is the folder containing `sampleSNPs.pl`. Supply your own path, as the default `Your/Script/Path` is a placeholder.

**match.categories** gives the column names of the variables to match on when sampling SNPs. `MAF` is computed from `FRQ`. Any number of categories can be used, but each must be present in `freqs.file` and `match.pop.file` unless it is `MAF`.

**match.bins** is a list of breaks for `cut()`, with one element per match category. For example, `seq(0, 0.5, 0.02)` splits the minor allele frequency into bins of 0.02, and `c(2)` splits a binary variable into two bins. Breaks must span all values, and each bin holding a GWAS SNP must hold at least one SNP in `match.pop.file`.

**cov.SNPs.per.cycle** is the minimum number of SNPs sampled in each batch to compute the covariance matrix. The number actually sampled is rounded up to a multiple of the number of GWAS SNPs.

**cov.cycles** is the number of batches to average the covariance matrix over.

**null.phenos.per.cycle** is the number of null genetic values computed per batch.

**null.cycles** is the number of batches of null genetic values. The paper used 10,000 replicates, that is, `null.phenos.per.cycle = 1000` and `null.cycles = 10`.

**load.cov.mat**: if `TRUE`, and `path/uncenteredcovmat.Rdata` exists, the matrix is loaded rather than computed. Note that Jeremy Berg's README calls this flag `load.cov.matrix` and looks in `path/Output`. A matrix computed with other software can be supplied this way, as an object named `uncentered.cov.mat` with populations in sorted order.

**PCProj**: if `TRUE`, partitions Qx across the principal components of the covariance matrix. This is applied to the unadjusted genetic values.

**sim.null**: if `FALSE`, only asymptotic p-values are computed, from the $\chi^2$ distribution for Qx and the standard normal distribution for the Z scores.

**check.allele.orientation** (default `TRUE`): compares `A1` and `A2` between `gwas.data.file` and `freqs.file`. Alleles polarized in opposite directions are flipped. The function stops if a strand issue is detected, but it does not detect one if this flag is `FALSE`.

## Output

`PolygenicAdaptationFunction` returns a list with the elements below, which are also saved in `path/Output`.

**stats** (`theStats.Robj`)

| Element | Description |
| --- | --- |
| `Qx` | Qx statistic |
| `Fst.comp`, `LD.comp` | Parts of Qx due to differences at single SNPs and to covariance between SNPs. They sum to `Qx` |
| `Qx.contrib` | Contribution of each population to `Qx`, named by population |
| `betas`, `pearson.rs`, `spearman.rhos` | Association with each environmental variable, one element per file |
| `reg.Z`, `ind.Z` | Regional and per-population Z scores |

With covariates, the following are also given:

| Element | Description |
| --- | --- |
| `Qx.res` | Covariate-adjusted Qx |
| `Fst.comp.res`, `LD.comp.res` | Components of `Qx.res` |
| `Qx.res.contrib` | Contribution of each population to `Qx.res` |
| `Qx.cov.contrib` | Contribution of each population to the part of Qx explained by the covariates. It sums to `Qx - Qx.res` |

**asymptotic.p.vals** (`asymptotic.pVals.Robj`) contains `Qx`, `reg.Z` and `ind.Z`, plus `Qx.res` with covariates.

**p.vals** (`pVals.Robj`) is returned if `sim.null = TRUE`, and contains the empirical p-values, including `Qx.res`, `Fst.comp.res` and `LD.comp.res` with covariates. A p-value of 0 means that no null replicate exceeded the observed value, so it should be reported as less than 1 / (number of replicates).

Other files are:

| File | Content |
| --- | --- |
| `path/Output/nullStats.Robj` | Null distributions |
| `path/Output/genetic.values.Robj` | Genetic values of the populations |
| `path/uncenteredcovmat.Rdata` | Matrix $U$ |
| `path/GWAS.file.RData` | GWAS data after allele alignment |
| `path/cov.SNPs*`, `path/cov.samples*`, `path/null.SNPs*`, `path/null.samples*` | Sampled SNPs and their frequencies. These can be deleted after a run |
| `path/figs/PCAQx.pdf`, `path/Output/PCVarProj.Rdata`, `path/Output/PCs.Rdata` | Principal component partition of Qx, if `PCProj = TRUE` |

## Differences from Jeremy Berg's code

- Added `covar.file`, `Qx.res` and its components, and the per-population contributions.
- $U$ is computed from frequencies mean-centered per SNP, rather than with `cov()` on standardized frequencies.
- `cov.SNPs.per.cycle` now sets the number of SNPs sampled. It was hard-coded to 5000.
- The sign of null effects is set with `FRQ <= 0.5` rather than `FRQ >= 0.5`, so that `MA.EFF` is applied to the minor allele.
- Match categories other than `MAF` are read from `freqs.file` rather than `match.pop.file`.
- Added `scripts_path`, and `fread` to read large files.
- Columns of the sampled data are selected by name rather than by position.

## Assortative mating simulations

`AM_simulations.R` tests whether assortative mating alone can inflate Qx. Each replicate simulates populations from a common ancestral frequency under assortative mating on a trait, with causal SNPs of positive effect and the same number of neutral SNPs. Qx is then computed using the neutral SNPs for the covariance matrix.

```
Rscript AM_simulations.R <path> <num.pops> <nSNP> <h2_eq> <spousal.corr> <n.reps> <num.cores>
```

| Argument | Description |
| --- | --- |
| `path` | Output folder |
| `num.pops` | Number of populations |
| `nSNP` | Number of causal SNPs, and of neutral SNPs |
| `h2_eq` | Heritability at equilibrium |
| `spousal.corr` | Phenotypic correlation between mates |
| `n.reps` | Number of replicates |
| `num.cores` | Number of cores |

For example, for 13 populations, 14,000 SNPs (half causal), $h^2 = 0.8$, $r = 0.2$ and 1,000 replicates on 16 cores:

```
Rscript AM_simulations.R ./results 13 7000 0.8 0.2 1000 16
```

The output is a tab-delimited file, `Simu_AM_<num.pops>Pops_10Gens_<nSNP>SNPs_<h2_eq>h2eq_r<spousal.corr>.txt`, with columns `Fst_component`, `LD_component` and `Qx`, and one line per replicate. The population size (`N = 1000`) and the number of generations (`ng = 10`) are set in the call to `simAM` near the end of the script. If you change `ng`, note that `10Gens` in the output file name is fixed.

## Citation

Please cite according to the code you use.

**Qx analysis (`functions_Qx_covar.R`).** This code extends the method and code of Berg & Coop, so please cite both papers:

- Hivert V, et al. Converging evidence of positive selection at height-associated loci in Europe. bioRxiv (2026). [doi:10.64898/2026.01.27.702172](https://doi.org/10.64898/2026.01.27.702172)
- Berg JJ, Coop G. A population genetic signal of polygenic adaptation. PLoS Genetics 10, e1004412 (2014). [doi:10.1371/journal.pgen.1004412](https://doi.org/10.1371/journal.pgen.1004412)

**Assortative mating simulations (`AM_simulations.R`).** This code is new to our study, so please cite:

- Hivert V, et al. Converging evidence of positive selection at height-associated loci in Europe. bioRxiv (2026). [doi:10.64898/2026.01.27.702172](https://doi.org/10.64898/2026.01.27.702172)

## Acknowledgments

The annotation of the scripts and this documentation were written with the help of Claude AI (Anthropic).

## License

This repository is distributed under the GNU General Public License v3.0. The code that `functions_Qx_covar.R` extends is copyright (C) 2014 by Jeremy J Berg; see [his repository](https://github.com/jjberg2/PolygenicAdaptationCode) for its license.
