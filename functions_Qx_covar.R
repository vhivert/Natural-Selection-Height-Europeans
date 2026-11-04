###############################################################################
#
# functions_Qx_covar.R
#
# Compute the Qx test of polygenic adaptation (Berg & Coop 2014, PLoS Genet,
# doi:10.1371/journal.pgen.1004412), with an optional adjustment for
# population-level covariates such as ancient admixture.
#
# Edited from Jeremy Berg's code:
#   https://github.com/jjberg2/PolygenicAdaptationCode
#
# Valentin Hivert
# 04 November 2026
#
# Main additions to the original code
#   - covar.file: GLS adjustment of genetic values for J population-level
#     covariates, giving the residual statistic Qx.res
#   - per-population contributions to Qx (Qx.contrib, Qx.res.contrib,
#     Qx.cov.contrib)
#   - covariance matrix computed from SNP-wise mean-centered frequencies
#   - scripts_path argument, and data.table::fread for the large files
#
# Dependencies
#   R packages : plyr, data.table, MASS
#   Perl script: sampleSNPs.pl, from Scripts/ in Jeremy Berg's repository.
#                It must be executable and is located with scripts_path.
#
# Notation used in the comments
#   M          number of populations
#   L          number of GWAS SNPs
#   J          number of covariates
#   p[l,i]     frequency of the effect allele of SNP l in population i
#   beta[l]    GWAS effect of SNP l
#   eps[l]     mean of p[l,] across populations
#   Z[i]       genetic value of population i, sum over l of beta[l] * p[l,i]
#   V          sum over l of beta[l]^2 * eps[l] * (1 - eps[l])
#
# Input file formats and arguments are described in README.md, which extends
# the documentation of Jeremy Berg's repository.
#
# The annotation of this script was written with the help of Claude AI
# (Anthropic).
#
###############################################################################



#' Run the Qx analysis
#'
#' Top-level function. Reads the input files, estimates the population
#' covariance matrix, computes Qx and its related statistics, and optionally
#' builds their null distributions from matched SNPs.
#'
#' @param gwas.data.file Path to the GWAS file. Required columns: SNP, EFF, FRQ;
#'   plus A1 and A2 if check.allele.orientation = T.
#' @param freqs.file Path to the frequencies of the GWAS SNPs in the test
#'   populations, one row per SNP-population pair. Required columns: SNP, CLST,
#'   FRQ; plus A1 and A2 if check.allele.orientation = T. Must also contain any
#'   match category other than MAF.
#' @param env.var.data.files List of paths to environmental variable files, with
#'   columns CLST, ENV and REG. At least one is required, even if the
#'   environmental statistics are not of interest. Must be a list when more than
#'   one file is supplied.
#' @param match.pop.file Path to the file of candidate SNPs for matching, from
#'   one population. Required columns: SNP, FRQ, plus any match category other
#'   than MAF.
#' @param full.dataset.file Path to the genome-wide frequencies in the test
#'   populations. Required columns: SNP, CLST, A1, A2, FRQ. Each SNP in
#'   match.pop.file must be present in every test population.
#' @param covar.file Optional path to the covariate file: a CLST column followed
#'   by J numeric columns of covariates. Each covariate must be mean-centered across
#'   the test populations before loading; this function does not center them.
#'   Use NULL or "NULL" for no adjustment.
#' @param path Output directory, created if needed.
#' @param scripts_path Directory that contains sampleSNPs.pl. Replace the default
#'   with your own path.
#' @param match.categories Column names to match on when sampling SNPs. MAF is
#'   computed from FRQ; other columns must be supplied as described above.
#' @param match.bins List of breaks passed to cut(), one element per match
#'   category. Breaks must span all the values in the data.
#' @param cov.SNPs.per.cycle Minimum number of matched SNPs sampled per cycle to
#'   estimate the covariance matrix.
#' @param cov.cycles Number of independent batches to average the covariance
#'   matrix over.
#' @param null.phenos.per.cycle Number of null genetic values per cycle.
#' @param null.cycles Number of null cycles. The total number of null replicates
#'   is null.phenos.per.cycle * null.cycles.
#' @param load.cov.mat If T, load <path>/uncenteredcovmat.Rdata when it exists
#'   rather than estimating the matrix again.
#' @param PCProj If T, project the genetic values on the principal components of
#'   the covariance matrix (see PCAProjections).
#' @param sim.null If T, compute empirical p-values from null simulations. If F,
#'   return asymptotic p-values only.
#' @param check.allele.orientation If T (default), align GWAS alleles to
#'   freqs.file with MatchAlleles().
#' @param log.file Unused; kept for the commented-out logging in MatchAlleles().
#'
#' @return A list with stats (output of CalcStats), asymptotic.p.vals and, when
#'   sim.null = T, p.vals. Results are also saved under <path>/Output/.
PolygenicAdaptationFunction <- function ( gwas.data.file , freqs.file , env.var.data.files , match.pop.file , full.dataset.file , covar.file=NULL, path , scripts_path="Your/Script/Path",
                                          match.categories , match.bins  , cov.SNPs.per.cycle = 5000 , cov.cycles = 4 , null.phenos.per.cycle = 1000 , null.cycles = 2 ,
                                          load.cov.mat = F , PCProj = F , sim.null = T , check.allele.orientation = T,log.file="log.txt") {
  library(plyr)
  library(data.table)
  library(MASS)
  ##### Read in GWAS file, frequencies file, match pop file, and environmental variables files; also sort environmental variable files alphanumerically, create output directory
  dir.create ( path )
  gwas.data <- read.table ( gwas.data.file , header = T , stringsAsFactors = F )
  freqs <- read.table ( freqs.file , header = T , stringsAsFactors = F )
  env.var.data <- lapply ( env.var.data.files , read.table , header = T , stringsAsFactors = F )
  env.var.data <- lapply ( env.var.data , function ( ENV.VAR ) ENV.VAR [ order ( ENV.VAR$CLST ) , ] )
  # Covariates are sorted by population, as are the environmental variables, so that their rows follow the same order as the frequency matrices
  covar.data=NULL
  if(!is.null(covar.file) && length(covar.file) > 0 && covar.file!="NULL"){
    print("include covariate file in analysis")
    covar.data <- read.table(covar.file,header=T,stringsAsFactors = F)
    covar.data <- covar.data[ order(covar.data$CLST), ]
    # Covariates are not centered here: they must be mean-centered across the populations before loading (see the note in CalcStats())
    #covar.data[,-c(1)] <- scale( covar.data[,-c(1)],center = T,scale = F ) #mean center the covariates
  }
  #match.pop <- read.table ( match.pop.file , header = T )
  # fread() replaces read.table() for speed on the large match file. Remove any comment lines from this file
  match.pop <- fread( match.pop.file , header = T, data.table=FALSE)
  dir.create ( paste ( path , "/Output" , sep = "" ) )
  
  # All populations and all SNPs must agree across the input files
  # sanity checks
  if ( any ( unlist ( lapply ( env.var.data , function ( y ) lapply ( env.var.data , function ( x ) !x$CLST %in% y$CLST ) ) ) ) ) stop ( "Environmental variable datasets differ in the number of populations." )
  if ( any ( !freqs$CLST %in% env.var.data [[ 1 ]]$CLST ) ) stop ( "Frequency dataset has data for populations not present in environmental variable datasets.")
  if ( any ( !env.var.data [[ 1 ]]$CLST %in% freqs$CLST ) ) stop ( "Environmental variable dataset has data for populations not present frequency dataset.")
  if ( any ( !gwas.data$SNP %in% freqs$SNP ) ) stop ( "GWAS dataset contains SNPs that are not in the frequency dataset.")
  if ( any ( !freqs$SNP %in% gwas.data$SNP ) ) stop ( "Freq dataset contains SNPs that are not in the gwas dataset." )
  
  # The covariate file must contain exactly the same populations
  if(!is.null(covar.data)){
    if ( any ( !covar.data$CLST %in% env.var.data [[ 1 ]]$CLST ) ) stop ( "Covariates dataset has data for populations not present in environmental variable datasets.")
    if ( any ( !env.var.data [[ 1 ]]$CLST %in% covar.data$CLST ) ) stop ( "Environmental dataset has data for populations not present in covariates datasets.")
  }
  
  if ( check.allele.orientation == T ) gwas.data <- MatchAlleles ( gwas.data , freqs  )
  save ( gwas.data , file = paste ( path , "/GWAS.file.RData" , sep = "" ) )
  
  
  # get minor allele frequency from effect allele frequency, and flip effects for minor alleles
  gwas.data$MAF <- ifelse ( gwas.data$FRQ > 0.5 , 1 - gwas.data$FRQ , gwas.data$FRQ )	
  gwas.data$MA.EFF <- ifelse ( gwas.data$FRQ > 0.5 , - gwas.data$EFF , gwas.data$EFF )	
  match.pop$MAF <- ifelse ( match.pop$FRQ > 0.5 , 1 - match.pop$FRQ , match.pop$FRQ )	
  

  # Match variables other than MAF are read from the first population in freqs.file, so they must be present in that file (Berg's code read them from match.pop)
  # add relevant matching columns that are not in gwas data 
  tmp=freqs[freqs$CLST==unique(freqs$CLST)[1],]
  gwas.data <- cbind ( gwas.data , tmp [ match ( gwas.data$SNP , tmp$SNP ) , match.categories [ !match.categories %in% names ( gwas.data ) ] ,drop = FALSE] )
  rm(tmp)
  #gwas.data <- cbind ( gwas.data , match.pop [ match ( gwas.data$SNP , match.pop$SNP ) , match.categories [ !match.categories %in% names ( gwas.data ) ] ,drop = FALSE] )
  
  # Each match category is cut into bins, separately for the GWAS SNPs (x) and the candidate SNPs (y). Together these define a multi-way contingency table
  # assign SNPs to appropriate bins in contingency table
  x <- mapply ( function ( CATS , BINS ) 
    cut ( gwas.data [ , CATS ] , breaks = BINS , include.lowest = T ) ,  
    CATS = match.categories , BINS = match.bins , SIMPLIFY = FALSE )
  names ( x ) <- bin.names <-  paste ( match.categories , ".BINS" ,sep = "" )
  gwas.data <- cbind ( gwas.data , x )
  
  y <- mapply ( function ( CATS , BINS ) 
    cut ( match.pop [ , CATS ] , breaks = BINS , include.lowest = T ) ,  
    CATS = match.categories , BINS = match.bins , SIMPLIFY = FALSE )
  names ( y ) <- paste ( match.categories , ".BINS" , sep = "" )
  match.pop <- cbind ( match.pop , y )
  rm ( x , y )
  
  # keep effects in the same order as the frequency matrix
  gwas.data <- gwas.data [ order ( gwas.data$SNP ) , ]

  # initialize some handy things
  pop.names <- env.var.data [[ 1 ]]$CLST
  num.pops <- length ( pop.names )
  
  # Sort by SNP then population. CalcStats() relies on this order when it reshapes FRQ into an M x L matrix
  freqs <- freqs [ with ( freqs , order ( SNP , CLST ) ) , ]
  
  # Estimate, or load, the M x M covariance matrix of standardized frequencies
  ## STEP 1 : Estimate variance-Covariance matrix of allele frequencies
  #recover()
  if ( file.exists ( paste ( path , "/uncenteredcovmat.Rdata" , sep = "" ) ) & load.cov.mat == T ) {
    load ( paste ( path , "/uncenteredcovmat.Rdata" , sep = "" ) )
  } else {
    uncentered.cov.mat <- SampleCovSNPs ( gwas.data , match.pop , pop.names , bin.names , SNPs.per.cycle = cov.SNPs.per.cycle , cycles = cov.cycles , path = path, scripts_path = scripts_path , full.dataset.file = full.dataset.file , env.var.data = env.var.data)
  }
  # T.mat: (M-1) x M mean-centering matrix. It centers on the mean of all populations and removes the last one (last in sorted order), which makes the centered covariance matrix invertible
  T.mat <- matrix ( rep ( c ( ( num.pops - 1 ) / num.pops , rep ( - 1 / ( num.pops ) , times = num.pops - 1 ) ) , times = num.pops - 1 ) , ncol = num.pops , nrow = num.pops - 1 )
  T.mat.full <- diag(num.pops) - matrix(1/num.pops, nrow=num.pops, ncol=num.pops)  # centering matrix
  # Off-center covariance and transformation matrices, used only for the regional and individual Z scores
  tmp1 <- lapply ( env.var.data , function ( x ) GetOffCovMats ( x , uncentered.cov.mat , gwas.data$EFF , freqs ) )
  regional.off.center.cov.mats <- lapply ( tmp1 , function ( x ) x [[ 1 ]] )
  regional.off.center.T.mats <- lapply ( tmp1 , function ( x ) x [[ 2 ]] )
  tmp2 <- GetOffCovMats ( data.frame ( env.var.data [[ 1 ]] [ , c ( 1 , 2 ) ] , REG = seq_len ( nrow ( env.var.data [[ 1 ]] ) ) ) , uncentered.cov.mat , gwas.data$EFF , freqs )
  individual.off.center.cov.mats <-  tmp2 [[ 1 ]]
  individual.off.center.T.mats <- tmp2 [[ 2 ]]
  # F.mat: centered covariance matrix of the M-1 retained populations. C.inv is the inverse of its Cholesky factor, used to rotate genetic values into independent standard units
  F.mat <- T.mat %*% uncentered.cov.mat %*% t ( T.mat )
  F.inv <- solve ( F.mat )
  C.inv <- solve ( t ( chol ( F.mat ) ) )
  #load(file = "CodeForRelease/HeightTest/uncentered.cov.mat.RData")
  # Additive genetic variance scaling term V. EFF is reordered here to follow the sorted frequencies
  epsilon.freqs <- tapply ( freqs$FRQ , freqs$SNP , mean )
  add.gen.var <- sum ( gwas.data$EFF [ order ( gwas.data$SNP ) ]^2  * epsilon.freqs * ( 1 - epsilon.freqs ) )
  
  
  ## STEP X : Compute the different statistics of interest
  the.stats <- CalcStats ( 
    freqs = matrix ( freqs$FRQ , nrow = num.pops ) ,
    effects = gwas.data$EFF ,
    env.var.data = env.var.data , 
    var = add.gen.var ,
    covar=covar.data,
    uncentered.cov.mat = uncentered.cov.mat ,
    T.mat = T.mat ,
    T.mat.full = T.mat.full,
    regional.off.center.cov.mats = regional.off.center.cov.mats , 
    regional.off.center.T.mats = regional.off.center.T.mats , 
    individual.off.center.cov.mats = individual.off.center.cov.mats , 
    individual.off.center.T.mats = individual.off.center.T.mats , 
    path = path ,
    PCProj = PCProj ,
    null = F
  )
  rownames ( the.stats$reg.Z ) <- sprintf("%s %d", "Region", seq_len ( nrow ( the.stats$reg.Z ) ) )
  colnames ( the.stats$reg.Z ) <- sprintf("%s %d", "Env File", seq_len ( ncol ( the.stats$reg.Z ) ) )
  names ( the.stats$Qx.contrib ) <- env.var.data[[1]]$CLST
  if(!is.null(covar.data)){
    names ( the.stats$Qx.res.contrib ) <- env.var.data[[1]]$CLST
    names ( the.stats$Qx.cov.contrib ) <- env.var.data[[1]]$CLST
  }
  names ( the.stats$ind.Z ) <- env.var.data[[1]]$CLST
  
  # Asymptotic p-values: chi-squared with M-1 degrees of freedom for Qx, M-J-1 for the covariate-adjusted statistic, and standard normal for the Z scores
  asymp.p.vals <- list()
  asymp.p.vals$Qx <- 1 - pchisq ( the.stats$Qx , nrow ( T.mat ))
  if(!is.null(covar.data)){
    asymp.p.vals$Qx.res <- 1 - pchisq ( the.stats$Qx.res , nrow ( T.mat ) - ( ncol ( covar.data ) - 1 ) ) # df = M - J - 1, as covar.data has J + 1 columns (CLST and J covariates)
  }
  
  asymp.p.vals$reg.Z <- 2 * (1-pnorm ( abs ( the.stats$reg.Z ) ) )
  asymp.p.vals$ind.Z <- 2 * (1-pnorm ( abs ( the.stats$ind.Z ) ) )
  rownames ( asymp.p.vals$reg.Z ) <- sprintf("%s %d", "Region", seq_len ( nrow ( asymp.p.vals$reg.Z ) ) )
  colnames ( asymp.p.vals$reg.Z ) <- sprintf("%s %d", "Env File", seq_len ( ncol ( asymp.p.vals$reg.Z ) ) )
  names ( asymp.p.vals$ind.Z ) <- env.var.data[[1]]$CLST
  save ( the.stats , file = paste ( path , "/Output/theStats.Robj" , sep = "" ) )
  save ( asymp.p.vals , file = paste ( path , "/Output/asymptotic.pVals.Robj" , sep = "" ) )
  if ( sim.null == F ){
    return ( list ( stats = the.stats , asymptotic.p.vals = asymp.p.vals ) )
  } else {
    # Empirical p-values: proportion of null replicates at least as extreme as the observed statistic. The p-value is exactly 0 if no replicate exceeds it, and should then be reported as < 1 / (number of replicates)
    p.vals <- list ()
    print("Start null simulations\n")
    null.stats <- NullStats ( gwas.data , match.pop , env.var.data , pop.names , uncentered.cov.mat , F.inv , C.inv , T.mat , T.mat.full , covar.data, regional.off.center.cov.mats , regional.off.center.T.mats , individual.off.center.cov.mats , individual.off.center.T.mats , bin.names , null.phenos.per.cycle , null.cycles , path , scripts_path = scripts_path,full.dataset.file )
    print("done")
    #recover()
    p.vals$Qx <- sum ( the.stats$Qx < null.stats$Qx ) / length ( null.stats$Qx )
    p.vals$Fst.comp <- sum ( the.stats$Fst.comp < null.stats$Fst.comp ) / length ( null.stats$Fst.comp )
    p.vals$LD.comp <- sum ( the.stats$LD.comp < null.stats$LD.comp ) / length ( null.stats$LD.comp )
    if(!is.null(covar.data)){
      p.vals$Qx.res <- sum ( the.stats$Qx.res < null.stats$Qx.res ) / length ( null.stats$Qx.res )
      p.vals$Fst.comp.res <- sum ( the.stats$Fst.comp.res < null.stats$Fst.comp.res ) / length ( null.stats$Fst.comp.res )
      p.vals$LD.comp.res <- sum ( the.stats$LD.comp.res < null.stats$LD.comp.res ) / length ( null.stats$LD.comp.res )
    }
    
    tmp.upper.tail <- rowSums ( null.stats$betas < unlist ( the.stats$betas ) ) / ncol ( null.stats$betas )
    p.vals$betas <- ifelse ( tmp.upper.tail > 0.5 , 2 * ( 1 - tmp.upper.tail ) , 2 * tmp.upper.tail )
    p.vals$pearson.rs <- rowSums ( null.stats$pearson.rs^2 > unlist ( the.stats$pearson.rs )^2 ) / ncol ( null.stats$betas )
    tmp.upper.tail <- rowSums ( null.stats$spearman.rhos < unlist ( the.stats$spearman.rhos ) ) / ncol ( null.stats$spearman.rhos )
    p.vals$spearman.rhos <- ifelse ( tmp.upper.tail > 0.5 , 2 * ( 1 - tmp.upper.tail ) , 2 * tmp.upper.tail )
    tmp.upper.tail <- mapply ( function ( x , y ) rowSums ( x > y ) / ncol ( x ) , x = null.stats$reg.Z , y = split ( t ( the.stats$reg.Z ) , 1:ncol ( the.stats$reg.Z ) ) )
    p.vals$reg.Z <- ifelse ( tmp.upper.tail > 0.5 , 2 * ( 1 - tmp.upper.tail ) , 2 * tmp.upper.tail )
    rownames ( p.vals$reg.Z ) <- sprintf("%s %d", "Region", seq_len ( nrow ( p.vals$reg.Z ) ) )
    colnames ( p.vals$reg.Z ) <- sprintf("%s %d", "Env File", seq_len ( ncol ( p.vals$reg.Z ) ) )
    tmp.upper.tail <- rowSums ( null.stats$ind.Z > the.stats$ind.Z ) / ncol ( null.stats$ind.Z )
    p.vals$ind.Z <- ifelse ( tmp.upper.tail > 0.5 , 2 * ( 1 - tmp.upper.tail ) , 2 * tmp.upper.tail )
    save ( p.vals , file = paste ( path , "/Output/pVals.Robj" , sep = "" ) )
    save ( null.stats , file = paste ( path , "/Output/nullStats.Robj" , sep = "" ) )
    return ( list ( stats = the.stats , p.vals = p.vals , asymptotic.p.vals = asymp.p.vals ) )
  }
}



#' Estimate the covariance matrix of standardized allele frequencies
#'
#' In each cycle, samples SNPs matched to the GWAS SNPs in every occupied bin of
#' the contingency table, extracts their frequencies from full.dataset.file and
#' computes the M x M matrix
#'   U = sum over l of s[l,] %*% t(s[l,]) / (nSNP - 1),
#' where s[l,i] = (p[l,i] - eps[l]) / sqrt(eps[l] * (1 - eps[l])). The matrix is
#' then averaged over cycles.
#'
#' Unlike Berg's code, which used cov() on p / sqrt(eps * (1 - eps)), frequencies
#' are mean-centered per SNP across populations, and not across SNPs. Rows and
#' columns of U therefore sum to 0: U is singular and must be used through T.mat
#' or through a generalized inverse. The name uncentered.cov.mat is kept from the
#' original code.
#'
#' @param gwas.data GWAS table with the *.BINS columns.
#' @param match.pop Candidate SNP table with the *.BINS columns.
#' @param pop.names Test population names.
#' @param bin.names Names of the *.BINS columns.
#' @param SNPs.per.cycle Minimum number of SNPs sampled per cycle. The number
#'   actually sampled is ceiling(SNPs.per.cycle / L) * L.
#' @param cycles Number of batches.
#' @param path Output directory.
#' @param full.dataset.file Path to the genome-wide frequency file.
#' @param env.var.data List of environmental tables, used for population names.
#' @param scripts_path Directory that contains sampleSNPs.pl.
#'
#' @return The M x M matrix U, also saved to <path>/uncenteredcovmat.Rdata.
SampleCovSNPs <- function ( gwas.data , match.pop , pop.names , bin.names , SNPs.per.cycle , cycles , path , full.dataset.file , env.var.data , scripts_path) {
  
  #recover()
  num.pops <- length ( pop.names )
  #T.mat <- matrix ( rep ( c ( ( num.pops - 1 ) / num.pops , rep ( - 1 / ( num.pops ) , times = num.pops ) ) , times = num.pops ) , ncol = num.pops , nrow = num.pops )
  gwas.cont.table <- table ( gwas.data [ , bin.names ] ) 
  my.cov.bins <- which ( gwas.cont.table != 0 )
  # Number of matched SNPs to sample per GWAS SNP, so that the bin proportions of the GWAS SNPs are preserved in every cycle
  total.cov.reps <- ceiling ( SNPs.per.cycle / nrow ( gwas.data ) ) ##WARNING: TO INSPECT
  this.many <- total.cov.reps * gwas.cont.table	
  
  snp.by.pop <- vector("list", cycles)
  epsilon.cov.snps <- vector("list", cycles)
  var.cov.snps <- vector("list", cycles)
  scaled.snp.by.pop <- vector("list", cycles)
  uncentered.cov.mat <- vector("list", cycles)
  
    for(k in 1:cycles){
      sampled.SNPs <- list ()
      j = 1 
      for ( BIN in my.cov.bins ) {
        this.many [ BIN ]
        this.bin <- mapply ( function ( x , y ) x [ y ] ,
                             x = dimnames ( gwas.cont.table ) , y = arrayInd ( BIN , dim ( this.many ) ) , SIMPLIFY = FALSE )
        in.this.bin <- list ()
        for ( i in 1 : length ( this.bin ) ) {
          in.this.bin [[ i ]] <- match.pop [ , bin.names [ i ] ] %in% this.bin [ i ] 
        }
        in.this.bin <- do.call ( cbind , in.this.bin )
        matched.SNPs <- match.pop [ rowSums ( in.this.bin ) == length ( bin.names ) , ]$SNP
        # Sampling is with replacement, so a SNP can be drawn several times; its frequencies are duplicated accordingly at cov.freqs below
        sampled.SNPs [[ j ]] <- as.character ( sample ( matched.SNPs , this.many [ BIN ] , replace = T ) )
        j = j + 1 	
      }
      sampled.SNPs.count <- table ( unlist ( sampled.SNPs ) )
      sampled.SNPs <- unlist ( sampled.SNPs )
      write ( unlist ( sampled.SNPs ) , file = paste ( path , "/cov.SNPs" , k , sep = "" ) , ncolumns = 1 )
      # Extract the sampled SNPs from the full dataset (sampleSNPs.pl keeps the header and the listed SNPs)
      system ( paste ( scripts_path,"/sampleSNPs.pl " , path , "/cov.SNPs" , k , " " , full.dataset.file , " > " , path , "/cov.samples" , k , sep = "" ) )
      sampled.cov.data <- fread ( paste ( path , "/cov.samples" , k , sep = "" ) , stringsAsFactors = FALSE , h = TRUE ,data.table = FALSE)
      sampled.cov.data <- sampled.cov.data [ , c("SNP" , "CLST" , "A1" , "A2" , "FRQ") ]
      colnames ( sampled.cov.data ) <- c ( "SNP" , "CLST" , "A1" , "A2" , "FRQ" )
      sampled.cov.data <- sampled.cov.data [ sampled.cov.data$CLST %in% env.var.data[[1]]$CLST , ]
      sampled.cov.data$FRQ <- as.numeric ( sampled.cov.data$FRQ )
      #sampled.SNPs <- read.table ( paste ( path , "/cov.SNPs" , k , sep = "" ) , stringsAsFactors = F ) 
      sampled.cov.data <- sampled.cov.data [ with ( sampled.cov.data , order ( SNP , CLST ) ) , ]
      split.sampled.cov.data <- split ( sampled.cov.data$FRQ , sampled.cov.data$SNP )
      cov.freqs <- mapply ( rep , x = split.sampled.cov.data , times = sampled.SNPs.count )
      
      snp.by.pop [[ k ]] <- t ( matrix ( unlist ( cov.freqs ) , nrow = num.pops ) )
      epsilon.cov.snps [[ k ]] <- apply ( snp.by.pop [[ k ]] , 1 , mean )
      
      # average of ratios
      # Standardize per SNP, then take the cross-product across SNPs
      var.cov.snps [[ k ]] <- epsilon.cov.snps [[ k ]] * ( 1 - epsilon.cov.snps [[ k ]] )
      scaled.snp.by.pop [[ k ]] <- (snp.by.pop [[ k ]] - epsilon.cov.snps [[ k ]]) / c ( sqrt ( var.cov.snps [[ k ]] ))
      uncentered.cov.mat [[ k ]] <- crossprod(scaled.snp.by.pop[[k]]) / (nrow(scaled.snp.by.pop[[k]])-1)
    }
    
  uncentered.cov.mat <- Reduce ( "+" , uncentered.cov.mat)/cycles
  save ( uncentered.cov.mat , file = paste ( path , "/uncenteredcovmat.Rdata" , sep = "" ) )
  return ( uncentered.cov.mat )
}


#' Build the off-center matrices used for the Z scores
#'
#' For each group in env.var.data$REG, builds a transformation that centers
#' populations outside the group on their own mean, and leaves populations in the
#' group uncentered. One row is dropped at random among populations outside the
#' group to make the covariance matrix invertible. These matrices are used only
#' by the regional and individual Z scores, which were not used in our study.
#'
#' @param env.var.data Table with columns CLST and REG.
#' @param uncentered.cov.mat The M x M matrix U.
#' @param effects GWAS effects.
#' @param freqs Frequency table.
#'
#' @return A list of two lists, one element per group: the off-center covariance
#'   matrices and the corresponding transformation matrices.
GetOffCovMats <- function ( env.var.data , uncentered.cov.mat , effects , freqs ) {
  #recover()
  mat.cols <- list ()
  add.vars.regional <- list ()
  add.vars.individual <- list ()
  k <- 1
  for ( i in sort ( unique ( env.var.data$REG ) ) ) {
    
    num.pops <- nrow ( env.var.data )
    num.center.pops <- sum ( env.var.data$REG != i )
    j <- 1
    mat.cols [[ k ]] <- lapply ( env.var.data$REG , function ( x ) 
      if ( x != i ) { 
        y <- rep ( - 1 / num.center.pops , num.pops ) ; 
        y [ j ] <- ( num.center.pops - 1 ) / num.center.pops ; 
        j <<- j + 1
        return ( y ) 
      } else { 
        y <-  rep ( 0 , num.pops ) ;
        y [ j ] <- 1 ;
        j <<- j + 1 ;
        return ( y ) 
      } 
    )
    epsilons <- with ( subset ( freqs , freqs$CLST %in% env.var.data$CLST [ env.var.data$REG != i ] ) , tapply ( FRQ , SNP , mean ) )
    add.vars.regional [[ k ]] <- 4 * sum ( epsilons * ( 1 - epsilons ) * effects ^2 )
    k <- k + 1
  }
  T.mats <- lapply ( mat.cols , function ( x ) do.call ( cbind , x ) )
  T.mats <- lapply ( T.mats , function ( x ) x [ - sample ( which ( x[1,] != 1 & x[1,] != 0 ) , 1 ) , ] )
  off.center.cov.mats <- lapply ( T.mats , function ( x ) x %*% uncentered.cov.mat %*% t ( x ) )
  return ( list ( off.center.cov.mats , T.mats ) )
}


#' Simulate null statistics from matched SNPs
#'
#' For each GWAS SNP, samples reps matched SNPs from the same joint bin and gives
#' them the effect of that GWAS SNP, signed so that it applies to the minor
#' allele. CalcStats() is then applied to each of the reps null genetic values.
#' Each cycle repeats this procedure, so the output contains reps * cycles
#' replicates.
#'
#' @param gwas.data GWAS table with MA.EFF and the *.BINS columns.
#' @param match.pop Candidate SNP table with the *.BINS columns.
#' @param env.var.data List of environmental tables.
#' @param pop.names Test population names.
#' @param uncentered.cov.mat The M x M matrix U.
#' @param F.inv,C.inv Unused; kept from the original code.
#' @param T.mat (M-1) x M mean-centering matrix.
#' @param T.mat.full M x M mean-centering matrix.
#' @param covar Covariate table, or NULL.
#' @param regional.off.center.cov.mats,regional.off.center.T.mats,individual.off.center.cov.mats,individual.off.center.T.mats
#'   Output of GetOffCovMats().
#' @param bin.names Names of the *.BINS columns.
#' @param reps Number of null genetic values per cycle.
#' @param cycles Number of cycles.
#' @param path Output directory.
#' @param scripts_path Directory that contains sampleSNPs.pl.
#' @param full.dataset.file Path to the genome-wide frequency file.
#'
#' @return A list of null statistics: vectors Qx, Fst.comp, LD.comp (and Qx.res,
#'   Fst.comp.res, LD.comp.res with covariates); matrices betas, pearson.rs,
#'   spearman.rhos and ind.Zs with one column per replicate; and the list reg.Zs,
#'   with one matrix per environmental file.
NullStats <- function ( gwas.data , match.pop , env.var.data , pop.names , uncentered.cov.mat , F.inv , C.inv , T.mat, T.mat.full , covar=NULL, regional.off.center.cov.mats , regional.off.center.T.mats , individual.off.center.cov.mats , individual.off.center.T.mats , bin.names , reps , cycles , path , scripts_path, full.dataset.file) {
  
  num.pops <- length ( pop.names )
  gwas.cont.table <- table ( gwas.data [ , bin.names ] )
  my.null.bins <- which ( gwas.cont.table != 0 )
  this.many <- reps * gwas.cont.table 
  null.stats <- list ()
  
    for(j in 1:cycles){
      sampled.SNPs <- list ()
      effects.list <- list ()
      for ( i in 1 : nrow ( gwas.data ) ) {
        
        # sample a set of <reps> SNPs to match gwas SNP <i> such that they match its properties jointly for all matching categories
        this.snp <- gwas.data [ i , bin.names ]
        in.this.bin <- list ()
        for ( k in 1 : length ( this.snp ) ) {
          in.this.bin [[ k ]] <-  match.pop [ , bin.names [ k ] ] %in% this.snp [[ k ]] 
        }
        in.this.bin <- Reduce ( "&" , in.this.bin )
        matched.SNPs <- match.pop$SNP[ in.this.bin ]
        sampled.SNPs [[ i ]] <- as.character ( sample ( matched.SNPs , reps , replace = T ) )
      # match.pop$FRQ is the frequency of the allele tracked in the full dataset. If it is the minor allele, it takes MA.EFF; otherwise it takes -MA.EFF
        effects.list [[ i ]] <- ifelse ( 
          match.pop [ match ( sampled.SNPs [[ i ]] , match.pop$SNP ) , ]$FRQ <= 0.5, # Berg's code used >=, which flips the sign of all null effects and leaves Qx unchanged
          gwas.data$MA.EFF [ i ] ,
          - gwas.data$MA.EFF [ i ] 
        )
      }
      all.sampled.SNPs <- unique ( unlist ( sampled.SNPs ) )
      write ( all.sampled.SNPs , file = paste ( path , "/null.SNPs" , j , sep = "" ) , ncolumns = 1 )
      system ( paste ( scripts_path,"/sampleSNPs.pl " , path , "/null.SNPs" , j , " " , full.dataset.file , " > " , path , "/null.samples" , j , sep = "" ) )
      sampled.null.data <- fread ( paste ( path , "/null.samples" , j , sep = "" ) , stringsAsFactors = F , h = T  ,data.table = F)
      sampled.null.data <- sampled.null.data [ , c ( "SNP" , "CLST" , "A1" , "A2" , "FRQ" ) ]
      colnames ( sampled.null.data ) <- c ( "SNP" , "CLST" , "A1" , "A2" , "FRQ" )
      sampled.null.data <- sampled.null.data [ sampled.null.data$CLST %in% env.var.data[[1]]$CLST , ]
      sampled.null.data <- sampled.null.data [ with ( sampled.null.data  , order ( SNP , CLST ) ) , ]
      split.sampled.null.data <- split ( sampled.null.data , sampled.null.data$CLST )
      null.SNP.mat <- do.call ( cbind , sampled.SNPs )
    # One reps x L frequency matrix per population, in the order of sampled.SNPs
      null.freqs <- lapply(split.sampled.null.data, function(this.pop) {
        frqs <- this.pop[match(null.SNP.mat, this.pop$SNP), "FRQ"]
        as.numeric(frqs)
      })
      null.freqs<- lapply(null.freqs, function(vec) {
        matrix(vec, nrow = reps, byrow = FALSE)
      })
      
      effects.mat <- do.call ( cbind , effects.list )
      new.effects.list <- alply ( effects.mat ,1 , function ( x ) x )
      null.gvs <- lapply ( null.freqs , function ( x ) rowSums ( x * effects.mat ) )
      mean.null.freqs <- Reduce ( "+" , null.freqs ) / length ( null.freqs )
      
      null.vars <- rowSums ( effects.mat^2 * mean.null.freqs * ( 1 - mean.null.freqs ) )
    # Reshape to a list of reps matrices of dimension M x L, as expected by CalcStats()
      null.freqs <- alply ( aperm ( laply ( null.freqs , function (x ) x ) , c ( 1, 3, 2 ) ) , 3 , function ( y ) y )
      
      null.stats[[j]] <- mapply(function(FREQ, VAR, EFFECTS) {
        CalcStats(FREQ, EFFECTS, env.var.data, VAR, uncentered.cov.mat, T.mat, T.mat.full, regional.off.center.cov.mats, regional.off.center.T.mats, individual.off.center.cov.mats, individual.off.center.T.mats, covar)
      }, FREQ = null.freqs, VAR = null.vars, EFFECTS = new.effects.list, SIMPLIFY = FALSE)
      
    }		
    
  # The output of CalcStats() has more elements with covariates, hence the two sets of indices below
  Qx <- unlist ( lapply ( null.stats , function ( y ) lapply ( y , function ( x ) x [[ 1 ]] ) ) )
  Fst.comp <- unlist ( lapply ( null.stats , function ( y ) lapply ( y , function ( x ) x [[ 2 ]] ) ) )
  LD.comp <- unlist ( lapply ( null.stats , function ( y ) lapply ( y , function ( x ) x [[ 3 ]] ) ) )
   if(!(is.null(covar))){
     Qx.res <- unlist ( lapply ( null.stats , function ( y ) lapply ( y , function ( x ) x [[ 5 ]] ) ) )
     Fst.comp.res <- unlist ( lapply ( null.stats , function ( y ) lapply ( y , function ( x ) x [[ 6 ]] ) ) )
     LD.comp.res <- unlist ( lapply ( null.stats , function ( y ) lapply ( y , function ( x ) x [[ 7 ]] ) ) )
     betas <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) unlist ( x [[ 10 ]] ) ) ) ) )
     pearson.rs <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) unlist ( x [[ 11 ]] ) ) ) ) )
     spearman.rhos <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) unlist ( x [[ 12 ]] ) ) ) ) )
     reg.Zs <- lapply (  1 : length ( env.var.data ) , function ( u ) do.call ( cbind , lapply ( lapply ( null.stats , function ( w ) lapply ( lapply ( 1 : length ( env.var.data ) , function ( y ) lapply ( w , function ( x ) x [[ 13 ]] [ , y ] ) ) , function ( z ) as.matrix ( do.call ( cbind , z  ) ) ) ) , function ( v ) v [[ u ]] ) ) )
     ind.Zs <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) x [[ 14 ]] ) ) ) )
     rownames ( ind.Zs ) <- env.var.data [[ 1 ]]$CLST
      output <- list ( Qx = Qx , Fst.comp = Fst.comp , LD.comp = LD.comp , Qx.res = Qx.res, Fst.comp.res = Fst.comp.res, LD.comp.res = LD.comp.res , betas = betas , pearson.rs = pearson.rs , spearman.rhos = spearman.rhos , reg.Zs = reg.Zs , ind.Zs = ind.Zs )
   }else{
     betas <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) unlist ( x [[ 5 ]] ) ) ) ) )
     pearson.rs <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) unlist ( x [[ 6 ]] ) ) ) ) )
     spearman.rhos <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) unlist ( x [[ 7 ]] ) ) ) ) )
     reg.Zs <- lapply (  1 : length ( env.var.data ) , function ( u ) do.call ( cbind , lapply ( lapply ( null.stats , function ( w ) lapply ( lapply ( 1 : length ( env.var.data ) , function ( y ) lapply ( w , function ( x ) x [[ 8 ]] [ , y ] ) ) , function ( z ) as.matrix ( do.call ( cbind , z  ) ) ) ) , function ( v ) v [[ u ]] ) ) )
     ind.Zs <- do.call ( cbind , lapply ( null.stats , function ( y ) do.call ( cbind , lapply ( y , function ( x ) x [[ 9 ]] ) ) ) )
     rownames ( ind.Zs ) <- env.var.data [[ 1 ]]$CLST
     output <- list ( Qx = Qx , Fst.comp = Fst.comp , LD.comp = LD.comp , betas = betas , pearson.rs = pearson.rs , spearman.rhos = spearman.rhos , reg.Zs = reg.Zs , ind.Zs = ind.Zs ) 
   }
  
  return ( output )
}


#' Partition Qx across principal components
#'
#' Computes, for each eigenvector of F.mat * var, the squared standardized
#' projection of the centered genetic values. These sum to Qx. The projection is
#' computed on the unadjusted genetic values, also when covariates are supplied.
#'
#' @param F.mat Centered (M-1) x (M-1) covariance matrix.
#' @param cent.gvs Centered genetic values, T.mat %*% Z.
#' @param var Scaling term V.
#' @param path Output directory.
#'
#' @return Nothing. Writes <path>/figs/PCAQx.pdf, <path>/Output/PCVarProj.Rdata
#'   and <path>/Output/PCs.Rdata.
PCAProjections <- function ( F.mat , cent.gvs , var , path) {
  #recover()
  my.eig <- eigen ( F.mat * var )
  var.proj <- numeric ( length ( my.eig$values ) )
  for ( i in 1 : length ( my.eig$values ) ) {
    var.proj [ i ] <- ( ( my.eig$vectors [ , i ] %*% cent.gvs ) / sqrt ( my.eig$values [ i ] ) )^2
  }
  dir.create ( paste ( path , "/figs" , sep = "" ) )
  pdf ( paste ( path , "/figs/PCAQx.pdf" , sep = "" ) )
  plot ( var.proj , pch = 20 , xlab = "Principal Component" , ylab = expression ( paste ( "PCA " , Q[X] , sep = "" ) ) )
  dev.off()
  save ( var.proj , file = paste ( path , "/Output/PCVarProj.Rdata" , sep = "" ) )
  save ( my.eig , file = paste ( path , "/Output/PCs.Rdata" , sep = "" ) )
}


#' Compute Qx and related statistics for one set of frequencies and effects
#'
#' Applied once to the GWAS SNPs and once per null replicate.
#'
#' Qx is computed from the M-1 centered contributions as
#'   Qx = t(T %*% Z) %*% solve(F.mat) %*% (T %*% Z) / V.
#' It is split into Fst.comp (diagonal of the SNP-by-SNP matrix, the signal in
#' single-SNP divergence) and LD.comp (off-diagonal, the signal in covariance
#' between SNPs), so that Qx = Fst.comp + LD.comp.
#'
#' With J covariates in the mean-centered matrix X, GLS regression of the
#' centered genetic values on X gives
#'   beta.hat = solve(t(W) %*% F.inv %*% W) %*% t(W) %*% F.inv %*% (T %*% Z),
#' where W is X without the dropped population. Qx.res is the same quadratic form
#' on the residuals, with chi-squared distribution on M-J-1 df under neutrality.
#'
#' @param freqs M x L frequency matrix, populations in rows and SNPs in columns,
#'   SNPs sorted as in effects.
#' @param effects Vector of L effects.
#' @param env.var.data List of environmental tables.
#' @param var Scaling term V.
#' @param uncentered.cov.mat The M x M matrix U.
#' @param T.mat (M-1) x M mean-centering matrix.
#' @param T.mat.full M x M mean-centering matrix, diag(M) - 1/M.
#' @param regional.off.center.cov.mats,regional.off.center.T.mats,individual.off.center.cov.mats,individual.off.center.T.mats
#'   Output of GetOffCovMats().
#' @param covar Covariate table, or NULL.
#' @param path Output directory, used when null = F or PCProj = T.
#' @param PCProj If T, call PCAProjections().
#' @param null If F, save the genetic values to Output/genetic.values.Robj.
#'
#' @return A list with Qx, Fst.comp, LD.comp, Qx.contrib, betas, pearson.rs,
#'   spearman.rhos, reg.Z and ind.Z. With covariates, also Qx.res, Fst.comp.res,
#'   LD.comp.res, Qx.res.contrib and Qx.cov.contrib. See README.md.
CalcStats <- function ( freqs , effects , env.var.data , var , uncentered.cov.mat , T.mat, T.mat.full , regional.off.center.cov.mats , regional.off.center.T.mats , individual.off.center.cov.mats , individual.off.center.T.mats, covar=NULL , path , PCProj = F , null = T ) {
  
  # mean center and drop a pop
  F.mat <- T.mat %*% ( uncentered.cov.mat ) %*% t ( T.mat )
  F.inv <- solve ( F.mat )
  C.inv <- solve ( t ( chol ( F.mat ) ) )
 
  # Qx and components
  contribs <- T.mat %*% t ( t ( freqs ) * effects )
  
  std.contribs <- C.inv  %*% contribs
  adj.cov.mat <- t ( std.contribs ) %*% std.contribs
  Fst.comp <- sum ( diag ( adj.cov.mat ) ) / var
  LD.comp <- sum ( adj.cov.mat [ row ( adj.cov.mat ) != col ( adj.cov.mat ) ] ) / var
  Qx <- Fst.comp + LD.comp 
  
  # Per-population contributions. With Zc = 2 * (Z - mean(Z)), Qx.contrib[i] = Zc[i] * (F+ %*% Zc)[i] / (4 * V), where F+ is the Moore-Penrose inverse of the M x M centered matrix. Contributions sum to Qx and can be negative
  gvs <- freqs %*% effects
  Zc <- (2*gvs) - mean(2*gvs)
  F.full.inv=ginv(T.mat.full %*% uncentered.cov.mat %*% t(T.mat.full))
  Qx.contrib <- as.numeric( (1/(4*var))*Zc * (F.full.inv %*% Zc))  
  
  if(!is.null(covar)){
    # NOTE: W is taken directly from the covariate file, so it equals T.mat %*% X only if each covariate is mean-centered across the M populations. Center covariates before running: otherwise Qx.res is wrong and no longer equals the sum of Qx.res.contrib
    drop <- which ( apply ( T.mat , 2 , function ( x ) !any ( x <= 1 & x >= 0  ) ) ) # line related to figuring out which population i dropped in earlier steps
    W=as.matrix(covar[-c(drop),-c(1)])
    #contribs=T.mat%*%contribs
    # Residualize each SNP's contributions. The GLS fit is linear, so residuals sum across SNPs to the residuals of the genetic values, and Qx.res splits into Fst and LD components as for Qx
    residuals_matrix <- apply(X=contribs,MARGIN=2, function(snp_contribs,W,F.inv) {
      # Fit a GLS model of SNP contributions on covariates
      Xt_Finv <- t(W) %*% F.inv
      beta_gls <- solve(Xt_Finv %*% W) %*% Xt_Finv %*% snp_contribs
      residuals= snp_contribs - W%*%beta_gls # Extract residuals
      return(residuals)
    },W=W,F.inv=F.inv)
    # Transpose the result to match the original dimensions
    if (all(dim(residuals_matrix) == dim(contribs))) {
      contribs <- residuals_matrix
    }else{contribs <- t(residuals_matrix)}
    
    # Same quadratic form and decomposition as above, on the residuals
    std.contribs <- C.inv  %*% contribs
    adj.cov.mat <- t ( std.contribs ) %*% std.contribs
    Fst.comp.res <- sum ( diag ( adj.cov.mat ) ) / var
    LD.comp.res <- sum ( adj.cov.mat [ row ( adj.cov.mat ) != col ( adj.cov.mat ) ] ) / var
    Qx.res <- Fst.comp.res + LD.comp.res
    
    # Repeat the GLS fit on all M populations to split Qx across populations: Qx.res.contrib sums to Qx.res, and Qx.cov.contrib to Qx - Qx.res (the part explained by the covariates)
    F.full.inv=ginv(T.mat.full %*% uncentered.cov.mat %*% t(T.mat.full))
    W=as.matrix(covar[,-c(1)])
    Xt_Finv <- t(W) %*% F.full.inv
    beta_gls <- solve(Xt_Finv %*% W) %*% Xt_Finv %*% Zc
    residuals= Zc - W%*%beta_gls # Extract residuals
    #t(W%*%beta_gls)%*%F.full.inv%*%(W%*%beta_gls)/(4*var) + t(residuals)%*%F.full.inv%*%residuals/(4*var) #Should be equal to Qx
    
    Qx.res.contrib <- as.numeric( (1/(4*var))*residuals * (F.full.inv %*% residuals))  # elementwise multiply
    Qx.cov.contrib <- as.numeric( (1/(4*var))* (W%*%beta_gls) * (F.full.inv %*% (W%*%beta_gls)))  # elementwise multiply
    #sum(c.cov) + sum(c.res)
    
  }
  
  # env var stats
  gvs <- freqs %*% effects
  cent.gvs <- T.mat %*% gvs
  std.gvs <- C.inv %*% cent.gvs
  #t(std.gvs)%*%std.gvs/var
  
  
  std.env <- lapply ( env.var.data , function ( x ) C.inv %*% T.mat %*% x$ENV )
  betas <- lapply ( std.env , function ( x ) sum ( std.gvs * x ) / sum ( x^2 ) )
  pearson.rs <- lapply ( std.env , function ( x ) sum ( std.gvs * x ) / sqrt ( sum ( x^2 )*sum ( std.gvs^2 ) ) )
  spearman.rhos <- lapply ( std.env , function ( x ) cor ( std.gvs , x , method = "spearman" ) )
  if ( null == F ) save ( gvs , file = paste ( path , "/Output/genetic.values.Robj" , sep = "" ) )
  if ( PCProj == T ) {
    PCAProjections ( F.mat , cent.gvs , var , path )
  }
  # regional Z scores (not used in our study)
  reg.Z.scores <- mapply ( function ( COV , TMAT ) 
    mapply ( function ( THIS.COV , THIS.TMAT ) 
      ZStats ( gvs , var , THIS.COV , THIS.TMAT ) , 
      THIS.COV = COV , THIS.TMAT = TMAT ) ,
    COV = regional.off.center.cov.mats , TMAT = regional.off.center.T.mats
  )
  ## individual Z scores (not used in our study)
  ind.Z.scores <- mapply ( function ( THIS.COV , THIS.TMAT ) 
    ZStats ( gvs , var , THIS.COV , THIS.TMAT ) , 
    THIS.COV = individual.off.center.cov.mats , THIS.TMAT = individual.off.center.T.mats 
  )
  
  if(!is.null(covar)){
    output=list ( Qx = Qx , Fst.comp = Fst.comp , LD.comp = LD.comp , Qx.contrib = Qx.contrib, Qx.res = Qx.res , Fst.comp.res = Fst.comp.res , LD.comp.res = LD.comp.res , Qx.res.contrib = Qx.res.contrib, Qx.cov.contrib = Qx.cov.contrib ,betas = betas , pearson.rs = pearson.rs , spearman.rhos = spearman.rhos , reg.Z = reg.Z.scores , ind.Z = ind.Z.scores )
  }else{
    output=list ( Qx = Qx , Fst.comp = Fst.comp , LD.comp = LD.comp , Qx.contrib = Qx.contrib , betas = betas , pearson.rs = pearson.rs , spearman.rhos = spearman.rhos , reg.Z = reg.Z.scores , ind.Z = ind.Z.scores )
    }
  return ( output ) 
}


#' Compute the Z score of one group of populations
#'
#' Compares the summed genetic value of the group with its expectation given the
#' other populations, using condNormal(). Not used in our study.
#'
#' @param gvs Vector of M genetic values.
#' @param var Scaling term V.
#' @param cov.mat Off-center covariance matrix.
#' @param T.mat Off-center transformation matrix.
#'
#' @return The Z score.
ZStats <- function ( gvs , var , cov.mat , T.mat ) {
  #recover()
  #cent.gvs <- T.mat %*% gvs
  drop <- which ( apply ( T.mat , 2 , function ( x ) !any ( x <= 1 & x >= 0  ) ) ) # line related to figuring out which population i dropped in earlier steps
  need <- which ( T.mat [ 1 , ] == 1 | T.mat [ 1 , ] == 0 )
  conditioned <- seq_along ( gvs ) [ - c ( need , drop ) ]
  shift.conditioned <- c ( conditioned [ conditioned < drop ] , conditioned [ conditioned > drop ] -1 )
  shift.need <- ifelse ( need > drop , need - 1 , need )
  cond.dist <- condNormal ( gvs [ conditioned ] , rep ( mean ( gvs [ - need ] ) , length ( gvs ) - 1 ) , var * cov.mat , shift.conditioned , shift.need )
  my.Z <- ( sum ( gvs [ need ] ) - sum ( cond.dist$condMean ) ) / sqrt ( sum ( cond.dist$condVar ) )
  return ( my.Z )
  
}

condNormal <- function(x.given, mu, sigma, given.ind, req.ind){
  # Returns conditional mean and variance of x[req.ind] 
  # Given x[given.ind] = x.given
  # where X is multivariate Normal with
  # mean = mu and covariance = sigma
  #recover()
  B <- sigma[req.ind, req.ind]
  C <- sigma[req.ind, given.ind, drop=FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% solve(D)
  cMu <- c(mu[req.ind] + CDinv %*% (x.given - mu[given.ind]))
  cVar <- B - CDinv %*% t(C)
  return ( list(condMean=cMu, condVar=cVar) )
}


#' Align GWAS alleles with the frequency file
#'
#' Orders the GWAS table by SNP. Where A1 and A2 are swapped relative to
#' assoc.loci.freqs, swaps them, and flips EFF and FRQ. Stops with "Strand Issue"
#' if alleles match in neither orientation, so strand-ambiguous errors are not
#' detected at A/T and C/G SNPs where both files have the same alleles.
#'
#' @param gwas.data GWAS table with columns SNP, A1, A2, EFF and FRQ.
#' @param assoc.loci.freqs Frequency table with columns SNP, CLST, A1 and A2.
#'
#' @return The GWAS table sorted by SNP and aligned to assoc.loci.freqs.
MatchAlleles <- function ( gwas.data , assoc.loci.freqs ) {
  
  #protect <- gwas.data
  #protect -> gwas.data
  correct.orientation <- subset ( assoc.loci.freqs , CLST == assoc.loci.freqs$CLST [ 1 ] )
  correct.orientation <- correct.orientation [ order ( correct.orientation$SNP ) , ]
  #gwas.data <- gwas.data [ match ( correct.orientation$SNP , gwas.data$SNP ) , ]
  gwas.data <- gwas.data [ order ( gwas.data$SNP ) , ]
  if ( nrow ( gwas.data ) != nrow ( correct.orientation ) ) {
    stop ( "GWAS dataset and Orientation Matching dataset contain differing numbers of SNPs" )
  }
  
  flip <- gwas.data$A1 == correct.orientation$A2 & gwas.data$A2 == correct.orientation$A1
  dont.flip <- gwas.data$A1 == correct.orientation$A1 & gwas.data$A2 == correct.orientation$A2
  for ( i in 1 : nrow ( gwas.data ) ) {
    if ( flip [ i ] ) {
      gwas.data$A1 [ i ] <- correct.orientation$A1 [ i ]
      gwas.data$A2 [ i ] <- correct.orientation$A2 [ i ]
      gwas.data$EFF [ i ] <- - gwas.data$EFF [ i ]
      gwas.data$FRQ [ i ] <- 1 - gwas.data$FRQ [ i ]
    } else if ( dont.flip [ i ] ) {
      #do nothing
    } else {
      stop ( "Strand Issue")		
    }
  }
  # if ( sum ( flip ) != 0 ) {
  # cat ( "The following allele orientations were flipped:\n" , file = log.file )
  # cat ( c ( gwas.data [ flip , "SNP" ] ) , sep = "\n" , file = log.file )	
  # }
  return ( gwas.data )
}


## Example command to run a Qx analysis. Replace the objects below with your own
## file paths. The example is inside if ( FALSE ) so that it does not run when
## this file is sourced: run its content directly.

if ( FALSE ) {
  library(plyr)
  library(data.table)
  library(MASS)
  PolygenicAdaptationFunction (
    gwas.data.file =  gwas.file,
    freqs.file = freq.file ,
    env.var.data.files = EnvFile, # Note: you can supply as many env.var.data.files concurrently as you want, supplied as a list
    match.pop.file = pop.file ,
    full.dataset.file = data.file,
    path = yourPath,
    scripts_path=path.to.function.script,
    match.categories = c ("MAF" , "ldscore") ,
    match.bins = list ( seq ( 0 , 0.5 , 0.02 ), c(0.4363, 13.4999, 20.9151, 32.0346, 305.1991002) ) ,
    cov.SNPs.per.cycle = 20000 ,
    cov.cycles = 5 ,
    null.phenos.per.cycle = 10000 ,
    null.cycles = 10 ,
    load.cov.mat = F ,
    sim.null = T ,
    check.allele.orientation = T,
    PCProj=T,
    covar.file=NULL #NULL if no covariates
  )
}
