#'The Superwise function
#'
#' Pairwise statistical tests! Wiser than Supermatrix. This function allows you to perform multiple pairwise association tests between two sets of data, containing numeric/categorical and categorical data, respectively, and to adjust reulting p-values by different methods. See \code{details} for more information.
#' @param x A data frame containing your data. See \code{details} for more information.
#' @param v1 A range of columns containing phenotypic data in x. Data can be either quantitative (numeric) or qualitative (factor). Example: 1:20, for 20 columns of trait data, located in the first 20 columns.
#' @param v2 A range of columns containing genotypic data in x. Example: 20:40, for 20 sites located in columns 20 to 40.
#' @param fisher.method Indicates which method should be used for Fisher's tests. \code{"exact"} is default. For big datasets with a high number of categories, other methods might be useful. "sim" simulates p-values using 20000 Monte Carlo iterations, and "hybrid" will use hybrid method (see \code{\link{fisher.test}} for more details on Fisher's methods).
#' @param correction In case that p-value correction for multiple tests is required, use this parameter to specify which correction must be used. See ?p.adjust for details of available correction methods. It is possible to have more than one correction in one go, using a vector with correction names in the \code{\link{p.adjust}} function fashion (e.g. c("bonferroni", "fdr", "hochberg")). Each correction will add a column to the output.
#' @param control Used as a internal control helper, to check that supermatrix performed the tests accordingly to variables' type.If TRUE, two columns indicating whether (1) the trait value was identified as numeric, and (2) the name of the test used for that column vs. genotypic data.
#'
#' @return A data frame with the results of each test, including P-values and adjusted P-values (if requested).
#'
#' @details Superwise function was designed to perform multiple pairwise association analyses between two sets of phenotypic vs genotypic data.
#' 		Set1 must ideally contain phenotypic data (in columns), which can be quantitative (e.g. number of offspring, weight, etc) or categorical unordered (e.g. diurnal/nocturnal, herbivore/carnivore/omnivore) data.
#'		Set2 must ideally contain categorical unordered data, such as which amino-acid in a specific site in a proteic sequence  Such data should be distributed in columns, using rows for each species to be included in the analyses.
#' 		So far, superwise performs only Kruskal-Wallis test for numeric vs. categorical variables, and Fisher's Exact Test for categorical vs. categorical variables.
#' 		The recommended format is a single data.frame containing the two sets of data to be compared pairwise.
#'		NOTE: Categorical data must be of factor class.
#' @export

spw <- function(x, v1, v2, fisher.method="exact", correction = NULL, control = TRUE){
  #We start by making sure that all character and logical columns are properly transformed to factors
  x[sapply(x, is.character)] <- lapply(x[sapply(x, is.character)], as.factor)
  x[sapply(x, is.logical)] <- lapply(x[sapply(x, is.logical)], as.factor)
  x[sapply(x[v2], is.numeric)] <- lapply(x[sapply(x[v2], is.numeric)], as.factor)
  #Then we our dataframe and creating empty vectors
  v1.tab <- x[,v1] #Dataframe with characteristics data
  v2.tab <- x[,v2] #Dataframe with mutation data
  Var1 <- character()
  Var2 <- character()
  numeric <- logical()
  Test <- character()
  FET.IC1 <- numeric()
  FET.IC2 <- numeric()
  FET.OR <- numeric()
  KW.statistic <- numeric()
  KW.df <- integer()
  P.value <- numeric()
  i <- 1 #Guide for vectors

for(z in 1:length(v1.tab[1,])){
    if (is.numeric(v1.tab[,z])==TRUE){
      for (y in 1:length(v2.tab[1,])){
        cc <- na.omit(cbind(v1.tab[,z], v2.tab[,y]))
        if (nlevels(as.factor(cc[,2])) >= 2L){
          KW <- kruskal.test(v1.tab[,z], v2.tab[,y])
          Var1[i] <- names(v1.tab[z])
          Var2[i] <- names(v2.tab[y])
          numeric[i] <- is.numeric(v1.tab[,z])
          Test[i] <- "Kruskal-Wallis"
          FET.IC1[i] <- NA
          FET.IC2[i] <- NA
          FET.OR[i] <- NA
          KW.statistic[i] <- KW$statistic[[1]]
          KW.df[i] <- KW$parameter[[1]]
          P.value[i] <- KW$p.value
          i <- i + 1
        }
        else {
          Var1[i] <- names(v1.tab[z])
          Var2[i] <- names(v2.tab[y])
          numeric[i] <- TRUE
          Test[i] <- "No test applied"
          FET.IC1[i] <- NA
          FET.IC2[i] <- NA
          FET.OR[i] <- NA
          KW.statistic[i] <- NA
          KW.df[i] <- NA
          P.value[i] <- NA
          i <- i + 1
              }
      }
    }
    else{ for (y in 1:length(v2.tab[1,])){
            cc <- na.omit(cbind(v1.tab[,z], v2.tab[,y]))
            if (nlevels(as.factor(cc[,2])) >=2L){
              if (fisher.method == "exact"){
                fet <- fisher.test(table(v1.tab[,z], v2.tab[,y]), workspace = 2e5)
                dins <- sum(dim(table(v1.tab[,z], v2.tab[,y])))
                Var1[i] <- names(v1.tab[z])
                Var2[i] <- names(v2.tab[y])
                numeric[i] <- is.numeric(v1.tab[,z])
                Test[i] <- "Fisher's Exact Test"
                if (dins == 4){
                FET.IC1[i] <- fet$conf.int[[1]]
                FET.IC2[i] <- fet$conf.int[[2]]
                FET.OR[i] <- fet$estimate[[1]]
                }
                else{
                FET.IC1[i] <- NA
                FET.IC2[i] <- NA
                FET.OR[i] <- NA
                }
                KW.statistic[i] <- NA
                KW.df[i] <- NA
                P.value[i] <- fet$p.value
                i <- i + 1
              }
              if (fisher.method == "sim"){
                fet <- fisher.test(table(v1.tab[,z], v2.tab[,y]), simulate.p.value = T, B = 20000)
                dins <- sum(dim(table(v1.tab[,z], v2.tab[,y])))
                Var1[i] <- names(v1.tab[z])
                Var2[i] <- names(v2.tab[y])
                numeric[i] <- is.numeric(v1.tab[,z])
                Test[i] <- "Fisher's Exact Test (Sim)"
                if (dins == 4){
                  FET.IC1[i] <- fet$conf.int[[1]]
                  FET.IC2[i] <- fet$conf.int[[2]]
                  FET.OR[i] <- fet$estimate[[1]]
                }
                else{
                  FET.IC1[i] <- NA
                  FET.IC2[i] <- NA
                  FET.OR[i] <- NA
                }
                KW.statistic[i] <- NA
                KW.df[i] <- NA
                P.value[i] <- fet$p.value
                i <- i + 1
              }
              if (fisher.method == "hybrid"){
                fet <- fisher.test(table(v1.tab[,z], v2.tab[,y]), hybrid = T, workspace = 2e8)
                dins <- sum(dim(table(v1.tab[,z], v2.tab[,y])))
                Var1[i] <- names(v1.tab[z])
                Var2[i] <- names(v2.tab[y])
                numeric[i] <- is.numeric(v1.tab[,z])
                Test[i] <- "Fisher's Exact Test (Hybrid)"
                if (dins == 4){
                  FET.IC1[i] <- fet$conf.int[[1]]
                  FET.IC2[i] <- fet$conf.int[[2]]
                  FET.OR[i] <- fet$estimate[[1]]
                }
                else{
                  FET.IC1[i] <- NA
                  FET.IC2[i] <- NA
                  FET.OR[i] <- NA
                }
                KW.statistic[i] <- NA
                KW.df[i] <- NA
                P.value[i] <- fet$p.value
                i <- i + 1
              }
            }
            else{
              Var1[i] <- names(v1.tab[z])
              Var2[i] <- names(v2.tab[y])
              numeric[i] <- is.numeric(v1.tab[,z])
              Test[i] <- "No test applied"
              FET.IC1[i] <- NA
              FET.IC2[i] <- NA
              FET.OR[i] <- NA
              KW.statistic[i] <- NA
              KW.df[i] <- NA
              P.value[i] <- NA
              i <- i + 1
            }
          }
        }
}
  df <- data.frame(Var1,Var2,numeric,Test, FET.IC1,FET.IC2, FET.OR, KW.statistic, KW.df, P.value)
#Shall we adjust those P.values?
  if (is.null(correction) == FALSE) {
      ps <- as.vector(df$P.value)
      for (e in 1:length(correction)){
        if (correction[e] == "bonferroni"){
          BF.P.value <- p.adjust(ps,method = correction[e])
          df <- cbind(df,BF.P.value)
        }
        if (correction[e] == "hochberg"){
          Hochberg.P.value <- p.adjust(ps,method = correction[e])
          df <- cbind(df,Hochberg.P.value)
        }
        if (correction[e] == "holm"){
          Holm.P.value <- p.adjust(ps,method = correction[e])
          df <- cbind(df,Holm.P.value)
        }
        if (correction[e] == "hommel"){
          Hommel.P.value <- p.adjust(ps,method = correction[e])
          df <- cbind(df,Hommel.P.value)
        }
        if (correction[e] == "fdr"){
          FDR.P.value <- p.adjust(ps,method = correction[e])
          df <- cbind(df,FDR.P.value)
        }
        if (correction[e] == "BH"){
          BH.P.value <- p.adjust(ps,method = correction[e])
          df <- cbind(df,BH.P.value)
        }
        if (correction[e] == "BY"){
          BY.P.value <- p.adjust(ps,method = correction[e])
          df <- cbind(df,BY.P.value)
        }
      }
  if (control == FALSE){
    df <- df[,-c(3:4)]
  }
  return(df)
  }
}
