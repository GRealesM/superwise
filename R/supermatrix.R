#' The Supermatrix function
#'
#' Not really a matrix! This function allows you to perform multiple pairwise association tests between numerical/categorical variables (phenotypes) and aligned amino-acid sites (categorical), create a matrix of p-values and adjust them by several p-value adjustment methods. For numeric vs. categorical data, Kruskal-Wallis test is used. For categorical data, Fisher Exact test is used. See details for more information.
#'
#' @param x A data frame containing your data. See details for more information.
#' @param var A range of columns containing trait data (Set 1) in x. Data can be either quantitative (numeric) or qualitative (factor). Example: 1:20, for 20 columns of trait data, located in the first 20 columns.
#' @param mut A range of columns containing genetypic data (Set 2) in x. Example: 20:40, for 20 sites located in columns 20 to 40.
#' @param fisher.method Indicates which method should be used for fisher's tests. "exact" is default. For big datasets, other methods might be useful. "sim" simulates p-values using 20000 Monte Carlo iterations, and "hybrid" will use hybrid method (see ?fisher.test for more details on Fisher's methods).
#' @param correction In case that p-value correction for multiple tests is required, use this parameter to specify which correction must be used. See \code{\link{?p.adjust}} for details of available correction methods.
#' @param control Used as a internal control helper, to check that supermatrix performed the tests accordingly to variables' type.If TRUE, two columns indicating whether (1) the trait value was identified as numeric, and (2) the name of the test used for that column vs. genotypic data.
#'
#' @return A matrix-like dataframe with P-values for each test performed.
#'
#' @details Supermatrix function was designed to perform multiple association analyses between phenotypic/life history traits, and proteic sites across multiple animal species. Such data should be distributed in columns, using rows for each species to be included in the analyses. The recommended format is a single data.frame containing the two sets of data to be compared pairwise.
#' 		Set1 must ideally contain phenotypic data (in columns), which can be quantitative (e.g. number of offspring, weight, etc) or categorical unordered (e.g. diurnal/nocturnal, herbivore/carnivore/omnivore) data.
#'		Set2 must ideally contain categorical unordered data, such as which amino-acid in a specific site in a proteic sequence.
#' 		IMPORTANT: Categorical data must be of factor class.
#' 		Since \code{supermatrix} only creates a data frame with P-values (which may be insufficient in most cases), we encourage using \code{\link{spw}} function instead.
#' @export

supermatrix <- function(x, var, mut, fisher.method = "exact", correction = NULL, control = TRUE){
  #We create an empty data frame
  rows <- length(var)
  cols <- length(mut) + 2 #+2 in order to leave space for "Categorical" and "Test" control columns
  df <- data.frame(matrix(nrow = rows, ncol = cols))
  names.cols <- c("numeric", "model", names(x[,mut]))
  names.rows <- names(x[,var])
  names(df) <- names.cols
  row.names(df) <- names.rows
  #Let's interrogate about numeric!
  for (i in 1:length(df[,1])) {
    df[i,1] <- is.numeric(x[,i])
  }
  #Now let's perform some statistical tests!
  n <- 1 #row guide variable for df
  m <- 3 #col guide variable for df
  #v1 <- 1 #col guide variable for var.tab
  m2 <- 1 #col guide variable for mut.tab
  var.tab <- x[,var] #Dataframe with characteristics data
  mut.tab <- x[,mut] #Dataframe with mutation data
  for (l in 1:length(var.tab[1,])) { #To run the test row-wise
    if (df[n,1] == TRUE) { #Numeric, not normal
      df[n,2] <- "Kruskal-Wallis"
      for (k in 1:length(mut.tab[1,])) { #Run column-wise
        cc <- na.omit(cbind(var.tab[,n], mut.tab[,m2])) #To check that at least levels are kept after NA correction
        if (nlevels(as.factor(cc[,2])) >= 2L) { #Run only if at least one mut is present after correction
          KW <- kruskal.test(var.tab[,n], mut.tab[,m2])
          df[n,m] <- KW$p.value
        }
        m <- m + 1
        m2 <- m2 + 1
      }
    }
    if (df[n,1] == FALSE) { #Categorical
      df[n,2] <- "Fisher's exact test"
      for (k in 1:length(mut.tab[1,])) {
        cc <- na.omit(cbind(var.tab[,n], mut.tab[,m2]))
        if (nlevels(as.factor(cc[,2])) >= 2L) {
          if (fisher.method == "exact") {
            #print(paste("char", n))
            #print(paste("mutation", m2))
            fet <- fisher.test(table(var.tab[,n], mut.tab[,m2]), workspace = 2e5)
            df[n,m] <- fet$p.value
          }

          if (fisher.method == "sim") {
            fet <- fisher.test(table(var.tab[,n], mut.tab[,m2]),simulate.p.value = T, B = 20000)
            df[n,m] <- fet$p.value
          }

          if (fisher.method == "hybrid") {
            fet <- fisher.test(table(var.tab[,n], mut.tab[,m2]), hybrid = T, workspace = 2e8)
            df[n,m] <- fet$p.value
          }
        }
        m <- m + 1
        m2 <- m2 + 1
      }
    }
    m <- 3
    m2 <- 1
    n <- n + 1
  }
  #Finally we remove annoying constant NA columns in the results
  df <- df[, !apply(is.na(df), 2, all)]
  if (is.null(correction) == FALSE) {
    ps <- as.vector(as.matrix(df[,3:ncol(df)]))
    pcor <- p.adjust(ps, method = correction)
    pcor <- matrix(pcor, nrow = nrow(df), ncol = ncol(df)-2)
    df[,3:ncol(df)] <- pcor
  }
  if (control == FALSE){
    df <- df[,-c(1:2)]
  }

  return(df)
}
