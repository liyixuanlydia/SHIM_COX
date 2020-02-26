#' Fit Strong Heredity Interaction Model For Cox Model
#'
#' @description function to fit the Strong Heredity Interaction Model for a
#'   sequence of tuning parameters. This is a penalized regression method that
#'   ensures the interaction term is non-zero only if its corresponding
#'   main-effects are non-zero.
#' @param x Design matrix of dimension \code{n x q}, where \code{n} is the
#'   number of subjects and q is the total number of variables; each row is an
#'   observation vector. This must include all main effects and interactions as
#'   well, with column names corresponding to the names of the main effects
#'   (e.g. \code{x1, x2, E}) and their interactions (e.g. \code{x1:E, x2:E}).
#'   All columns should be scaled to have mean 0 and variance 1; this is done
#'   internally by the \code{\link{shim}} function.
#' @param y a Survival dataset containing either 2 columns with one denotes possible 
#'   censoning time and one for censored status
#' @param main.effect.names character vector of main effects names. MUST be
#'   ordered in the same way as the column names of \code{x}. e.g. if the column
#'   names of \code{x} are \code{"x1","x2"} then \code{main.effect.names =
#'   c("x1","x2")}
#' @param interaction.names character vector of interaction names. MUST be
#'   separated by a colon (e.g. x1:x2), AND MUST be ordered in the same way as
#'   the column names of \code{x}
#' @param lambda sequence of tuning parameters for the penalty. 
#' @param alpha sequence of tuning parameters for the adaption different level 
#'   in selection between main terms and interaction terms
#' @param nlambda total number of tuning parameters. Now setting the lambda and 
#'   alpha has same length and used as pairs
#' @param threshold Convergence threshold for coordinate descent. Each
#'   coordinate-descent loop continues until the change in the objective
#'   function after all coefficient updates is less than threshold. Default
#'   value is \code{1e-4}.
#' @param max.iter Maximum number of passes over the data for all tuning
#'   parameter values; default is 100.
#' @param center Should \code{x}  be centered. Default is
#'   \code{TRUE}. Centering \code{y} applies to \code{family="gaussian"} only.
#' @param normalize Should \code{x} be scaled to have unit variance. Default is
#'   \code{TRUE}




shim_cox <- function(x, y, main.effect.names, interaction.names,
                     lambda, alpha, nlambda,
                     threshold = 1e-4, max.iter = 100,
                     center, normalize) {
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Get the list of the lambda and creat matrix to store results
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
  # Get the list of lambda_beta and lambda_gamma from lambda and alpha
  lambda_beta_list <- lambda * alpha
  ones <- rep(1,times=nlambda)
  lambda_gamma_list <- lambda * (ones- alpha)
  # Name all the tuning parameters
  tuning_params_mat <- matrix(c(lambda_gamma_list, lambda_beta_list),
                              nrow = 2, ncol = nlambda, byrow = T)
  dimnames(tuning_params_mat)[[1]] <- c("lambda.gamma","lambda.beta")
  dimnames(tuning_params_mat)[[2]] <- paste0("s",seq_len(nlambda))
  lambdaNames <- dimnames(tuning_params_mat)[[2]]
  
  # matrix to store results of betas and alphas on standardized scale
  betaMat <- matrix(nrow = length(main.effect.names), ncol = nlambda,
                    dimnames = list(c(main.effect.names),
                                    lambdaNames))
  
  gammaMat <- matrix(nrow = length(interaction.names), ncol = nlambda,
                     dimnames = list(c(interaction.names),
                                     lambdaNames))
  
  coefficientMat <- matrix(nrow = length(c(main.effect.names, interaction.names)),
                           ncol = nlambda,
                           dimnames = list(c(main.effect.names, interaction.names),
                                           lambdaNames))
  
  # matrix to store the outprint results
  outPrint <- matrix(NA, nrow = nlambda, ncol = 6,
                     dimnames = list(lambdaNames,
                                     c("dfBeta","dfAlpha","deviance",
                                       "percentDev",
                                       "lambdaBeta", "lambdaGamma")))
  
  # Other Initial setings for iterations
  
  # for all the x_tilde in zero_x_tilde, return the following matrix with 0 for each coefficient
  # this is like a place holder.
  # use in the iteration for gamma_next when all the x_tilde are 0
  coef_zero_gamma_matrix <- matrix(data = 0,
                                   nrow = length(interaction.names),
                                   ncol = 1,
                                   dimnames = list(interaction.names))
  
  # index data.frame to figure out which j < j'
  # use in the iteration for beta_next 
  index <- data.frame(main.effect.names, seq_along(main.effect.names),
                      stringsAsFactors = F)
  colnames(index) <- c("main.effect.names","index")
  
  
  m <- 1 # iteration counter
  delta <- 1 # threshold initialization
  converged <- 0
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # standardize  $step 1 of algorithm
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  # only standardize x matrix
  obj <- standardize(x = x, center = center, normalize = normalize)
  x <- obj$x
  # x means
  bx <- obj$bx
  # x standard deviation
  sx <- obj$sx
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Initialization  $step 2 of algorithm
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  # get Initial value for betas and alphas
  betas_and_alphas <- as.matrix(
    coef(glmnet::glmnet(x = x, y = y, 
                        family = "cox",
                        alpha = 0,
                        ##?
                        lambda = 0.005,
                        standardize = F))[,,drop = F])
  
  # Get Initial value for betas and gammas
  # Use function convert transfer value from alpha to gamma
  uni_start <- convert(betas_and_alphas, main.effect.names = main.effect.names,
                       interaction.names = interaction.names)
  
  beta_hat_previous <- uni_start[main.effect.names, , drop = F]
  
  gamma_hat_previous <- uni_start[interaction.names, , drop = F]
  
  # Get adaptive weights
  adaptive.weights <- betas_and_alphas
  for (j in main.effect.names) {
    adaptive.weights[j,] <- abs(1/betas_and_alphas[j,])
  }
  for (k in interaction.names) {
    # Get names of main effects corresponding to interaction and setting interaction weights
    main <- strsplit(rownames(betas_and_alphas[k, , drop = F]),":")[[1]]
    adaptive.weights[k,] <- abs(prod(betas_and_alphas[main,])/betas_and_alphas[k,])
  }
  
  
  for (LAMBDA in lambdaNames) {
    
    lambdaIndex <- which(LAMBDA==lambdaNames)
    lambda_beta <- tuning_params_mat["lambda.beta",LAMBDA][[1]]
    lambda_gamma <- tuning_params_mat["lambda.gamma",LAMBDA][[1]]
    
    
    while (converged == 0 && m < max.iter){
      ##### One iteration begins     
      
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # update gamma (interaction parameter)  $step 3 of algorithm
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      y_tilde <- x[,main.effect.names,drop = F] %*% beta_hat_previous
      
      x_tilde <- xtilde(interaction.names = interaction.names,
                        data.main.effects = x[,main.effect.names, drop = F],
                        beta.main.effects = beta_hat_previous)
      
      # indices of the x_tilde matrices that have all 0 columns
      zero_x_tilde <- is.null(colnames(check_col_0(x_tilde)))
      
      # this will store the results but will be shorter than nlambda
      gamma_hat_next <- if (zero_x_tilde) coef_zero_gamma_matrix else
        as.matrix(coef(glmnet::glmnet(
          x = x_tilde,
          y = y,
          family = "cox", 
          offset = y_tilde,
          lambda = lambda_gamma,
          penalty.factor = adaptive.weights[interaction.names,,drop=F],
          standardize = F))[,,drop = F])
      
      
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # update beta (main effect parameter)  $step 4 of algorithm
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      beta_hat_next <- beta_hat_previous
      
      for (j in main.effect.names) {
        
        # determine the main effects not in j
        j_prime_not_in_j <- setdiff(main.effect.names,j)
        
        y_tilde_2 <- x[,j_prime_not_in_j, drop = F] %*%
          beta_hat_next[j_prime_not_in_j, , drop = F] +
          rowSums(
            xtilde_mod(beta.main.effects = beta_hat_next[j_prime_not_in_j, , drop = F],
                       gamma.interaction.effects = gamma_hat_next,
                       interaction.names = interaction.names[-grep(paste0("\\b",j,"\\b"), interaction.names)],
                       data.main.effects = x[,j_prime_not_in_j, drop = F]))
        
        # j' less than j 
        j.prime.less <- index[which(index[,"index"] < index[which(index$main.effect.names == j),2]),
                              "main.effect.names"]
        
        # need to make sure paste(j.prime.less,j,sep=":") are variables in x matrix
        # this is to get around situations where there is only interactions with the E variable
        j.prime.less.interaction <- intersect(paste(j.prime.less,j, sep = ":"), colnames(x))
        
        # need to get the main effects that are in j.prime.greater.interaction
        j.prime.less <- gsub("\\:(.*)", "", j.prime.less.interaction)
        
        
        # the if conditions in term1 and term2 are to check if there are
        # any variables greater or less than j
        # lapply is faster than mclapply here
        term_1 <- if (length(j.prime.less.interaction) != 0 ) {
          x[,j.prime.less.interaction] %*%
            (gamma_hat_next[j.prime.less.interaction,, drop = F] *
               beta_hat_next[j.prime.less, , drop = F])} else 0
        
        # j' greater than j
        j.prime.greater <- index[which(index[,"index"] >
                                         index[which(index$main.effect.names == j),2]),
                                 "main.effect.names"]
        
        # need to make sure j.prime.greater is a variable in x matrix
        # this is to get around situations where there is only interactions with the E variable
        j.prime.greater.interaction <- intersect(paste(j,j.prime.greater,sep = ":"), colnames(x))
        
        # need to get the main effects that are in j.prime.greater.interaction
        j.prime.greater <- if (all(gsub("\\:(.*)", "", j.prime.greater.interaction) == j))
          gsub("(.*)\\:", "", j.prime.greater.interaction) else gsub("\\:(.*)", "", j.prime.greater.interaction)
        
        term_2 <- if (length(j.prime.greater) != 0) {
          x[,j.prime.greater.interaction] %*%
            (gamma_hat_next[j.prime.greater.interaction,, drop = F] *
               beta_hat_next[j.prime.greater,,drop = F]) } else 0
        
        x_tilde_2 <- x[,j, drop = F] +  term_1 + term_2
        
        # glmnet refuse to deal with 1 column x conditions so write soft function to calculate beta
        beta_hat_next_j <- soft(x = x_tilde_2,
                                y = y,
                                offset = y_tilde_2,
                                weight = adaptive.weights[j,,drop=F],
                                lambda = lambda_beta)
        
        beta_hat_next[j,] <- beta_hat_next_j
        
      }
      
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # solution test  $step 5 of algortihm
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      delta <- abs(beta_hat_previous - beta_hat_next) + abs(gamma_hat_next - gamma_hat_previous)
      converged <- as.numeric(delta <= threshold)
      
      m <- m + 1
      
      beta_hat_previous <- beta_hat_next
      gamma_hat_previous <- gamma_hat_next
      
      ##### One iteration ends
    }
    
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # store the results for one lambda & alpha
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # Convert Gamma to Alpha
    Betas_and_Alphas <- convert2(beta_hat_next, gamma_hat_next,
                                 main.effect.names = main.effect.names,
                                 interaction.names = interaction.names)
    
    betaMat[,lambdaIndex] <- beta_hat_next
    gammaMat[,lambdaIndex] <- gamma_hat_next
    coefficientMat[,lambdaIndex] <- Betas_and_Alphas
  }
  
  lambda.beta <- unlist(lambda_beta_list)
  names(lambda.beta) <- paste0("s",1:nlambda,".beta")
  lambda.gamma <- unlist(lambda_gamma_list)
  names(lambda.gamma) <- paste0("s",1:nlambda, ".gamma")
  
  out <- list(coefficient = coefficientMat,
              gamma = gammaMat,
              lambda.beta = lambda.beta,
              lambda.gamma = lambda.gamma,
              tuning.parameters = tuning_params_mat,
              converged = converged, x = x, y = y, bx = bx, by = by, sx = sx,
              center = center, normalize = normalize,
              nlambda = nlambda,
              interaction.names = interaction.names,
              main.effect.names = main.effect.names)
  
  class(out) <- "lspath"
  return(out)
}
