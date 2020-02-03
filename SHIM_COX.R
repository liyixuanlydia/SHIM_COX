
shim_cox <- function(x, y, main.effect.names, interaction.names,
                   lambda.beta, lambda.gamma,
                   threshold, max.iter,
                   center, normalize) {

   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   # standardize step 1 of algortihm
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   # only standardize x matrix
   obj <- standardize(x = x, center = center, normalize = normalize)
   x <- obj$x
   # x means
   bx <- obj$bx
   # x standard deviation
   sx <- obj$sx


   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   # Initialization step 2 of algortihm
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   # get Initial value for betas and alphas
   betas_and_alphas <- as.matrix(
   coef(glmnet::glmnet(x = x, y = y, 
                       family = "cox",
                       alpha = 0,
                       standardize = F))[-1,,drop = F])
   
   ## get adaptive weights
      # create output matrix
      adaptive.weights <- matrix(nrow = nrow(betas.and.alphas)) %>%
      magrittr::set_rownames(rownames(betas.and.alphas))
      # main effects weights
      for (j in main.effect.names) {
        adaptive.weights[j,] <- abs(1/betas.and.alphas[j,])
           }
      for (k in interaction.names) {
      # get names of main effects corresponding to interaction
      main <- strsplit(rownames(betas.and.alphas[k, , drop = F]),":")[[1]]
      adaptive.weights[k,] <- abs(prod(betas.and.alphas[main,])/betas.and.alphas[k,])
           }
  
   lambda_gamma <- lambda.gamma
   lambda_beta <- lambda.beta
   
   # this converts the alphas to gammas 
   uni_start <- convert(betas_and_alphas, main.effect.names = main.effect.names,
                        interaction.names = interaction.names)

   beta_hat_previous <- uni_start[main.effect.names, , drop = F]

   gamma_hat_previous <- uni_start[interaction.names, , drop = F]

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


# Iteration Begin 
  while (converged == 0 && m < max.iter){
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # update gamma (interaction parameter) step 3 of algortihm
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        standardize = F, intercept = F))[-1,,drop = F])
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # update beta (main effect parameter) step 4 of algortihm
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    beta_hat_next <- beta_hat_previous
    
    for (j in main.effect.names) {
      
      #j="x1"
      # determine the main effects not in j
      j_prime_not_in_j <- setdiff(main.effect.names,j)
      
      # for (notconverged in not_converged) {
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
      
      beta_hat_next_j <- as.matrix(coef(glmnet::glmnet(
          x = x_tilde_2,
          y = y,
          family = "cox", 
          offset = y_tilde_2,
          lambda = lambda_beta,
          penalty.factor = adaptive.weights[j,,drop=F],
          standardize = F, intercept = F))[-1,,drop = F])
  
      
      beta_hat_next[j,] <- beta_hat_next_j
      
    }
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # solution test step 5 of algortihm
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    delta <- abs(beta_hat_previous - beta_hat_next) + abs(gamma_hat_next - gamma_hat_previous)
    converged <- as.numeric(delta <= threshold)
    
    m <- m + 1
    
    beta_hat_previous <- beta_hat_next
    gamma_hat_previous <- gamma_hat_next
  }
    # Convert Gamma to Alpha
    Betas_and_Alphas <- convert2(beta_hat_next, gamma_hat_next,
                                 main.effect.names = main.effect.names,
                                 interaction.names = interaction.names)
    out <- list(
                beta = beta_hat_next,
                alpha = Betas_and_Alphas[interaction.names,],
                gamma = gamma_hat_next,
                converged = converged, 
                iteration.num = m,
                interaction.names = interaction.names,
                main.effect.names = main.effect.names)
    return(out)
}












#' Standardize Data

standardize <- function(x, center = TRUE, normalize = TRUE) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  if (center) {
    bx <- colMeans(x)
    x <- scale(x,bx,FALSE)
  } else {
    bx <- rep(0,p)
  }
  if (normalize) {
    sx <- sqrt(colSums(x ^ 2) / n)
    x <- scale(x,FALSE,sx)
  } else {
    sx <- rep(1,p)
  }
  
  return(list(x = x, bx = bx, sx = sx))
}

#' Convert alphas to gammas
#'
#' @description function that takes a vector of betas (which are the main
#'   effects) and alphas (which are the interaction effects) and converts the
#'   alphas to gammas.
#' @param betas.and.alphas q x 1 data.frame or matrix of main effects and
#'   interaction estimates. For example the output from the \code{uni_fun}
#'   function. The rownames must be appropriately labelled because these labels
#'   will be used in other functions
#' @param main.effect.names character vector of main effects names. MUST be
#'   ordered in the same way as the column names of \code{x}. e.g. if the column
#'   names of \code{x} are \code{\"x1\",\"x2\"} then \code{main.effect.names =
#'   c("x1","x2")}
#' @param interaction.names character vector of interaction names. MUST be
#'   separated by a colon (e.g. x1:x2), AND MUST be
#'   ordered in the same way as the column names of \code{x}
#' @param epsilon threshold to avoid division by a very small beta e.g. if any
#'   of the main effects are less than epsilon, set gamma to zero. This should
#'   not really be an important parameter because this function is only used in
#'   the initialization step, where the intial estimates are from OLS or ridge
#'   regression and therefor should not be very close to 0
#' @details note that \deqn{y = \beta_0 + \beta_1 x_1 + ... + \beta_p x_p +
#'   \alpha_{12} x_1 x_2 + ... + \alpha_{p-1,p} x_p x_{p-1} } and
#'   \deqn{\alpha_{ij} = \gamma_{ij} * \beta_i*\beta_j , i < j}
#'
#'   This function is used because the fitting algorithm estimates the gammas,
#'   and furthermore, the L1 penalty is placed on the gammas. It is used only in
#'   the initialization step in the \code{\link{shim}} function
#'
#' @seealso \code{\link{shim}}, \code{\link{Q_theta}}
#' @return a labelled q x 1 data.frame of betas and gammas

convert <- function(betas.and.alphas, main.effect.names, interaction.names,
                    epsilon = 1e-5) {
  
  betas_and_gammas <- matrix(nrow = nrow(betas.and.alphas)) %>%
    magrittr::set_rownames(rownames(betas.and.alphas))
  
  for (k in interaction.names) {
    
    # get names of main effects corresponding to interaction
    main <- strsplit(rownames(betas.and.alphas[k, , drop = F]),":")[[1]]
    
    # convert alpha to gamma BUT NEED TO CHECK IF BETAS ARE 0
    betas_and_gammas[k,] <- if (any(abs(betas.and.alphas[main,]) < epsilon )) 0 else
      betas.and.alphas[k,]/prod(betas.and.alphas[main,])
  }
  
  # add back the main effects which dont need to be transformed
  for (j in main.effect.names) {
    betas_and_gammas[j,] <- betas.and.alphas[j,]
  }
  
  return(betas_and_gammas)
}




#' Calculate working X's to update Gammas.
#'
#' @description function used to calculate working X's (xtilde) in step 3 of
#'   algorithm
#' @param interaction.names character vector of interaction names. must be
#'   separated by a ':' (e.g. x1:x2)
#' @param data.main.effects data frame or matrix containing the main effects
#'   data
#' @param beta.main.effects data frame or matrix containing the coefficients of
#'   main effects
#' @param nlambda number of tuning parameters
#' @return matrix of working X's (xtilde)

xtilde <- function(interaction.names, data.main.effects, beta.main.effects){
  
  # create output matrix
  xtildas <- matrix(ncol = length(interaction.names),
                    nrow = nrow(data.main.effects))
  colnames(xtildas) <- interaction.names
  
  for (k in interaction.names) {
    
    # get names of main effects corresponding to interaction
    main <- strsplit(k,":")[[1]]
    
    # step 3 to calculate x tilda
    xtildas[,k] <- prod(beta.main.effects[main,]) *
      data.main.effects[,main[1],drop = F] *
      data.main.effects[,main[2],drop = F]
  }
  
  return(xtildas)
}

#' Calculate working X's to update Betas
#'
#' @description function used to calculate working X's (xtilde) in step 4 of
#'   algorithm
#' @param interaction.names character vector of interaction names. must be
#'   separated by a ':' (e.g. x1:x2)
#' @param data.main.effects data frame or matrix containing the main effects
#'   data
#' @param beta.main.effects data frame or matrix containing the coefficients of
#'   main effects
#' @param gamma.interaction.effects data frame or matrix containing the gamma
#'   parameters
#' @return matrix of working X's (xtilde) of dimension n x (p*(p-1)/2)
#' @note this function is a modified x_tilde for step 4 because we thought maybe
#'   there was a typo. Math and results suggests that there is a typo in the
#'   original paper.

xtilde_mod <- function(interaction.names, data.main.effects, beta.main.effects,
                       gamma.interaction.effects){
  
  # create output matrix. no pipe is faster
  xtildas <- matrix(ncol = length(interaction.names),
                    nrow = nrow(data.main.effects))
  colnames(xtildas) <- interaction.names
  
  for (k in interaction.names) {
    
    # get names of main effects corresponding to interaction
    main <- strsplit(k, ":")[[1]]
    
    # step 4 to calculate x tilda
    xtildas[,k] <- prod(beta.main.effects[main,]) * gamma.interaction.effects[k,] *
      data.main.effects[,main[1],drop = F] *
      data.main.effects[,main[2],drop = F]
  }
  
  return(xtildas)
}


#' check_col_0 check number of zero cols
check_col_0 <- function(M) {
  M[, colSums(abs(M)) != 0, drop = F]
}

#' Convert gammas to alphas
#'
#' @description function that takes a vector of betas (which are the main
#'   effects) and gammas and converts the alphas to gammas. This function is
#'   used to calculate the linear predictor of the likelihood function (the Q
#'   function in the fitting algorithm)
#' @param betas.and.gammas q x 1 data.frame or matrix of betas and gamma
#'   estimates. For example the output from the \code{convert} function. The
#'   rownames must be appropriately labelled because these labels will be used
#'   in other functions
#' @inheritParams convert
#' @return a labelled q x 1 data.frame of betas and alphas

convert2 <- function(beta, gamma, main.effect.names, interaction.names,
                     intercept = NULL) {
  
  betas.and.gammas <- rbind2(beta,gamma)
  
  # create output matrix
  betas.and.alphas <- matrix(nrow = nrow(betas.and.gammas)) %>%
    magrittr::set_rownames(rownames(betas.and.gammas))
  
  for (k in interaction.names) {
    
    # get names of main effects corresponding to interaction
    main <- strsplit(rownames(betas.and.gammas[k, , drop = F]), ":")[[1]]
    
    # convert gamma to alpha
    betas.and.alphas[k,] <- betas.and.gammas[k,]*prod(betas.and.gammas[main,])
  }
  
  # add back the main effects which dont need to be transformed
  for (j in main.effect.names) {
    betas.and.alphas[j,] <- betas.and.gammas[j,]
  }
  
  # add back intercept if it is non-NULL
  if (!is.null(intercept)) betas.and.alphas["(Intercept)",] <- intercept
  
  return(betas.and.alphas)
  
}
