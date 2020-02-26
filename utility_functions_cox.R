#' #'Standardize Data
#'
#' @description function that standardize data
#' @param x Design matrix of dimension \code{n x q}, where \code{n} is the
#'   number of subjects and q is the total number of variables; each row is an
#'   observation vector. This must include all main effects and interactions as
#'   well, with column names corresponding to the names of the main effects
#'   (e.g. \code{x1, x2, E}) and their interactions (e.g. \code{x1:E, x2:E}).
#' @param center Should \code{x} be centered. Default is
#'   \code{TRUE}
#' @param normalize Should \code{x} be scaled to have unit variance. Default is
#'   \code{TRUE}


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
#'   interaction estimates. 
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
#'
#'   This function is used because the fitting algorithm estimates the gammas,
#'   and furthermore, the L1 penalty is placed on the gammas. It is used only in
#'   the initialization step of the algorithm


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


#' Soft Thresholding Function
#' @param x vector of the x_tilda for beta
#' @param y a Survival dataset containing either 2 columns with one denotes possible 
#'   censoning time and one for censored status
#' @param beta vector of regression coefficients to be thresholded
#' @param lambda tuning parameters
#' @param weight vector of weights for each beta
#' @param offset vector of offset for each beta
#' @return matrix of thresholded regression coefficients
#' @note user must supply \code{x} AND \code{y}, or \code{beta}, but not both. I
#'   set it up this way because to get the sequence of lambdas, I use the
#'   \code{beta} argument so that I only compute this once. I use the \code{x,
#'   y} argument for the CV folds. \code{lambda} can be a vector and this
#'   functions will return each thresholded beta for each lambda

soft <- function(x, y, beta, lambda, weight, offset) {
  
  if (missing(x) & missing(y) & missing(beta)) stop("user must supply x AND y, or beta but not both")
  if (missing(x) & missing(y)) return(list("beta" = sign(beta) * pmax(0, abs(beta) - lambda * weight)))
  if (missing(beta)) {
    time = y[,1]
    status = y[,2]
    
    
    beta <- coef(coxph.fit(x = x, y = Surv(time,status), strata = NULL, 
                           control = coxph.control(), 
                           method = "efron", 
                           rownames = NULL,
                           offset = offset))[1]
    
    b_lasso <- sign(beta)* pmax(0, abs(beta) - lambda*weight)
    
    # need to return a matrix, because this is used in the step to
    # calculate y_tilde in the shim function
    return(matrix(b_lasso, ncol = 1))
  }
  
}



