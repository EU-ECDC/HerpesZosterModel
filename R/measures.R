#' Matthew's correlation coefficient
#' 
#' Matthew's correlation coefficient (\eqn{MCC}) is given by
#' \deqn{MCC = (TP \cdot TN - FP \cdot FN) / (\sqrt((TP + FP)(TP + FN)(TN + FP)(TN + FN)))}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' It is considered to be better than F1 and Rand index (accuracy) for binary classification
#' \insertCite{@see @Chicco2017 for details}{zostmod}. 
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [f1() rand()]
#' @keywords classif
#' @family classification scores
#' @references 
#' \insertRef{Chicco2017}{zostmod} 
#' @export
#' @examples 
#' mcc(tp = 45, fp = 15, fn = 25, tn = 15)
mcc <- function(tp, tn, fp, fn, ...){
  n <- tp + tn + fp + fn
  s <- (tp + fn) / n
  p <- (tp + fp) / n
  score <- (tp / n - s * p) / sqrt(p * s * (1 - s) * (1 - p))
  return(score)
}

#' F score
#' 
#' The F score is given by
#' \deqn{F = ((1 + \beta^2) \cdot TP) / ((1 + \beta^2) \cdot TP + \beta^2 \cdot FP + FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, \eqn{FN} denotes false negatives, and
#' \eqn{\beta} is a weighting parameter. In our implementation \eqn{\beta}
#' has a default value of \eqn{1}, meaning our F score defaults to the F1 score 
#' (hence the function name).
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @param beta The weighting parameter, defaults to 1
#' @seealso [mcc()]
#' @keywords classif
#' @family classification scores
#' @references 
#' \insertRef{Soerensen1948}{zostmod}
#' \insertRef{Tversky1977}{zostmod}
#' @export
#' @examples 
#' f1(tp = 45, fp = 15, fn = 25, tn = 15)
f1 <- function(tp, tn, fp, fn, beta = 1, ...){
  score <- ((1 + beta ^ 2) * tp) / ((1 + beta ^ 2) * tp + beta ^ 2 * fp + fn)
  return(score)
}

#' Rand index/Accuracy
#' 
#' The Rand index \insertCite{Rand1971}{zostmod} (also known as accuracy) is given by
#' \deqn{R = (TP + TN) / (TP + TN + FP + FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#'
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [mcc(), err()]
#' @keywords classif
#' @family classification scores
#' @references
#' \insertRef{Rand1971}{zostmod}
#' @export
#' @examples 
#' rand(tp = 45, fp = 15, fn = 25, tn = 15)
#' @name rand
#' @aliases acc
rand <- function(tp, tn, fp, fn, ...){ 
  n <- tp + tn + fp + fn
  score <- (tp + tn) / n
  return(score)
}
acc <- rand

#' Optimsation precision
#' 
#' Optimisation precision is given by
#' \deqn{OP = R - |TP / (TP + FN) - TN / (TN + FP)| / (TP / (TP + FN) + TN / (TN + FP))}
#' where \eqn{R} denotes the Rand index, \eqn{TP} denotes true positives,
#' \eqn{TN} denotes true negatives, \eqn{FP} denotes false positives, and 
#' \eqn{FN} denotes false negatives.
#'
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [rand()]
#' @keywords classif
#' @family classification scores
#' @export
#' @examples 
#' op(tp = 45, fp = 15, fn = 25, tn = 15)
op <- function(tp, tn, fp, fn, ...){ 
  r <- rand(tp, tn, fp, fn)
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  score <- r - abs(sens - spec) / (sens + spec)
  return(score)
}

#' Error rate
#' 
#' The error rate is given by
#' \deqn{R = (FP + FN) / N}
#' where \eqn{N = TP + TN + FP + FN}, and \eqn{TP} denotes true positives,
#' \eqn{TN} denotes true negatives, \eqn{FP} denotes false positives, \eqn{FN} denotes false negatives.
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [rand()]
#' @keywords classif
#' @family classification scores
#' @export
#' @examples 
#' err(tp = 45, fp = 15, fn = 25, tn = 15)
err <- function(tp, tn, fp, fn, ...){
  n <- tp + tn + fp + fn
  score <- (fp + fn) / n
  return(score)
}

#' Cohen's kappa
#' 
#' Cohen's kappa (\eqn{\kappa}), also known as Heidke's skill score, is given by
#' \deqn{\kappa = 1 - (1 - p_o) / (1 - p_e)}
#' where \eqn{p_o} is the observed frequency and \eqn{p_e} is the 
#' expected frequency. We calculate \eqn{p_o} and \eqn{p_e} based on true positives
#' (\eqn{TP}), true negatives (\eqn{NP}), false positives (\eqn{FP}), and
#' false negatives (\eqn{FN}).
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [tss()]
#' @keywords classif
#' @family classification scores
#' @references
#' \insertRef{Cohen1960}{zostmod}
#' \insertRef{Heidke1926}{zostmod}
#' @export
#' @examples
#' kappa(tp = 45, fp = 15, fn = 25, tn = 15)
kappa <- function(tp, tn, fp, fn, ...){
  po <- rand(tp, tn, fp, fn)
  n <- tp + tn + fp + fn
  pe <-  ((tp + fp) / n) * ((tp + fn)/ n) +
    ((tn + fn) / n) * ((tn + fp) / n)
  score <- 1 - (1 - po) / (1 - pe)
  return(score)
}

#' True skill statistic
#' 
#' The true skill statistic (\eqn{TSS}) also known as Youden's index or markedness
#' is a measure often used in ecology for assessing species distribution models.
#' It is given by
#' \deqn{TSS = (TP / (TP + FN)) + (TN / (TN + FP)) - 1}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [kappa()]
#' @keywords classif
#' @family classification scores
#' @references
#' \insertRef{Youden1950}{zostmod}
#' @export
#' @examples 
#' tss(tp = 45, fp = 15, fn = 25, tn = 15)
tss <- function(tp, tn, fp, fn, ...){
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  score <- sens + spec - 1
  return(score)
}
mark <- tss

#' Jaccard index
#' 
#' The Jaccard index (also known as Tanimoto similarity) is given by
#' \deqn{Jaccard = TP / (TP + FP + FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#'  
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @keywords classif
#' @family classification scores
#' @references 
#' \insertRef{Jaccard1908}{zostmod}
#' \insertRef{Tversky1977}{zostmod}
#' @export
#' @examples 
#' jacc(tp = 45, fp = 15, fn = 25, tn = 15)
jacc <- function(tp, tn, fp, fn, ...){
  score <- tp / (tp + fp + fn)
  return(score)
}

#' Fowlkes-Mallows index
#' 
#' Fowlkes-Mallows' index is given by
#' \deqn{FM = \sqrt((TP / (TP + FP)) \cdot (TP / (TP + FN)))}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @keywords classif
#' @family classification scores
#' @family classification scores
#' @references 
#' \insertRef{FowlkesMallows1983}{zostmod}
#' @export
#' @examples
#' fm(tp = 45, fp = 15, fn = 25, tn = 15)
fm <- function(tp, tn, fp, fn, ...){
  score <- sqrt((tp / (tp + fp)) * (tp / (tp + fn)))
  return(score)
}

#' Positive likelihood ratio
#' 
#' The positive likelhihood ratio is given by
#' \deqn{LR_+ = (TP / (TP + FN)) / (1 - (TN / (TN + FP)))}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [dp() dor() gain()]
#' @keywords classif
#' @family classification scores
#' @export
#' @examples
#' lr(tp = 45, fp = 15, fn = 25, tn = 15)
lr <- function(tp, tn, fp, fn, ...){
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  score <- sens / (1 - spec)
  return(score)
}

#' Diagnostic odds ratio
#' 
#' The diagnostic odds ratio is given by
#' \deqn{DOR_+ = (TP \cdot TN) / (FP \cdot FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [lr()]
#' @keywords classif
#' @family classification scores
#' @export
#' @examples
#' dor(tp = 45, fp = 15, fn = 25, tn = 15)
dor <- function(tp, tn, fp, fn, ...){
  score <- (tp * tn) / (fp * fn)
  return(score)
}

#' Kvamme's gain statistic
#' 
#' Kvamme's gain statistic is given by
#' \deqn{gain = 1 - ((1 - (TN / (TN + FP)) / (TP / (TP + FN))}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [dp() lr()]
#' @keywords classif
#' @family classification scores
#' @references
#' \insertRef{Kvamme1988}{zostmod}
#' @export
#' @examples 
#' gain(tp = 45, fp = 15, fn = 25, tn = 15)
gain <- function(tp, tn, fp, fn, ...){
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  score <- 1 - (1 / lr(tp, tn, fp, fn))
  return(score)
}

#' Discriminant power
#' 
#' Discriminant power is given by
#' \deqn{DP = \sqrt(3) / pi \cdot (log((TP (FP + TN)) / (FP (FN + TP))) + 
#' log((TN (FN + TP)) / (FN (FP + TN))))}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [gain() lr()]
#' @keywords classif
#' @family classification scores
#' @export
#' @examples 
#' dp(tp = 45, fp = 15, fn = 25, tn = 15)
dp <- function(tp, tn, fp, fn, ...){
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  score <- (sqrt(3) / pi) * (log(sens / (1 - spec)) + log(spec / (1 - sens)))
  return(score)
}

#' Balanced classification rate
#' 
#' Balanced classification rate is given by
#' \deqn{BCR = (TP / (TP + FN) + TN / (TN + FP)) / 2}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @keywords classif
#' @family classification scores
#' @export
#' @examples 
#' bcr(tp = 45, fp = 15, fn = 25, tn = 15)
bcr <- function(tp, tn, fp, fn, ...){
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  score <- (sens + spec) / 2
  return(score)
}

#' Geometric mean
#' 
#' The geometric mean is given by
#' \deqn{GM = \sqrt(TP / (TP + FN) \cdot TN / (TN + FP))}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' The adjusted geometric mean is given by
#' \deqn{AGM = (GM + (TN / (TN + FP)) + (FP + TN)) / (1 + FP + TN) if TP / (TP + FN) > 0}
#' and \eqn{0} otherwise, where \eqn{TP} again denotes true positives, \eqn{TN} 
#' true negatives, \eqn{FP} false positives, and \eqn{FN} false negatives.
#'  
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @param adj If TRUE calculate the adjusted mean, defaults to FALSE
#' @keywords classif
#' @family classification scores
#' @export
#' @examples
#' gm(tp = 45, fp = 15, fn = 25, tn = 15)
gm <- function(tp, tn, fp, fn, adj = FALSE, ...){
  if(!is.logical(adj))
    stop("adj must be TRUE or FALSE")
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  m <- sqrt(sens * spec)
  if(adj == TRUE)
    score <- m
  else
    score <- ifelse(sens > 0, (m + spec * (fp + tn)) / (1 - fp + tn),0)
  return(score)
}

#' Brier score
#' 
#' The Brier score originates from weather forecasting. It is given by
#' \deqn{Brier = 1 / N \sum_{i = 1}^N (p_i - o_i)^2}
#' and is the mean squared error between the predicted probability, \eqn{p},
#' and the actual occurance, \eqn{o}.
#' 
#' @param rev whether to return 1 - Brier (TRUE) or not, defaults to FALSE
#' @family classification scores
#' @references
#' \insertRef{Brier1971}{zostmod}
#' @export
brier <- function(rev = FALSE, ...){
  if(!is.logical(rev))
    stop("rev must be TRUE or FALSE")
  score <- "TO BE ADDED"
  #if(adj == TRUE)
  #  score <- 1 - score
  return(score)
}

#' Area under ROC curve
#' 
#' The area under the receiver operating characteristic (ROC) curve for a
#' 2 x 2 confusion matrix is given by
#' \deqn{AUC = ((TP / (TP + FN) - (TN / (TN + FP)) + 1) / 2}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @keywords classif
#' @family classification scores
#' @references
#' \insertRef{Powers2007}{zostmod}
#' @export
#' @examples
#' auc(tp = 45, fp = 15, fn = 25, tn = 15)
auc <- function(tp, tn, fp, fn, ...){
  sens <- tp / (tp + fn) #TPR
  spec <- tn / (tn + fp)
  score <- (sens - spec + 1) / 2
  return(score)
}

#' Ignorance score
#' 
#' The ignorance score, also referred to as the logarithmic score
#' 
#' @keywords classif
#' @family classification scores
#' @references
#' \insertRef{Good1952}{zostmod}
#' \insertRef{RoulstonSmith2002}{zostmod}
#' @export
ignr <- function(...){
  score <- "TO BE ADDED"
  return(score)
}

#' Equitable threat score
#' 
#' The equitable threat score is given by
#' \deqn{ET = (TP - (TP + FP)(TP + FN) / N) / (TP - (TP + FP)(TP + FN) / N + FP + FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives,
#' and \eqn{N} is the sum of them.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @keywords classif
#' @family classification scores
#' @references 
#' \insertRef{Gilbert1884}{zostmod}
#' @export
#' @examples 
#' et(tp = 45, fp = 15, fn = 25, tn = 15)
et <- function(tp, tn, fp, fn, ...){
  n <- tp + tn + fp + fn
  score <- (tp - ((tp + fp) * (tp + fn)) / n) / (tp - ((tp + fp) * (tp + fn)) / n + fp + fn)
  return(score)
}

#' Gower and Legendre
#' 
#' The Gower and Legandre coefficient is given by
#' \deqn{GL = (TP + TN) / (TP + TN + \theta (FP + FN))}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, \eqn{FN} denotes false negatives, and
#' \eqn{\theta} is a weighting parameter. In our implementation \eqn{\theta}
#' has a default value of \eqn{1}, making it simple matching.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [sm() ss() rt()]
#' @param theta The weighting parameter, defaults to 1
#' @keywords classif
#' @family classification scores
#' @references 
#' \insertRef{GowerLegendre1986}{zostmod}
#' @export
#' @examples 
#' g1(tp = 45, fp = 15, fn = 25, tn = 15)
gl <- function(tp, tn, fp, fn, theta = 1, ...){
  t <- tp + tn
  f <- fp + fn
  score <- t / (t + theta * f)
  return(score)
}

#' Simple matching
#' 
#' The simple matching coefficient is given by
#' \deqn{SM = (TP + TN) / (TP + TN + FP + FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [gl() rt() ss()]
#' @keywords classif
#' @family classification scores
#' @references 
#' \insertRef{SokalMichener1958}{zostmod}
#' @export
#' @examples 
#' sm(tp = 45, fp = 15, fn = 25, tn = 15)
sm <- function(tp, tn, fp, fn, ...){
  t <- tp + tn
  n <- tp + tn + fp + fn
  score <- t / n
  return(score)
}

#' Rogers-Tanimoto
#' 
#' The Rogers-Tanimoto coefficient is given by
#' \deqn{RT = (TP + TN) / (TP + TN + 2FP + 2FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [gl() sm() ss()]
#' @keywords classif
#' @family classification scores
#' @export
#' @examples 
#' rt(tp = 45, fp = 15, fn = 25, tn = 15)
rt <- function(tp, tn, fp, fn, ...){
  t <- tp + tn
  f2 <- 2 * fp + 2 * fn
  score <- t / (t + f2)
  return(score)
}


#' Sokal and Sneath
#' 
#' The Sokal and Sneath coefficient is given by
#' \deqn{SS = (2TP + 2TN) / (2TP + 2TN + FP + FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [gl() sm() ss()]
#' @keywords classif
#' @family classification scores
#' @references 
#' \insertRef{SokalSneath1963}{zostmod}
#' @export
#' @examples 
#' ss(tp = 45, fp = 15, fn = 25, tn = 15)
ss <- function(tp, tn, fp, fn, ...){
  t2 <- 2 * tp + 2 * tn
  f <- fp + fn
  score <- t2 / (t2 + f)
  return(score)
}

#' Faith measure
#' 
#' The Faith measure is given by
#' \deqn{(TP - FP - FN) / (TP + TN + FP + FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [rus()]
#' @keywords classif
#' @family classification scores
#' @references 
#' \insertRef{Faith1983}{zostmod}
#' @export
#' @examples 
#' fai(tp = 45, fp = 15, fn = 25, tn = 15)
fai <- function(tp, tn, fp, fn, ...){
  n <- tp + tn + fp + fn
  score <- (tp - fp - fn) / n
  return(score)
}

#' Russell-Rao
#' 
#' The Russell-Rao coefficient is given by
#' \deqn{TP / (TP + TN + FP + FN)}
#' where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
#' \eqn{FP} denotes false positives, \eqn{FN} denotes false negatives.
#' 
#' @param tp True positives
#' @param tn True negatives
#' @param fp False positives
#' @param fn False negatives
#' @seealso [fai()]
#' @keywords classif
#' @family classification scores
#' @export
#' @examples 
#' rus(tp = 45, fp = 15, fn = 25, tn = 15)
rus <- function(tp, tn, fp, fn, ...){
  n <- tp + tn + fp + fn
  score <- tp / n
  return(score)
}

#' Canberra distance
#'
#' Calculates the Canberra distance of two vectors x and y.
#'
#' @param x, y
#' @seealso [braycurtis()]
#' @examples
#' x <- c(1, 2, -3)
#' y <- c(6, 4, 5)
#' canberra(x, y)
#' @export
#'
canberra_dist <- function(x, y){
  sum(abs(x - y) / (abs(x) + abs(y)))  
}

#' Bray-Curtis dissimilarity
#'
#' Calculates the Bray-Curtis dissimilarity of two vectors x and y.
#'
#' @param x, y
#' @seealso [canberra()]
#' @examples
#' x <- c()
#' y <- c()
#' braycurtis(x, y)
#' @export
#'
braycurtis <- function(x, y){
  if(any(x + y == 0))
    stop("Division by 0")
  sum(abs(x - y) / (x + y))
}