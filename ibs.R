#----------------------------------------------------------------------#
#         A function to compute IBS ... sbrier fxn in the ipred package (v.‘0.9.14’) 
#         with small change; class(pred) --> class(pred)[1]   #
#----------------------------------------------------------------------#

my_sbrier <- function (obj, pred, btime = range(obj[, 1])) 
{
  if (!inherits(obj, "Surv")) 
    stop("obj is not of class Surv")
  class(obj) <- NULL
  if (attr(obj, "type") != "right") 
    stop("only right-censoring allowed")
  N <- nrow(obj)
  time <- obj[, 1]
  ot <- order(time)
  cens <- obj[ot, 2]
  time <- time[ot]
  if (is.null(btime)) 
    stop("btime not given")
  if (length(btime) < 1) 
    stop("btime not given")
  if (length(btime) == 2) {
    if (btime[1] < min(time)) 
      warning("btime[1] is smaller than min(time)")
    if (btime[2] > max(time)) 
      warning("btime[2] is larger than max(time)")
    btime <- time[time >= btime[1] & time <= btime[2]]
  }
  ptype <- class(pred)[1] # changed this only!
  if (is.null(ptype)) {
    if (is.vector(pred)) 
      ptype <- "vector"
    if (is.list(pred)) 
      ptype <- "list"
  }
  if (ptype == "numeric" && is.vector(pred)) 
    ptype <- "vector"
  survs <- NULL
  switch(ptype, survfit = {
    survs <- getsurv(pred, btime)
    survs <- matrix(rep(survs, N), nrow = length(btime))
  }, list = {
    if (!inherits(pred[[1]], "survfit")) stop("pred is not a list of survfit objects")
    if (length(pred) != N) stop("pred must be of length(time)")
    pred <- pred[ot]
    survs <- matrix(unlist(lapply(pred, getsurv, times = btime)), 
                    nrow = length(btime), ncol = N)
  }, vector = {
    if (length(pred) != N) stop("pred must be of length(time)")
    if (length(btime) != 1) stop("cannot compute integrated Brier score with pred")
    survs <- pred[ot]
  }, matrix = {
    if (all(dim(pred) == c(length(btime), N))) survs <- pred[, 
                                                             ot] else stop("wrong dimensions of pred")
  })
  if (is.null(survs)) 
    stop("unknown type of pred")
  hatcdist <- prodlim(Surv(time, cens) ~ 1, reverse = TRUE)
  csurv <- predict(hatcdist, times = time, type = "surv")
  csurv[csurv == 0] <- Inf
  csurv_btime <- predict(hatcdist, times = btime, type = "surv")
  csurv_btime[is.na(csurv_btime)] <- min(csurv_btime, na.rm = TRUE)
  csurv_btime[csurv_btime == 0] <- Inf
  bsc <- rep(0, length(btime))
  if (length(btime) > 1) {
    for (j in 1:length(btime)) {
      help1 <- as.integer(time <= btime[j] & cens == 1)
      help2 <- as.integer(time > btime[j])
      bsc[j] <- mean((0 - survs[j, ])^2 * help1 * (1/csurv) + 
                       (1 - survs[j, ])^2 * help2 * (1/csurv_btime[j]))
    }
    idx <- 2:length(btime)
    RET <- diff(btime) %*% ((bsc[idx - 1] + bsc[idx])/2)
    RET <- RET/diff(range(btime))
    names(RET) <- "integrated Brier score"
    attr(RET, "time") <- range(btime)
  }
  else {
    help1 <- as.integer(time <= btime & cens == 1)
    help2 <- as.integer(time > btime)
    cs <- predict(hatcdist, times = btime, type = "surv")
    if (cs == 0) 
      cs <- Inf
    RET <- mean((0 - survs)^2 * help1 * (1/csurv) + (1 - 
                                                       survs)^2 * help2 * (1/cs))
    names(RET) <- "Brier score"
    attr(RET, "time") <- btime
  }
  RET
}