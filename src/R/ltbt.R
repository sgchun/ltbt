library(pROC)

#' likelihood of an individual being a case given being a carrier 
#' 
#' @param gamma.seq assumed rare-variant effect sizes (gamma)
#' 
#' @param prs polygenic risk score of carrier
#' 
#' @param b adjustment for effect size from the baseline gamma (1 if no adjustment neeed)
#' 
#' @param Kprev.x prevalence of case in population
#' 
#' @param h2.prs.x heritability explained by PRS (h2L[x])
#' 
#' @return log likelihood values
#'
lik.case <- function (gamma.seq, prs, b, Kprev.x, h2.prs.x) {
    tau.std <- qnorm(1 - Kprev.x, mean=0, sd=1)

    gamma <- gamma.seq * b    

    tau.mut.std <- tau.std - gamma
    
    i.mut <- dnorm(tau.mut.std, mean=0, sd=1) /
        pnorm(tau.mut.std, mean=0, sd=1, lower=FALSE)
    
    prs.mean.case.exomepos <- i.mut * h2.prs.x

    prs.sd.case.exomepos <- 
        sqrt(h2.prs.x * (1 - h2.prs.x * i.mut * (i.mut - tau.mut.std)))

    lik <-
        dnorm(prs, 
              mean=prs.mean.case.exomepos,
              sd=prs.sd.case.exomepos, log=TRUE) 

    lik
}

#' likelihood of being a carrier individual being a case 
#' 
#' @param gamma.seq assumed rare-variant effect sizes (gamma)
#' 
#' @param genofreq frequency of carrier in populations
#' 
#' @param b adjustment for effect size from the baseline gamma (1 if no adjustment neeed)
#' 
#' @param Kprev.x prevalence of case in population
#' 
#' @return likelihood values
#'
case.x <- function (gamma.seq, genofreq, b, Kprev.x) {
    tau.std <- qnorm(1 - Kprev, mean=0, sd=1)

    tau.mut.std <- tau.std - gamma.seq * b

    x <- pnorm(tau.mut.std, mean=0, sd=1, lower=FALSE) 
    
    x * genofreq
}

#' likelihood of an individual being a control given being a carrier
#'
#' @param gamma.seq assumed rare-variant effect sizes (gamma)
#'
#' @param prs polygenic risk score of carrier
#'
#' @param b adjustment for effect size from the baseline gamma (1 if no adjustment neeed)
#'
#' @param Kprev.x prevalence of case in population
#'
#' @param h2.prs.x heritability explained by PRS (h2L[x])
#'
#' @return log likelihood values
#'
lik.ctrl <- function (gamma.seq, prs, b, Kprev.x, h2.prs.x) {
    tau.std <- qnorm(1 - Kprev.x, mean=0, sd=1)

    gamma <- gamma.seq * b    

    tau.mut.std <- tau.std - gamma

    v.mut <- -dnorm(tau.mut.std, mean=0, sd=1) /
        (1 - pnorm(tau.mut.std, mean=0, sd=1, lower=FALSE))

    prs.mean.ctrl.exomepos <- v.mut * h2.prs.x

    prs.sd.ctrl.exomepos <- 
        sqrt(h2.prs.x * (1 - h2.prs.x * v.mut * (v.mut - tau.mut.std))) 
        
    lik <-
        dnorm(prs, 
              mean=prs.mean.ctrl.exomepos,
              sd=prs.sd.ctrl.exomepos, log=TRUE)

    lik
}

#' likelihood of being a carrier individual being a control 
#'
#' @param gamma.seq assumed rare-variant effect sizes (gamma)
#'
#' @param frequency of carriers in population
#'
#' @param b adjustment for effect size from the baseline gamma (1 if no adjustment neeed)
#'
#' @param Kprev.x prevalence of case in population
#'
#' @return likelihood values
#'
ctrl.x <- function (gamma.seq, genofreq, b, Kprev.x) {
    tau.std <- qnorm(1 - Kprev.x, mean=0, sd=1)

    tau.mut.std <- tau.std - gamma.seq * b

    x <- pnorm(tau.mut.std, mean=0, sd=1, lower=TRUE) 

    x * genofreq
}

#' Frequency of carriers of all functional rare variants. Maximum likelihood estimator was used to account for the case/control ascertainment
#' 
#' @param gamma.seq assumed rare-variant effect sizes (gamma)
#'
#' @param Kprev.x prevalence of case in population
#'
#' @param p0 ascertainment ratio for control, i.e., Pcontrol / (1 - K) 
#'
#' @param p1 ascertainment ratio for case, i.e., Pcase / K
#' 
#' @param n.mut number of carriers in case/control
#' 
#' @param n number of all subjects in case/control
#'
#' @return most likely estimate of frequency of carriers in population
#'
caf.mle <- function (gamma.seq, Kprev.x, p0, p1, n.mut, n) {
    tau.std <- qnorm(1 - Kprev.x, mean=0, sd=1)

    n.nomut <- n - n.mut    
    expr1 <- -(p0 - p1) * Kprev.x + p0
    expr2 <- (p0 - p1) *
        (1 - Kprev.x - pnorm(tau.std - gamma.seq, lower=TRUE))
    
    caf.hat <- n.mut * expr1 / (n * expr1 - n.nomut * expr2)
    caf.hat
}


#' Likelihood Ratio Test of liability threshold burden test (LTBT)
#' 
#' @param gamma.seq assumed rare-variant effect sizes (gamma)
#' 
#' @param gt genotype vector of cohort (0: non-carrier, 1: carrier)
#' 
#' @param prs polygenic risk score of carrier
#' 
#' @param outcome outcome vector of cohort (0: control, 1: case)
#' 
#' @param Kprev.x prevalence of case in population
#' 
#' @param h2.prs.x heritability explained by PRS (h2L[x]). Estimated from cohort if it is left as NULL.
#' 
#' @return LTBT test results
#' @export  
lrt.fe.2df <- function (gamma.seq, gt, prs, outcome, Kprev.x, h2.prs.x=NULL) {

    stopifnot("matrix" %in% class(gt))
    stopifnot(nrow(gt) == length(outcome))
    stopifnot(length(prs) == length(outcome))

    gamma0.idx <- which(gamma.seq == 0)
    
    case.idx <- which(outcome > 0)
    ctrl.idx <- which(outcome == 0)

    ## case/control ascertainment parameter
    case.p <- length(case.idx) / length(outcome)
    p1 <- case.p / Kprev.x
    p0 <- (1 - case.p) / (1 - Kprev.x)

    ## burden collapse by individual
    burden <- apply(gt, 1, sum)

    ## Estimate h2.prs
    if (is.null(h2.prs.x)) {
        ## Estimate AUC First
        auc.est <- auc(roc(outcome, prs))

        ## Wray et al 2010
        Q <- qnorm(auc.est)

        tau0 <- qnorm(1 - Kprev.x, mean=0, sd=1)
        
        i0 <- ( dnorm(tau0, mean=0, sd=1) /
                pnorm(tau0, mean=0, sd=1, lower=FALSE) )
        
        v0 <- ( -dnorm(tau0, mean=0, sd=1) /
                (1 - pnorm(tau0, mean=0, sd=1, lower=FALSE)) )

        h2.prs.x <- (
             (2 * Q ** 2) /
              ((v0 - i0)**2 + Q**2 * i0 * (i0 - tau0) + v0 * (v0 - tau0)) )
    }

    ## No rare variant
    if (all(burden == 0)) {
        return(list(gamma=gamma.seq,
                    llr=0,
                    pvalue=1, 
                    h2Lx=h2.prs.x,
                    gamma.mle=0,
                    caf.mle=0))
    }


    ## rescale so that var(prs) = h2.prs.x
    prs <- prs * sqrt(h2.prs.x)

    logL <- rep(0, length(gamma.seq))

    ## dominant
    carriers <- which(burden > 0)

    for (I in intersect(carriers, case.idx)) {
        logL <- logL +
            lik.case(gamma.seq, prs[I], 1, Kprev.x, h2.prs.x)
    }
    
    for (I in intersect(carriers, ctrl.idx)) {
        logL <- logL +
            lik.ctrl(gamma.seq, prs[I], 1, Kprev.x, h2.prs.x)
    }

    caf <- caf.mle(gamma.seq, Kprev.x, p0, p1, sum(burden > 0),
                   length(outcome))

    x.y1 <- case.x(gamma.seq, caf, 1, Kprev.x)

    logL <- logL + log(p1 * x.y1) * sum(case.idx %in% carriers) +
        log(p1 * (1 - caf) * Kprev.x) * sum(!(case.idx %in% carriers))
    
    x.y0 <- ctrl.x(gamma.seq, caf, 1, Kprev.x)

    logL <- logL + log(p0 * x.y0) * sum(ctrl.idx %in% carriers) +
        log(p0 * (1 - caf) * (1 - Kprev.x)) *
            sum(!(ctrl.idx %in% carriers))

    ## normalize by denominator
    logL <-
        logL - (log(p1 * x.y1 + p0 * x.y0 +
                    (1 - caf) * (p1 * Kprev.x + p0 * (1 - Kprev.x))) *
                length(outcome))

    maxL.idx <- which.max(logL)

    llr <- -2 * (logL[gamma0.idx] - max(logL, na.rm=TRUE))

    list(gamma=gamma.seq,
         llr=llr,
         pvalue=pchisq(llr, df=1, lower.tail=FALSE),
         h2Lx=h2.prs.x,
         gamma.mle=gamma.seq[maxL.idx],
         caf.mle=caf[maxL.idx],
         logL=logL)
}
