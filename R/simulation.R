#' Title helper function: betapar
#'
#' @param mu mean of VAF
#' @param var variance of VAF
#'
#' @return
#' @export
#'
#' @keywords internal
betapar <- function(mu, var) {
  if(var==0|var>=0.25){
    return(params = c(alpha = NA, beta = NA))
  }else{
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = c(alpha = alpha, beta = beta))
  }
}

#' Title helper function of simdat
#'
#' @param mus means of VAF
#' @param vars variances of VAF
#' @param n readdepth
#' @param N sample size
#'
#' @return
#' @export
#'
#' @keywords internal
simdat <- function(mus, vars, n, N){
  M <- length(mus)
  Ns <- rep(N, M)
  pars <- mapply(betapar, mus, vars, SIMPLIFY=T)


  ps <- data.frame(mapply(rbeta,
               n=Ns,
               shape1=pars["alpha",],
               shape2=pars["beta",]
               ))
  NAindex <- which(is.na(pars["alpha",]))
  if(length(NAindex)>0){
    for(i in 1:length(NAindex)){
      ps[,NAindex[i]] <- mus[NAindex[i]]
    }
  }
  ps <- reshape(ps, direction = "long", varying = list(colnames(ps)),
                v.names = "p",
                timevar = "group")
  ps$vcs <- rbinom(M*N, n, ps$p)
  ps$n <- n
  vcf <- ps[,c("p", "vcs", "n", "group")]
  names(vcf)[2:3] <- c("vc", "dp")
  vcf$group <- as.factor(vcf$group)
  return(vcf)
}

#' Title
#'
#' @param mus means of groups's VAF
#' @param vars variance of VAF for different groups
#' @param n read depth
#' @param N sample size
#' @param nsim number of simulations
#' @param alpha.level
#'
#' @return simulation result lists including "prob.reject.fisher", "prob.reject.pearson",  "prob.reject.betabinom", and so on
#' @export
#'
#' @examples power.vaf(c(0.2, 0.2), c(0.01, 0.01), 50, 3, 100, 0.05)
#'
power.vaf <- function(mus, vars, n, N, nsim, alpha.level){
  #defaultW <- getOption("warn")
  #options(warn = -1)
  rslt <- data.frame(matrix(NA,    # Create empty data frame
                            nrow = nsim,
                            ncol = 8))
  names(rslt) <- c("P.value.fisher",
                   "Reject.fisher",
                   "P.value.pearson",
                   "Reject.pearson",
                   "P.value.betabinom",
                   "Reject.betabinom",
                   "Converge.betabinom",
                   "Vars.betabinom"
                   )

  for(s in 1:nsim){
    temp <- simdat(mus, vars, n, N)
    temp %>%
      group_by(group) %>%
      summarise(dp = sum(dp), vc= sum(vc), vaf=sum(vc)/sum(dp)) -> vaf

      tab <- as.matrix(cbind(vaf$vc, vaf$dp-vaf$vc))
      fisherrslt <- fisher.test(tab, alternative = "two.sided")
      pvalue.fisher <- fisherrslt$p.value
      #test <- "Fisher's exact"
      reject.fisher <- ifelse(pvalue.fisher<alpha.level, 1, 0)

      chirslt <- prop.test(vaf$vc, vaf$dp, alternative = "two.sided")
      pvalue.pearson <- chirslt$p.value
      #test <- "Chi-squared"
      reject.pearson <- ifelse(pvalue.pearson<alpha.level, 1, 0)

      warn <- 0
      if(length(unique(vars))>1){
        rand <- as.formula("~ group")
      }else{
        rand <- as.formula("~ 1")
      }

      tryCatch({m1 <- betabin(cbind(vc, dp - vc) ~ group, rand, temp, hessian=F)
                m0 <- betabin(cbind(vc, dp - vc) ~ 1, rand, temp, hessian=F)},
               warning=function(w) warn <<- 1)
      #betarslt <- wald.test(b = coef(m1), Sigma = vcov(m1), Terms=2:length(levels(temp$group)))
      betarslt <- anova(m1,m0, test="LRT")
      pvalue.betabinom <- betarslt@anova.table[["P(> Chi2)"]][2]

      #test <- "Beta-Binomial LRT"
      reject.betabinom <- ifelse(pvalue.betabinom<alpha.level, 1, 0)
      converge.betabinom <- warn
      newdata = data.frame(group=levels(temp$group))
      newdata$group <- factor(newdata$group)
      X <- model.matrix(~group, newdata)
      ps <- inv.logit(X%*%coef(m1))
      phis <- m1@random.param
      a2 <- (1-phis)/((phis*ps)/(1-ps)+phis)
      a1 <- ps*a2/(1-ps)
      vars.m1 <- (a1*a2)/((a1+a2)^2*(a1+a2+1))
      vars.m1[which(vars >= 0.25)] <- NA

      newdata = data.frame(group=1)
      newdata$group <- factor(newdata$group)
      X <- model.matrix(~1, newdata)
      ps <- rep(inv.logit(X%*%coef(m0)), length(levels(temp$group)))
      phis <- m0@random.param
      a2 <- (1-phis)/((phis*ps)/(1-ps)+phis)
      a1 <- ps*a2/(1-ps)
      vars.m0 <- (a1*a2)/((a1+a2)^2*(a1+a2+1))
      vars.m0[which(vars >= 0.25)] <- NA
      vars.m0m1 <- round(c(vars.m0, vars.m1), digits = 2)
      vars.betabinom <- paste("(reduced:",vars.m0m1[1],", ", vars.m0m1[2],"); ",
                              "(full:",vars.m0m1[3],", ", vars.m0m1[4],")",
                              collapse="")

    rslt[s,] <- c(pvalue.fisher, reject.fisher,
                  pvalue.pearson, reject.pearson,
                  pvalue.betabinom,reject.betabinom,
                  converge.betabinom, vars.betabinom)
  }
for(i in 1:7){
  rslt[,i] <- as.numeric(rslt[,i])
}
N.fisher <- nsim-sum(is.na(rslt$Reject.fisher))
n.fisher <- sum(rslt$Reject.fisher, na.rm = T)
N.pearson <- nsim-sum(is.na(rslt$Reject.pearson))
n.pearson <- sum(rslt$Reject.pearson, na.rm = T)
betabinom.rslt <- subset(rslt, Converge.betabinom==0)
N.betabinom <- dim(betabinom.rslt)[1]-sum(is.na(betabinom.rslt$Reject.betabinom))
n.betabinom <- sum(betabinom.rslt$Reject.betabinom, na.rm = T)
prob.reject.fisher <-n.fisher/N.fisher
prob.reject.pearson <-n.pearson/N.pearson
prob.reject.betabinom <-n.betabinom/N.betabinom

return(list(sim.result=rslt,
            N.fisher=N.fisher,
            n.fisher=n.fisher,
            prob.reject.fisher=prob.reject.fisher,
            N.pearson=N.pearson,
            n.pearson=n.pearson,
            prob.reject.pearson=prob.reject.pearson,
            N.betabinom=N.betabinom,
            n.betabinom=n.betabinom,
            prob.reject.betabinom=prob.reject.betabinom
            ))
#options(warn = defaultW)
}

#sim.result <- power.vaf(mus, vars, n, N, nsim, alpha.level)
#
#view(sim.result$sim.result)
#
#sim.result$prob.reject.fisher
#sim.result$N.fisher
#sim.result$prob.reject.pearson
#sim.result$N.pearson
#sim.result$prob.reject.betabinom
#sim.result$N.betabinom
#
##Scenarios
#
##two groups
#
##Mus
##equal mutation rates: 0.25, 0.75
##different mutation rates: (0.25, 0.5), (0.25, 0.75)
#
##Vars
##(0, 0.1, 0.2)
#
##n
##(50, 100, 150)
#
##N
##(3,5,10)
#
##Mus, Vars, n, N, prob of reject fisher, prob of reject pearson, prob of reject betabinomial,
#
