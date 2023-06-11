#' Estimate parameters from Beta Binomial model
#'
#' @param data VAF data: contains Locus,chrom, vc, dp, smapleID
#'
#' @return rslt	Beta Binomail parameter estimates for each locus each group
#' @export
#'
#' @examples data(VAF)
#' @examples temp <- BetaBinompars(VAF)
#' @examples summary(temp)
#'
#' @references{
#'
#'   \insertRef{aod}{easyVAF}
#'
#' }
#'
#' @importFrom Rdpack reprompt
#'

BetaBinompars <- function(VAF){
  ID <- unique(VAF$Locus)
  groups <- unique(VAF$group)

  rslt <- data.frame(matrix(NA,    # Create empty data frame
                            nrow = length(ID)*length(groups),
                            ncol = 8))

  names(rslt) <- c("ID", "Group", "p", "phi", "a1", "a2", "var", "converge.issue")


  for(i in 1:length(ID)){
    temp <- subset(VAF, ID==ID[i])
    temp <- subset(temp, dp!=0)
    warn <- 0
    tryCatch({m1 <- betabin(cbind(vc, dp - vc) ~ group, ~ group, temp, hessian=F)},
             warning=function(w) warn <<- 1)
    warns <- rep(warn, 4)
    group=c(1:length(groups))
    id <- rep(ID[i], length(groups))
    newdata = data.frame(group=c(1:length(groups)))
    newdata$group <- factor(newdata$group)

    X <- model.matrix(~group, newdata)
    ps <- inv.logit(X%*%coef(m1))
    phis <- m1@random.param
    a2 <- (1-phis)/((phis*ps)/(1-ps)+phis)
    a1 <- ps*a2/(1-ps)
    vars <- (a1*a2)/((a1+a2)^2*(a1+a2+1))
    vars[which(vars >= 0.25)] <- NA
    rslt[(i*length(groups)-(length(groups)-1)):(i*length(groups)),] <- cbind(id, group, ps, phis, a1, a2, vars, warns
    )

  }
  return(rslt)
}

#temp <- BetaBinompars(VAF)
#
#
##plot
#temp$var <- as.numeric(temp$var)
#summary(temp$var)
#quantile(temp$var[which(temp$converge.issue==0)], probs=c(0.05, 0.25, 0.75, 0.95))
#
#library (dplyr)
#p1 <- temp %>%
#  summarize(lower = quantile(var, probs = .025),
#            upper = quantile(var, probs = .975))
#library(ggplot2)
#ggplot(temp, aes(x = var)) +
#  geom_histogram(aes(y = ..density..)) +
#  geom_density() +
#  geom_vline(data = p1, aes(xintercept = lower)) +
#  geom_vline(data = p1, aes(xintercept = upper)) + labs(x = "Estimated variance")
#
#tempmain <- VAFmain(VAF, groups=c(1:4))
#table(tempmain$Overdispersion)
#
#
