#' Main comparison of VAFs among N groups
#'
#' We recommend VAF analysis work flow as following: 1). Start with exploratory plots for Variant Allele Count, Read Depth, and VAF for quality checking (i.e., unexpected biological variability, batch effect, technical effect); 2). Conduct statistical test to assess the variability of overall VAF distribution among the experiment samples (i.e., test if the heterogeneity of experiment samples is significant within each treatment group); 3).The main comparison of VAFs will be conducted as described below: a). For each locus, the goodness of fit test for binomial distribution (overdispersion) is conducted first; b). Appropriate method (model-based or non-parametric) will be selected to perform the VAF comparison among treatment groups; c).	The raw and adjusted p-values will be reported for each locus, accompanied with the estimated VAFs, difference in VAFs and the corresponding confidence intervals (only available for two group comparisons).
#' @param data VAF data: contains Locus,chrom, vc, dp, smapleID
#' @param method Choose from "betabinom", "Exact", "Pearson", the default is NULL, if NULL then conduct main comparison of VAFs following step 3).
#' @param groups vector of group names, at least two groups
#' @param digits digit place in result table
#'
#' @return rslt	Analysis result table
#' @export
#'
#' @examples data(VAF)
#' @examples groups <- unique(VAF$group)[c(1:4)]
#' @examples rslt <- VAFmain(data=VAF, groups=groups)
#'
#' @references{
#'   \insertRef{fisher1970ra}{easyVAF}
#'
#'   \insertRef{pearson1900x}{easyVAF}
#'
#'   \insertRef{prentice1986binary}{easyVAF}
#'
#'   \insertRef{tarone1979testing}{easyVAF}
#'
#'   \insertRef{aod}{easyVAF}
#'
#' }
#'
#' @importFrom Rdpack reprompt
#'
VAFmain <- function(data, method=NULL, groups, digits=3){
  defaultW <- getOption("warn")
  options(warn = -1)
  ID <- unique(data$Locus)
  def <- is.null(method)
  #initialization
  rslt <- data.frame(matrix(NA,    # Create empty data frame
                            nrow = length(ID),
                            ncol = 3*length(groups)+5+1))
  names(rslt) <- c("ID",
                   as.vector(outer(c("Read.Depth", "Variant.Count", "VAF"),groups, paste, sep=".")),
                   "P.value",
                   "Test",
                   "Effect.size",
                   "95% CI",
                   "Overdispersion"
  )

  for(i in 1:length(ID)){
    temp <- subset(data, ID==ID[i]&group %in% groups)
    temp <- subset(temp, dp!=0)
    #Tarone.test
    dispersion <- rep(NA, length(groups))
    for(j in 1:length(groups)){
      NMtemp <- subset(temp, group==groups[j], select=c("dp", "vc"))
      testemp <- Tarone.test(NMtemp$dp, NMtemp$vc)
      dispersion[j] <- ifelse(testemp$p.value < 0.05, T, F)
    }
    overdispersion <- sum(dispersion, na.rm = T)
    Odispersion <- ifelse(overdispersion>0, "Yes", "No")

    if(def){
      method <- ifelse(overdispersion>0, "betabinom", "Exact")
    }



      temp %>%
        group_by(group) %>%
        summarise(dp = sum(dp), vc= sum(vc), vaf=sum(vc)/sum(dp)) -> vaf

    if(method=="Exact"){
      tab <- as.matrix(cbind(vaf$vc, vaf$dp-vaf$vc))
      fisherrslt <- fisher.test(tab, alternative = "two.sided")
      pvalue <- fisherrslt$p.value
      test <- "Fisher's exact"
    }else if(method=="Pearson"){
      chirslt <- prop.test(vaf$vc, vaf$dp, alternative = "two.sided")
      pvalue <- chirslt$p.value
      test <- "Chi-squared"
    }else if(method=="betabinom"){
      m1 <- betabin(cbind(vc, dp - vc) ~ group, ~ group, temp, hessian=F)
      m0 <- betabin(cbind(vc, dp - vc) ~ 1, ~ group, temp, hessian=F)
      betarslt <- anova(m0, m1, "LRT")
      pvalue <- betarslt@anova.table[["P(> Chi2)"]][2]
      test <- "Beta-Binomial LRT"
    }
      chirslt <- prop.test(vaf$vc, vaf$dp, alternative = "two.sided")
      if(length(groups)==2){
        effsize <- paste("Diff in prop = ", round(chirslt$estimate[1]-chirslt$estimate[2],digits), sep="")
        CI <- paste("(",paste(round(chirslt$conf.int,digits), collapse=", "),")", sep="")
      }else{
        effsize <- "-"
        CI <- "-"
      }


    rslt[i,] <- c(ID[i],
                  round(as.vector(t(vaf[,-1])), digits),
                  pvalue,
                  test,
                  effsize,
                  CI,
                  Odispersion
    )

  }

  rslt$p.adjust <- p.adjust(rslt$P.value, "BH")
  rslt$sig.diff <- ifelse(rslt$P.value < 0.05, "Difference", "No difference")
  rslt$sig.diff.fdr <- ifelse(rslt$p.adjust < 0.05, "Difference", "No difference")

  if(length(groups)==2){
    rslt$sig.change.20 <- ifelse(abs(as.numeric(rslt[,1+3])-as.numeric(rslt[,1+2*3])) > 0.2, "Difference", "No difference")

    rslt$Change.direction <- ifelse(as.numeric(rslt[,1+3])-as.numeric(rslt[,1+2*3]) > 0, "Group1 > Group2",
                                    ifelse(as.numeric(rslt[,1+3])-as.numeric(rslt[,1+2*3]) < 0, "Group1 < Group2",
                                           NA
                                    ))
  }else{
    rslt$sig.change.20 <- "-"
    rslt$Change.direction <- "-"
  }
  options(warn = defaultW)
  return(rslt)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Method for checking biological variability among samples (mice or human)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Quality checking for biological variability among samples
#'
#' @param data VAF data: contains locus,chrom, vc, dp, smapleID
#' @param method choose from "lm" or "lmer"
#'
#' @return rslt QC result table
#' @export
#'
#' @examples data(VAF)
#' @examples rslt <- QCchecking(data=VAF, method="lm")
#'
#' @references{
#'   \insertRef{ggplot2}{easyVAF}
#'
#'   \insertRef{lme4}{easyVAF}
#'
#' }
QCchecking <- function(data, method="lm"){
  data$vaf <- data$vc/data$dp
  if(method=="lmer"){
    m0 <- lmer(vaf ~ group+(1|chrom), dat=data, REML=F)
    m1 <- lmer(vaf ~ group/sample+(1|chrom), dat=data, REML=F)
    rslt <- anova(m1, m0)
  }else if(method=="lm"){
    m0 <- lm(vaf ~ group, dat=data)
    m1 <- lm(vaf ~ group/sample, dat=data)
    rslt <- anova(m1, m0)
  }
  p <- ggplot(data, aes(x=vaf, color=sample, fill=sample)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
    labs(title="VAF histogram plot",x="VAF", y = "Density")+
    theme_classic() + facet_grid(rows = vars(group))
  print(p)
  return(rslt)
}

#temp <- QCchecking(dat, "lmer")

