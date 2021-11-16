library("survival")

cox_data <- read.csv("E:/sirebrowser/STAD/mRNA/gdac.broadinstitute.org_STAD.mRNAseq_Preprocess.Level_3.2016012800.0.0/cox_data.csv", header=T,row.names= 1)#存活1 死亡2

fcox <- function(x){
  FML <- as.formula(paste0('Surv(cox_data$time, cox_data$status) ~ ', x))
  res.cox <- coxph(FML, data = cox_data)
  GSum <- summary(res.cox)
  HR <- round(GSum$coefficients[, 2], 2)
  PValue <- round(GSum$coefficients[, 5], 3)
  CI <- paste0(round(GSum$conf.int[, 3:4], 2), collapse = "-")
  fcox <- data.frame('Characteristics' = x,
                       'Hazard Ratio' = HR,
                       'CI95' = CI,
                       'P Value' = PValue)
  return (fcox)
}

mirna <- colnames(cox_data)[c(2:length(cox_data[1, ])-2)]
univar_analysis <- lapply(mirna, fcox)

n <- 0 #the number of survival-related mirnas
diff_list <- list() #survival-related mirnas
HR<-list()
CI<-list()
P_value<-list()
for(i in mirna){
  a <- fcox(i)
  HR[[i]]<-a$Hazard.Ratio
  CI[[i]]<-a$CI95
  P_value[[i]]<-a$P.Value
  if(a$P.Value < 0.05){
    n <- n + 1
    diff_list[n] <- as.character(a$Characteristics)
  } 
}
n 

# output
write.csv(HR, "E:/sirebrowser/STAD/mRNA/分析/All_HR_CESC-COXPH.csv", row.names = F, quote = F)
write.csv(CI, "E:/sirebrowser/STAD/mRNA/分析/All_CI_CESC-COXPH.csv",  row.names = F, quote = F)
write.csv(P_value, "E:/sirebrowser/STAD/mRNA/分析/All_P_value_CESC-COXPH.csv", row.names = F, quote = F)
write.csv(diff_list, "E:/sirebrowser/STAD/mRNA/分析/New-CESC-COX-diff_list.csv", row.names = F, quote = F)



