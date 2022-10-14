set.seed(100)

ANOVA_by_group <- function(var1){
  temp_db <- data.frame(apply(MOMcomb_sub[,c(103:106)],2,function(x)
    data.frame(y = MOMcomb_sub[[var1]],x,
               Var = var1) %>%
      filter(y!="NA"&complete.cases(x)) %>% 
      mutate(p.value = round(kruskal.test(x~y)$p.value,4)) %>%
      group_by(y) %>% 
      summarise(Median = median(x),
                Q1 = summary(x)[[2]],
                Q3 = summary(x)[[5]],
                Var = rep(unique(Var),length(Q1)),
                p.val = rep(unique(p.value),length(Q1)))))
  return(temp_db)
}



spearcor <- function (mat, ...)
{
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(x = mat[, i], y = mat[, j],method = "spearman")
      tmp1 <- tryCatch(print(spearmanCI(mat[, i],mat[, j],level = 0.95, method = "Euclidean", plot = FALSE)),
                       error = function(e) {NA})
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- (1-tmp1)/2
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp1+(1-tmp1)/2
    }
  }
  list(p = p.mat, lowCI = lowCI.mat, uppCI = uppCI.mat)
}


test <- function(x,y, data = ""){
  db <- data.frame(data)
  uni_c <- data.frame()
  if(is.numeric(db[[x]])){
    temp_a <-glm(y ~ x, data = db, familiy = "binomial")
    temp_b <- cbind(summary(temp_a)$coefficients[,c(1,4)],
                    exp(cbind(OR = coef(temp_a), 
                              confint(temp_a))),
                    Outcome = j,
                    Var = i)
    uni_c <- rbind(uni_c,temp_b)
  }else{is.character(db[[x]])&length(unique(db[[x]]))>1}{
    temp_a <- wilcox.test(x ~ y, data = db, 
                          conf.int = T, conf.level=0.95)
    temp_b <- data.frame(OR = 1-temp_a$estimate,
                         Lower = 1-temp_a$conf.int[2],
                         Upper = 1-temp_a$conf.int[1],
                         p.value = temp_a$p.value,
                         Outcome = j,
                         Var = i)
    uni_c <- rbind(uni_c,temp_b)
  }
  return(uni_c)
}


uni_c <- data.frame()
for(j in colnames(FS_FA_Del)[c(136:137)]){
  for (i in colnames(FS_FA_Del)[c(4:8,84,85,45,44,9:22,24:26,
                                  139:151)]){
    if (all(FS_FA_Del[,i]=="No"|FS_FA_Del[,i]=="1No"|
            FS_FA_Del[,i]=="Fresh")){
      print(i)
    }
    else {
      print(i)
      temp_a <- glm(as.formula(paste0(j, "~",i)),
                    FS_FA_Del,family ="gaussian")
      temp_b <- cbind(summary(temp_a)$coefficients[,c(1)],
                      confint(temp_a,level = 0.95),
                      summary(temp_a)$coefficients[,c(4)],
                      Outcome = j,
                      Var = i)
      uni_c <- rbind(uni_c,temp_b)
    }
  }
}

for(j in colnames(FS_FA_Del)[c(138,153)]){
  for (i in colnames(FS_FA_Del)[c(4:8,84,85,45,44,9:22,24:26,
                                  139:151)]){
    if (all(FS_FA_Del[,i]=="No"|FS_FA_Del[,i]=="1No"|
            FS_FA_Del[,i]=="Fresh")){
      print(i)
    }
    else {
      print(i)
      temp_a <- glm(as.formula(paste0(j, " ~",i)),
                    FS_FA_Del,family ="binomial")
      temp_b <- cbind(summary(temp_a)$coefficients[,c(1,4)],
                      exp(cbind(OR = coef(temp_a), 
                                confint(temp_a))),
                      Outcome = j,
                      Var = i)
      temp_b <- temp_b[,c(3:5,2,6,7)]
      colnames(temp_b) <- colnames(uni_c)
      uni_c <- rbind(uni_c,temp_b)
    }
  }
}


uni_c <- uni_c %>% filter(!grepl("Intercep",rownames(uni_c)))
# rownames(uni_c) <- colnames(HLMdb)[c(4:7,13,14,18,22:34,36)]

uni_c[,c(1:4)] <- apply(uni_c[,c(1:4)], 2, function(x) round(as.numeric(x),4))
uni_c <- data.frame(uni_c)
# uni_c$Pr...z.. <- ifelse(uni_c$Pr...z..<0.001, paste("<0.001"),paste(uni_c$Pr...z..))

colnames(uni_c) <- c("Est/OR","Lower","Upper","p.value","Outcome","Var")
uni_c$p.value <- ifelse(uni_c$p.value<0.001, paste("<0.001"),paste(uni_c$p.value))
uni_c$Var_strat <- rownames(uni_c)
writexl::write_xlsx(uni_c,"Results/FA_FS/FS_Mother_Gcomp_Uni_21062022_all.xlsx")
uni_c <- uni_c %>% filter(p.value<0.06)
writexl::write_xlsx(uni_c,"Results/FA_FS/FS_Mother_Gcomp_Uni_21062022_sig.xlsx")


mglmtb <- data.frame()

for(j in colnames(FS_FA_Del)[c(136:137)]){
  for (i in colnames(FS_FA_Del)[31:43]){
    if (paste0(i,"_cut") %in% uni_c$Var[uni_c$Outcome==j]){
      print(i)
      temp_a <- glm(as.formula(paste0(j, "~",i, "_cut", "+ Age + Pre_Weight + SexOB +", 
                                      paste0(uni_c$Var[!grepl("cut",
                                                              uni_c$Var)&uni_c$Outcome==j],
                                             collapse = "+"))),
                    FS_FA_Del,family ="gaussian")
      temp_b <- data.frame(cbind(summary(temp_a)$coefficients[,c(1)],
                                 confint(temp_a,level = 0.95),
                                 summary(temp_a)$coefficients[,c(4)],
                                 Outcome = j,
                                 Var = i))
      mglmtb <- rbind(mglmtb,temp_b[2,])
    }
    else {
      next
    }
  }
}

for(j in colnames(FS_FA_Del)[c(138,153)]){
  for (i in colnames(FS_FA_Del)[31:43]){
    if (paste0(i,"_cut") %in% uni_c$Var[uni_c$Outcome==j]){
      print(i)
      temp_a <- glm(as.formula(paste0(j, "~",i, "_cut", "+ Age + Pre_Weight +", 
                                      paste0(uni_c$Var[!grepl("cut",uni_c$Var)&uni_c$Outcome==j],
                                             collapse = "+"))),
                    FS_FA_Del,family ="binomial")
      temp_b <- data.frame(cbind(summary(temp_a)$coefficients[,c(1,4)],
                                 exp(cbind(OR = coef(temp_a), 
                                           confint(temp_a))),
                                 Outcome = j,
                                 Var = i))
      temp_b <- temp_b[,c(3:5,2,6,7)]
      colnames(temp_b) <- colnames(mglmtb)
      mglmtb <- rbind(mglmtb,temp_b[2,])
    }
    else {
      next
    }
  }
}

colnames(mglmtb) <-  c("Est/OR","Lower","Upper","p.value","Outcome","Var")

mglmtb[,1:4] <- apply(mglmtb[,1:4], 2, function(x) round(as.numeric(x),4))


mglmtb$p.value <- ifelse(mglmtb$p.value<0.001, paste("<0.001"),
                         paste(mglmtb$p.value))

mglmtb$Var_Strat <- rownames(mglmtb)
writexl::write_xlsx(mglmtb,"Results/FA_FS/FS_Mother_Gcomp_Mul_23062022_all.xlsx")
