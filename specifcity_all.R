library(janitor)
library(reshape)
library (ggplot2)
library(dplyr)
library(Hmisc)
library(tidyverse)
library(data.table)
library(scales)
library(ggpubr)
library(psych)
library(readr)

#Figure 2E and 2F

specie=read.table("table_initial_specie.csv",  header=TRUE,sep=",")
specie=adorn_percentages(specie, denominator = "row", na.rm = TRUE)

uniformis = specie[,c(32,197)]
colnames(uniformis) <- c("abundance", "Donor")
uniformis$abundance=as.numeric(uniformis$abundance)*100 
uniformis_summary <- uniformis %>% 
  group_by(Donor) %>%   
  summarise(mean = mean(abundance),  
            sd = sd(abundance), 
            n = n(),  
            SE = sd(abundance)/sqrt(n())) 
uniformisplot=ggplot(uniformis_summary, aes(Donor, mean)) + 
  geom_col() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2)
uniformisplot+labs(title="B uniformis in the initial inocula", y = "Relative abundance (%)", x="Donor")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))+
  scale_x_discrete(labels= c("D1", "D2", "D3", "D4", "D5", "D6", "D7","D8","D9","D10"))+
  scale_y_continuous(breaks=seq(0, 10, by = 1))

anaerostipes = specie[,c(88,197)]
colnames(anaerostipes) <- c("abundance", "Donor")
anaerostipes$abundance=as.numeric(anaerostipes$abundance)*100 
anaerostipes_summary <- anaerostipes %>% 
  group_by(Donor) %>%   
  summarise(mean = mean(abundance),  
            sd = sd(abundance), 
            n = n(),  
            SE = sd(abundance)/sqrt(n())) 
anaerosplot=ggplot(anaerostipes_summary, aes(Donor, mean)) + 
  geom_col() +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2)
anaerosplot+labs(title="Anaerostipes sp. in the initial inocula", y = "Relative abundance (%)", x="Donor")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))+
  scale_x_discrete(labels= c("D1", "D2", "D3", "D4", "D5", "D6", "D7","D8","D9","D10"))+
  scale_y_continuous(breaks=seq(0, 1, by = 0.1))

#Figure 3  

taxtable<-read_csv("level-7_gg.csv")
metadata<-read_csv("metadata.csv")
taxtable = left_join(taxtable,metadata) 
taxtable <- taxtable[, -c(193:198)]
taxtable=aggregate(taxtable, by=list(taxtable$Donor, taxtable$Treatment, taxtable$Time), FUN=mean)
taxtable <- taxtable[, -4]
setnames(taxtable, old=c("Group.1","Group.2","Group.3"), new=c("Donor", "Treatment", "Time"))
taxtable$Time <- factor(taxtable$Time,levels = c("Initial", "Final"))
taxtable$Donor <- factor(taxtable$Donor,levels = c("D1", "D2", "D3", "D4", "D5", "D6", "D7","D8","D9","D10"))
taxtable$Treatment <- factor(taxtable$Treatment,levels = c("Blank", "FOS", "RS", "Pectin", "Glucan"))
taxtable <- taxtable[ -c(195:197) ]
taxtable=adorn_percentages(taxtable, denominator = "row", na.rm = TRUE)

bifidobacterium = taxtable[,c(1,2,3,9)]
colnames(bifidobacterium) <- c("Donor","Treatment", "Time", "Bifidobacterium")
ggplot(bifidobacterium, aes(x=Time, y=Bifidobacterium, group=Donor, color=Donor)) +
  scale_y_log10(oob = scales::squish_infinite)+
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Relative abundance (log transformed)") +
  ggtitle("Bifidobacterium relative abundance before and after fiber treatment")

rbromii = taxtable[,c(1,2,3,145)]
colnames(rbromii) <- c("Donor","Treatment", "Time", "rbromii")
ggplot(rbromii, aes(x=Time, y=rbromii, group=Donor, color=Donor)) +
  scale_y_log10(oob = scales::squish_infinite)+
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Relative abundance (log transformed)") +
  ggtitle("R. bromii relative abundance before and after fiber treatment")

lachnospira = taxtable[,c(1,2,3,111)]
colnames(lachnospira) <- c("Donor","Treatment", "Time", "lachnospira")
ggplot(lachnospira, aes(x=Time, y=lachnospira, group=Donor, color=Donor)) +
  scale_y_log10(oob = scales::squish_infinite)+
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Relative abundance (log transformed)") +
  ggtitle("Lachnospira sp. relative abundance before and after fiber treatment")

anaerostipes = taxtable[,c(1,2,3,90)]
colnames(anaerostipes) <- c("Donor","Treatment", "Time", "Anaerostipes")
ggplot(anaerostipes, aes(x=Time, y=Anaerostipes, group=Donor, color=Donor)) +
  scale_y_log10(oob = scales::squish_infinite)+
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Relative abundance (log transformed)") +
  ggtitle("Anaerostipes log of relative abundance before and after fiber treatment")

buniformis = taxtable[,c(1,2,3,34)]
colnames(buniformis) <- c("Donor","Treatment", "Time", "buniformis")
ggplot(buniformis, aes(x=Time, y=buniformis, group=Donor, color=Donor)) +
  scale_y_log10(oob = scales::squish_infinite) +
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Relative abundance (log transformed)") +
  ggtitle("B. uniformis relative abundance before and after fiber treatment")

cramosum = taxtable[,c(1,2,3,166)]
colnames(cramosum) <- c("Donor","Treatment", "Time", "cramosum")
ggplot(cramosum, aes(x=Time, y=cramosum, group=Donor, color=Donor)) +
  scale_y_log10(oob = scales::squish_infinite)+
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Relative abundance (log transformed)") +
  ggtitle("Clostridium ramosum relative abundance before and after fiber treatment")

#Figure 6A,B and C 

SCFA<-read_csv("SCFA_forR_ correlation.csv")
SCFA <- SCFA[c(4,5,1,2,3)]
SCFA=adorn_percentages(SCFA, denominator = "row", na.rm = TRUE)
SCFA=aggregate(SCFA[, 3:5], list(SCFA$Donor,SCFA$Treatment), mean)
names(SCFA)[1]<-paste("Donor") 
names(SCFA)[2]<-paste("Treatment") 
SCFA$Donor = as.character(SCFA$Donor)
ggplot(SCFA,aes(x=Treatment, y=Butyrate, color=as.factor(Donor))) + 
  geom_point(position=position_jitterdodge(dodge.width=0.3),size = 2.5)+
  theme_bw() +
  geom_boxplot(color="black")+
  theme (legend.title = element_blank(), axis.title.x = element_blank()) 
ggplot(SCFA,aes(x=Treatment, y=Propionate, color=as.factor(Donor))) + 
  geom_point(position=position_jitterdodge(dodge.width=0.3),size = 2.5)+
  theme_bw() +
  geom_boxplot(color="black")+
  theme (legend.title = element_blank(), axis.title.x = element_blank()) 
ggplot(SCFA,aes(x=Treatment, y=Acetate, color=Donor)) + 
  geom_point(position=position_jitterdodge(dodge.width=0.3),size = 2.5)+
  theme_bw() +  
  geom_boxplot(color="black")+
  theme (legend.title = element_blank(), axis.title.x = element_blank()) 

##To make boxplots excluding outliers (D7, D9 and D10) - used for glucan only
SCFA_noout = SCFA[!(SCFA$Donor=="10" | SCFA$Donor=="9"| SCFA$Donor=="7"),]
ggplot(SCFA_noout,aes(x=Treatment, y=Butyrate, color=as.factor(Donor))) + 
  geom_point(position=position_jitterdodge(dodge.width=0.3),size = 2.5)+
  theme_bw() +
  geom_boxplot(color="black")+
  theme (legend.title = element_blank(), axis.title.x = element_blank()) 
ggplot(SCFA_noout,aes(x=Treatment, y=Propionate, color=as.factor(Donor))) + 
  geom_point(position=position_jitterdodge(dodge.width=0.3),size = 2.5)+
  theme_bw() +
  geom_boxplot(color="black")+
  theme (legend.title = element_blank(), axis.title.x = element_blank()) 
ggplot(SCFA_noout,aes(x=Treatment, y=Acetate, color=Donor)) + 
  geom_point(position=position_jitterdodge(dodge.width=0.3),size = 2.5)+
  theme_bw() +  
  geom_boxplot(color="black")+
  theme (legend.title = element_blank(), axis.title.x = element_blank()) 


#Figure 6D

bacteria<-read_csv("hightha1perc_forR.csv")
SCFA<-read_csv("SCFA_forR_ correlation.csv")
SCFA = SCFA %>% 
  group_by(SCFA$Donor, SCFA$Treatment) %>% 
  summarise_at(.vars = names(.)[1:3],
               .funs = c(mean))%>%
  ungroup()
bacteria = bacteria %>% 
  group_by(bacteria$donor, bacteria$treatment_II_B) %>% 
  summarise_at(.vars = names(.)[1:26],
               .funs = c(mean)) %>%
  ungroup()
setnames(SCFA, old=c("SCFA$Donor","SCFA$Treatment"), new=c("Donor", "Treatment"))
setnames(bacteria, old=c("bacteria$treatment_II_B","bacteria$donor"), new=c("Treatment", "Donor"))
bacteria=adorn_percentages(bacteria, denominator = "row", na.rm = TRUE)
SCFA=adorn_percentages(SCFA, denominator = "row", na.rm = TRUE)
SCFAbac = left_join(SCFA,bacteria) 
SCFAbac = as.matrix(SCFAbac[,3:31])
correlation = corr.test(SCFAbac,
                        use = "pairwise",
                        method="spearman",
                        adjust="fdr",     
                        alpha=.05)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
correlation= flattenCorrMatrix(correlation$r, correlation$p)
corbut = correlation[grep("Butyrate", correlation$row),]
corprop = correlation[grep("Propionate", correlation$row),]
corprop = corprop[-1,]
corace = correlation[grep("Acetate", correlation$row),]
corace = corace[-c(1,2),]
corunique = rbind(corbut, corprop, corace[, colnames(corbut)]) #junta as tabelas de acetato, propionato e butira
corunique$stars <- cut(corunique$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  
p=ggplot(aes(x=row, y=column, fill=cor), data=corunique)
p + geom_tile() + 
  scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6",
                       midpoint=0,
                       limits=c(-1,1)) + 
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="Spearsman's rho") + 
  theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))

#Figure 6E and SCFA distance matrix further imported into Qiime and mothur

SCFA_dm<-read_csv("SCFA_forR_ distancematrix_FN2.csv")
SCFA_dm <- SCFA_dm[c(4,1,2,3)]
SCFA_dm=adorn_percentages(SCFA_dm, denominator = "row", na.rm = TRUE)
D <-  data.frame(SCFA_dm$D)
SCFA_dm <- SCFA_dm[ -c(1) ]
p <- prcomp(SCFA_dm, scale=TRUE)
s <- summary(p)
pch.group <- c(rep(21,29), rep(22,30), rep(23,30), rep(24,30), rep(25,30))
col.group <- c(rep("skyblue2", times=29), rep("gold", times=30), rep("green2", times=30), rep("#6A3D9A", times=30), rep("#E31A1C", times=30))
plot(p$x[,1], p$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=3, las=1, asp=1) + 
  abline(v=0, lty=2, col="grey50") + 
  abline(h=0, lty=2, col="grey50") + 
  text(p$x[,1], p$x[,2],labels=D$SCFA_dm.D, font=2) 
s$importance
l.x <- p$rotation[,1]*3.5
l.y <- p$rotation[,2]*3.5
arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=0.15, lwd=1.5)
l.pos <- l.y 
lo <- which(l.y < 0) 
hi <- which(l.y > 0) 
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")
text(l.x, l.y, labels=row.names(p$rotation), col="red", pos=l.pos)
legend("topleft", legend=c("Blank","FOS","Glucan","Pectin","RS"), col="black", pt.bg=c( "skyblue2","gold","green2","#6A3D9A","#E31A1C"), pch=c(21, 22, 23, 24, 25), pt.cex=1.5)

#SCFA distance matrix further imported into Qiime and mothur 
library(tidyverse)
SCFA_dm<-read_csv("SCFA_forR_ distancematrix_noout.csv")
SCFA_dm <- SCFA_dm[c(4,1,2,3)]
library(janitor)
SCFA_dm=adorn_percentages(SCFA_dm, denominator = "row", na.rm = TRUE)
SCFA_dm=column_to_rownames(SCFA_dm,var = "index")
distance=dist(SCFA_dm)
distance <- as.matrix(distance)
write.csv(distance, file = "SCFA_distmat_noout.csv")

#Figure S3

# Functions 
#Data Pre-Processing
feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var = NULL, 
                                     out_cut = 0.05, zero_cut = 0.90, lib_cut = 1000, neg_lb){
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  # Match sample IDs between metadata and feature table
  sample_ID = intersect(meta_data[, sample_var], colnames(feature_table))
  feature_table = feature_table[, sample_ID]
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]
  
  # 1. Identify outliers within each taxon
  if (!is.null(group_var)) {
    group = meta_data[, group_var]
    z = feature_table + 1 # Add pseudo-count (1) 
    f = log(z); f[f == 0] = NA; f = colMeans(f, na.rm = T)
    f_fit = lm(f ~ group)
    e = residuals(f_fit)
    y = t(t(z) - e)
    
    outlier_check = function(x){
      # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
      mu1 = quantile(x, 0.25, na.rm = T); mu2 = quantile(x, 0.75, na.rm = T)
      sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T); sigma2 = sigma1
      pi = 0.75
      n = length(x)
      epsilon = 100
      tol = 1e-5
      score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
      while (epsilon > tol) {
        grp1_ind = (score >= 1)
        mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
        sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
        sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
        pi_new = sum(grp1_ind)/n
        
        para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
        if(any(is.na(para))) break
        
        score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
          ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))
        
        epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 + 
                         (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
        mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
      }
      
      if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
        if(pi < out_cut){
          out_ind = grp1_ind
        }else if(pi > 1 - out_cut){
          out_ind = (!grp1_ind)
        }else{
          out_ind = rep(FALSE, n)
        }
      }else{
        out_ind = rep(FALSE, n)
      }
      return(out_ind)
    }
    out_ind = t(apply(y, 1, function(i) unlist(tapply(i, group, function(j) outlier_check(j)))))
    feature_table[out_ind] = NA
  }
  
  # 2. Discard taxa with zeros  >=  zero_cut
  zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  taxa_del = which(zero_prop >= zero_cut)
  if(length(taxa_del) > 0){
    feature_table = feature_table[- taxa_del, ]
  }
  
  # 3. Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = T)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }
  
  # 4. Identify taxa with structure zeros
  if (!is.null(group_var)) {
    group = factor(meta_data[, group_var])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1
    
    p_hat = t(apply(present_table, 1, function(x)
      unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
    samp_size = t(apply(feature_table, 1, function(x)
      unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)
    
    struc_zero = (p_hat == 0) * 1
    # Whether we need to classify a taxon into structural zero by its negative lower bound?
    if(neg_lb) struc_zero[p_hat_lo <= 0] = 1
    
    # Entries considered to be structural zeros are set to be 0s
    struc_ind = struc_zero[, group]
    feature_table = feature_table * (1 - struc_ind)
    
    colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
  }else{
    struc_zero = NULL
  }
  
  # 5. Return results
  res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero)
  return(res)
}

# ANCOM main function
ANCOM = function(feature_table, meta_data, struc_zero = NULL, main_var, p_adj_method = "BH", 
                 alpha = 0.05, adj_formula = NULL, rand_formula = NULL){
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  comp_table = log(as.matrix(comp_table) + 1)
  n_taxa = dim(comp_table)[1]
  taxa_id = rownames(comp_table)
  n_samp = dim(comp_table)[2]
  
  # Determine the type of statistical test and its formula.
  if (is.null(rand_formula) & is.null(adj_formula)) {
    # Basic model
    # Whether the main variable of interest has two levels or more?
    if (length(unique(meta_data%>%pull(main_var))) == 2) {
      # Two levels: Wilcoxon rank-sum test
      tfun = exactRankTests::wilcox.exact
    } else{
      # More than two levels: Kruskal-Wallis test
      tfun = stats::kruskal.test
    }
    # Formula
    tformula = formula(paste("x ~", main_var, sep = " "))
  }else if (is.null(rand_formula) & !is.null(adj_formula)) {
    # Model: ANOVA
    tfun = stats::aov
    # Formula
    tformula = formula(paste("x ~", main_var, "+", adj_formula, sep = " "))
  }else if (!is.null(rand_formula)) {
    # Model: Mixed-effects model
    tfun = nlme::lme
    # Formula
    if (is.null(adj_formula)) {
      # Random intercept model
      tformula = formula(paste("x ~", main_var))
    }else {
      # Random coefficients/slope model
      tformula = formula(paste("x ~", main_var, "+", adj_formula))
    }
  }
  
  # Calculate the p-value for each pairwise comparison of taxa.
  p_data = matrix(NA, nrow = n_taxa, ncol = n_taxa)
  colnames(p_data) = taxa_id
  rownames(p_data) = taxa_id
  for (i in 1:(n_taxa - 1)) {
    # Loop through each taxon.
    # For each taxon i, additive log ratio (alr) transform the OTU table using taxon i as the reference.
    # e.g. the first alr matrix will be the log abundance data (comp_table) recursively subtracted 
    # by the log abundance of 1st taxon (1st column) column-wisely, and remove the first i columns since:
    # the first (i - 1) columns were calculated by previous iterations, and
    # the i^th column contains all zeros.
    alr_data = apply(comp_table, 1, function(x) x - comp_table[i, ]) 
    # apply(...) allows crossing the data in a number of ways and avoid explicit use of loop constructs.
    # Here, we basically want to iteratively subtract each column of the comp_table by its i^th column.
    alr_data = alr_data[, - (1:i), drop = FALSE]
    n_lr = dim(alr_data)[2] # number of log-ratios (lr)
    alr_data = cbind(alr_data, meta_data) # merge with the metadata
    
    # P-values
    if (is.null(rand_formula) & is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        tfun(tformula, data = data.frame(x, alr_data, check.names = FALSE))$p.value
      }
      ) 
    }else if (is.null(rand_formula) & !is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE), 
                   na.action = na.omit)
        summary(fit)[[1]][main_var, "Pr(>F)"]
      }
      )
    }else if (!is.null(rand_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(fixed = tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE),
                   random = formula(rand_formula),
                   na.action = na.omit)
        anova(fit)[main_var, "p-value"]
      }
      ) 
    }
  }
  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
  diag(p_data) = 1 # let p-values on diagonal equal to 1
  
  # Multiple comparisons correction.
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = p_adj_method))
  
  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W = apply(q_data, 2, function(x) sum(x < alpha))
  
  # Organize outputs
  out = data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
  # Declare a taxon to be differentially abundant based on the quantile of W statistic.
  # We perform (n_taxa - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_taxa - 1).
  out = out%>%mutate(detected_0.9 = ifelse(W > 0.9 * (n_taxa -1), TRUE, FALSE),
                     detected_0.8 = ifelse(W > 0.8 * (n_taxa -1), TRUE, FALSE),
                     detected_0.7 = ifelse(W > 0.7 * (n_taxa -1), TRUE, FALSE),
                     detected_0.6 = ifelse(W > 0.6 * (n_taxa -1), TRUE, FALSE))
  
  res = out
}


#DATA ANALYSIS 
otu_data = read_csv("taxtable_ancom.csv")
otu_id = otu_data$`#SampleID`
otu_data = data.frame(otu_data[, -1], check.names = FALSE)
row.names(otu_data) = otu_id
meta_data1 = read_csv("metadata.csv")
meta_data1 = meta_data1 %>% dplyr::rename(Sample.ID = 'index')

#PECTIN
meta_data <- subset(meta_data1, Treatment == c("Pectin", "Blank"))
# Step 1: Data preprocessing
feature_table = otu_data; sample_var = "Sample.ID"; group_var = "Treatment"
out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info
# Step 2: ANCOM
main_var = "Treatment"; p_adj_method = "holm"; alpha = 0.05
adj_formula = "Time"; rand_formula = "~ 1 | Donor"
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start
resPECT = res

#GLUCAN
meta_data <- subset(meta_data1, Treatment == c("Glucan", "Blank"))
# Step 1: Data preprocessing
feature_table = otu_data; sample_var = "Sample.ID"; group_var = "Treatment"
out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info
# Step 2: ANCOM
main_var = "Treatment"; p_adj_method = "holm"; alpha = 0.05
adj_formula = "Time"; rand_formula = "~ 1 | Donor"
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start
resGLUCAN = res

#FOS
meta_data <- subset(meta_data1, Treatment == c("FOS", "Blank"))
# Step 1: Data preprocessing
feature_table = otu_data; sample_var = "Sample.ID"; group_var = "Treatment"
out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info
# Step 2: ANCOM
main_var = "Treatment"; p_adj_method = "holm"; alpha = 0.05
adj_formula = "Time"; rand_formula = "~ 1 | Donor"
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start
resFOS = res

#RS
meta_data <- subset(meta_data1, Treatment == c("RS", "Blank"))
# Step 1: Data preprocessing
feature_table = otu_data; sample_var = "Sample.ID"; group_var = "Treatment"
out_cut = 0.05; zero_cut = 0.90; lib_cut = 0; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info
# Step 2: ANCOM
main_var = "Treatment"; p_adj_method = "holm"; alpha = 0.05
adj_formula = "Time"; rand_formula = "~ 1 | Donor"
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start
resRS = res

#Figure S4

taxtable<-read_csv("level-7_gg.csv")
taxtable$donor = as.character(taxtable$donor)
taxtable$time_0_1_2 = as.character(taxtable$time_0_1_2)
taxtable=adorn_percentages(taxtable, denominator = "row", na.rm = TRUE)
#Anaerostipes
anaero_glucan = taxtable[,c(88,193,195,197)]
colnames(anaero_glucan) <- c("Anaerostipes","Treatment", "Time","Donor")
anaero_glucan = subset(anaero_glucan, Treatment== 'Glucan' | Treatment== 'II')
anaero_glucan$Donor <- factor(anaero_glucan$Donor,levels = c("1", "2", "3", "4", "5", "6", "7","8","9","10"))
anaero_glucan$Time <- factor(anaero_glucan$Time,levels = c("Initial","Final"))
anaero_glucan <- anaero_glucan %>% 
  group_by(Donor, Time) %>%   
  dplyr::summarise(mean = mean(Anaerostipes),  # calculates the mean of each group
                   sd = sd(Anaerostipes), # calculates the standard deviation of each group
                   n = n(),  # calculates the sample size per group
                   SE = sd(Anaerostipes)/sqrt(n())) # calculates the standard error of each group
ggplot(data=anaero_glucan, aes(x=Donor, y=mean, fill=Time)) +
  geom_bar(stat="identity", position=position_dodge(0.5))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2, position=position_dodge(0.5))
#Buniformis
Buniformis_glucan = taxtable[,c(32,193,195,197)]
colnames(Buniformis_glucan) <- c("Buniformis","Treatment", "Time","Donor")
Buniformis_glucan = subset(Buniformis_glucan, Treatment== 'Glucan' | Treatment== 'II')
Buniformis_glucan$Donor <- factor(Buniformis_glucan$Donor,levels = c("1", "2", "3", "4", "5", "6", "7","8","9","10"))
Buniformis_glucan$Time <- factor(Buniformis_glucan$Time,levels = c("Initial","Final"))
Buniformis_glucan <- Buniformis_glucan %>% # 
  group_by(Donor, Time) %>%   
  dplyr::summarise(mean = mean(Buniformis),  # calculates the mean of each group
                   sd = sd(Buniformis), # calculates the standard deviation of each group
                   n = n(),  # calculates the sample size per group
                   SE = sd(Buniformis)/sqrt(n())) # calculates the standard error of each group
ggplot(data=Buniformis_glucan, aes(x=Donor, y=mean, fill=Time)) +
  geom_bar(stat="identity", position=position_dodge(0.5))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2, position=position_dodge(0.5))
cramosum_glucan = taxtable[,c(164,193,195,197)]
colnames(cramosum_glucan) <- c("cramosum","Treatment", "Time","Donor")
cramosum_glucan = subset(cramosum_glucan, Treatment== 'Glucan' | Treatment== 'II')
cramosum_glucan$Donor <- factor(cramosum_glucan$Donor,levels = c("1", "2", "3", "4", "5", "6", "7","8","9","10"))
cramosum_glucan$Time <- factor(cramosum_glucan$Time,levels = c("Initial","Final"))
cramosum_glucan <- cramosum_glucan %>% # 
  group_by(Donor, Time) %>%  
  dplyr::summarise(mean = mean(cramosum),  # calculates the mean of each group
                   sd = sd(cramosum), # calculates the standard deviation of each group
                   n = n(),  # calculates the sample size per group
                   SE = sd(cramosum)/sqrt(n())) # calculates the standard error of each group
ggplot(data=cramosum_glucan, aes(x=Donor, y=mean, fill=Time)) +
  geom_bar(stat="identity", position=position_dodge(0.5))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2, position=position_dodge(0.5))


