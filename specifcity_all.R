library(janitor)
library (ggplot2)
library(dplyr)
library(Hmisc)
library(tidyverse)
library(data.table)
library(scales)
library(ggpubr)
library(psych)

#Figure 1E and 1F

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
anaerostipes$abundance=as.numeric(anaerostipes$abundance)*100 #transforma em percentagem
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

#Figure 3 and S2 

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

blautia = taxtable[,c(1,2,3,91)]
colnames(blautia) <- c("Donor","Treatment", "Time", "blautia")
ggplot(blautia, aes(x=Time, y=blautia, group=Donor, color=Donor)) +
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Blautia relative abundance(%)") 

bacteroidessp = taxtable[,c(1,2,3,27)]
colnames(bacteroidessp) <- c("Donor","Treatment", "Time", "bacteroidessp")
ggplot(bacteroidessp, aes(x=Time, y=bacteroidessp, group=Donor, color=Donor)) +
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Bacteroides relative abundance (%)") 

lachnospira = taxtable[,c(1,2,3,111)]
colnames(lachnospira) <- c("Donor","Treatment", "Time", "lachnospira")
ggplot(lachnospira, aes(x=Time, y=lachnospira, group=Donor, color=Donor)) +
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Lachnospira relative abundance (%)")

anaerostipes = taxtable[,c(1,2,3,90)]
colnames(anaerostipes) <- c("Donor","Treatment", "Time", "Anaerostipes")
ggplot(anaerostipes, aes(x=Time, y=Anaerostipes, group=Donor, color=Donor)) +
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("Anaerostipes relative abundance (%)")

buniformis = taxtable[,c(1,2,3,34)]
colnames(buniformis) <- c("Donor","Treatment", "Time", "buniformis")
ggplot(buniformis, aes(x=Time, y=buniformis, group=Donor, color=Donor)) +
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("B. uniformis relative abundance (%)") 

rbromii = taxtable[,c(1,2,3,145)]
colnames(rbromii) <- c("Donor","Treatment", "Time", "rbromii")
ggplot(rbromii, aes(x=Time, y=rbromii, group=Donor, color=Donor)) +
  geom_point(size=1) +
  geom_line() +
  facet_grid(~Treatment) + 
  theme_bw() +
  theme (legend.title = element_blank(), axis.title.x = element_blank()) +
  ylab("R. brommi relative abundance (%)")

#Figure 5A 

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

#Figure 5C and Figure S3

SCFA<-read_csv("2(1).csv")
ggplot(SCFA,aes(x=Treatment, y=perc, color=SCFA))  + 
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white",outlier.colour = NA, 
               position = position_dodge(width=0.9))

ggplot(SCFA, aes(x=Treatment, y=perc, color=Treatment)) + 
  geom_jitter(position=position_jitter(0.2))+
  geom_boxplot() 
