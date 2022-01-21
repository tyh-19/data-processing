pico_sample_lingshi <- read.csv("C:/Users/Tao/Desktop/远程办公文件/lulab 2 (1)/lulab/Liquid Biopsy/01.sample info/2020.03.27-pico-sample/pico-sample.csv",header=TRUE)
head(pico_sample,10)
class(pico_sample)
pico_sample$RNA_mass_divide5 <- as.numeric(pico_sample$RNA_mass.ng./5)
pico_sample$RNA_mass_log2 <- log(as.numeric(pico_sample$RNA_mass.ng.),2)
pico_sample$RNA_mass_log2
library(ggplot2)
##肿瘤种类柱状???
pico_sample$Disease <- factor(pico_sample$Disease, levels=c("ESCA","LUAD","HCC","STAD","CRC","Health"), ordered=TRUE)
ggplot(pico_sample)+geom_bar(aes(x=Disease,fill=Stage_summary),position="fill")

##肿瘤种类散点???,不太???
pico_sample$Library_yield <- as.numeric(as.character(pico_sample$Library_yield))
pico_sample$Library_yield_log2 <- as.numeric(as.character(pico_sample$Library_yield_log2))
pico_sample$RNA_mass.ng. <- as.numeric(as.character(pico_sample$RNA_mass.ng.))
pico_sample$PCR_cycles <- as.factor(pico_sample_lingshi$PCR_cycles)
ggplot(pico_sample,aes(x=RNA_mass.ng.,y=Library_yield))+geom_point(aes(colour=Disease))+facet_wrap(~Source,scales="free")

##线性回归
#一元
plot(pico_sample$RNA_mass.ng.,pico_sample$Library_yield,main = "Library_yield ~ RNA_mass")
lm.reg <- lm(pico_sample$Library_yield~pico_sample$RNA_mass.ng.)
lm.reg
abline(lm.reg)
anova(lm.reg)

#多元
library(psych)
pico_sample$PCR_cycles <- as.numeric(pico_sample$PCR_cycles)+13
pairs(pico_sample[c("Library_yield_log2","RNA_mass_log2","PCR_cycles")])
pairs.panels(pico_sample[c("Library_yield_log2","RNA_mass_log2","PCR_cycles")])
model <- lm(Library_yield_log2~RNA_mass_log2+PCR_cycles,data=pico_sample)
summary(model)
model
##直方图，信息量和bar差不???
pico_sample_nohealth <- read.csv("C:/Users/Tao/Desktop/pico-sample -withoutHealth.csv",header=TRUE)
ggplot(pico_sample_nohealth)+geom_histogram(aes(x=Disease,fill=Stage_summary),stat="count",position="dodge")+theme_bw()

##风琴玫瑰???
pico_sample$Disease <- as.character(pico_sample$Disease)
head(pico_sample$Disease)
ggplot(pico_sample)+
geom_bar(aes(x=Disease,fill=Stage_summary))+
coord_polar(direction=1)+
theme_bw()

##密度???
ggplot(pico_sample)+geom_density(aes(x=Age,colour=Disease),size=1.5)+theme_light()

ggplot(pico_sample)+geom_boxplot(aes(x=Disease,y=RNA_mass.ng.,colour = Source),position = "dodge")

## boxplot
compare_means(value ~ Treament, data = pico_sample, group.by = "variable")
ggplot(pico_sample)+geom_boxplot(aes(x=PCR_cycles,y=Library_yield,colour=Disease))+theme_light()

## ggpubr kruskal-wallis秩和检验
install.packages("ggpubr")
library(ggpubr)
my_comparisons <- list( c("14", "15"), c("14", "16"), c("15", "16") )
options(repr.plot.width=4, repr.plot.height=4)
ggplot(pico_sample, aes(x=as.character(PCR_cycles), y=Library_yield,fill=PCR_cycles)) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA) +
  #geom_jitter(width = 0.3, size=0.01) +# , aes(color=supp) +
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 2600, label.x = 1.6) +
  theme_bw()
