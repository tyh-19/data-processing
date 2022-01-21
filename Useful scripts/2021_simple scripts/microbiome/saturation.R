setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/level_3_features/microbe/0420_subsample_DNA RNA comparison/")
saturation <- read.csv("Comparison_0420.csv")
saturation$sample_ID <- factor(saturation$sample_ID,levels=unique(saturation$sample_ID))

#pre-test
{
#read pairs
ggplot(saturation,aes(x=saturation$total_read_pairs,y=saturation$microbe_read,color = saturation$type, shape = saturation$sample_ID))+geom_point()+
  stat_smooth(formula = y ~log(x))+
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=20),
    axis.text.y = element_text(face="bold",  color="black", size=20),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  xlab("Total read pairs")+
  ylab("microbe read pairs")

#genus number
saturation$x <- as.factor(saturation$total_read_pairs)
saturation$type <- factor(saturation$type,levels = c("RNA","DNA-PRJEB28329","DNA-GSE81314"))
ggplot(saturation,aes(x=saturation$total_read_pairs,y=saturation$genus_num,color = saturation$type))+
  geom_point()+
  #geom_boxplot(aes(shape = saturation$x))+
  stat_smooth(formula = y ~log(x),span = 0.9)+
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=20),
    axis.text.y = element_text(face="bold",  color="black", size=20),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  xlab("Total read pairs")+
  ylab("Genus number")+
  scale_x_continuous(breaks = c(0,10000000,20000000,30000000),labels = c("0","10 M","20 M","30 M"))+
  scale_color_aaas()

#species number
ggplot(saturation,aes(x=saturation$total_read_pairs,y=saturation$specie_num.read.6.,color = saturation$type, shape = saturation$sample_ID))+geom_point()+
  stat_smooth(formula = y ~log(x))+
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=20),
    axis.text.y = element_text(face="bold",  color="black", size=20),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  xlab("Total read pairs")+
  ylab("Species number")

#genus number between cancer and health
RNA_summary <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/level_3_features/microbe/pico/stat_kraken_report_summary.txt",sep = "\t")
View(RNA_summary)
RNA_summary$type <- lapply(strsplit(as.character(RNA_summary$sample),"-",fixed=TRUE), function(x) x[1])
RNA_summary$type <- as.character(lapply(strsplit(as.character(RNA_summary$type),"_",fixed=TRUE), function(x) x[1]))

RNA_summary$type <- factor(RNA_summary$type,levels=c("NC","CRC","STAD","HCC","LUAD","ESCA"))
ggplot(RNA_summary,aes(x=RNA_summary$type,y=RNA_summary$genus_num,fill=type))+geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
  #geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
  geom_point(size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
  #geom_point(size = 3, position = position_dodge(width = 1))+
  theme_bw()+
  theme(#legend.position="right",
    #plot.margin = margin(2,2,2,50,unit="pt"),
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    legend.text= element_text(color="black",family = "Arial", size=24),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_text(color="black", size=24, vjust = 0),
    #axis.ticks.x = element_line(hjust=0),
    axis.ticks.length.x = unit(-0.2,"cm"),
    axis.text.y = element_text(color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  stat_compare_means(comparisons = list(c("NC","CRC"),c("NC","STAD"),c("NC","HCC"),c("NC","LUAD"),c("NC","ESCA")),
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     label = "p.signif")+
  xlab("")+
  ylab("Genus number")
}

#20210420 for mbRNA
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/level_3_features/microbe/0420_subsample_DNA RNA comparison/")
saturation <- read.csv("Comparison_0420.csv")
#saturation <- read.csv("./filtered/comparison_0518_2.csv")
#saturation <- read.csv("./filtered/comparison_0518_10.csv")
#saturation <- read.csv("./filtered/comparison_0518_50.csv")
saturation$sample_ID <- factor(saturation$sample_ID,levels=unique(saturation$sample_ID))

#genus number saturation plot
saturation$x <- as.factor(saturation$total_read_pairs)
saturation$type <- factor(saturation$type,levels = c("RNA","DNA-PRJEB28329","DNA-GSE81314"))
ggplot(saturation,aes(x=saturation$total_read_pairs,y=saturation$genus_num,color = saturation$type))+
  geom_point()+
  #geom_boxplot(aes(shape = saturation$x))+
  stat_smooth(formula = y ~log(x),span = 0.9)+
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=20),
    axis.text.y = element_text(face="bold",  color="black", size=20),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=18))+
  ylim(0,1600)+
  xlab("Total read pairs")+
  ylab("Genus number (read pair > 0)")+
  scale_x_continuous(breaks = c(0,10000000,20000000,30000000),labels = c("0","10 M","20 M","30 M"))+
  scale_color_jco()
}

#phylum composition
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/level_3_features/microbe/0420_phylum composition/RNA 263")
phylum <- read.csv("phylum-used_relative-abundance_no-chordata_summary_forplot.csv",header = TRUE)

stackplot <- melt(phylum)

stackplot$taxo <- factor(stackplot$taxo,levels = c("Proteobacteria","Firmicutes","Actinobacteria","Bacteroidetes","Cyanobacteria","Others","Virus"))
stackplot$variable <- factor(stackplot$variable,levels = c("CRC","STAD","HCC","LUAD","ESCA","NC"))

ggplot(stackplot,aes(x=stackplot$variable,y=stackplot$value,fill = stackplot$taxo)) + geom_bar(stat = "identity", width=0.5, col='black') +
  theme_bw()+
  theme(#legend.position="bottom",
    legend.position="right",
    panel.grid=element_blank(),
    legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=20,angle = 45,hjust = 1),
    axis.text.y = element_text(face="bold",  color="black", size=20),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  ylab("Microbe Fraction(%)")+
  xlab("")+
  #geom_vline(aes(xintercept=6.5))+
  scale_y_continuous(breaks = c(0,25,50,75,100),labels = c("0","25","50","75","100"),expand = c(0,0),limits = c(0,103))+
  scale_fill_jco(alpha = 0.8)


##diversity and evenness
library(vegan)
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/level_3_features/microbe/0420_subsample_DNA RNA comparison/diversity/")
DNA <- read.csv("30M_G_level_PRJEB28329_remove-homo_relative-abundance.txt",sep = "\t",header = TRUE, row.names = 1)
RNA <- read.csv("30M_G_level_RNA_subsample-remove_homo-relative_abundance.txt",sep = "\t",header = TRUE, row.names = 1)
RNA_pico <- read.csv("genus-remove_homo-relative_abudance.txt",sep = "\t",header = TRUE, row.names = 1)
#alpha
{
  DNA_diversity <- t(DNA)
  RNA_diversity <- t(RNA_pico)
  DNA_diversity[is.na(DNA_diversity)] <- 0
  RNA_diversity[is.na(RNA_diversity)] <- 0
  shannon_DNA <- diversity(DNA_diversity, index = "shannon", MARGIN = 1, base = exp(1))
  shannon_RNA <- diversity(RNA_diversity, index = "shannon", MARGIN = 1, base = exp(1))
  
  DNA_boxplot <- as.data.frame(shannon_DNA)
  DNA_boxplot$type <- "DNA"
  colnames(DNA_boxplot) <- c("Shannon diversity","type")
  RNA_boxplot <- as.data.frame(shannon_RNA)
  RNA_boxplot$type <- "RNA"
  colnames(RNA_boxplot) <- c("Shannon diversity","type")
  diversity_boxplot <- rbind(DNA_boxplot,RNA_boxplot)
  
  ggplot(diversity_boxplot,aes(x=type,y=`Shannon diversity`,fill = type))+geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    #geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    geom_point(size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
    #geom_point(size = 3, position = position_dodge(width = 1))+
    theme_bw()+
    theme(#legend.position="right",
      #plot.margin = margin(2,2,2,50,unit="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(color="black", size=24, vjust = 0),
      #axis.ticks.x = element_line(hjust=0),
      axis.ticks.length.x = unit(-0.2,"cm"),
      axis.text.y = element_text(color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    stat_compare_means(comparisons = list(c("DNA","RNA")),
                       method = "wilcox.test",
                       method.args = list(alternative = "greater"),
                       label = "p.signif")
  ##compere between cancer types
  diversity_boxplot$cancertype <- lapply(strsplit(as.character(rownames(diversity_boxplot)),".",fixed=TRUE), function(x) x[1])
  diversity_boxplot$cancertype <- as.character(lapply(strsplit(as.character(diversity_boxplot$cancertype),"_",fixed=TRUE), function(x) x[1]))
  diversity_boxplot$cancertype <- factor(diversity_boxplot$cancertype,levels=c("NC","CRC","STAD","HCC","LUAD","ESCA"))
  ggplot(diversity_boxplot[-which(diversity_boxplot$type=="DNA"),],aes(x=cancertype,y=`Shannon diversity`,fill=cancertype))+geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    #facet_wrap(~type)+
    #geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    geom_point(size = 1, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
    #geom_point(size = 3, position = position_dodge(width = 1))+
    theme_bw()+
    theme(#legend.position="right",
      #plot.margin = margin(2,2,2,50,unit="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_text(face="bold", color="transparent",family = "Arial", size=24),
      legend.text= element_text(color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(face="bold",color="black", size=24, hjust = 1,vjust  = 0.7,angle = 45),
      #axis.ticks.x = element_line(hjust=0),
      axis.ticks.length.x = unit(-0.2,"cm"),
      axis.text.y = element_text(face="bold",color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    stat_compare_means(comparisons = list(c("NC","CRC"),c("NC","STAD"),c("NC","HCC"),c("NC","LUAD"),c("NC","ESCA")),
                       method = "wilcox.test",
                       method.args = list(alternative = "greater"),
                       label = "p.signif")+
    xlab("")+
    ylab("Shannon diversity")+
    scale_fill_jco()
}

#specie number
{
  DNA_diversity <- t(DNA)
  RNA_diversity <- t(RNA_pico)
  DNA_diversity[is.na(DNA_diversity)] <- 0
  RNA_diversity[is.na(RNA_diversity)] <- 0
  DNA_S <- specnumber(DNA_diversity) ## rowSums(BCI > 0) does the same...
  RNA_S <- specnumber(RNA_diversity) ## rowSums(BCI > 0) does the same...
  DNA_boxplot_S <- as.data.frame(DNA_S)
  DNA_boxplot_S$type <- "DNA"
  colnames(DNA_boxplot_S) <- c("Genus Number","type")
  RNA_boxplot_S <- as.data.frame(RNA_S)
  RNA_boxplot_S$type <- "RNA"
  colnames(RNA_boxplot_S) <- c("Genus Number","type")
  genus_boxplot <- rbind(DNA_boxplot_S,RNA_boxplot_S)
  
  ggplot(genus_boxplot,aes(x=type,y=`Genus Number`,fill = type))+geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    #geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    geom_point(size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
    #geom_point(size = 3, position = position_dodge(width = 1))+
    theme_bw()+
    theme(#legend.position="right",
      #plot.margin = margin(2,2,2,50,unit="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(color="black", size=24, vjust = 0),
      #axis.ticks.x = element_line(hjust=0),
      axis.ticks.length.x = unit(-0.2,"cm"),
      axis.text.y = element_text(color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    stat_compare_means(comparisons = list(c("DNA","RNA")),
                       method = "wilcox.test",
                       method.args = list(alternative = "less"),
                       label = "p.signif")+
    scale_fill_jco(alpha = 0.8)
  
  ##compere between cancer types
  genus_boxplot$cancertype <- lapply(strsplit(as.character(rownames(genus_boxplot)),".",fixed=TRUE), function(x) x[1])
  genus_boxplot$cancertype <- as.character(lapply(strsplit(as.character(genus_boxplot$cancertype),"_",fixed=TRUE), function(x) x[1]))
  genus_boxplot$cancertype <- factor(genus_boxplot$cancertype,levels=c("NC","CRC","STAD","HCC","LUAD","ESCA"))
  ggplot(genus_boxplot[-which(genus_boxplot$type=="DNA"),],aes(x=cancertype,y=`Genus Number`,fill=cancertype))+geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    #facet_wrap(~type)+
    #geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    geom_point(size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
    #geom_point(size = 3, position = position_dodge(width = 1))+
    theme_bw()+
    theme(#legend.position="right",
      #plot.margin = margin(2,2,2,50,unit="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(color="black", size=24, vjust = 0),
      #axis.ticks.x = element_line(hjust=0),
      axis.ticks.length.x = unit(-0.2,"cm"),
      axis.text.y = element_text(color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    stat_compare_means(comparisons = list(c("NC","CRC"),c("NC","STAD"),c("NC","HCC"),c("NC","LUAD"),c("NC","ESCA")),
                       method = "wilcox.test",
                       method.args = list(alternative = "greater"),
                       label = "p.signif")+
    xlab("")+
    ylab("Genus Number")+
    scale_fill_jco(alpha = 0.8)
}

#eveness index
{
  #compare between DNA and RNA
  DNA_S <- specnumber(DNA_diversity) ## rowSums(BCI > 0) does the same...
  DNA_H <- diversity(DNA_diversity)
  DNA_J <- DNA_H/log(DNA_S)
  RNA_S <- specnumber(RNA_diversity) ## rowSums(BCI > 0) does the same...
  RNA_H <- diversity(RNA_diversity)
  RNA_J <- RNA_H/log(RNA_S)
  
  DNA_boxplot_J <- as.data.frame(DNA_J)
  DNA_boxplot_J$type <- "DNA"
  colnames(DNA_boxplot_J) <- c("Pielou's evenness (J)","type")
  RNA_boxplot_J <- as.data.frame(RNA_J)
  RNA_boxplot_J$type <- "RNA"
  colnames(RNA_boxplot_J) <- c("Pielou's evenness (J)","type")
  eveness_boxplot <- rbind(DNA_boxplot_J,RNA_boxplot_J)
  
  ggplot(eveness_boxplot,aes(x=type,y=`Pielou's evenness (J)`,fill = type))+geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    #geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    geom_point(size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
    #geom_point(size = 3, position = position_dodge(width = 1))+
    theme_bw()+
    theme(#legend.position="right",
      #plot.margin = margin(2,2,2,50,unit="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(color="black", size=24, vjust = 0),
      #axis.ticks.x = element_line(hjust=0),
      axis.ticks.length.x = unit(-0.2,"cm"),
      axis.text.y = element_text(color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    stat_compare_means(comparisons = list(c("DNA","RNA")),
                       method = "wilcox.test",
                       method.args = list(alternative = "greater"),
                       label = "p.signif")
  
  ##compere between cancer types
  eveness_boxplot$cancertype <- lapply(strsplit(as.character(rownames(eveness_boxplot)),".",fixed=TRUE), function(x) x[1])
  eveness_boxplot$cancertype <- as.character(lapply(strsplit(as.character(eveness_boxplot$cancertype),"_",fixed=TRUE), function(x) x[1]))
  eveness_boxplot$cancertype <- factor(eveness_boxplot$cancertype,levels=c("NC","CRC","STAD","HCC","LUAD","ESCA"))
  ggplot(eveness_boxplot[-which(eveness_boxplot$type=="DNA"),],aes(x=cancertype,y=`Pielou's evenness (J)`,fill=cancertype))+geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    #facet_wrap(~type)+
    #geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
    geom_point(size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
    #geom_point(size = 3, position = position_dodge(width = 1))+
    theme_bw()+
    theme(#legend.position="right",
      #plot.margin = margin(2,2,2,50,unit="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
      legend.text= element_text(color="black",family = "Arial", size=24),
      plot.title = element_text(hjust = 0.5,size=36,face="bold"),
      axis.text.x = element_text(color="black", size=24, vjust = 0),
      #axis.ticks.x = element_line(hjust=0),
      axis.ticks.length.x = unit(-0.2,"cm"),
      axis.text.y = element_text(color="black", size=24),
      axis.title.x = element_text(face="bold", color="black", size=24),
      axis.title.y = element_text(face="bold",color="black", size=24))+
    stat_compare_means(comparisons = list(c("NC","CRC"),c("NC","STAD"),c("NC","HCC"),c("NC","LUAD"),c("NC","ESCA")),
                       method = "wilcox.test",
                       method.args = list(alternative = "greater"),
                       label = "p.signif")+
    xlab("")+
    ylab("Pielou's evenness (J)")
}

#subsample alpha
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/level_3_features/microbe/0420_subsample_DNA RNA comparison/diversity/")
DNA <- read.csv("microbe_G_level_PRJEB28329_remove-homo_relative-abundance.txt",sep = "\t",header = TRUE, row.names = 1)
RNA <- read.csv("microbe_G_level_RNA_subsample-remove_homo-relative_abundance.txt",sep = "\t",header = TRUE, row.names = 1)

DNA_diversity <- t(DNA)
RNA_diversity <- t(RNA)
DNA_diversity[is.na(DNA_diversity)] <- 0
RNA_diversity[is.na(RNA_diversity)] <- 0
shannon_DNA <- diversity(DNA_diversity, index = "shannon", MARGIN = 1, base = exp(1))
shannon_RNA <- diversity(RNA_diversity, index = "shannon", MARGIN = 1, base = exp(1))

DNA_boxplot <- as.data.frame(shannon_DNA)
DNA_boxplot$type <- "DNA"
colnames(DNA_boxplot) <- c("Shannon diversity","type")
RNA_boxplot <- as.data.frame(shannon_RNA)
RNA_boxplot$type <- "RNA"
colnames(RNA_boxplot) <- c("Shannon diversity","type")
diversity_boxplot <- rbind(DNA_boxplot,RNA_boxplot)

ggplot(diversity_boxplot,aes(x=type,y=`Shannon diversity`,fill = type))+geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
  #geom_boxplot(alpha = 0, size = 1, position = position_dodge(1.1),outlier.size=0)+
  geom_point(size = 3, position = position_jitterdodge(dodge.width=1,jitter.width = 0.5))+
  #geom_point(size = 3, position = position_dodge(width = 1))+
  theme_bw()+
  theme(#legend.position="right",
    #plot.margin = margin(2,2,2,50,unit="pt"),
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(face="bold", color="black",family = "Arial", size=24),
    legend.text= element_text(color="black",family = "Arial", size=24),
    plot.title = element_text(hjust = 0.5,size=36,face="bold"),
    axis.text.x = element_text(color="black", size=24, vjust = 0),
    #axis.ticks.x = element_line(hjust=0),
    axis.ticks.length.x = unit(-0.2,"cm"),
    axis.text.y = element_text(color="black", size=24),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))

#alpha saturation plot
diversity_boxplot$x <- gsub("X","",as.character(lapply(strsplit(as.character(rownames(diversity_boxplot)),".",fixed=TRUE), function(x) x[1])))
diversity_boxplot$total_read_pairs <- as.numeric(diversity_boxplot$x)

diversity_boxplot$type <- factor(diversity_boxplot$type,levels = c("RNA","DNA"))
ggplot(diversity_boxplot,aes(x=diversity_boxplot$total_read_pairs,y=diversity_boxplot$`Shannon diversity`,color = diversity_boxplot$type))+
  #geom_point()+
  #geom_boxplot(aes(shape = diversity_boxplot$x))+
  stat_smooth(formula = y ~log(x),span = 100)+
  theme_bw()+
  theme(#legend.position="right",
    legend.position="right",
    panel.grid=element_blank(),
    legend.title = element_text(face="bold", color="transparent",family = "Arial", size=20),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.text.x = element_text(face="bold", color="black", size=20),
    axis.text.y = element_text(face="bold",  color="black", size=20),
    axis.title.x = element_text(face="bold", color="black", size=24),
    axis.title.y = element_text(face="bold",color="black", size=24))+
  xlab("Total read pairs")+
  ylab("Shannon diversity")+
  scale_x_continuous(breaks = c(0,5000000,20000000,30000000),labels = c("0","5 M","20 M","30 M"),limits = c(0,30000000))

##20210424 subtype barplot for mbRNA
##barplot

#filped_barplot
library(ggplot2)
{
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/level_3_features/microbe/0424_BARPLOT/results/")
bar <- read.csv("CRC_right_vs_left.txt",sep = "\t",header=T,row.names = 1)
head(bar)
bar$cancertype <- c(rep("CRC", nrow(bar)))
bar$color <- c(rep("#f2be45",nrow(bar)))
#bar$cancertype <- c(rep("STAD", nrow(bar)))
#bar$color <- c(rep("#1685a9",nrow(bar)))
#bar$cancertype <- c(rep("HCC", nrow(bar)))
#bar$color <- c(rep("#be002f",nrow(bar)))
#bar$cancertype <- c(rep("LUAD", nrow(bar)))
#bar$color <- c(rep("#006633",nrow(bar)))

bar <- bar[bar$padj<0.05,] # for CRC colon vs rectum. HCC CHB positive vs negative, AFP+ vs AFP-, LUAD smoker vs never smoking
#bar <- bar[bar$pvalue < 0.05,] #only for CRC right vs left 
positive_top10 <- bar[bar$log2FoldChange>1,]
positive_top10 <- head(positive_top10[order(positive_top10$log2FoldChange,decreasing = TRUE),],10)

negative_top10 <- bar[bar$log2FoldChange < -1,]
negative_top10 <- head(negative_top10[order(negative_top10$log2FoldChange,decreasing = TRUE),],10)

bar_forplot <- rbind(negative_top10,positive_top10)

bar_forplot$genus_name <- as.character(lapply(strsplit(rownames(bar_forplot),"|",fixed = TRUE),function(x) x[2]))

bar_forplot$positive_label <- bar_forplot$genus_name
bar_forplot$negative_label <- bar_forplot$genus_name
bar_forplot[which(bar_forplot$log2FoldChange< -1),]$positive_label <- NA
bar_forplot[which(bar_forplot$log2FoldChange> 1),]$negative_label <- NA

#level change
bar_forplot$genus_name <- factor(bar_forplot$genus_name, levels = bar_forplot$genus_name)
bar_plot <- ggplot(data=bar_forplot,aes(x=genus_name,y=log2FoldChange))+
  geom_bar(stat="identity",width=0.8,fill=bar_forplot$color)+
  #geom_bar(stat="identity",width=0.8,fill=-1*log10(bar_forplot$padj))+
  #scale_fill_jco()+
  #scale_fill_continuous(low = "white", high = "#1685a9")+
  labs(x=NULL,
       y="Log2Foldchange",
       fill = expression(-log[10](P.value)))+
  #geom_vline(aes(xintercept=10.5))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position = "topright",
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black", size = 14 ,face = "bold"),
        axis.text.y = element_blank(),
        axis.title = element_text(color = "black", size = 18, face = "bold"))+
  coord_flip()
bar_plot
bar_plot+
  geom_text(aes(bar_forplot$genus_name,y=-0.1,label=bar_forplot$positive_label),hjust=1,vjust=0.4,colour="black",size=6)+
  geom_text(aes(bar_forplot$genus_name,y=0.1,label=bar_forplot$negative_label),hjust=0,vjust=0.4,colour="black",size=6)+
  ##CRC colon vs rectum
  geom_label(aes(x=6,y=-14,label="CRC "),hjust=0,vjust=0,colour="#f2be45",family="Arial",size=8)+
  geom_label(aes(x=1,y=-14,label="Left sided"),hjust=0,vjust=0,colour="black",family="Arial",size=6)+
  geom_label(aes(x=7,y=15,label="Right sided"),hjust=1,vjust=0,colour="transparent",family="Arial",size=6)
  
  ##STAD HER2 positive vs negative
  #geom_label(aes(x=16,y=-14,label="STAD "),hjust=0,vjust=0,colour="#1685a9",family="Arial",size=8)+
  #geom_label(aes(x=1,y=-14,label="HER2 negative"),hjust=0,vjust=0,colour="black",family="Arial",size=6)+
  #geom_label(aes(x=16,y=10,label="HER2 positive"),hjust=1,vjust=0,colour="black",family="Arial",size=6)
  
  ##HCC AFP high vs low
  #geom_label(aes(x=16,y=-14,label="HCC "),hjust=0,vjust=0,colour="#be002f",family="Arial",size=8)+
  #geom_label(aes(x=1,y=-14,label="AFP-"),hjust=0,vjust=0,colour="black",family="Arial",size=6)+
  #geom_label(aes(x=16,y=10,label="AFP+"),hjust=1,vjust=0,colour="black",family="Arial",size=6)

  ##HCC HBV+ vs HBV-
  #geom_label(aes(x=19,y=-14,label="HCC "),hjust=0,vjust=0,colour="#be002f",family="Arial",size=8)+
  #geom_label(aes(x=1,y=-14,label="HBV-"),hjust=0,vjust=0,colour="black",family="Arial",size=6)+
  #geom_label(aes(x=19,y=10,label="HBV+"),hjust=1,vjust=0,colour="black",family="Arial",size=6)

  ##LUAD Smoker vs Non smoker
  #geom_label(aes(x=19,y=-14,label="LUAD "),hjust=0,vjust=0,colour="#006633",family="Arial",size=8)+
  #geom_label(aes(x=1,y=-14,label="Non smoking"),hjust=0,vjust=0,colour="black",family="Arial",size=6)+
  #geom_label(aes(x=19,y=10,label="Smoker"),hjust=1,vjust=0,colour="black",family="Arial",size=6)
} 
  