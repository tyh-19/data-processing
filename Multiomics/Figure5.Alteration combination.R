
setwd("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20211114_yunfan/")

#alteration number
{
  titration <- read.csv("cv-feature-number.txt",header = TRUE,sep = "\t")
  RNA_color <- c("Alternative promoter"="#E9C2A6",
                 "RNA Expression"="#A5435C",
                 "RNA Splicing"="#C5E3BF",
                 "Alternative polyadenyltion"="#003F87",
                 "Chimeric RNA"="#FF3D0D",
                 "RNA editing"="#324F17",
                 "Allele specific expression"="#800080",
                 "RNA SNP"="#333333",
                 "DNA copy number"="#FFE234",
                 "DNA nucleosome occupancy"="#44A5F9",
                 "DNA methylation"="#C7C2C2",
                 "Combined"="#FFFFFF",
                 "miRNA"="#000000")
  
  
  titration$data <- gsub("expression","RNA Expression",titration$data)
  titration$data <- gsub("alt-promoter","Alternative promoter",titration$data)
  titration$data <- gsub("APA","Alternative polyadenyltion",titration$data)
  titration$data <- gsub("ASE","Allele specific expression",titration$data)
  titration$data <- gsub("chimera","Chimeric RNA",titration$data)
  titration$data <- gsub("CNV","DNA copy number",titration$data)
  titration$data <- gsub("editing","RNA editing",titration$data)
  titration$data <- gsub("MeDIP","DNA methylation",titration$data)
  titration$data <- gsub("nucleasome","DNA nucleosome occupancy",titration$data)
  titration$data <- gsub("SNP","RNA SNP",titration$data)
  titration$data <- gsub("splicing","RNA Splicing",titration$data)
  titration$data <- gsub("combined","Combined",titration$data)

  titration$X. <- factor(as.character(titration$X.),levels = c("10","50","100","150","200"))
  ggplot(titration[(titration$comparison=="STAD-NC")&(titration$metric=="AUROC"),],aes(x=X.,y=value,color=data))+
    geom_point(position = position_dodge(width = 0.2),size = 3)+
    geom_line(aes(group = data),position = position_dodge(width = 0.2))+
    #geom_bar(aes(x=data.type,y=mean_early-0.3),stat = "identity",colour = "black",fill=alpha("white",0.4))+
    #geom_text(aes(x=data.type,y=mean-0.3+0.08,label=round(mean,digits = 3)),size = 4,angle = 0)+
    #scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels = c("0.3","0.4","0.5","0.6","0.7","0.8","0.9","1 "),expand = c(0,0),limits = c(0,0.7))+
    scale_color_manual(values = RNA_color)+
    xlab("")+
    ylab("AUROC")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(20,20,10,80),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.ticks.x = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))
  }

#all stage ROC barplot(without clinical stratify)
{
  plot <- read.csv("cv-individual-and-all.txt",header = TRUE,sep = "\t")
  sparseness <- read.csv("alteration_sparseness.csv",header = TRUE,row.names = 1)
  {
  sparseness$alteration <- gsub("Expression","RNA Expression",sparseness$alteration)
  sparseness$alteration <- gsub("Altpromoter","Alternative promoter",sparseness$alteration)
  sparseness$alteration <- gsub("APA","Alternative polyadenyltion",sparseness$alteration)
  sparseness$alteration <- gsub("ASE","Allele specific expression",sparseness$alteration)
  sparseness$alteration <- gsub("chimeric","Chimeric RNA",sparseness$alteration)
  sparseness$alteration <- gsub("CNV-CPM","DNA copy number",sparseness$alteration)
  sparseness$alteration <- gsub("EDIT","RNA editing",sparseness$alteration)
  sparseness$alteration <- gsub("Medip","DNA methylation",sparseness$alteration)
  sparseness$alteration <- gsub("Nucleosome","DNA nucleosome occupancy",sparseness$alteration)
  sparseness$alteration <- gsub("SNP","RNA SNP",sparseness$alteration)
  sparseness$alteration <- gsub("Splicing","RNA Splicing",sparseness$alteration)
  sparseness$alteration <- gsub("miRNA","miRNA",sparseness$alteration)
  }
  
  RNA_color <- c("Alternative promoter"="#E9C2A6",
                 "RNA Expression"="#A5435C",
                 "RNA Splicing"="#C5E3BF",
                 "Alternative polyadenyltion"="#003F87",
                 "Chimeric RNA"="#FF3D0D",
                 "RNA editing"="#324F17",
                 "Allele specific expression"="#800080",
                 "RNA SNP"="#333333",
                 "DNA copy number"="#FFE234",
                 "DNA nucleosome occupancy"="#44A5F9",
                 "DNA methylation"="#C7C2C2",
                 "Combined-all"="#FFFFFF",
                 "Combined-DNA"="#FFFFFF",
                 "Combined-RNA"="#FFFFFF",
                 "miRNA"="#000000")
  
  {
  plot$data.type <- gsub("expression","RNA Expression",plot$data.type)
  plot$data.type <- gsub("alt-promoter","Alternative promoter",plot$data.type)
  plot$data.type <- gsub("APA","Alternative polyadenyltion",plot$data.type)
  plot$data.type <- gsub("ASE","Allele specific expression",plot$data.type)
  plot$data.type <- gsub("chimera","Chimeric RNA",plot$data.type)
  plot$data.type <- gsub("CNV","DNA copy number",plot$data.type)
  plot$data.type <- gsub("editing","RNA editing",plot$data.type)
  plot$data.type <- gsub("MeDIP","DNA methylation",plot$data.type)
  plot$data.type <- gsub("nucleasome","DNA nucleosome occupancy",plot$data.type)
  plot$data.type <- gsub("SNP","RNA SNP",plot$data.type)
  plot$data.type <- gsub("splicing","RNA Splicing",plot$data.type)
  plot$data.type <- gsub("combined-all","Combined-all",plot$data.type)
  plot$data.type <- gsub("combined-DNA","Combined-DNA",plot$data.type)
  plot$data.type <- gsub("combined-RNA","Combined-RNA",plot$data.type)
  }
  
  {
  #plot$value <- as.numeric(as.character(plot$value))
  comparison<- "STAD-NC"
  #ROC_STAD
  {
  rm(plot_at50feature)
  plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="AUROC"),] 
  plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
  colnames(plot_at50feature_mean) <- c("data.type","mean")
  plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
  colnames(plot_at50feature_sd) <- c("data.type","sd")
  plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
  plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
  plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
  plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
  plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
  plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "miRNA"),]
  plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-all"),]
  plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-DNA"),]
  plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-RNA"),]
  
  #AUC
  order <- plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"]
  #order <- order[-which(order=="Combined-all")]
  #order <- order[-which(order=="Combined-DNA")]
  #order <- order[-which(order=="Combined-RNA")]
  #order <- order[-which(order=="miRNA")]
  #order <- append("Combined-RNA",order)
  #order <- append("Combined-DNA",order)
  #order <- append("Combined-all",order)
  #order <- append(order,"miRNA")
  plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
  bar <- ggplot(plot_at50feature,aes(x=data.type,y=mean,fill=data.type))+
    geom_bar(stat = "identity",colour = "black")+
    geom_text(aes(x=data.type,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
    geom_errorbar(aes(x=data.type, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
    xlab("")+
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
    scale_fill_manual(values = RNA_color)+
    #geom_vline(xintercept=1.5,linetype = 2)+
    #geom_vline(xintercept=3.5,linetype = 2)+
    xlab("")+
    ylab("AUROC")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(15,5,-10,45),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.ticks.x = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))
  
  sparseness_bar <- 
    ggplot(plot_at50feature,aes(x=data.type,y=0,fill=X1.sparseness))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "black")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-10,10,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  ggarrange(bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.56),legend = "none")
  ROC_STAD <- plot_at50feature
  ROC_STAD$comparison <- "STAD"
  }
  #recall at 95% specificity STAD
  {
  rm(plot_at50feature)
  plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="recall"),] 
  plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
  colnames(plot_at50feature_mean) <- c("data.type","mean")
  plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
  colnames(plot_at50feature_sd) <- c("data.type","sd")
  plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
  plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
  plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
  plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
  plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
  plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "miRNA"),]
  plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-all"),]
  plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-DNA"),]
  plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-RNA"),]
  
  order <- plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"]
  #order <- order[-which(order=="Combined-all")]
  #order <- order[-which(order=="Combined-DNA")]
  #order <- order[-which(order=="Combined-RNA")]
  #order <- order[-which(order=="miRNA")]
  #order <- append("Combined-RNA",order)
  #order <- append("Combined-DNA",order)
  #order <- append("Combined-all",order)
  #order <- append(order,"miRNA")
  plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
  recall_bar <- ggplot(plot_at50feature,aes(x=data.type,y=mean,fill=data.type))+
    geom_bar(stat = "identity",colour = "black")+
    geom_text(aes(x=data.type,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
    geom_errorbar(aes(x=data.type, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
    xlab("")+
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
    scale_fill_manual(values = RNA_color)+
    #geom_vline(xintercept=1.5,linetype = 2)+
    #geom_vline(xintercept=3.5,linetype = 2)+
    xlab("")+
    ylab("Sensitivity at 95% specificity")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(15,5,-10,45),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.ticks.x = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_text(face="bold", color="black", size=20),
      axis.title.y = element_text(face="bold",color="black", size=20))
  
  sparseness_bar <- 
    ggplot(plot_at50feature,aes(x=data.type,y=0,fill=X1.sparseness))+
    geom_tile()+
    scale_fill_gradient(low="white",high = "black")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-10,10,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  ggarrange(recall_bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.56),legend = "none")
  recall_STAD <- plot_at50feature
  recall_STAD$comparison <- "STAD"
  }
  }
  {
    #plot$value <- as.numeric(as.character(plot$value))
    comparison<- "CRC-NC"
    #ROC_STAD
    {
      rm(plot_at50feature)
      plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="AUROC"),] 
      plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
      colnames(plot_at50feature_mean) <- c("data.type","mean")
      plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
      colnames(plot_at50feature_sd) <- c("data.type","sd")
      plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
      plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
      plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "miRNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-all"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-DNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-RNA"),]
      
      #AUC
      order <- plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"]
      #order <- order[-which(order=="Combined-all")]
      #order <- order[-which(order=="Combined-DNA")]
      #order <- order[-which(order=="Combined-RNA")]
      #order <- order[-which(order=="miRNA")]
      #order <- append("Combined-RNA",order)
      #order <- append("Combined-DNA",order)
      #order <- append("Combined-all",order)
      #order <- append(order,"miRNA")
      plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
      bar <- ggplot(plot_at50feature,aes(x=data.type,y=mean,fill=data.type))+
        geom_bar(stat = "identity",colour = "black")+
        geom_text(aes(x=data.type,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
        geom_errorbar(aes(x=data.type, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
        xlab("")+
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
        scale_fill_manual(values = RNA_color)+
        #geom_vline(xintercept=1.5,linetype = 2)+
        #geom_vline(xintercept=3.5,linetype = 2)+
        xlab("")+
        ylab("AUROC")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(15,5,-10,45),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
      
      sparseness_bar <- 
        ggplot(plot_at50feature,aes(x=data.type,y=0,fill=X1.sparseness))+
        geom_tile()+
        scale_fill_gradient(low="white",high = "black")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(-10,10,10,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      ggarrange(bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.56),legend = "none")
      ROC_CRC <- plot_at50feature
      ROC_CRC$comparison <- "CRC"
    }
    #recall at 95% specificity STAD
    {
      rm(plot_at50feature)
      plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="recall"),] 
      plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
      colnames(plot_at50feature_mean) <- c("data.type","mean")
      plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
      colnames(plot_at50feature_sd) <- c("data.type","sd")
      plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
      plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
      plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "miRNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-all"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-DNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-RNA"),]
      
      order <- plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"]
      #order <- order[-which(order=="Combined-all")]
      #order <- order[-which(order=="Combined-DNA")]
      #order <- order[-which(order=="Combined-RNA")]
      #order <- order[-which(order=="miRNA")]
      #order <- append("Combined-RNA",order)
      #order <- append("Combined-DNA",order)
      #order <- append("Combined-all",order)
      #order <- append(order,"miRNA")
      plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
      recall_bar <- ggplot(plot_at50feature,aes(x=data.type,y=mean,fill=data.type))+
        geom_bar(stat = "identity",colour = "black")+
        geom_text(aes(x=data.type,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
        geom_errorbar(aes(x=data.type, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
        xlab("")+
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
        scale_fill_manual(values = RNA_color)+
        #geom_vline(xintercept=1.5,linetype = 2)+
        #geom_vline(xintercept=3.5,linetype = 2)+
        xlab("")+
        ylab("Sensitivity at 95% specificity")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(15,5,-10,45),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
      
      sparseness_bar <- 
        ggplot(plot_at50feature,aes(x=data.type,y=0,fill=X1.sparseness))+
        geom_tile()+
        scale_fill_gradient(low="white",high = "black")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(-10,10,10,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      ggarrange(recall_bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.56),legend = "none")
      recall_CRC <- plot_at50feature
      recall_CRC$comparison <- "CRC"
    }
  }
  {
    #plot$value <- as.numeric(as.character(plot$value))
    comparison<- "cancer-NC"
    #ROC_STAD
    {
      rm(plot_at50feature)
      plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="AUROC"),] 
      plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
      colnames(plot_at50feature_mean) <- c("data.type","mean")
      plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
      colnames(plot_at50feature_sd) <- c("data.type","sd")
      plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
      plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
      plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "miRNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-all"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-DNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-RNA"),]
      
      #AUC
      order <- plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"]
      #order <- order[-which(order=="Combined-all")]
      #order <- order[-which(order=="Combined-DNA")]
      #order <- order[-which(order=="Combined-RNA")]
      #order <- order[-which(order=="miRNA")]
      #order <- append("Combined-RNA",order)
      #order <- append("Combined-DNA",order)
      #order <- append("Combined-all",order)
      #order <- append(order,"miRNA")
      plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
      bar <- ggplot(plot_at50feature,aes(x=data.type,y=mean,fill=data.type))+
        geom_bar(stat = "identity",colour = "black")+
        geom_text(aes(x=data.type,y=mean+0.08,label=round(mean,digits = 3)),size = 6,angle = 0)+
        geom_errorbar(aes(x=data.type, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
        xlab("")+
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
        scale_fill_manual(values = RNA_color)+
        #geom_vline(xintercept=1.5,linetype = 2)+
        #geom_vline(xintercept=3.5,linetype = 2)+
        xlab("")+
        ylab("AUROC")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(15,5,-10,78),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
      
      sparseness_bar <- 
        ggplot(plot_at50feature,aes(x=data.type,y=0,fill=X1.sparseness))+
        geom_tile()+
        scale_fill_gradient(low="white",high = "black")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(-10,10,10,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      ggarrange(bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.56),legend = "none")
      ROC_cancer <- plot_at50feature
      ROC_cancer$comparison <- "GIC"
    }
    #recall at 95% specificity STAD
    {
      rm(plot_at50feature)
      plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="recall"),] 
      plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
      colnames(plot_at50feature_mean) <- c("data.type","mean")
      plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
      colnames(plot_at50feature_sd) <- c("data.type","sd")
      plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
      plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
      plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "miRNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-all"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-DNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-RNA"),]
      
      order <- plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"]
      #order <- order[-which(order=="Combined-all")]
      #order <- order[-which(order=="Combined-DNA")]
      #order <- order[-which(order=="Combined-RNA")]
      #order <- order[-which(order=="miRNA")]
      #order <- append("Combined-RNA",order)
      #order <- append("Combined-DNA",order)
      #order <- append("Combined-all",order)
      #order <- append(order,"miRNA")
      plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
      recall_bar <- ggplot(plot_at50feature,aes(x=data.type,y=mean,fill=data.type))+
        geom_bar(stat = "identity",colour = "black")+
        geom_text(aes(x=data.type,y=mean+0.05,label=round(mean,digits = 3)),size = 6,angle = 0)+
        geom_errorbar(aes(x=data.type, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
        xlab("")+
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
        scale_fill_manual(values = RNA_color)+
        #geom_vline(xintercept=1.5,linetype = 2)+
        #geom_vline(xintercept=3.5,linetype = 2)+
        xlab("")+
        ylab("Sensitivity at 95% specificity")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(15,5,-10,78),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
      
      sparseness_bar <- 
        ggplot(plot_at50feature,aes(x=data.type,y=0,fill=X1.sparseness))+
        geom_tile()+
        scale_fill_gradient(low="white",high = "black")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(-10,10,10,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      ggarrange(recall_bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.56),legend = "none")
      recall_cancer <- plot_at50feature
      recall_cancer$comparison <- "GIC"
    }
  }
  {
    #plot$value <- as.numeric(as.character(plot$value))
    comparison<- "CRC-STAD"
    #ROC_STAD
    {
      rm(plot_at50feature)
      plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="AUROC"),] 
      plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
      colnames(plot_at50feature_mean) <- c("data.type","mean")
      plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
      colnames(plot_at50feature_sd) <- c("data.type","sd")
      plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
      plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
      plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "miRNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-all"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-DNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-RNA"),]
      
      #AUC
      order <- plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"]
      #order <- order[-which(order=="Combined-all")]
      #order <- order[-which(order=="Combined-DNA")]
      #order <- order[-which(order=="Combined-RNA")]
      #order <- order[-which(order=="miRNA")]
      #order <- append("Combined-RNA",order)
      #order <- append("Combined-DNA",order)
      #order <- append("Combined-all",order)
      #order <- append(order,"miRNA")
      plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
      bar <- ggplot(plot_at50feature,aes(x=data.type,y=mean,fill=data.type))+
        geom_bar(stat = "identity",colour = "black")+
        geom_text(aes(x=data.type,y=mean+0.08,label=round(mean,digits = 3)),size = 6,angle = 0)+
        geom_errorbar(aes(x=data.type, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
        xlab("")+
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
        scale_fill_manual(values = RNA_color)+
        #geom_vline(xintercept=1.5,linetype = 2)+
        #geom_vline(xintercept=3.5,linetype = 2)+
        xlab("")+
        ylab("AUROC")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(15,5,-10,78),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
      
      sparseness_bar <- 
        ggplot(plot_at50feature,aes(x=data.type,y=0,fill=X1.sparseness))+
        geom_tile()+
        scale_fill_gradient(low="white",high = "black")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(-10,10,10,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      ggarrange(bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.56),legend = "none")
    }
    #recall at 95% specificity STAD
    {
      rm(plot_at50feature)
      plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="recall"),] 
      plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
      colnames(plot_at50feature_mean) <- c("data.type","mean")
      plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
      colnames(plot_at50feature_sd) <- c("data.type","sd")
      plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
      plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
      plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "miRNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-all"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-DNA"),]
      plot_at50feature <- plot_at50feature[-which(plot_at50feature$data.type == "Combined-RNA"),]
      
      order <- plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"]
      #order <- order[-which(order=="Combined-all")]
      #order <- order[-which(order=="Combined-DNA")]
      #order <- order[-which(order=="Combined-RNA")]
      #order <- order[-which(order=="miRNA")]
      #order <- append("Combined-RNA",order)
      #order <- append("Combined-DNA",order)
      #order <- append("Combined-all",order)
      #order <- append(order,"miRNA")
      plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
      recall_bar <- ggplot(plot_at50feature,aes(x=data.type,y=mean,fill=data.type))+
        geom_bar(stat = "identity",colour = "black")+
        geom_text(aes(x=data.type,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
        geom_errorbar(aes(x=data.type, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
        xlab("")+
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
        scale_fill_manual(values = RNA_color)+
        #geom_vline(xintercept=1.5,linetype = 2)+
        #geom_vline(xintercept=3.5,linetype = 2)+
        xlab("")+
        ylab("Sensitivity at 95% specificity")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(15,5,-10,45),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
      
      sparseness_bar <- 
        ggplot(plot_at50feature,aes(x=data.type,y=0,fill=X1.sparseness))+
        geom_tile()+
        scale_fill_gradient(low="white",high = "black")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(-10,10,10,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      ggarrange(recall_bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.56),legend = "none")
    }
  }
  
  ROC_final <- rbind(ROC_STAD,ROC_CRC,ROC_cancer)
  recall_final <- rbind(recall_STAD,recall_CRC,recall_cancer)
  #ROC
  {
    plot_at50feature <- ROC_final
    #AUC
    order <- plot_at50feature[order(plot_at50feature[which(plot_at50feature$comparison=="GIC"),]$mean,decreasing = TRUE),"data.type"]
    #order <- order[-which(order=="Combined-all")]
    #order <- order[-which(order=="Combined-DNA")]
    #order <- order[-which(order=="Combined-RNA")]
    #order <- order[-which(order=="miRNA")]
    #order <- append("Combined-RNA",order)
    #order <- append("Combined-DNA",order)
    #order <- append("Combined-all",order)
    #order <- append(order,"miRNA")
    plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
    plot_at50feature$x <- paste(plot_at50feature$data.type,plot_at50feature$comparison,sep = "-")
    
    i=1
    order_x <- {}
    while(i <= length(order)){
      a <- paste(order[i],"GIC",sep = "-")
      b <- paste(order[i],"STAD",sep = "-")
      c <- paste(order[i],"CRC",sep = "-")
      order_x <- c(order_x,a,b,c)
      i=i+1
    }
    
    plot_at50feature$x <- factor(plot_at50feature$x,levels = order_x)
    bar <- ggplot(plot_at50feature,aes(x=x,y=mean,fill=data.type))+
      geom_bar(stat = "identity",colour = "black")+
      geom_text(aes(x=x,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
      geom_errorbar(aes(x=x, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
      xlab("")+
      scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
      scale_fill_manual(values = RNA_color)+
      #geom_vline(xintercept=1.5,linetype = 2)+
      #geom_vline(xintercept=3.5,linetype = 2)+
      xlab("")+
      ylab("AUROC")+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(15,5,-10,70),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))
    
    sparseness_bar <- 
      ggplot(plot_at50feature,aes(x=x,y=0,fill=X1.sparseness))+
      geom_tile()+
      scale_fill_gradient(low="white",high = "black")+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(-10,10,10,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    ggarrange(bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.7),legend = "none")
  }
  #recall
  {
    plot_at50feature <- recall_final
    #AUC
    order <- plot_at50feature[order(plot_at50feature[which(plot_at50feature$comparison=="GIC"),]$mean,decreasing = TRUE),"data.type"]
    #order <- order[-which(order=="Combined-all")]
    #order <- order[-which(order=="Combined-DNA")]
    #order <- order[-which(order=="Combined-RNA")]
    #order <- order[-which(order=="miRNA")]
    #order <- append("Combined-RNA",order)
    #order <- append("Combined-DNA",order)
    #order <- append("Combined-all",order)
    #order <- append(order,"miRNA")
    plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
    plot_at50feature$x <- paste(plot_at50feature$data.type,plot_at50feature$comparison,sep = "-")
    
    i=1
    order_x <- {}
    while(i <= length(order)){
      a <- paste(order[i],"GIC",sep = "-")
      b <- paste(order[i],"STAD",sep = "-")
      c <- paste(order[i],"CRC",sep = "-")
      order_x <- c(order_x,a,b,c)
      i=i+1
    }
    
    plot_at50feature$x <- factor(plot_at50feature$x,levels = order_x)
    bar <- ggplot(plot_at50feature,aes(x=x,y=mean,fill=data.type))+
      geom_bar(stat = "identity",colour = "black")+
      geom_text(aes(x=x,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
      geom_errorbar(aes(x=x, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
      xlab("")+
      scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
      scale_fill_manual(values = RNA_color)+
      #geom_vline(xintercept=1.5,linetype = 2)+
      #geom_vline(xintercept=3.5,linetype = 2)+
      xlab("")+
      ylab("Sensitivity at 95% specificity")+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(15,5,-10,70),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
        plot.title = element_text(hjust = 0.5,size=24,face="bold"),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_text(face="bold", color="black", size=20),
        axis.title.y = element_text(face="bold",color="black", size=20))
    
    sparseness_bar <- 
      ggplot(plot_at50feature,aes(x=x,y=0,fill=X1.sparseness))+
      geom_tile()+
      scale_fill_gradient(low="white",high = "black")+
      theme_bw()+
      theme(#legend.position="right",
        plot.margin = unit(x=c(-10,10,10,10),units="pt"),
        legend.position="right",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
        legend.title = element_blank(),
        legend.text= element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
    
    ggarrange(bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.7),legend = "none")
  }
  
  
  #miRNA
  {
    {
    comparison <- "STAD-NC"
    #ROC_STAD
    {
      rm(plot_at50feature)
      plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="AUROC"),] 
      plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
      colnames(plot_at50feature_mean) <- c("data.type","mean")
      plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
      colnames(plot_at50feature_sd) <- c("data.type","sd")
      plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
      plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
      plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
      plot_at50feature <- plot_at50feature[which(plot_at50feature$data.type == "miRNA"),]
      ROC_STAD <- plot_at50feature
      ROC_STAD$comparison <- "STAD"
    }
    #recall at 95% specificity
    {
      rm(plot_at50feature)
      plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="recall"),] 
      plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
      colnames(plot_at50feature_mean) <- c("data.type","mean")
      plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
      colnames(plot_at50feature_sd) <- c("data.type","sd")
      plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
      plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
      plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
      plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
      plot_at50feature <- plot_at50feature[which(plot_at50feature$data.type == "miRNA"),]
      recall_STAD <- plot_at50feature
      recall_STAD$comparison <- "STAD"
    }
    }
    {
      comparison <- "CRC-NC"
      #ROC_CRC
      {
        rm(plot_at50feature)
        plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="AUROC"),] 
        plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
        colnames(plot_at50feature_mean) <- c("data.type","mean")
        plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
        colnames(plot_at50feature_sd) <- c("data.type","sd")
        plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
        plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
        plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
        plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
        plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
        plot_at50feature <- plot_at50feature[which(plot_at50feature$data.type == "miRNA"),]
        ROC_CRC <- plot_at50feature
        ROC_CRC$comparison <- "CRC"
      }
      #recall at 95% specificity
      {
        rm(plot_at50feature)
        plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="recall"),] 
        plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
        colnames(plot_at50feature_mean) <- c("data.type","mean")
        plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
        colnames(plot_at50feature_sd) <- c("data.type","sd")
        plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
        plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
        plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
        plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
        plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
        plot_at50feature <- plot_at50feature[which(plot_at50feature$data.type == "miRNA"),]
        recall_CRC <- plot_at50feature
        recall_CRC$comparison <- "CRC"
      }
    }
    {
      comparison <- "cancer-NC"
      #ROC_cancer
      {
        rm(plot_at50feature)
        plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="AUROC"),] 
        plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
        colnames(plot_at50feature_mean) <- c("data.type","mean")
        plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
        colnames(plot_at50feature_sd) <- c("data.type","sd")
        plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
        plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
        plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
        plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
        plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
        plot_at50feature <- plot_at50feature[which(plot_at50feature$data.type == "miRNA"),]
        ROC_cancer <- plot_at50feature
        ROC_cancer$comparison <- "GIC"
      }
      #recall at 95% specificity
      {
        rm(plot_at50feature)
        plot_at50feature_tmp <- plot[(plot$cancer==comparison)&(plot$metric=="recall"),] 
        plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
        colnames(plot_at50feature_mean) <- c("data.type","mean")
        plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
        colnames(plot_at50feature_sd) <- c("data.type","sd")
        plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
        plot_at50feature <- left_join(plot_at50feature,sparseness,by = c("data.type"="alteration"))
        plot_at50feature[plot_at50feature$data.type=="Combined-all","X1.sparseness"] <- 0
        plot_at50feature[plot_at50feature$data.type=="Combined-RNA","X1.sparseness"] <- 0
        plot_at50feature[plot_at50feature$data.type=="Combined-DNA","X1.sparseness"] <- 0
        plot_at50feature <- plot_at50feature[which(plot_at50feature$data.type == "miRNA"),]
        recall_cancer <- plot_at50feature
        recall_cancer$comparison <- "GIC"
      }
    }
    ROC_miRNA <- rbind(ROC_STAD,ROC_CRC,ROC_cancer)
    recall_miRNA <- rbind(recall_STAD,recall_CRC,recall_cancer)
    
    #ROC
    {
      plot_at50feature <- ROC_miRNA
      #AUC
      order <- plot_at50feature[order(plot_at50feature[which(plot_at50feature$comparison=="GIC"),]$mean,decreasing = TRUE),"data.type"]
      #order <- order[-which(order=="Combined-all")]
      #order <- order[-which(order=="Combined-DNA")]
      #order <- order[-which(order=="Combined-RNA")]
      #order <- order[-which(order=="miRNA")]
      #order <- append("Combined-RNA",order)
      #order <- append("Combined-DNA",order)
      #order <- append("Combined-all",order)
      #order <- append(order,"miRNA")
      plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
      plot_at50feature$x <- paste(plot_at50feature$data.type,plot_at50feature$comparison,sep = "-")
      
      i=1
      order_x <- {}
      while(i <= length(order)){
        a <- paste(order[i],"GIC",sep = "-")
        b <- paste(order[i],"STAD",sep = "-")
        c <- paste(order[i],"CRC",sep = "-")
        order_x <- c(order_x,a,b,c)
        i=i+1
      }
      
      plot_at50feature$x <- factor(plot_at50feature$x,levels = order_x)
      bar <- ggplot(plot_at50feature,aes(x=x,y=mean,fill=data.type))+
        geom_bar(stat = "identity",colour = "black")+
        geom_text(aes(x=x,y=mean+0.08,label=round(mean,digits = 3)),size = 4,angle = 0)+
        geom_errorbar(aes(x=x, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
        xlab("")+
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
        scale_fill_manual(values = RNA_color)+
        #geom_vline(xintercept=1.5,linetype = 2)+
        #geom_vline(xintercept=3.5,linetype = 2)+
        xlab("")+
        ylab("AUROC")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(15,5,-10,70),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
      
      sparseness_bar <- 
        ggplot(plot_at50feature,aes(x=x,y=0,fill=X1.sparseness))+
        geom_tile()+
        scale_fill_gradient(low="white",high = "black")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(-10,10,10,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      ggarrange(bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.35),legend = "none")
    }
    #recall
    {
      plot_at50feature <- recall_miRNA
      #AUC
      order <- plot_at50feature[order(plot_at50feature[which(plot_at50feature$comparison=="GIC"),]$mean,decreasing = TRUE),"data.type"]
      #order <- order[-which(order=="Combined-all")]
      #order <- order[-which(order=="Combined-DNA")]
      #order <- order[-which(order=="Combined-RNA")]
      #order <- order[-which(order=="miRNA")]
      #order <- append("Combined-RNA",order)
      #order <- append("Combined-DNA",order)
      #order <- append("Combined-all",order)
      #order <- append(order,"miRNA")
      plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
      plot_at50feature$x <- paste(plot_at50feature$data.type,plot_at50feature$comparison,sep = "-")
      
      i=1
      order_x <- {}
      while(i <= length(order)){
        a <- paste(order[i],"GIC",sep = "-")
        b <- paste(order[i],"STAD",sep = "-")
        c <- paste(order[i],"CRC",sep = "-")
        order_x <- c(order_x,a,b,c)
        i=i+1
      }
      
      plot_at50feature$x <- factor(plot_at50feature$x,levels = order_x)
      bar <- ggplot(plot_at50feature,aes(x=x,y=mean,fill=data.type))+
        geom_bar(stat = "identity",colour = "black")+
        geom_text(aes(x=x,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
        geom_errorbar(aes(x=x, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
        xlab("")+
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
        scale_fill_manual(values = RNA_color)+
        #geom_vline(xintercept=1.5,linetype = 2)+
        #geom_vline(xintercept=3.5,linetype = 2)+
        xlab("")+
        ylab("Sensitivity at 95% specificity")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(15,5,-10,70),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
          plot.title = element_text(hjust = 0.5,size=24,face="bold"),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold",color="black", size=20))
      
      sparseness_bar <- 
        ggplot(plot_at50feature,aes(x=x,y=0,fill=X1.sparseness))+
        geom_tile()+
        scale_fill_gradient(low="white",high = "black")+
        theme_bw()+
        theme(#legend.position="right",
          plot.margin = unit(x=c(-10,10,10,10),units="pt"),
          legend.position="right",
          panel.grid=element_blank(),
          panel.border=element_blank(),
          #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
          legend.title = element_blank(),
          legend.text= element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_blank(),
          #axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(face="bold", color="black", size=0, angle = 0,hjust=1,vjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
      
      ggarrange(bar,sparseness_bar,ncol = 1,nrow = 2,align = "v",heights = c(1,0.35),legend = "none")
    }
  }
  
}

#performance with clinical stratify
{
#all stage ROC barplot
{
plot <- read.csv("2021-11-05-cv-performance-individual-and-all.txt",header = TRUE,sep = "\t")
plot_early <- read.csv("early-2021-11-05-cv-performance-individual-and-all.txt",header = TRUE,sep = "\t")

plot <- rbind(plot,plot_early)

RNA_color <- c("Alternative promoter"="#E9C2A6",
               "RNA Expression"="#A5435C",
               "RNA Splicing"="#C5E3BF",
               "Alternative polyadenyltion"="#003F87",
               "Chimeric RNA"="#FF3D0D",
               "RNA editing"="#324F17",
               "Allele specific expression"="#800080",
               "RNA SNP"="#333333",
               "DNA copy number"="#FFE234",
               "DNA nucleosome occupancy"="#44A5F9",
               "DNA methylation"="#C7C2C2",
               "Combined"="#FFFFFF",
               "miRNA"="#000000")


plot$data.type <- gsub("expression","RNA Expression",plot$data.type)
plot$data.type <- gsub("alt-promoter","Alternative promoter",plot$data.type)
plot$data.type <- gsub("APA","Alternative polyadenyltion",plot$data.type)
plot$data.type <- gsub("ASE","Allele specific expression",plot$data.type)
plot$data.type <- gsub("chimera","Chimeric RNA",plot$data.type)
plot$data.type <- gsub("CNV","DNA copy number",plot$data.type)
plot$data.type <- gsub("editing","RNA editing",plot$data.type)
plot$data.type <- gsub("MeDIP","DNA methylation",plot$data.type)
plot$data.type <- gsub("nucleasome","DNA nucleosome occupancy",plot$data.type)
plot$data.type <- gsub("SNP","RNA SNP",plot$data.type)
plot$data.type <- gsub("splicing","RNA Splicing",plot$data.type)
plot$data.type <- gsub("combined","Combined",plot$data.type)

rm(plot_at50feature)
plot_at50feature_tmp <- plot[(plot$cancer=="CRC-NC")&(plot$metric=="AUROC"),] 
plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
colnames(plot_at50feature_mean) <- c("data.type","mean")
plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
colnames(plot_at50feature_sd) <- c("data.type","sd")
plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))

plot_at50feature_tmp <- plot[(plot$cancer=="CRC-NC")&(plot$metric=="AUROC_early"),] 
plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, min)
colnames(plot_at50feature_mean) <- c("data.type","mean")
plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
colnames(plot_at50feature_sd) <- c("data.type","sd")
plot_at50feature_early <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))
colnames(plot_at50feature_early) <- c("data.type","mean_early","sd_early")

plot_at50feature <- left_join(plot_at50feature,plot_at50feature_early,by=c("data.type"="data.type"))

#AUC
order <- plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"]
order <- order[-which(order=="miRNA")]
order <- append(order,"miRNA")
plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = order)
ggplot(plot_at50feature)+
  geom_bar(aes(x=data.type,y=mean-0.3,fill=data.type),stat = "identity",colour = "black")+
  geom_bar(aes(x=data.type,y=mean_early-0.3),stat = "identity",colour = "black",fill=alpha("white",0.4))+
  geom_text(aes(x=data.type,y=mean-0.3+0.08,label=round(mean,digits = 3)),size = 4,angle = 0)+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels = c("0.3","0.4","0.5","0.6","0.7","0.8","0.9","1 "),expand = c(0,0),limits = c(0,0.7))+
  scale_fill_manual(values = RNA_color)+
  geom_errorbar(aes(x=data.type, ymin=mean-0.3-sd, ymax=mean-0.3+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
  xlab("")+
  ylab("AUROC")+
  theme_bw()+
  geom_vline(xintercept=1.5,linetype = 2)+
  theme(#legend.position="right",
    plot.margin = unit(x=c(20,20,10,80),units="pt"),
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
    legend.title = element_blank(),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
    axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
    axis.title.x = element_text(face="bold", color="black", size=20),
    axis.title.y = element_text(face="bold",color="black", size=20))

#recall at 95% specificity
plot_at50feature_tmp <- plot[(plot$cancer=="STAD-NC")&(plot$metric=="recall"),] 

plot_at50feature_mean <- aggregate(value ~ data.type, plot_at50feature_tmp, mean)
colnames(plot_at50feature_mean) <- c("data.type","mean")
plot_at50feature_sd <- aggregate(value ~ data.type, plot_at50feature_tmp, sd)
colnames(plot_at50feature_sd) <- c("data.type","sd")
plot_at50feature <- left_join(plot_at50feature_mean,plot_at50feature_sd, by =c("data.type"="data.type"))


plot_at50feature$data.type <- factor(plot_at50feature$data.type,levels = plot_at50feature[order(plot_at50feature$mean,decreasing = TRUE),"data.type"])
ggplot(plot_at50feature,aes(x=data.type,y=mean,fill=data.type))+
  geom_bar(stat = "identity",colour = "black")+
  geom_text(aes(x=data.type,y=mean+0.05,label=round(mean,digits = 3)),size = 4,angle = 0)+
  geom_errorbar(aes(x=data.type, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.5, size=1.5)+
  xlab("")+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),labels = c("0","0.2","0.4","0.6","0.8"),expand = c(0,0),limits = c(0,1))+
  scale_fill_manual(values = RNA_color)+
  xlab("")+
  ylab("Sensitivity at 95% specificity")+
  theme_bw()+
  theme(#legend.position="right",
    plot.margin = unit(x=c(15,5,10,45),units="pt"),
    legend.position="right",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
    legend.title = element_blank(),
    legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
    plot.title = element_text(hjust = 0.5,size=24,face="bold"),
    axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_text(face="bold", color="black", size=16, angle = 45,hjust = 1,vjust = 1),
    axis.text.y = element_text(face="bold", color="black", size=18, angle = 0,hjust=1,vjust = 0.5),
    axis.title.x = element_text(face="bold", color="black", size=20),
    axis.title.y = element_text(face="bold",color="black", size=20))
}

#early stage ROC plot
{
ML_samples <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/ML_samples.csv",header = TRUE)
performance_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20211105_yunfan/performance"
early_roc <- data.frame(matrix(ncol=4))
colnames(early_roc) <- c("data-type","metric","cancer","value")

n=1
a=1
cancers <- dir(performance_dir)
while(a<= length(cancers)){
  b=1
  datatype <- dir(paste0(performance_dir,"/",cancers[a]))
  while(b<=length(datatype)){
    c=1
    while(c<=10){
    predicted_all <- read.csv(paste0(performance_dir,"/",cancers[a],"/",datatype[b],"/",c,"/","predictions.txt"),sep = "\t")

    predicted_meta <- left_join(predicted_all,ML_samples,by = c("sample_id"="sample_id"))
    predicted_early <- predicted_meta[predicted_meta$Stage_numeric<2,c(1,2,3)]

    roc_early <- roc(predicted_early$y_true,predicted_early$y_pred,levels=c(0,1),direction=c("<"))
  
    early_roc[n,"cancer"]<- cancers[a]
    early_roc[n,"data-type"]<- datatype[b]
    early_roc[n,"metric"]<- "AUROC_early"
    early_roc[n,"value"]<- as.numeric(roc_early$auc)
    n=n+1
    c=c+1
    }
  b=b+1
  }
a=a+1
}
write.table(early_roc,"early-2021-11-05-cv-performance-individual-and-all.txt",sep= "\t",quote = FALSE,row.names = FALSE)
}
}

#combination
{
  library(ggbeeswarm)
  plot <- read.csv("combination-performance.txt",header = TRUE,sep = "\t")
  
  plot <- plot[(plot$metric=="AUROC")&(plot$cancer=="cancer-NC"),]
  
  plot %>%
      group_by(n) %>%
    summarise(median=median(value)) -> plot_median
  
  #plot <- plot[-which(plot$n==10),]
  plot$n <- factor(as.character(plot$n),levels = as.character(seq(from=1,to=10)))
  plot$red <- plot$value
  plot[plot$red <=0.9,"red"] <- NA
  #plot[plot$red > 0.9,"red"] <- "red"
 
  
  ggplot(plot,aes(x=plot$n,y=plot$value))+
    scale_y_continuous(breaks = c(0,0.25,0.50,0.75,0.9,1),labels = c("0","0.25","0.50","0.75","0.90",""),expand = c(0,0),limits = c(0.25,1))+
    #geom_boxplot(alpha = 1, aes(fill = as.numeric(as.character(plot$n))))+
    geom_hline(yintercept=0.9,linetype = 2, color = "red")+
    geom_quasirandom(method = "smiley",size = 1,width = 0.3,shape = 21,stroke = 0.6,aes(fill = as.numeric(as.character(plot$n)),color = as.numeric(as.character(plot$n))))+
    geom_quasirandom(aes(x=plot$n,y=plot$red),method = "smiley",size = 1.5,width = 0,shape = 21,stroke = 0.8,color = "red",fill="red")+
    #geom_line(data=plot[plot$max=="max",],aes(x=plot[plot$max=="max",]$n,y=plot[plot$max=="max",]$value,group = max),color = alpha("purple", alpha = 0.5), size = 1) +
    #geom_line(data=plot_median,aes(x=plot_median$n,y=plot_median$median,group = 1),color = alpha("purple", alpha = 0.5), size = 1) +
    #geom_smooth(aes(group = as.numeric(as.character(plot$n)))) +
    #geom_violin(alpha = 1, aes(linetype = NA,fill = as.numeric(as.character(plot$n))))+
    #geom_jitter(alpha = 0.4, size = 0.5, shape = 21, aes(linetype = NA,fill = as.numeric(as.character(plot$n))),position = position_jitter(width = 0.2))+
    scale_fill_gradient2()+
    scale_color_gradient2()+
    xlab("")+
    ylab("AUROC")+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(15,5,10,10),units="pt"),
      legend.position= "none",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.text= element_text(face="bold", color="black",family = "Arial", size=20),
      plot.title = element_text(hjust = 0.5,size=24,face="bold"),
      axis.ticks.x = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.text.x = element_text(color="black", size=23, angle = 0,hjust = 0.5,vjust = 0.5),
      axis.text.y = element_text(color="black", size=16, angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_text(color="black", size=23),
      axis.title.y = element_text(color="black", size=23))
}

#plot ROC-curve
{
library(pROC)
#report ROC
files <- grep("-",dir("cross-validation-combination/cancer-NC/5/6/"),value = TRUE)
i=1
while(i<=length(files)){
  predicted <- read.csv(paste0("cross-validation-combination/cancer-NC/5/6/",files[i]),sep = "\t")
  roc.curve <- roc(predicted$y_true,predicted$y_pred,levels=c(0,1),direction=c("<"))
  ci.auc(roc.curve,conf.level = 0.95)
  record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
  print(paste0(files[i],":",roc.curve$auc))
  i=i+1
}

#test
{
outdir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20211114_yunfan"
predicted <- read.csv("cross-validation-combination/cancer-NC/3/alt-promoter-nucleasome-splicing.txt",sep = "\t")
roc.curve <- roc(predicted$y_true,predicted$y_pred,levels=c(0,1),direction=c("<"))
ci.auc(roc.curve,conf.level = 0.95)
record <- ci.coords(roc.curve,x="best",conf.level = 0.95,ret = c("recall","specificity","precision"),best.method="youden",best.policy="random")
print(paste0(files[i],":",roc.curve$auc))
write.table(roc.curve$auc,paste0(outdir,"AUC.txt"),quote = FALSE,sep="\t")
#write.table(record,paste0(outdir,"matrices_record.txt"),quote = FALSE,sep="\t")

#plot
pdf(paste0("testAUROC.pdf"),width=7.8,height=6)
roc_curve_prediction_df <- PRROC::roc.curve(scores.class0=predicted$y_pred, weights.class0 = predicted$y_true,
                                            curve=TRUE,rand.compute = TRUE) # rand.compute 是0.5的随机分类器，为了画对角线的

par(cex.axis=2, cex.lab =1.6,font.axis =1.6, mar = c(5,5,2,3))
plot(roc_curve_prediction_df,rand.plot = TRUE,auc.main=FALSE,legend =4,lwd=5)
text(1.3, 0.2, paste('AUROC =', ' ',round(roc_curve_prediction_df$auc, digits = 4), sep =""), cex=2, font=2, col='black')
dev.off()
}

predicted_all <- read.csv("cross-validation-combination/cancer-NC/3/MeDIP-nucleasome-CNV.txt",sep = "\t")
roc(predicted_all$y_true,predicted_all$y_pred,plot=TRUE,legacy.axes=TRUE,
    xlab="False Positive Percentage",ylab="True Positive Percentage",
    col="red", lwd=4,print.auc=TRUE,print.auc.y= 0.25,print.auc.x= 0.35)

predicted_all <- read.csv("cross-validation-combination/cancer-NC/7/alt-promoter-ASE-editing-expression-SNP-chimera-splicing.txt",sep = "\t")
roc(predicted_all$y_true,predicted_all$y_pred,plot=TRUE,legacy.axes=TRUE,
    xlab="False Positive Percentage",ylab="True Positive Percentage",
    col="red", lwd=4,print.auc=TRUE,print.auc.y= 0.25,print.auc.x= 0.35)

predicted_all <- read.csv("cross-validation-combination/cancer-NC/9/MeDIP-alt-promoter-ASE-editing-expression-nucleasome-SNP-chimera-splicing.txt",sep = "\t")
roc(predicted_all$y_true,predicted_all$y_pred,plot=TRUE,legacy.axes=TRUE,
    xlab="False Positive Percentage",ylab="True Positive Percentage",
    col="red", lwd=4,print.auc=TRUE,print.auc.y= 0.25,print.auc.x= 0.35)


ML_samples <- read.csv("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/ML_samples.csv",header = TRUE)
predicted_all <- read.csv("cross-validation-combination/cancer-NC/5/6/MeDIP-alt-promoter-editing-nucleasome-SNP-chimera.txt",sep = "\t")

predicted_meta <- left_join(predicted_all,ML_samples,by = c("sample_id"="sample_id"))

#alteration
{
  predicted_medip <- read.csv("cross-validation-combination/cancer-NC/5/1/MeDIP.txt",sep = "\t")
  predicted_altpromoter <- read.csv("cross-validation-combination/cancer-NC/5/1/alt-promoter.txt",sep = "\t")
  predicted_editing <- read.csv("cross-validation-combination/cancer-NC/5/1/editing.txt",sep = "\t")
  predicted_nucleasome <- read.csv("cross-validation-combination/cancer-NC/5/1/nucleasome.txt",sep = "\t")
  predicted_SNP <- read.csv("cross-validation-combination/cancer-NC/5/1/SNP.txt",sep = "\t")
  predicted_chimeric <- read.csv("cross-validation-combination/cancer-NC/5/1/chimera.txt",sep = "\t")
  
  RNA_color <- c("Alternative promoter"="#E9C2A6",
                 "RNA Expression"="#A5435C",
                 "RNA Splicing"="#C5E3BF",
                 "Alternative polyadenyltion"="#003F87",
                 "Chimeric RNA"="#FF3D0D",
                 "RNA editing"="#324F17",
                 "Allele specific expression"="#800080",
                 "RNA SNP"="#333333",
                 "DNA copy number"="#FFE234",
                 "DNA nucleosome occupancy"="#44A5F9",
                 "DNA methylation"="#C7C2C2",
                 "Combined"="white")
  
  par(pty="s")
  roc(predicted_all$y_true,predicted_all$y_pred,plot=TRUE,legacy.axes=TRUE,
      xlab="False Positive Percentage",ylab="True Positive Percentage",
      col="purple", lwd=4,print.auc=TRUE,print.auc.y= 0.35,print.auc.x= 0.45,cex.axis = 1.5, cex.lab = 2)
  plot.roc(predicted_medip$y_true,predicted_medip$y_pred,
           col="#C7C2C2", lwd=4, print.auc=TRUE, print.auc.y=0.30, print.auc.x= 0.45, add=TRUE)
  plot.roc(predicted_altpromoter$y_true,predicted_altpromoter$y_pred,
           col="#E9C2A6", lwd=4, print.auc=TRUE, print.auc.y=0.25, print.auc.x= 0.45, add=TRUE)
  plot.roc(predicted_editing$y_true,predicted_editing$y_pred,
           col="#49E20E", lwd=4, print.auc=TRUE, print.auc.y=0.20, print.auc.x= 0.45, add=TRUE)
  plot.roc(predicted_nucleasome$y_true,predicted_nucleasome$y_pred,
           col="#44A5F9", lwd=4, print.auc=TRUE, print.auc.y=0.15, print.auc.x= 0.45, add=TRUE)
  plot.roc(predicted_SNP$y_true,predicted_SNP$y_pred,
           col="#333333", lwd=4, print.auc=TRUE, print.auc.y=0.10, print.auc.x= 0.45, add=TRUE)
  plot.roc(predicted_chimeric$y_true,predicted_chimeric$y_pred,
           col="#FF3D0D", lwd=4, print.auc=TRUE, print.auc.y=0.05, print.auc.x= 0.45, add=TRUE)
  legend("bottomright",legend=c("6-Combined: 0.914","RNA altpromoter: 0.835","RNA editing: 0.811","Nucleasome occupancy: 0.750","RNA chimeric: 0.710","RNA SNP: 0.674","DNA methylation: 0.552"),
         col=c("purple","#E9C2A6","#49E20E","#44A5F9","#FF3D0D","#333333","#C7C2C2"),lwd=4,
         cex = 1.2)
}

#stage
{
predicted_1 <- predicted_meta[(predicted_meta$Stage_numeric==1)|(predicted_meta$Stage_numeric==0),c(1,2,3)]
predicted_2 <- predicted_meta[(predicted_meta$Stage_numeric==2)|(predicted_meta$Stage_numeric==0),c(1,2,3)]
predicted_3 <- predicted_meta[(predicted_meta$Stage_numeric==3)|(predicted_meta$Stage_numeric==0),c(1,2,3)]
predicted_4 <- predicted_meta[(predicted_meta$Stage_numeric==4)|(predicted_meta$Stage_numeric==0),c(1,2,3)]

par(pty="s")
roc(predicted_1$y_true,predicted_1$y_pred,plot=TRUE,legacy.axes=TRUE,
    xlab="False Positive Percentage",ylab="True Positive Percentage",
    col="red", lwd=4,print.auc=TRUE,print.auc.y= 0.20,print.auc.x= 0.35,cex.axis = 1.5, cex.lab = 2)
plot.roc(predicted_2$y_true,predicted_2$y_pred,
         col="#EE9A00", lwd=4, print.auc=TRUE, print.auc.y=0.15, print.auc.x= 0.35, add=TRUE)
plot.roc(predicted_3$y_true,predicted_3$y_pred,
         col="grey", lwd=4, print.auc=TRUE, print.auc.y=0.10, print.auc.x= 0.35, add=TRUE)
plot.roc(predicted_4$y_true,predicted_4$y_pred,
         col="black", lwd=4, print.auc=TRUE, print.auc.y=0.05, print.auc.x= 0.35, add=TRUE)
legend("bottomright",legend=c("Stage I: 0.952","Stage II: 0.875","Stage III: 0.922","Stage IV: 0.838"),
       col=c("red","#EE9A00","grey","black"),lwd=4,cex=1.6)
}

#gender
{
  predicted_male <- predicted_meta[predicted_meta$Gender=="M",c(1,2,3)]
  predicted_female <- predicted_meta[(predicted_meta$Gender=="F"),c(1,2,3)]
  
  par(pty="s")
  roc(predicted_all$y_true,predicted_all$y_pred,plot=TRUE,legacy.axes=TRUE,
      xlab="False Positive Percentage",ylab="True Positive Percentage",
      col="purple", lwd=4,print.auc=TRUE,print.auc.y= 0.20,print.auc.x= 0.35,cex.axis = 1.5, cex.lab = 2)
  plot.roc(predicted_male$y_true,predicted_male$y_pred,
           col="dark blue", lwd=4, print.auc=TRUE, print.auc.y=0.15, print.auc.x= 0.35, add=TRUE)
  plot.roc(predicted_female$y_true,predicted_female$y_pred,
           col="pink", lwd=4, print.auc=TRUE, print.auc.y=0.10, print.auc.x= 0.35, add=TRUE)
  legend("bottomright",legend=c("All gender: 0.914","Male: 0.994","Female: 0.754"),
         col=c("purple","dark blue","pink"),lwd=4,cex=1.6)
}

#Age
{
  predicted_older <- predicted_meta[predicted_meta$Age>=60,c(1,2,3)]
  predicted_younger <- predicted_meta[(predicted_meta$Age<60),c(1,2,3)]
  
  par(pty="s")
  roc(predicted_all$y_true,predicted_all$y_pred,plot=TRUE,legacy.axes=TRUE,
      xlab="False Positive Percentage",ylab="True Positive Percentage",
      col="purple", lwd=4,print.auc=TRUE,print.auc.y= 0.20,print.auc.x= 0.35,cex.axis = 1.5, cex.lab = 2)
  plot.roc(predicted_older$y_true,predicted_older$y_pred,
           col="brown", lwd=4, print.auc=TRUE, print.auc.y=0.15, print.auc.x= 0.35, add=TRUE)
  plot.roc(predicted_younger$y_true,predicted_younger$y_pred,
           col="green", lwd=4, print.auc=TRUE, print.auc.y=0.10, print.auc.x= 0.35, add=TRUE)
  legend("bottomright",legend=c("All: 0.914","Older(>=60): 0.861","Younger(<60): 0.940"),
         col=c("purple","brown","green"),lwd=4,cex=1.6)
}
}

#plot features of a specific combination(50 features)
{
  RNA_color <- c("Alternative promoter"="#E9C2A6",
                 "RNA Expression"="#A5435C",
                 "RNA Splicing"="#C5E3BF",
                 "Alternative polyadenyltion"="#003F87",
                 "Chimeric RNA"="#FF3D0D",
                 "RNA editing"="#49E20E",
                 "Allele specific expression"="#800080",
                 "RNA SNP"="#333333",
                 "DNA copy number"="#FFE234",
                 "DNA nucleosome occupancy"="#44A5F9",
                 "DNA methylation"="#C7C2C2",
                 "Combined"="#FFFFFF",
                 "miRNA"="#000000")
  features_dir <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20211114_yunfan/features"
  group <- "cancer-NC"
  alterations <- c("MeDIP","alt-promoter","editing","nucleasome","SNP","chimera")
  number <- length(alterations)
  seed <- "5"
  
  #summarize features for a specific combination
  i=1
  features <- {}
  while(i<=length(alterations)){
  features_LOOCV <- dir(paste0(features_dir,"/",group,"/",alterations[i],"/",seed))
  j=1
  single_features <- {}
  while(j<=length(features_LOOCV)){
  single_alteration_features <- read.csv(paste0(features_dir,"/",group,"/",alterations[i],"/",seed,"/",features_LOOCV[j]),sep = "\t",header = FALSE)
  single_alteration_features <- head(single_alteration_features,floor(50/number))
  single_alteration_features$type <- alterations[i]
  single_features <- rbind(single_features,single_alteration_features)
  j=j+1
  }
  features <- rbind(features,single_features)
  i=i+1
  }
  
  
  #
  features$count <- 1
  colnames(features) <- c("alteration","importance","type","count")
  features$alteration <- paste(features$type,features$alteration,sep = "|")
  feature_importance <- aggregate(importance ~ alteration, features, mean)
  feature_reccurence <- aggregate(count ~ alteration, features, sum)
  
  feature_plot <- left_join(feature_reccurence,feature_importance, by = c("alteration"="alteration"))
  
  final_plot <- feature_plot[feature_plot$count>=(length(features_LOOCV))*0.5,]
  
  final_plot <- head(final_plot[order(final_plot$importance,decreasing = TRUE),],20)
  #write.csv(final_plot,"demo_feature.csv")
  final_plot <- read.csv("demo_feature.csv")
  final_plot <- head(final_plot[order(final_plot$importance,decreasing = TRUE),],15)
  final_plot$type <- gsub("Nucleosome occupancy","DNA nucleosome occupancy",fixed = TRUE, final_plot$type)
  final_plot$name <- factor(final_plot$name,levels = final_plot$name)
  final_plot$name2 <- factor(final_plot$name2,levels = final_plot$name2)
  feature_bar <- ggplot(final_plot,aes(x=importance,y=name,fill=-log10(pvalue)))+
    geom_bar(stat='identity', color = "black")+
    scale_fill_gradient2(low="#B0E2FF",high="#162252",mid="#26466D",midpoint = 5,na.value = "white")+
    scale_x_continuous(expand = c(0,0.0001))+
    #geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
    #geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
    labs(color="Importance",
         fill=expression(-log[10]("P value")),
         x="")+
    theme_bw()+
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.text = element_text(size = rel(1.3),face="bold",colour = "black"),
          axis.text.y = element_text(size = 0,colour = "black"),
          axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black",angle = 45,vjust=1,hjust=1),
          axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
          axis.title.y = element_blank())
  
  type_bar <- 
    ggplot(final_plot,aes(x=0,y=name2,fill=type))+
    geom_tile()+
    scale_fill_manual(values = RNA_color)+
    theme_bw()+
    theme(#legend.position="right",
      plot.margin = unit(x=c(-10,-10,10,10),units="pt"),
      legend.position="right",
      panel.grid=element_blank(),
      panel.border=element_blank(),
      #legend.title = element_text(face="bold", color="black",family = "Arial", size=20),
      legend.title = element_blank(),
      legend.text= element_blank(),
      plot.title = element_blank(),
      axis.ticks = element_blank(),
      #axis.text.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(face="bold", color="black", size=0, angle = 45,hjust = 1,vjust = 1),
      axis.text.y = element_text(face="bold", color="black", size=rel(1.3), angle = 0,hjust=1,vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  ggarrange(type_bar,feature_bar,ncol = 2,nrow = 1,align = "h",heights = c(2,1),legend = "none")
}

#sparseness
{
path <- "/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20211105_yunfan/ML_matrix"

alterations <- gsub("_ML.txt","",dir(path),fixed=TRUE)
alteration_sparseness <- data.frame(matrix(ncol=3))
colnames(alteration_sparseness) <- c("alteration","sparseness","alteration_number")
i=1
while(i<=length(alterations)){
  matrix <- read.csv(paste0(path,"/",alterations[i],"_ML.txt"),sep = "\t",header = TRUE)
  matrix[is.na(matrix)] <- 0
  total <- (nrow(matrix)*ncol(matrix))
  zero <- sum(matrix == 0)
  sparseness <- zero/total
  alteration_sparseness[i,"alteration"] <- alterations[i]
  alteration_sparseness[i,"sparseness"] <- sparseness
  alteration_sparseness[i,"alteration_number"] <- nrow(matrix)
  i=i+1
}
alteration_sparseness$`1-sparseness` <- 1-alteration_sparseness$sparseness
write.csv(alteration_sparseness,"alteration_sparseness.csv")
}


#plot feature of 10 combination(250 features)
{
  features <- {}
  cancers <- c("CRC-NC","CRC-STAD","STAD-NC","cancer-NC")
  k=1
  while(k<=length(cancers)){
    l=1
    while(l<=10){
      feature_tmp <- read.csv(paste0("/Users/yuhuan/Desktop/Seafile/Lulab/2021/02.exOmics/multiomics_RNA/08.MachineLearning/20211114_yunfan/feature-importance/",cancers[k],"-",l,".txt"),header = FALSE,sep = "\t")
      feature_tmp$cancer <- cancers[k]
      features <- rbind(features,feature_tmp)
      l=l+1
    }
    k=k+1
  }
  features$count <- 1
  colnames(features) <- c("alteration","importance","cancer","count")
  feature_STADvsNC <- features[which(features$cancer=="cancer-NC"),]
  feature_importance <- aggregate(importance ~ alteration, feature_STADvsNC, mean)
  feature_reccurence <- aggregate(count ~ alteration, feature_STADvsNC, sum)
  
  feature_plot <- left_join(feature_reccurence,feature_importance, by = c("alteration"="alteration"))
  
  final_plot <- feature_plot[feature_plot$count=="10",]
  
  final_plot <- head(final_plot[order(final_plot$importance,decreasing = TRUE),],20)
  write.csv(final_plot,"demo_feature.csv")
  final_plot <- read.csv("demo_feature.csv")
  ggplot(final_plot,aes(x=importance,y=name,fill=-log10(pvalue)))+
    geom_bar(stat='identity')+
    scale_colour_gradient2(low="#162252",high="#3B5E2B",mid="#008B00",midpoint = 0)+
    #geom_vline(xintercept = c(6.5,10.5,16.5),linetype = "dashed",color="black")+
    #geom_vline(xintercept = 25.5,linetype = "dashed",color="grey")+
    labs(color="Importance",
         #size=expression("Reccurence"),
         x="")+
    theme_bw()+
    theme(strip.text = element_text(size = rel(1.3),face="bold",colour = "black"),
          axis.text.y = element_text(size = rel(1.3),face="bold",colour = "black"),
          axis.text.x = element_text(size=rel(1.3),face="bold",colour = "black",angle = 45,vjust=1,hjust=1),
          axis.title.x = element_text(size=rel(1.3),face="bold",colour = "black"),
          axis.title.y = element_blank())
}
