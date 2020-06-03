library(ggplot2)

p_value <- t(as.matrix(read.table("2020-05-29_Pvalue_HQ.txt",sep="\t",header=T,row.names=1)))
table <- t(log(as.matrix(read.table("2020-05-26_table2_HQ.txt",sep="\t",header=T,row.names=1))))[-1,]
colnames(p_value)[7] <-">20"
p_value_melt <- reshape::melt(p_value)
table_melt <- reshape::melt(table)

names(p_value_melt) <- c("genotype","mut_type","value")
names(table_melt) <- c("genotype","mut_type","value")

Data <- merge(table_melt,p_value_melt,by=c("genotype","mut_type"))

names(Data) <- c("genotype","mut_type","fold_change","p_value")
Data1 <- data.frame(Var1=rep(unique(Data$genotype),time=length(unique(Data$mut_type))),
                    Var2=rep(unique(Data$mut_type),each=length(unique(Data$genotype))))
Data$p_value <- ifelse(Data$p_value < 0.0001,4,ifelse(
  Data$p_value < 0.001,3,ifelse(
    Data$p_value < 0.01,2,ifelse(
      Data$p_value < 0.05,1,0 ))))


p1 <- ggplot(Data,aes(x=mut_type,y=genotype))+
     geom_point(aes(color=fold_change,size=p_value))+
     geom_tile(data=Data1,aes(x=Var2,y=Var1),fill=NA,color="gray")+
    scale_color_gradient2(limit=c(-max(abs(Data$fold_change)),max(abs(Data$fold_change))),
                                                    low="blue",mid="white",high = "red",name="Fold change")+
     scale_size_continuous(range=c(1,10),breaks=0:4,labels=c("0-0.05","0.05-0.01","0.001-0.01","0.0001-0.001","<0.0001"),name="p value")+
     scale_x_discrete(position="top",limit=c("mut freq","A-T mut","C:G tansversion","C:G transition"))+
scale_y_discrete(limit=rev(c("IgHWT.12_MSH2...","IgHWT.12_MLH1...","IgHWT.12_Pms2...","IgHWT.12_Polh...","IgHWT.12_Exo1...","IgHWT.12_Ape2...","IgHWT.12_UNG...","IgHWT.12_53BP1...","IgHWT.12_ATM...","IgHWT.12_H2AX...","IgHWT.12_XLF..")),labels=rev(c("MSH2","MLH1","Pms2","Polh","Exo1","Ape2","UNG","53BP1","ATM","H2AX","XLF")))+
theme(panel.background = element_rect(fill="white",color=NA),
                     panel.grid = element_blank(),axis.ticks = element_blank())+
     geom_point(data=Data[Data$fold_change>1.5 | Data$fold_change < -1.5,],aes(size=p_value),shape=1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0))+theme(legend.position='none')




p2 <- ggplot(Data,aes(x=mut_type,y=genotype))+
     geom_point(aes(color=fold_change,size=p_value))+
     geom_tile(data=Data1,aes(x=Var2,y=Var1),fill=NA,color="gray")+
    scale_color_gradient2(limit=c(-max(abs(Data$fold_change)),max(abs(Data$fold_change))),
                                                    low="blue",mid="white",high = "red",name="Fold change")+
     scale_size_continuous(range=c(1,10),breaks=0:4,labels=c("0-0.05","0.05-0.01","0.001-0.01","0.0001-0.001","<0.0001"),name="p value")+
     scale_x_discrete(position="top",limit=c("1","2_20",">20"),labels=c("1","2_20",">20                   "))+
    theme(panel.background = element_rect(fill="white",color=NA),
                     panel.grid = element_blank(),axis.ticks = element_blank())+
     geom_point(data=Data[Data$fold_change>1.5 | Data$fold_change < -1.5,],aes(size=p_value),shape=1)+
theme(axis.text.x = element_text(angle = 45, hjust = 0))+
scale_y_discrete(limit=rev(c("IgHWT.12_MSH2...","IgHWT.12_MLH1...","IgHWT.12_Pms2...","IgHWT.12_Polh...","IgHWT.12_Exo1...","IgHWT.12_Ape2...","IgHWT.12_UNG...","IgHWT.12_53BP1...","IgHWT.12_ATM...","IgHWT.12_H2AX...","IgHWT.12_XLF..")),labels=rev(c("MSH2","MLH1","Pms2","Polh","Exo1","Ape2","UNG","53BP1","ATM","H2AX","XLF")))+theme(legend.position='none')


p3 <- ggplot(Data,aes(x=mut_type,y=genotype))+
     geom_point(aes(color=fold_change,size=p_value))+
     geom_tile(data=Data1,aes(x=Var2,y=Var1),fill=NA,color="gray")+
    scale_color_gradient2(limit=c(-max(abs(Data$fold_change)),max(abs(Data$fold_change))),
                                                    low="blue",mid="white",high = "red",name="Fold change")+
     scale_size_continuous(range=c(1,10),breaks=0:4,labels=c("0-0.05","0.05-0.01","0.001-0.01","0.0001-0.001","<0.0001"),name="p value")+
     scale_x_discrete(position="top",limit=c("CSR"),labels=c("% IgG1              "))+
	theme(panel.background = element_rect(fill="white",color=NA),
                     panel.grid = element_blank(),axis.ticks = element_blank())+
     geom_point(data=Data[Data$fold_change>1.5 | Data$fold_change < -1.5,],aes(size=p_value),shape=1)+
theme(axis.text.x = element_text(angle = 45, hjust = 0))+
scale_y_discrete(limit=rev(c("IgHWT.12_MSH2...","IgHWT.12_MLH1...","IgHWT.12_Pms2...","IgHWT.12_Polh...","IgHWT.12_Exo1...","IgHWT.12_Ape2...","IgHWT.12_UNG...","IgHWT.12_53BP1...","IgHWT.12_ATM...","IgHWT.12_H2AX...","IgHWT.12_XLF..")),labels=rev(c("MSH2","MLH1","Pms2","Polh","Exo1","Ape2","UNG","53BP1","ATM","H2AX","XLF")))


pdf("ggplot.heatmap.1.pdf", width = 9, height = 6)

ggarrange(p1,p2,p3,ncol=3,widths=c(2.3,1.7,2))
dev.off()

#cowplot::plot_grid(plotlist = list(p1,p2,p3),nrows=1)
