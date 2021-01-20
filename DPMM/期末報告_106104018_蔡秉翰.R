
options(scipen=999)
#read data and package
dta<-read.table('data.txt',header=T,sep="")
library(ggplot2)
library(gridExtra)


#normality test
shapiro.test(dta$Valence)
shapiro.test(dta$Arousal)
shapiro.test(dta$Control)

# outlier analysis
ans<-median(dta$Valence)
M<-mad(dta$Valence)
(dta$Valence-ans)/M>2.5
(dta$Valence-ans)/M < -2.5


M2<-mad(dta$Arousal)
ans2<-median(dta$Arousal)
(dta$Arousal-ans2)/M>2.5
(dta$Arousal-ans2)/M< -2.5



M3<-mad(dta$Control)
ans3<-median(dta$Control)
(dta$Control-ans3)/M>2.5
(dta$Control-ans3)/M< -2.5


#draw scatter part of Arousal and valence
ggplot(dta,aes(Arousal,Valence))+geom_point(col='SkyBlue',size=3)+labs(x='Arousal',y='Valence')+scale_y_continuous(breaks=seq(0,9,1))+
  theme_classic()+theme(axis.title.y=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5,vjust=0.6))+
  theme(axis.text=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5))+
  theme(legend.text=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5))+
  theme(legend.title=element_text(family="Times New Roman",size=24,colour='Black'))+
  theme(axis.line = element_line(colour = 'black', size = 2),axis.ticks = element_line(colour = "black", size = 2))+
  theme(axis.title.x=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5,vjust=0.6))


#coorelation test
for (i in c(1:360)){
  dta$Category[i]<-ifelse(dta$Valence[i]>5,'Positive','Negative' )
}

dta$Category<-as.factor(dta$Category)

dta_p<-subset(dta,Category=='Positive')
dta_n<-subset(dta,Category=='Negative')

cor.test(dta_p$Valence,dta_p$Arousal,alternative='greater',method='kendall')
cor.test(dta_n$Valence,dta_n$Arousal,alternative='less',method='kendall')


#read data with cluster
dta_clus<-read.table('data_cluster.txt',header=T,sep='')
dta_clus$Cluster<-as.factor(dta_clus$Cluster)
dta_clus

#draw boxplot 
a<-ggplot(dta_clus,aes(Cluster,Valence,fill=Cluster))+geom_boxplot()+labs(x='Cluster',y='Valence')+theme_classic()+theme(legend.position='none')+scale_y_continuous(breaks=seq(0,9,1))+
  theme(axis.title.y=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5,vjust=0.6))+
  theme(axis.text=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5))+
  theme(axis.line = element_line(colour = 'black', size = 2),axis.ticks = element_line(colour = "black", size = 2))+
  theme(axis.title.x=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5,vjust=0.6))
b<-ggplot(dta_clus,aes(Cluster,Arousal,fill=Cluster))+geom_boxplot()+labs(x='Cluster',y='Arousal')+theme_classic()+theme(legend.position='none')+scale_y_continuous(breaks=seq(0,9,1))+
  theme(axis.title.y=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5,vjust=0.6))+
  theme(axis.text=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5))+
  theme(axis.line = element_line(colour = 'black', size = 2),axis.ticks = element_line(colour = "black", size = 2))+
  theme(axis.title.x=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5,vjust=0.6))
c<-ggplot(dta_clus,aes(Cluster,Dominance,fill=Cluster))+geom_boxplot()+labs(x='Cluster',y='Dominance')+theme_classic()+theme(legend.position='none')+scale_y_continuous(breaks=seq(0,9,1))+
  theme(axis.title.y=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5,vjust=0.6))+
  theme(axis.text=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5))+
  theme(axis.line = element_line(colour = 'black', size = 2),axis.ticks = element_line(colour = "black", size = 2))+
  theme(axis.title.x=element_text(family="Times New Roman",size=24,colour='Black',hjust=0.5,vjust=0.6))
grid.arrange(a,b,c,ncol=3)




