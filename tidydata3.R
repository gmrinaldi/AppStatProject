setwd("C:/Users/nb2/Desktop/Politecnico/StatApp/AppStatProject")
library(tidyverse)
library(reshape2)
library(magrittr)

# Importo i dati e li pulisco
# Nota: per vedere meglio cosa succede, usate il comando head() dopo ogni passaggio
mydata<-read.csv("Data2.csv", header=T, sep=";")%>%select(1:43)

colnames(mydata)[1]<-"Time"
# "Sciolgo" il dataset, utile ad esempio per i grafici
reshaped<-melt(mydata, id=c("Time"), variable.name="Treatment",value.name="Impedance")

# Creo una "lookup table", in modo da associare nuove variabili ai vari trattamenti (GF sì/no, 
# Inibitore, Concentrazione, Numero della misura)
treatNames<-unique(reshaped$Treatment)
treatGF<-c(rep(TRUE,39),rep(FALSE,3))
treatInhib<-c(rep("no",4),rep("A",11),rep("B",12),rep("A+B",12),rep("no",3))
treatConc<-c(rep("C0",4), rep(c(rep("C4",3),rep("C1",3),rep("C025",3),rep("C006",3)),3),rep("C0",3))%>%.[-16]
treatMeas<-c(1:4,rep(1:3,13))%>%.[-16]%>%as.factor()

lookUp<-data.frame(treatNames,treatInhib,treatConc,treatGF,treatMeas) %>%
    plyr::rename(c("treatNames"="Treatment","treatInhib"="Inhibitor","treatConc"="Concentration","treatMeas"="Measure"))

# creo nuove colonne con le nuove informazioni usando la lookup table (ed elimino la colonna Treatment ormai inutile)
processed<-inner_join(reshaped,lookUp,by="Treatment")%>%
  .[c("Treatment","treatGF","Inhibitor","Concentration","Measure","Time","Impedance")]%>%
  select(-Treatment)

# metto il tempo nel formato giusto
processed$Time<-floor(processed$Time)+(processed$Time-floor(processed$Time))*100/60

# Grafici!
library(wesanderson)
library(RColorBrewer)

# processed_plot<-processed
processed_plot<-processed%>%group_by(Time)%>%mutate(base=mean(Impedance[treatGF==F]),baseGF=mean(Impedance[treatGF==T & Inhibitor=="no"]))%>%ungroup()


p1 <- ggplot(processed_plot%>%filter(Inhibitor!="no"),aes(x = Time, y=Impedance, group=interaction(Inhibitor,Concentration,Measure)))
p1 + geom_line(aes(color = Inhibitor),size=1)+geom_line(aes(x=Time,y=baseGF),linetype="longdash",size=1.25)+geom_line(aes(x=Time,y=base),linetype="longdash",size=1.25)+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

(p2 <- p1 + geom_line(aes(x=Time,y=base),size=.7,linetype="dashed",colour="lemonchiffon3")+geom_line(aes(x=Time,y=baseGF),size=.7,linetype="dashed",colour="lemonchiffon3")+
    geom_line(aes(color=Measure),size=.85) + scale_color_manual(values=wes_palette(n=3, name="Zissou1"))+
    facet_grid(factor(Inhibitor,levels=c("A","B","A+B"))~Concentration))+theme_minimal()

p3<- ggplot(filter(processed_plot,Inhibitor!="no"),aes(x = Concentration, y=Impedance))
p3 + geom_boxplot(aes(fill=Inhibitor))+facet_wrap(~Inhibitor)+theme_minimal()

p4 <- ggplot(filter(processed_plot,Inhibitor!="no"),aes(x = Inhibitor, y=Impedance))
p4 +geom_boxplot(aes(fill=Inhibitor))+facet_wrap(~Concentration, ncol=4) +theme_minimal()

# PCA

library(broom)
library(ggfortify)

# Preparo il dataset
processed_ext<-processed%>%
  spread(key=Time,value=Impedance)

# Faccio la PCA
performedPCA<-processed_ext%>%
      select(-(1:4)) %>% prcomp(scale=FALSE) 


# Per vedere la variabilità spiegata dalle varie componenti principali
performedPCA%>%.$x%>%as.data.frame()%>%
  summarize_all(funs(var)) %>%
    gather(key = pc, value = variance) %>% 
      mutate(var_exp = variance/sum(variance))

# Grafico dei loadings

load.data <- performedPCA$rotation %>% as.data.frame()

p5<-ggplot(load.data%>%rownames_to_column(var="Time")%>%select(c(Time,PC1,PC2,PC3))%>%melt(id=c("Time"),value.name = "Loadings"), aes(x=Time,y=Loadings))
p5+geom_bar(stat="identity",width=1,color="black",fill="snow1")+facet_wrap(~variable,nrow=3)+theme_minimal()+theme(axis.text.x = element_blank())

# Grafico delle prime due PC
autoplot(performedPCA,data=processed_ext,
         colour="Inhibitor",label=T, label.label="Concentration", label.repel=T)+theme_minimal()

# Clusters (ancora da fare, bisogna decidere che cosa usare)
# Nota: possiamo usare i cluster per aggiungere il dato mancante! 

# Calcolo dei livelli e anova

Levels <- processed%>%group_by(treatGF,Inhibitor,Concentration,Measure)%>%
  summarize(Level=mean(Impedance))

model<-aov(Level~Inhibitor+Concentration,data=Levels%>%filter(Inhibitor!="no",Inhibitor!="B")%>%droplevels())
anova(model)

# Analisi del fattore forma
MeanProfiles<- processed_ext %>% group_by(treatGF, Inhibitor, Concentration) %>%
        summarize_at(vars(-Measure),funs(mean))%>%ungroup()

mean_across_treats<-MeanProfiles%>%summarize_at(vars(-(1:3)),funs(mean))
mean_across_time<-MeanProfiles%>% gather(-(1:3),key="Time",value="Impedance")%>%
      group_by(treatGF,Inhibitor,Concentration)%>%summarize(mean=mean(Impedance))
grand_mean<-MeanProfiles%>% gather(-(1:3),key="Time",value="Impedance")%>%.$Impedance%>%mean()

# Matrice dei residui
ResidualMatr<- MeanProfiles%>%select(-(1:3))%>%as.matrix() %>% -grand_mean-
        (matrix(rep(mean_across_time$mean,111),nrow=14)-grand_mean)-
        (t(matrix(rep(as.numeric(mean_across_treats),14),ncol=14))-grand_mean)

# ResidualMatr_Plot<-MeanProfiles
# ResidualMatr_Plot[,4:114]<-ResidualMatr
# ResidualMatr_Plot<-gather(ResidualMatr_Plot,-(1:3),key="Time",value="Impedance")
# p6 <- ggplot(ResidualMatr_Plot,aes(x = Time, y=Impedance, group=interaction(Inhibitor,Concentration,treatGF)))
# p6 + geom_line(aes(color = Inhibitor),size=1)+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

modPCA<-prcomp(ResidualMatr,scale=F)

modPCA%>%.$x%>%as.data.frame()%>%
  summarize_all(funs(var)) %>%
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance))

load.data <- modPCA$rotation %>% as.data.frame()

p7<-ggplot(load.data%>%rownames_to_column(var="Time")%>%select(c(Time,PC1,PC2,PC3))%>%melt(id=c("Time"),value.name = "Loadings"), aes(x=Time,y=Loadings))
p7+geom_bar(stat="identity",width=1,color="black",fill="snow1")+facet_wrap(~variable,nrow=3)+theme_minimal()+theme(axis.text.x = element_blank())

autoplot(modPCA,data=processed_ext,
         colour="Inhibitor",label=T,label.label="Concentration", label.repel=T)+theme_minimal()

# Posso tenere solo la prima componente principale ~97%
shape_score<-processed_ext%>%mutate(Score=processed_ext%>%select(-(1:4))%>%
                as.matrix()%>% multiply_by_matrix(modPCA$rotation[,1])%>%as.vector())%>%select((1:4),Score)

model<-aov(Score~Inhibitor+Concentration,data=shape_score%>%filter(Inhibitor!="no")%>%droplevels())
anova(model)

# Proviamo a visualizzare la forma
shape_as_factor<-shape_score%>%mutate(Score_factor=round(Score)%>%as.factor())%>%
    select(-Score)%>%inner_join(processed)%>%plyr::rename(c("Score_factor"="Score"))

p8 <- ggplot(shape_as_factor,aes(x = Time, y=Impedance, group=interaction(Inhibitor,Concentration,Measure,treatGF)))
p8 + geom_line(aes(color =Inhibitor),size=.75)+facet_wrap(~Score)+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))+theme_minimal()

# Consideriamo levels e shape insieme
joinedLevelShape<-inner_join(Levels, shape_score)

p9<-ggplot(joinedLevelShape,aes(x=Level,y=Score))
p9 + geom_point(aes(color=Inhibitor))+theme_minimal()

# Fittiamo un modello lineare
model<-lm(Score~poly(Level,2),data=joinedLevelShape)
summary(model)
plot(model)

p10<-ggplot(joinedLevelShape,aes(x=Level,y=Score))
p10 + geom_point(aes(color=Inhibitor))+geom_smooth(se=T,method='lm',formula=y~poly(x,2))+theme_minimal()


# Prova di LME
library(lme4)
processed_lme<-processed%>%group_by(Time)%>%
                            mutate(Sample=1:42)%>%ungroup()

basemodel<-lmer(Impedance~1+Time+(1+Time|Sample),data=processed_lme,REML=F)
interceptmodel<-lmer(Impedance~1+Time+interaction(Inhibitor,Concentration,treatGF)+(1+Time|Sample),data=processed_lme,REML=F)
lineartermmodel<-lmer(Impedance~1+Time+Inhibitor:Time+(1+Time|Sample),data=processed_lme,REML=F)
secondordermodel<-lmer(Impedance~1+poly(Time,2)+interaction(Inhibitor,Concentration,treatGF)*Time+(1+Time|Sample),data=processed_lme,REML=F)


ggplot(processed_lme, aes(Time, Impedance, color=interaction(Inhibitor,Concentration,treatGF))) +
   # stat_summary(fun.data=mean_se, geom="pointrange") +
  stat_summary(aes(y=fitted(basemodel), linetype=Inhibitor),
               fun.y=mean, geom="line")


# Prove con la FDA
library(fda)


mytimes<-processed$Time%>%unique()
mytimes<-(mytimes-mytimes[1])*6
imped_basis<-create.bspline.basis(range(mytimes),114,5,mytimes)
imped_par<-fdPar(imped_basis,3,1e-5)
imped_data<-processed$Impedance%>%matrix(nrow=111)

result<-smooth.monotone(mytimes,imped_data,imped_par)

Wfd<-result$Wfdobj
beta<-result$beta

imped_hat<-t(beta[1,]+beta[2,]*t(eval.monfd(mytimes,Wfd)))
Dimped_hat<-t(beta[2,]*t(eval.monfd(mytimes,Wfd,1)))
D2imped_hat<-t(beta[2]*t(eval.monfd(mytimes,Wfd,2)))

# plotfit.fd(imped_data,mytimes,result$yhatfd,residual=T)

smoothed_processed<-processed%>%mutate(Smoothed=imped_hat%>%as.vector(),FirstDerivative=Dimped_hat%>%as.vector(), SecondDerivative=D2imped_hat%>%as.vector())
smoothed_processed$Time<-(smoothed_processed$Time-smoothed_processed$Time[1])*6
p11 <- ggplot(smoothed_processed,aes(x = Time, y=Impedance))
p11 + geom_point(size=.75)+geom_line(aes(x=Time,y=Smoothed,color=Inhibitor,group=interaction(treatGF,Inhibitor,Concentration,Measure)))+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

p12 <- ggplot(smoothed_processed,aes(x = Time, y = FirstDerivative))
p12 + geom_line(aes(x=Time,y=log10(FirstDerivative),color=Inhibitor,group=interaction(treatGF,Inhibitor,Concentration,Measure)))+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

p13 <- ggplot(smoothed_processed,aes(x = Time, y = SecondDerivative))
p13 + geom_line(aes(x=Time,y=SecondDerivative,color=Inhibitor,group=interaction(treatGF,Inhibitor,Concentration,Measure)))+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

# # Scegliamo lambda (più o meno)
# loglam<-seq(-6,1,.25)
# nlam<- length(loglam)
# dfsave<-rep(NA,nlam)
# gcvsave<-rep(NA,nlam)
# 
# for (ilam in 1:nlam) {
#   cat(paste('log10 lambda =', loglam[ilam],'\n'))
#   lambda <- 10^loglam[ilam]
#   fdParobj <- fdPar(imped_basis,3,lambda)
#   smoothlist<-smooth.basis(mytimes,imped_data,fdParobj)
#   dfsave[ilam]<-smoothlist$df
#   gcvsave[ilam]<-sum(smoothlist$gcv)
# 
# }
# plot(loglam,gcvsave)

# Functional pca
Dimped.pca<-pca.fd(result$yhatfd,2)
Dimped.pca$varprop
plot(Dimped.pca)

Dimped.pca<-pca.fd(result$Wfdobj,4)
Dimped.pca$varprop
plot(Dimped.pca)

# Dimped.rot.pca<-varmx.pca.fd(Dimped.pca)
# plot(Dimped.rot.pca)

# Proviamo a eliminare la parte lineare?

# Clustering
library(fdakma)
kma_result<-kma(x=mytimes,y0=t(imped_hat),y1=t(Dimped_hat),n.clust=1)
kma.show.results(kma_result)
#meglio
library(funHDDC)
clust_result<-funHDDC(result$yhatfd,4)




