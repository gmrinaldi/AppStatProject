setwd("C:/Users/nb2/Desktop/Politecnico/StatApp/AppStatProject")
library(tidyverse)
library(reshape2)
library(magrittr)
library(lme4)
library(fda)

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

# Manova
fit<-manova(cbind(joinedLevelShape$Level,joinedLevelShape$Score)~Inhibitor*Concentration,data=joinedLevelShape)
summary.manova(fit)
summary.aov(fit)

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

# plotfit.fd(imped_data,mytimes,result$yhatfd,residual=F)

smoothed_processed<-processed%>%mutate(Smoothed=imped_hat%>%as.vector(),FirstDerivative=Dimped_hat%>%as.vector(), SecondDerivative=D2imped_hat%>%as.vector())
smoothed_processed$Time<-(smoothed_processed$Time-smoothed_processed$Time[1])*6
p11 <- ggplot(smoothed_processed,aes(x = Time, y=Impedance))
p11 + geom_point(size=.75)+geom_line(aes(x=Time,y=Smoothed,color=Inhibitor,group=interaction(treatGF,Inhibitor,Concentration,Measure)))+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

p12 <- ggplot(smoothed_processed,aes(x = Time, y = FirstDerivative))
p12 + geom_line(aes(x=Time,y=log10(FirstDerivative),color=Inhibitor,group=interaction(treatGF,Inhibitor,Concentration,Measure)))+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))+labs(y="Log10 of First Derivative")

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
imped.pca<-pca.fd(result$yhatfd,2)
imped.pca$varprop
plot(imped.pca)

Dimped.pca<-pca.fd(result$Wfdobj,4)
Dimped.pca$varprop
plot(Dimped.pca)

# Hierarchical model
require(lme4)
require(lattice)
scores.pca1<-imped.pca$scores[,1]
scores.pca2<-imped.pca$scores[,2]

processed_lme<-processed_ext%>%select((1:4))%>%arrange(desc(treatGF),desc(Inhibitor))
processed_lme[5:42,]<-arrange(processed_lme[5:42,],Inhibitor,desc(Concentration))
processed_lme[16:39,]<-arrange(processed_lme[16:39,],desc(Inhibitor))
processed_lme<-processed_lme%>%mutate(Scores=scores.pca1,Scores2=scores.pca2)%>%filter(Inhibitor!="no")%>%select(-treatGF)
processed_lme$Concentration<-plyr::revalue(processed_lme$Concentration, c(C0="0",C006="0.06",C025="0.25",C1="1",C4="4"))%>%as.character()%>%as.numeric()
processed_lme$Concentration<-log(processed_lme$Concentration)
# processed_lme<-processed_lme%>%mutate(ID=c(1,1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12))

basemodel<-lmer(Scores~1+poly(Concentration,2)*Inhibitor+(1+poly(Concentration,2)|Measure),data=processed_lme,REML=F)

ggplot(processed_lme, aes(Concentration, Scores, color=Inhibitor)) +
  stat_summary(aes(y=fitted(basemodel)), fun.y=mean, geom="line")+
    geom_point(aes(x=Concentration,y=Scores))+theme_minimal()

lmodel<-lm(Scores~poly(Concentration,2)*Inhibitor,data=processed_lme)

ggplot(processed_lme, aes(Concentration, Scores, color=Inhibitor)) +
  stat_summary(aes(y=fitted(lmodel)), fun.y=mean, geom="line")+
  geom_point(aes(x=Concentration,y=Scores))+theme_minimal()

x11()
par(mfrow=c(2,2))
plot(lmodel)

plot(basemodel)
plot(basemodel, type = c("p", "smooth"))
plot(basemodel, sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))
qqmath(basemodel,id=.05)
shapiro.test(residuals(basemodel))
iqrvec<-sapply(simulate(basemodel,1000),IQR)
obsval<-IQR(processed_lme$Scores)
post.pred.p<-mean(obsval>=c(obsval,iqrvec))

# lme_var<-processed_lme%>%group_by(Inhibitor,Concentration)%>%summarise(var=sd(Scores)^2)
# lme_var<-inner_join(lme_var,processed_lme)
# weightedmodel<-lmer(Scores~1+Concentration*Inhibitor+(1+Concentration|Measure),weights=1/var,data=lme_var,REML=F)
# 
# ggplot(lme_var, aes(Concentration, Scores, color=Inhibitor)) +
#   stat_summary(aes(y=fitted(weightedmodel)), fun.y=mean, geom="line")+
#     geom_point(aes(x=Concentration,y=Scores))+theme_minimal()
# 
# # weightedlinearmodel<-lm(Scores~1+Concentration*Inhibitor,weights=1/var,data=lme_var)
# # 
# # ggplot(lme_var, aes(Concentration, Scores, color=Inhibitor)) +
# #   stat_summary(aes(y=fitted(weightedlinearmodel)), fun.y=mean, geom="line")+
# #   geom_point(aes(x=Concentration,y=Scores))+theme_minimal()
# 
# plot(weightedmodel)
# plot(weightedmodel, type = c("p", "smooth"))
# plot(weightedmodel, sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))
# qqmath(weightedmodel,id=.05)
# 
# iqrvec<-sapply(simulate(weightedmodel,10000),median)
# obsval<-median(lme_var$Scores)
# post.pred.p<-mean(obsval>=c(obsval,iqrvec))


basemodel2<-lmer(Scores2~1+poly(Concentration,2)*Inhibitor+(1|Measure),data=processed_lme,REML=F)

ggplot(processed_lme, aes(Concentration, Scores2, color=Inhibitor)) +
  stat_summary(aes(y=fitted(basemodel2)), fun.y=mean, geom="line")+
  geom_point(aes(x=Concentration,y=Scores2))+theme_minimal()


plot(basemodel2)
plot(basemodel2, type = c("p", "smooth"))
plot(basemodel2, sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))
qqmath(basemodel2,id=.05)

# lme_var2<-processed_lme%>%group_by(Inhibitor,Concentration)%>%summarise(var=sd(Scores2)^2)
# lme_var2<-inner_join(lme_var2,processed_lme)
# weightedmodel2<-lmer(Scores2~1+poly(Concentration,2)+(1+poly(Concentration,2)|Inhibitor),weights=1/var,data=lme_var2%>%filter(Inhibitor!="B"),REML=F)
# 
# ggplot(lme_var2%>%filter(Inhibitor!="B"), aes(Concentration, Scores2, color=Inhibitor)) +
#   stat_summary(aes(y=fitted(weightedmodel2)), fun.y=mean, geom="line")+
#   geom_point(aes(x=Concentration,y=Scores2))+theme_minimal()
# 
# plot(weightedmodel2)
# plot(weightedmodel2, type = c("p", "smooth"))
# plot(weightedmodel2, sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))
# qqmath(weightedmodel2,id=.05)

mean_curve<-mean(result$yhatfd)

load("myFunction.Rdata")
d_p<-downscale_processed(processed,"A",c(.125,.5,2),basemodel,basemodel2,mean_curve,imped.pca)

lookUp2<-data.frame(Concentration=c(unique(d_p$Concentration)),Estimated=c(rep("FALSE",4),rep("TRUE",3)))
d_plot<-inner_join(d_p,lookUp2,by="Concentration")

ggplot(d_plot,aes(x=Time,y=Impedance,color=Concentration%>%as.factor(),group=Concentration))+geom_line(aes(linetype=Estimated))+theme_minimal()+guides(linetype=FALSE)+scale_color_discrete(name="Concentration")

# Dimped.rot.pca<-varmx.pca.fd(Dimped.pca)
# plot(Dimped.rot.pca)

# Proviamo a eliminare la parte lineare?

# Clustering
library(fdakma)
kma_result<-kma(x=mytimes,y0=t(imped_hat),y1=t(Dimped_hat),n.clust=1)
kma.show.results(kma_result)

kma_result_der<-kma(x=mytimes,y0=t(log(Dimped_hat)),y1=t(D2imped_hat/Dimped_hat),n.clust=2)
kma.show.results(kma_result_der)

#meglio
library(funHDDC)
clust_result<-funHDDC(result$yhatfd,5)
plot(result$yhatfd,col=clust_result$class)

clustered_processed<-smoothed_processed%>%mutate(Cluster=rep(clust_result$class,each=111))
p14 <- ggplot(clustered_processed%>%filter(Inhibitor=="A",Concentration=="C006"),aes(x = Time, y=Impedance))
p14 +geom_line(aes(colour=Cluster%>%as.factor(),group=interaction(treatGF,Inhibitor,Concentration,Measure)))+theme_minimal()

p15<- ggplot(processed_lme,aes(x=Scores,y=Scores2))
p15+geom_point(aes(color=Inhibitor))

# Prova di FLM
Treatments<-processed$Inhibitor%>%unique()
p<-length(Treatments)+1
TreatmentsList<-vector("list",p)
TreatmentsList[[1]]<-c(rep(1,42),0)
for(j in 2:p) {
  xj=processed_ext$Inhibitor==Treatments[j-1]
  TreatmentsList[[j]]<-c(xj,1)
}

coef<-result$yhatfd$coefs
coef43<-cbind(coef,matrix(0,114,1))
Imped43fd<-fd(coef43,imped_basis)

betabasis<-create.bspline.basis(range(mytimes),114,5,mytimes)
betafdPar<-fdPar(betabasis)
betaList<-vector("list",p)
for(j in 1:p) betaList[[j]]<-betafdPar

fRegressList<-fRegress(Imped43fd,TreatmentsList,betaList)

betaestList<-fRegressList$betaestlist
treatmentsFit<-fRegressList$yhatfdobj
treatments<-c("Mean",Treatments)
par(mfrow=c(2,3),cex=1)
for (j in 1:p) plot(betaestList[[j]]$fd,lwd=2,main=Treatments[j])
plot(treatmentsFit,lwd=2,col=1,lty=1)

# Missing value
load("missingA_C006.RData")
processed_nomissing<-rbind(processed,data.frame(treatGF=T,Inhibitor="A",Concentration="C006",
                                                Measure=3,Time=unique(processed$Time),Impedance=missingobs))

processed_plot2<-processed_nomissing%>%group_by(Time)%>%mutate(base=mean(Impedance[treatGF==F]),baseGF=mean(Impedance[treatGF==T & Inhibitor=="no"]))%>%ungroup()


p16 <- ggplot(processed_plot2%>%filter(Inhibitor!="no"),aes(x = Time, y=Impedance, group=interaction(Inhibitor,Concentration,Measure)))
p16 + geom_line(aes(color = Inhibitor),size=1)+geom_line(aes(x=Time,y=baseGF),linetype="longdash",size=1.25)+geom_line(aes(x=Time,y=base),linetype="longdash",size=1.25)+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

(p17 <- p16 + geom_line(aes(x=Time,y=base),size=.7,linetype="dashed",colour="lemonchiffon3")+geom_line(aes(x=Time,y=baseGF),size=.7,linetype="dashed",colour="lemonchiffon3")+
    geom_line(aes(color=Measure),size=.85) + scale_color_manual(values=wes_palette(n=3, name="Zissou1"))+
    facet_grid(factor(Inhibitor,levels=c("A","B","A+B"))~Concentration))+theme_minimal()

mytimes<-processed$Time%>%unique()
mytimes<-(mytimes-mytimes[1])*6
imped_basis<-create.bspline.basis(range(mytimes),114,5,mytimes)
imped_par<-fdPar(imped_basis,3,1e-5)
imped_data<-processed_nomissing$Impedance%>%matrix(nrow=111)

# result<-smooth.monotone(mytimes,imped_data,imped_par)
load("fitted_spline.RData")

Wfd<-result$Wfdobj
beta<-result$beta

imped_hat<-t(beta[1,]+beta[2,]*t(eval.monfd(mytimes,Wfd)))
Dimped_hat<-t(beta[2,]*t(eval.monfd(mytimes,Wfd,1)))
D2imped_hat<-t(beta[2]*t(eval.monfd(mytimes,Wfd,2)))

# plotfit.fd(imped_data,mytimes,result$yhatfd,residual=F)

smoothed_processed<-processed_nomissing%>%mutate(Smoothed=imped_hat%>%as.vector(),FirstDerivative=Dimped_hat%>%as.vector(), SecondDerivative=D2imped_hat%>%as.vector())
smoothed_processed$Time<-(smoothed_processed$Time-smoothed_processed$Time[1])*6
p18 <- ggplot(smoothed_processed,aes(x = Time, y=Impedance))
p18 + geom_point(size=.75)+geom_line(aes(x=Time,y=Smoothed,color=Inhibitor,group=interaction(treatGF,Inhibitor,Concentration,Measure)))+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

p19 <- ggplot(smoothed_processed,aes(x = Time, y = FirstDerivative))
p19 + geom_line(aes(x=Time,y=log10(FirstDerivative),color=Inhibitor,group=interaction(treatGF,Inhibitor,Concentration,Measure)))+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))+labs(y="Log10 of First Derivative")

p20 <- ggplot(smoothed_processed,aes(x = Time, y = SecondDerivative))
p20 + geom_line(aes(x=Time,y=SecondDerivative,color=Inhibitor,group=interaction(treatGF,Inhibitor,Concentration,Measure)))+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

imped.pca<-pca.fd(result$yhatfd,2)
imped.pca$varprop
plot(imped.pca)
plotscores(imped.pca)

Dimped.pca<-pca.fd(result$Wfdobj,4)
Dimped.pca$varprop
plot(Dimped.pca)



# Hierarchical model
scores.pca1<-imped.pca$scores[,1]
scores.pca2<-imped.pca$scores[,2]

processed_lme<-processed_nomissing%>%spread(key=Time,value=Impedance)%>%select((1:4))%>%arrange(desc(treatGF),desc(Inhibitor))
processed_lme[5:43,]<-arrange(processed_lme[5:43,],Inhibitor,desc(Concentration))
processed_lme[17:40,]<-arrange(processed_lme[17:40,],desc(Inhibitor))
processed_lme<-processed_lme%>%mutate(Scores=scores.pca1,Scores2=scores.pca2)%>%select(-treatGF)
temprow<-processed_lme[43,4:5]
processed_lme[17:43,4:5]<-processed_lme[16:42,4:5]
processed_lme[16,4:5]<-temprow
processed_lme<-processed_lme%>%filter(Inhibitor!="no",Inhibitor!="untreated")
processed_lme$Concentration<-plyr::revalue(processed_lme$Concentration, c(C0="0",C006="0.06",C025="0.25",C1="1",C4="4"))%>%as.character()%>%as.numeric()
processed_lme$Concentration<-log(processed_lme$Concentration)
# processed_lme<-processed_lme%>%mutate(ID=c(1,1,1,2,2,2,3,3,3,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12))

basemodel<-lmer(Scores~1+poly(Concentration,2)*Inhibitor+(1|Measure),data=processed_lme,REML=F)

ggplot(processed_lme, aes(Concentration, Scores, color=Inhibitor)) +
  stat_summary(aes(y=fitted(basemodel)), fun.y=mean, geom="line")+
  geom_point(aes(x=Concentration,y=Scores))+theme_minimal()

shapiro.test(residuals(basemodel))

plot(basemodel)
plot(basemodel, type = c("p", "smooth"))
plot(basemodel, sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))
qqmath(basemodel,id=.05)

basemodel2<-lmer(Scores2~1+poly(Concentration,2)*Inhibitor+(1|Measure),data=processed_lme,REML=F)

ggplot(processed_lme, aes(Concentration, Scores2, color=Inhibitor)) +
  stat_summary(aes(y=fitted(basemodel2)), fun.y=mean, geom="line")+
  geom_point(aes(x=Concentration,y=Scores2))+theme_minimal()


plot(basemodel2)
plot(basemodel2, type = c("p", "smooth"))
plot(basemodel2, sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))
qqmath(basemodel2,id=.05)

mean_curve<-mean(result$yhatfd)

load("myFunction.Rdata")
d_p<-downscale_processed(processed_nomissing,"A+B",c(.125,.5,2),basemodel,basemodel2,mean_curve,imped.pca)

lookUp2<-data.frame(Concentration=c(unique(d_p$Concentration)),Estimated=c(rep("FALSE",4),rep("TRUE",3)))
d_plot<-inner_join(d_p,lookUp2,by="Concentration")

ggplot(d_plot,aes(x=Time,y=Impedance,color=Concentration%>%as.factor(),group=Concentration))+geom_line(aes(linetype=Estimated))+theme_minimal()+guides(linetype=FALSE)+scale_color_discrete(name="Concentration")

# Grafico dato mancante
pmissing<-ggplot(processed_nomissing%>%filter(Inhibitor=="A",Concentration=="C006")%>%mutate(Estimated=Measure==3),aes(x=Time,y=Impedance,group=interaction(Inhibitor,Concentration,Measure)))+guides(Color="none")
pmissing+geom_line(aes(linetype=Estimated,color=Measure))+theme_minimal()
# Prova di FLM
levels(processed_nomissing$Inhibitor)<-c(levels(processed_nomissing$Inhibitor),"untreated")
processed_nomissing[processed_nomissing$treatGF==F,]$Inhibitor<-"untreated"

Treatments<-processed_nomissing$Inhibitor%>%unique()
p<-length(Treatments)+1
TreatmentsList<-vector("list",p)
TreatmentsList[[1]]<-c(rep(1,43),0)
Tnames<-c(rep("no",4),rep("A",12),rep("B",12),rep("A+B",12),rep("untreated",3))

for(j in 2:p) {
  xj=Tnames==Treatments[j-1]
  TreatmentsList[[j]]<-c(xj,1)
}

coef<-result$yhatfd$coefs
coef44<-cbind(coef,matrix(0,114,1))
Imped44fd<-fd(coef44,imped_basis)

betabasis<-create.bspline.basis(range(mytimes),114,5,mytimes)
betafdPar<-fdPar(betabasis)
betaList<-vector("list",p)
for(j in 1:p) {betaList[[j]]<-betafdPar
  betaList[[j]]$lambda<-10}


fRegressList<-fRegress(Imped44fd,TreatmentsList,betaList)

betaestList<-fRegressList$betaestlist
treatmentsFit<-fRegressList$yhatfdobj
treatments<-c("Mean","GF","A","B","A+B","untreated")
par(mfrow=c(2,3),cex=1)
for (j in 1:p) plot(betaestList[[j]]$fd,lwd=2,main=treatments[j])
plot(treatmentsFit,lwd=2,col=1,lty=1)
plot(betaestList[[1]])
plot(betaestList[[1]]$fd+betaestList[[2]]$fd,ylim=c(0,2.5))
plot(betaestList[[1]]$fd+betaestList[[3]]$fd,ylim=c(0,2.5))
plot(betaestList[[1]]$fd+betaestList[[4]]$fd,ylim=c(0,2.5))
plot(betaestList[[1]]$fd+betaestList[[5]]$fd,ylim=c(0,2.5))
plot(betaestList[[1]]$fd+betaestList[[6]]$fd,ylim=c(0,2.5))

fittedval_flm<-NULL
for(j in 2:6)
  fittedval_flm<-c(fittedval_flm,eval.fd(betaestList[[1]]$fd+betaestList[[j]]$fd,mytimes))

plot_flm<-data.frame(Treatment=rep(treatments[-1],each=111),Time=rep(mytimes,5),Impedance=fittedval_flm)
myplot<-ggplot(plot_flm,aes(x=Time,y=Impedance))
myplot+geom_line(size=.75,aes(color=Treatment))+theme_minimal()

yhatmat  <- predict(fRegressList$yhatfdobj,mytimes)
ymat     <- eval.fd(mytimes,Imped44fd)
temprmat <- ymat[,1:43] - yhatmat[,1:43]
SigmaE <- var(t(temprmat))

par(mfrow=c(1,1))
contour(SigmaE, xlab="Day", ylab="Day")
lines(c(0, 365), c(0, 365),lty=4)

par(mfrow=c(1,1), mar=c(5,5,3,2), pty="m")
stddevE <- sqrt(diag(SigmaE))
plot(mytimes, stddevE, type="l",
     xlab="Day", ylab="Standard error (deg C)")

#  Repeat regression, this time outputting results for
#  confidence intervals
ciao<-smooth.basis(mytimes,imped_data,imped_par)


stderrList <- fRegress.stderr(fRegressList, ciao$y2cMap, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  plot regression function standard errors

op <- par(mfrow=c(2,3), pty="s")
for (j in 1:p) {
  betastderrj <- eval.fd(mytimes, betastderrlist[[j]])
  plot(mytimes, betastderrj,
       type="l",lty=1, xlab="Day", ylab="Reg. Coeff.",
       main=treatments[j])
}
par(op)

#  plot regression functions with confidence limits

op <- par(mfrow=c(2,3))
for (j in 1:p) {
  betafdParj  <- betaestList[[j]]
  betafdj     <- betafdParj$fd
  betaj       <- eval.fd(mytimes, betafdj)
  betastderrj <- eval.fd(mytimes, betastderrlist[[j]])
  matplot(mytimes, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
          type="l",lty=c(1,4,4), xlab="Day", ylab="Reg. Coeff.",
          main=treatments[j])
}
par(op)

# Test
F.res = myFperm.fd(Imped44fd, TreatmentsList, betaList,cex.axis=1.5,cex.lab=1.5)

t.res = mytperm.fd(Imped44fd[5:16],Imped44fd[17:28],nperm=1000)
t.res = mytperm.fd(Imped44fd[5:16],Imped44fd[29:40],nperm=1000)
t.res = mytperm.fd(Imped44fd[17:28],Imped44fd[29:40],nperm=1000)

myplot<-ggplot(plot_flm%>%filter(Treatment %in% c("A","B")),aes(x=Time,y=Impedance))
myplot+geom_line(size=.75,aes(color=Treatment))+theme_minimal()+scale_color_discrete(h=c(0,240))
myplot<-ggplot(plot_flm%>%filter(Treatment %in% c("A","A+B")),aes(x=Time,y=Impedance))
myplot+geom_line(size=.75,aes(color=Treatment))+theme_minimal()+scale_color_discrete(h=c(0,120))
myplot<-ggplot(plot_flm%>%filter(Treatment %in% c("B","A+B")),aes(x=Time,y=Impedance))
myplot+geom_line(size=.75,aes(color=Treatment))+theme_minimal()+scale_color_discrete(h=c(120,240))

# Proviamo a normalizzare
processed_norm<-processed_nomissing%>%group_by(Time)%>%mutate(base=mean(Impedance[treatGF==F]),baseGF=mean(Impedance[treatGF==T & Inhibitor=="no"]))%>%
  mutate(normalized=(Impedance-base)/(baseGF-base))%>%ungroup()%>%filter(Inhibitor!="no")

imped_data_norm<-processed_norm$normalized%>%matrix(nrow=111)

imped_par_norm<-fdPar(imped_basis,3,10^-1.5)
result_norm<-smooth.basis(mytimes,imped_data_norm,imped_par_norm)

Treatments<-processed_norm$Inhibitor%>%unique()
p<-length(Treatments)+1
TreatmentsList<-vector("list",p)
TreatmentsList[[1]]<-c(rep(1,36),0)
for(j in 2:p) {
  xj=processed_norm%>%select(-base,-baseGF,-Impedance)%>%spread(key=Time,value=normalized)%>%.$Inhibitor==Treatments[j-1]
  TreatmentsList[[j]]<-c(xj,1)
}

coef<-result_norm$fd$coefs
coef44<-cbind(coef,matrix(0,114,1))
Imped44fd<-fd(coef44,imped_basis)

betabasis<-create.bspline.basis(range(mytimes),114,5,mytimes)
betafdPar<-fdPar(betabasis)
betaList<-vector("list",p)
for(j in 1:p) betaList[[j]]<-betafdPar

for (j in 1:p)
  betaList[[j]]$lambda<-10^4

fRegressList<-fRegress(Imped44fd,TreatmentsList,betaList)

betaestList<-fRegressList$betaestlist
treatmentsFit<-fRegressList$yhatfdobj
treatments<-c("Mean","A","B","A+B")
par(mfrow=c(2,3),cex=1)
for (j in 1:p) plot(betaestList[[j]]$fd,lwd=2,main=treatments[j])

plot(treatmentsFit,lwd=2,col=1,lty=1,main="Prediction")


mytperm.fd <- function(x1fd,x2fd,nperm=200,q=0.05,argvals=NULL,plotres=TRUE,...) # first and second 
{                                                                          # groups of data,
  if( !is.fd(x1fd) | !is.fd(x2fd) ){                                     # number permuts
    stop("x1fd and x2fd must both be functional data objects")         # quantile
  }                                                                      # where to evaluate
  # do I plot
  rangeobs = x1fd$basis$range
  rangehat = x2fd$basis$range
  
  
  if( !prod(rangeobs == rangehat) ){
    stop("x1fd and x2fd do not have the same range.")
  }
  
  if(is.null(argvals)){
    argvals = seq(rangeobs[1],rangeobs[2],length.out=101)
  }
  
  q = 1-q
  
  x1mat = eval.fd(argvals,x1fd)
  x2mat = eval.fd(argvals,x2fd)
  
  n1 = ncol(x1mat)
  n2 = ncol(x2mat)
  
  Xmat = cbind(x1mat,x2mat)
  
  Tnull = rep(0,nperm)
  
  Tnullvals = matrix(0,length(argvals),nperm)
  
  for(i in 1:nperm){
    tXmat = Xmat[,sample(n1+n2)]
    
    tmean1 = apply(tXmat[,1:n1],1,mean)
    tmean2 = apply(tXmat[,n1+(1:n2)],1,mean)
    
    tvar1 = apply(tXmat[,1:n1],1,var)/n1
    tvar2 = apply(tXmat[,n1+(1:n2)],1,var)/n2
    
    Tnullvals[,i] = abs(tmean1-tmean2)/sqrt(tvar1+tvar2)
    Tnull[i] = max(Tnullvals[,i])
  }
  
  mean1 = apply(Xmat[,1:n1],1,mean)
  mean2 = apply(Xmat[,n1+(1:n2)],1,mean)
  
  var1 = apply(Xmat[,1:n1],1,var)/n1
  var2 = apply(Xmat[,n1+(1:n2)],1,var)/n2
  
  Tvals = abs(mean1-mean2)/sqrt(var1+var2)
  Tobs = max(Tvals)
  
  pval = mean( Tobs < Tnull )
  qval = quantile(Tnull,q)
  
  pvals.pts = apply(Tvals<Tnullvals,1,mean)
  qvals.pts = apply(Tnullvals,1,quantile,q)
  
  if(plotres){
    
    if( is.null(names(x1fd$fdnames)) | is.null(names(x2fd$fdnames)) ){
      xlab='argvals'
    }	
    else if( prod(names(x1fd$fdnames)[1] == names(x2fd$fdnames)[1]) ){
      xlab = names(x1fd$fdnames)[1]
    }
    else{ xlab = 'argvals' }
    
    ylims = c( min(Tvals,qvals.pts),max(Tobs,qval))
    
    plot(argvals,Tvals,type='l',col=2,ylim=ylims,lwd=2,
         xlab=xlab,ylab='t-statistic',...)
    lines(argvals,qvals.pts,lty=3,col=4,lwd=2)
    abline(h=qval,lty=2,col=4,lwd=2)
    
    legendstr = c('Observed Statistic',
                  paste('pointwise',1-q,'critical value'),
                  paste('maximum',1-q,'critical value'))

    legend("bottomright",legend=legendstr,col=c(2,4,4),
           lty=c(1,3,2),lwd=c(2,2,2))
  }
  
  
  return( list(pval=pval,qval=qval,Tobs=Tobs,Tnull=Tnull,
               Tvals=Tvals,Tnullvals=Tnullvals,qvals.pts=qvals.pts,
               pvals.pts=pvals.pts,argvals=argvals) )
}

myFperm.fd <- function(yfdPar, xfdlist, betalist, wt=NULL,
                     nperm=200,argvals=NULL,q=0.05,plotres=TRUE,...)
{
  # yfdpar, xfdList, betalist, wt = standard inputs to fRegress
  # nperm = number of permutations,
  # argvals = where to evaluate functional responses,
  # q =  quantile to compare
  # plotres:  Do we plot the results?
  
  Fnull     = rep(0,nperm)
  Fnullvals = c()
  
  q = 1-q
  
  begin <- proc.time()
  fRegressList <- fRegress(yfdPar, xfdlist, betalist)
  elapsed.time <- max(proc.time()-begin,na.rm=TRUE)
  
  if( elapsed.time > 30/nperm ){
    print(paste('Estimated Computing time =',
                round(nperm*elapsed.time),'seconds.'))
  }
  
  yhat <- fRegressList$yhatfdobj
  if(is.list(yhat) && ('fd' %in% names(yhat))) yhat <- yhat$fd
  
  tFstat <- Fstat.fd(yfdPar,yhat,argvals)
  
  Fvals <- tFstat$F
  Fobs = max(Fvals)
  
  argvals = tFstat$argvals
  
  if(is.vector(yfdPar)){ 
    n = length(yfdPar) 
  }  else { 
    n = ncol(yfdPar$coefs) 
  }
  
  for(i in 1:nperm){
    
    tyfdPar = yfdPar[sample(n)]
    
    fRegressList <- fRegress(tyfdPar, xfdlist, betalist)
    
    yhat <- fRegressList$yhatfdobj
    if(is.list(yhat) && ('fd' %in% names(yhat))) yhat <- yhat$fd
    
    tFstat = Fstat.fd(tyfdPar,yhat,argvals)
    
    Fnullvals <- cbind(Fnullvals,tFstat$F)
    
    Fnull[i] = max(Fnullvals[,i])
  }
  
  pval = mean( Fobs < Fnull )
  qval = quantile(Fnull,q)
  
  pvals.pts = apply(Fvals<Fnullvals,1,mean)
  qvals.pts = apply(Fnullvals,1,quantile,q)
  
  if(plotres){
    if(is.fd(yfdPar)){
      ylims = c(min(c(Fvals,qval,qvals.pts)),max(c(Fobs,qval)))
      
      if( is.null(names(yhat$fdnames)) ){ xlab = 'argvals' }
      else{ xlab = names(yhat$fdnames)[1] }
      
      plot(argvals,Fvals,type="l",ylim=ylims,col=2,lwd=2,
           xlab=xlab,ylab='F-statistic',main='Permutation F-Test',...)
      lines(argvals,qvals.pts,lty=3,col=4,lwd=2)
      abline(h=qval,lty=2,col=4,lwd=2)
      
      legendstr = c('Observed Statistic',
                    paste('pointwise',1-q,'critical value'),
                    paste('maximum',1-q,'critical value'))
      
      legend("bottomright",legend=legendstr,col=c(2,4,4),
             lty=c(1,3,2),lwd=c(2,2,2))
    }
    else{
      xlims = c(min(c(Fnull,Fobs)),max(c(Fnull,Fobs)))
      hstat = hist(Fnull,xlim=xlims,lwd=2,xlab='F-value',
                   main = 'Permutation F-Test')
      abline(v = Fobs,col=2,lwd=2)
      abline(v = qval,col=4,lty=2,lwd=2)
      
      legendstr = c('Observed Statistic',
                    paste('Permutation',1-q,'critical value'))
      
      legend("bottomright",legend=legendstr,col=c(2,4),
             lty=c(1,2),lwd=c(2,2))
    }
  }
  
  return(list(pval=pval,qval=qval,Fobs=Fobs,Fnull=Fnull,
              Fvals=Fvals,Fnullvals=Fnullvals,pvals.pts=pvals.pts,qvals.pts=qvals.pts,
              fRegressList=fRegressList,argvals=argvals))
}

