setwd("C:/Users/nb2/Desktop/Politecnico/StatApp/Progetto")
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

# Grafici!
# Ignorate questo sotto per ora
# processed_plot<-processed%>%group_by(Time)%>%mutate(base=Impedance[treatGF==F],baseGF=Impedance[treatGF==T & Inhibitor=="no"])%>%ungroup()
library(wesanderson)
library(RColorBrewer)

processed_plot<-processed

p1 <- ggplot(processed_plot,aes(x = Time, y=Impedance, group=interaction(Inhibitor,Concentration,Measure,treatGF)))
p1 + geom_line(aes(color = Inhibitor),size=1)+theme_minimal()+scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))

(p2 <- p1 + geom_line(aes(color=Measure),size=1) + scale_color_brewer(palette = "Accent")+
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

p3<-ggplot(load.data%>%rownames_to_column(var="Time")%>%select(c(Time,PC1,PC2,PC3))%>%melt(id=c("Time"),value.name = "Loadings"), aes(x=Time,y=Loadings))
p3+geom_bar(stat="identity",width=1,color="black",fill="snow1")+facet_wrap(~variable,nrow=3)+theme_minimal()+theme(axis.text.x = element_blank())

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
modPCA<-prcomp(ResidualMatr,scale=F)

modPCA%>%.$x%>%as.data.frame()%>%
  summarize_all(funs(var)) %>%
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance))

autoplot(modPCA,data=processed_ext,
         colour="Inhibitor",label=T,label.label="Concentration", label.repel=T)+theme_minimal()

# Posso tenere solo la prima componente principale ~97%
shape_score<-processed_ext%>%mutate(Score=processed_ext%>%select(-(1:4))%>%
                as.matrix()%>% multiply_by_matrix(modPCA$rotation[,1]))%>%select((1:4),Score)

model<-aov(Score~Inhibitor+Concentration,data=shape_score%>%filter(Inhibitor!="no")%>%droplevels())
anova(model)
