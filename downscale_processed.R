downscale_processed<-function(processed,target_Inh,target_Conc,Fitted_Model,Fitted_Model2,mean_curve,imped.pca){
      require(fda)
      newdata<-data.frame(log(target_Conc), rep(target_Inh,length(target_Conc))) 
      colnames(newdata)<-c("Concentration","Inhibitor")
      predicted_values<-predict(Fitted_Model,newdata)
      predicted_values2<-predict(Fitted_Model2,newdata)
      
      Mean_Profiles<- processed%>%filter(Inhibitor==target_Inh)%>%spread(key=Time,value=Impedance)%>%select(-treatGF)%>%
            group_by(Inhibitor, Concentration) %>%
            summarize_at(vars(-Measure),funs(mean))%>%ungroup()%>%gather(-c(Concentration,Inhibitor),key="Time",value="Impedance")
      
      Mean_Profiles$Time<-as.numeric(Mean_Profiles$Time)
      Mean_Profiles$Time<-(Mean_Profiles$Time-Mean_Profiles$Time[1])*6
      mytimes<-Mean_Profiles$Time%>%unique()
      Mean_Profiles$Concentration<-plyr::revalue(Mean_Profiles$Concentration, c(C0="0",C006="0.06",C025="0.25",C1="1",C4="4"))%>%as.character()%>%as.numeric()
      
      for (i in 1:length(target_Conc)){
            estimated_curve<-imped.pca$harmonics[1]*predicted_values[i]+imped.pca$harmonics[2]*predicted_values2[i]+mean_curve
            Mean_Profiles<-Mean_Profiles%>%
                add_row(Inhibitor=rep(target_Inh,length(mytimes)), Concentration=rep(target_Conc[i],length(mytimes)), Time=mytimes,Impedance=eval.fd(mytimes,estimated_curve)%>%as.vector())
      }
      return(Mean_Profiles)
}