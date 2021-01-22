# End of skate_pred3.R. loads first half of that script to save time processing data that doesn't change with each climate model.

library(mgcv)
load("C:/Users/brian.grieve/Documents/thornyskate/Data/sub2/predictionworkspace.RData")
load("C:/Users/brian.grieve/Documents/thornyskate/rmodels.RData")

yearnum = 146; 

filename = 'D:/Backup/2018_6_01/Thorny_skate/Conditions_2D/DataIntoR/skate_var_exp_all_1.csv';
climdat = read.csv(filename,header=FALSE); names(climdat) = names(data)
gam.short.pred = matrix(data=NA,nrow = nrow(climdat), ncol = yearnum);
gampa.short.pred = matrix(data=NA,nrow = nrow(climdat), ncol = yearnum);
gam.long.pred = matrix(data=NA,nrow = nrow(climdat), ncol = yearnum);
gampa.long.pred = matrix(data=NA,nrow = nrow(climdat), ncol = yearnum);
gam.ll.pred = matrix(data=NA,nrow = nrow(climdat), ncol = yearnum);

for (yr in 1:yearnum) {
  if (is.element(yr,seq(1,yearnum,10))){
    print(yr)
  }
  filename = sprintf('D:/Backup/2018_6_01/Thorny_skate/Conditions_2D/DataIntoR/skate_var_exp_all_%d.csv',yr);
  climdat = read.csv(filename,header=FALSE); names(climdat) = names(data[1:ncol(climdat)])
  gam.short.pred[,yr] = predict(gam.short,climdat,type="response");  
  gampa.short.pred[,yr] = predict(gampa.short,climdat,type="response");
  gam.long.pred[,yr] = predict(gam.long,climdat,type="response");
  gampa.long.pred[,yr] = predict(gampa.long,climdat,type="response");
  gam.ll.pred[,yr] = predict(gam.ll,climdat,type="response")
}


# Export
gam.short.pred[which(is.na(gam.short.pred))]=-999; 
gampa.short.pred[which(is.na(gampa.short.pred))]=-999; 
gam.long.pred[which(is.na(gam.long.pred))]=-999; 
gampa.long.pred[which(is.na(gampa.long.pred))]=-999; 
gam.ll.pred[which(is.na(gam.ll.pred))]=-999; 
write.csv(gam.short.pred,"D:/Backup/2018_6_01/Thorny_Skate/Conditions_2D/ExportFromR2/gamshort_pred.csv",row.names=FALSE) 
write.csv(gam.long.pred,"D:/Backup/2018_6_01/Thorny_Skate/Conditions_2D/ExportFromR2/gamlong_pred.csv",row.names=FALSE)  
write.csv(gampa.short.pred,"D:/Backup/2018_6_01/Thorny_Skate/Conditions_2D/ExportFromR2/gamshort_pa.csv",row.names=FALSE) 
write.csv(gampa.long.pred,"D:/Backup/2018_6_01/Thorny_Skate/Conditions_2D/ExportFromR2/gamlong_pa.csv",row.names=FALSE)
write.csv(gam.ll.pred,"D:/Backup/2018_6_01/Thorny_Skate/Conditions_2D/ExportFromR2/gamll_pred.csv",row.names=FALSE) 


