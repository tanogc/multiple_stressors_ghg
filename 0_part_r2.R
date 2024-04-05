# Function to calculate variance partition for models with different predictors

hier.r2.res<-function(mod.set, sel.pred) {
  
  hier.r2.mod<-data.frame(matrix(0, length(mod.set), length(sel.pred)))
  names(hier.r2.mod)<-sel.pred
  rownames(hier.r2.mod)<-names(mod.set)
  
  for (j in 1:length(mod.set)){
    calcVarPart(mod.set[[j]])->part_var
    
    part_var[-length(part_var)]->part_var # removing residuals
    part_var[order(names(part_var))]->hier.r2.mod[j, which(names(hier.r2.mod) %in% names(part_var))]
  }
  
  return(hier.r2.mod)
} 