### from similarity matrices write a function that generates a full matrix

### sim combine

library(reshape)

sim_combine <- function(c_inputs, p_inputs, p_cutoff=0.05){
  
  n1=gsub("_vs_.*", "", basename(c_inputs))
  n2=gsub("\\.txt", "", gsub(".*_vs_", "", basename(c_inputs)))
  name_pairs=data.frame(n1=n1,n2=n2)
  orderer=data.frame(table(name_pairs$n1))
  orderer=orderer[order(orderer$Freq, decreasing=T),]
  
  c_df = data.frame(row_name=NA, col_name=NA, value=NA)
  p_df = data.frame(row_name=NA, col_name=NA, value=NA)
  type_df=data.frame(feat=NA, dataType=NA)
  for(i in 1:length(c_inputs)){
    c_mats1 = read.table(c_inputs[i], T, "\t", check.names = F)
    rownames(c_mats1) <- c_mats1[,1]
    c_mats1=c_mats1[,-1]
    p_mats1 = read.table(p_inputs[i], T, "\t")
    rownames(p_mats1) <- p_mats1[,1]
    p_mats1=p_mats1[,-1]
    
    doublecheck_r=intersect(rownames(p_mats1), rownames(c_mats1))
    doublecheck_c=intersect(colnames(p_mats1), colnames(c_mats1))
    
    c_mats1=c_mats1[doublecheck_r,doublecheck_c]
    p_mats1=p_mats1[doublecheck_r,doublecheck_c]
    
    colnames(c_mats1) <- gsub("CD45\\.", "CD45\\+", colnames(c_mats1))
    colnames(p_mats1) <- gsub("CD45\\.", "CD45\\+", colnames(p_mats1))
    colnames(c_mats1) <- gsub("CD25\\.", "CD25\\+", colnames(c_mats1))
    colnames(p_mats1) <- gsub("CD25\\.", "CD25\\+", colnames(p_mats1))
    
    if(grepl("X", colnames(c_mats1)[1])){
      colnames(c_mats1) <- gsub("X", "", colnames(c_mats1))
      colnames(p_mats1) <- gsub("X", "", colnames(p_mats1))
      
    }
    
    
    if(name_pairs$n1[i]==name_pairs$n2[i]){
      if(grepl("D_0", colnames(c_mats1)[1])){
        rownames(c_mats1) <- colnames(c_mats1)
        colnames(p_mats1) <- colnames(p_mats1)
      }
      #colnames(c_mats1) <- rownames(c_mats1)
      #colnames(p_mats1) <- rownames(p_mats1)
      type=data.frame(feat=rownames(c_mats1), dataType=rep(name_pairs$n1[i], nrow(c_mats1)))
      type_df=rbind(type_df, type)
    }
    
    row_name=rep(rownames(c_mats1), ncol(c_mats1))
    col_name=rep(colnames(c_mats1),each=nrow(c_mats1))
    
    c_nums = c(unlist(c_mats1))
    p_nums = c(unlist(p_mats1))
    c_df = rbind(c_df,data.frame(row_name=row_name, col_name=col_name, value=c_nums))
    p_df = rbind(p_df,data.frame(row_name=row_name, col_name=col_name, value=p_nums))
    
    if(name_pairs$n1[i]!=name_pairs$n2[i]){
      c_nums=c(unlist(t(c_mats1)))
      c_df = rbind(c_df,data.frame(row_name=col_name, col_name=row_name, value=c_nums))
      p_nums=c(unlist(t(p_mats1)))
      p_df = rbind(p_df,data.frame(row_name=col_name, col_name=row_name, value=p_nums))
      
    }
    
  }
  c_df=c_df[-1,]
  p_df=p_df[-1,]
  type_df <- type_df[-1,]
  rownames(type_df) <- type_df$feat
  
  c_df$p.value <- p_df$value
  c_df$p.value[is.na(c_df$p.value)]<-1
  
  c_df_sub <- subset(c_df, p.value < p_cutoff)
  c_df_sub$row_name <- factor(c_df_sub$row_name)
  c_df_sub$col_name <- factor(c_df_sub$col_name)
  
  nameVals <- sort(unique(unlist(c_df_sub[1:2])))
  myMat1 <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
  myMat1[as.matrix(c_df_sub[c("col_name", "row_name")])] <- c_df_sub[["value"]]
  myMat2 <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
  myMat2[as.matrix(c_df_sub[c("row_name", "col_name")])] <- c_df_sub[["value"]]
  
  myMat12 <- data.frame(myMat1+myMat2)
  
  fixers_up <- which(myMat12 > 1, arr.ind = T)
  fixers_down <- which(myMat12 < -1, arr.ind = T)
  
  fixers <- rbind(fixers_up, fixers_down)
  myMat_fix <- myMat12
  for( i in 1:nrow(fixers)){
    coord=fixers[i,]
    myMat_fix[coord[1],coord[2]] <- c(myMat_fix[coord[1],coord[2]]/2)
  }
  
  colnames(myMat_fix) <- rownames(myMat_fix)
  
  # c_full <- cast(row_name ~ col_name, data=c_df)
  # c_full <- cast(c_df, row_name ~ col_name, mean)
  # p_full <- cast(row_name ~ col_name, data=p_df)
  # 
  # huh=data.frame(colnames(c_full)[-1], c_full[,1])
  # x1=unlist(c(c_full[200,-1]))
  # x2=unlist(c(c_full[,201]))
  # huh2=data.frame(x1,x2)
  # 
  # c_flat = c(unlist(data.frame(c_full[,-1])))
  # p_flat = c(unlist(data.frame(p_full[,-1])))
  # c_flat[p_flat>p_cutoff]<-0
  # c_mats2 = matrix(c_flat, nrow(c_full), c(ncol(c_full)-1))
  # rownames(c_mats2) <- c_full[,1]
  # colnames(c_mats2) <- colnames(c_full)[-1]
  # 
  # c_mats2[is.na(c_mats2)] <- 0
  # idx=diag(ncol(c_mats2))
  # c_mats2[idx==1]<-0
  # c_mats2<-c_mats2[as.character(type_df$feat), as.character(type_df$feat)]
  # 
  # good_cols <- which(rowSums(abs(c_mats2))>0)
  
  print(paste0(nrow(myMat_fix), " non-empty features found & kept"))
  
  return(list(myMat_fix, type_df[rownames(myMat_fix),]))
}



