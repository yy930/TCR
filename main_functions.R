library(stringr)

#look for clonecounts for matched seq (ONLY merged folder)!!
DQ2_HC <- c("CD1386", "CD1408", "CD1409", "CD1431","CD1450")#dropped:"CD1331","CD1440"
DQ8_HC <- c("CD1370", "CD1390", "CD1453")#to drop:CD1443,CD1465
DQ2_CD <- c("CD1357", "CD1358", "CD1364", "CD1368", "CD1393", "CD1422", "CD1426", "CD1451")#to drop:CD1422,CD1451
CD_list<-list(DQ2_HC,DQ8_HC,DQ2_CD)
names(CD_list) <- c("DQ2_HC", "DQ8_HC", "DQ2_CD")

#function to create "clonecount_(setting)" data frame for all the matched sequences 
create_df_clonecount <-function(folder_report_matching_sequence_onesetting){
  files<-list.files(folder_report_matching_sequence_onesetting,pattern="CD*")[1:16]#matching seq files
  df_clonecount<-data.frame()
  mixcr_folder<-"/Users/yingy_adm/Documents/TCRdata/mixcr_umi_no_filter"
  for (icd in files){
    donor = substr(basename(icd),1,6)
    x = read.csv(paste(folder_report_matching_sequence_onesetting,"/",icd,sep = ""),sep = '\t')
    x$donor<-donor
    x = x[!is.na(x$matching_sequences),]
    x = x[x$matching_sequences!="",]
    x$matching_sequences = unlist(lapply(x$matching_sequences,FUN = function(x) str_split_fixed(x,",",2)[1]))
    #x$flag = paste("",x$matching_sequences,"")
    if (nrow(x)>0){
      cdgroup<-names(CD_list)[unlist(lapply(CD_list, function(x) donor%in%x))]
      x$cloneCount<-NA
      xA<-subset(x,chain=="A")
      xB<-subset(x,chain=="B")
      if(nrow(xB)>0){
        mixcr_fileB_path<-paste(mixcr_folder,"/",cdgroup,"/",donor,"/TRB/",donor,"_clones_TRB.csv",sep = "")
        mixcr_fileB <- read.csv(mixcr_fileB_path,sep = '\t')
        for (i in 1:nrow(xB)){
          cdr3<-xB$matching_sequences[i]
          v_gene<-xB$v_gene[i]
          j_gene<-xB$j_gene[i]
          logical_inds_match <- (grepl(cdr3,xB$sequence) & grepl(v_gene,xB$v_gene) &
                                   grepl(j_gene,xB$j_gene))
          logical_inds_mixcr <- (grepl(cdr3,mixcr_fileB$aaSeqCDR3) & grepl(v_gene,mixcr_fileB$allVHitsWithScore) &
                                   grepl(j_gene,mixcr_fileB$allJHitsWithScore))
          if(sum(logical_inds_match)==sum(logical_inds_mixcr)){
            xB$cloneCount[logical_inds_match]<-mixcr_fileB$cloneCount[logical_inds_mixcr]
          } else {print(c(icd,i))}}
        df_clonecount<-rbind(df_clonecount,xB)
        #print(c(donor, "beta", as.character(sum(mixcr_fileB$cloneCount))))
      }
      if(nrow(xA)>0){
        mixcr_fileA_path<-paste(mixcr_folder,"/",cdgroup,"/",donor,"/TRA/",donor,"_clones_TRA.csv",sep = "")
        mixcr_fileA <- read.csv(mixcr_fileA_path,sep = '\t')
        for (j in 1:nrow(xA)){
          cdr3<-xA$matching_sequences[j]
          v_gene<-xA$v_gene[j]
          j_gene<-xA$j_gene[j]
          logical_inds_match <- (grepl(cdr3,xA$sequence) & grepl(v_gene,xA$v_gene) &
                                   grepl(j_gene,xA$j_gene))
          logical_inds_mixcr <- (grepl(cdr3,mixcr_fileA$aaSeqCDR3) & grepl(v_gene,mixcr_fileA$allVHitsWithScore) &
                                   grepl(j_gene,mixcr_fileA$allJHitsWithScore))
          if(sum(logical_inds_match)==sum(logical_inds_mixcr)){
            xA$cloneCount[logical_inds_match]<-mixcr_fileA$cloneCount[logical_inds_mixcr]
          } else {
            xA$cloneCount[j]<-mixcr_fileA$cloneCount[grepl(cdr3,mixcr_fileA$aaSeqCDR3)]}}
        df_clonecount<-rbind(df_clonecount,xA)
        #print(c(donor, "alpha", as.character(sum(mixcr_fileA$cloneCount))))
      }
    }
  }
  filename_clonecount <- paste("clonecount",paste(str_split(basename(folder_report_matching_sequence_onesetting),"_")[[1]]
                                                  [7:length(str_split(basename(folder_report_matching_sequence_onesetting),"_")[[1]])],
                                                  collapse = '_'),sep = "_")
  df_clonecount<-subset(df_clonecount,donor!="CD1426")
  assign(filename_clonecount,df_clonecount,envir = parent.frame())#use "envir" for "assign" to access the environment outside the function!
}

#function to merge matchNum with repertoire_size_clonecounts
# merge_matchNum_with_repertoire_size_clonecounts<-function(matchNum){
#   return (merge(matchNum,repertoire_size_clonecounts,by.x=c("patient","chain"),by.y=c("patient","chain")))
# }
merge_matchNum_with_repertoire_size_clonecounts<-function(matchNum){
  return (merge(matchNum,repertoire_size_clonecounts,by.x=c("repertoire_identifier","chain"),by.y=c("repertoire_identifier","chain")))
}

#function to calculate proportions and add the columns to matchNum
get_proportion<-function(matchNum){
  matchNum$proportion_sequences_matched<-matchNum$count_of_sequences_matched/matchNum$repertoire_size
  matchNum$proportion_clones_matched<-matchNum$count_of_clones_matched/matchNum$repertoire_size_clonecounts
  return (matchNum)
}

rename_CD<-function(x) {
  colnames(x)[6]<-"CD_status"
  return(x)}

library(caret)
loocv_fit<-function(dataset){
  train_control <- trainControl(method="LOOCV",savePredictions=T)# define training control
  model_sequences <- train(CD_status~proportion_sequences_matched, data=dataset, method="glm",family=binomial(),trControl=train_control)
  dataset$pred_cd_proportion_sequences<-model_sequences$pred$pred
  model_clones <- train(CD_status~proportion_clones_matched, data=dataset, method="glm",family=binomial(),trControl=train_control)
  dataset$pred_cd_proportion_clones<-model_clones$pred$pred
  #return(dataset)
  return(list(dataset,model_sequences,model_clones))
}

#function to evalueate prediction using sequences
library(dplyr)
evaluate_pred_seq<-function(pred_df){
  pred_df <- data.frame(lapply(pred_df, as.character), stringsAsFactors=FALSE)
  pred_df <- mutate_all(pred_df, .funs=toupper)
  accuracy=sum(pred_df$pred_cd_proportion_sequences==pred_df$CD_status)/nrow(pred_df)
  sensitivity=sum(pred_df$pred_cd_proportion_sequences=="TRUE" & pred_df$CD_status=="TRUE")/sum(pred_df$CD_status=="TRUE")
  specificity=sum(pred_df$pred_cd_proportion_sequences=="FALSE" & pred_df$CD_status=="FALSE")/sum(pred_df$CD_status=="FALSE")
  precision=sum(pred_df$pred_cd_proportion_sequences=="TRUE" & pred_df$CD_status=="TRUE")/sum(pred_df$pred_cd_proportion_sequences=="TRUE")
  return(c(accuracy,precision,sensitivity,specificity))
}

#function to evalueate prediction using clonecounts
evaluate_pred_clone<-function(pred_df){
  pred_df <- data.frame(lapply(pred_df, as.character), stringsAsFactors=FALSE)
  pred_df <- mutate_all(pred_df, .funs=toupper)
  accuracy=sum(pred_df$pred_cd_proportion_clones==pred_df$CD_status)/nrow(pred_df)
  sensitivity=sum(pred_df$pred_cd_proportion_clones=="TRUE" & pred_df$CD_status=="TRUE")/sum(pred_df$CD_status=="TRUE")
  specificity=sum(pred_df$pred_cd_proportion_clones=="FALSE" & pred_df$CD_status=="FALSE")/sum(pred_df$CD_status=="FALSE")
  precision=sum(pred_df$pred_cd_proportion_clones=="TRUE" & pred_df$CD_status=="TRUE")/sum(pred_df$pred_cd_proportion_clones=="TRUE")
  return(c(accuracy,precision,sensitivity,specificity))
}


