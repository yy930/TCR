library(dplyr)
#folder="/Users/yingy_adm/Documents/TCRdata/all_results/clonal_percentage/report_matching_sequences_percentage_mixcr_umi_no_filter_only_beta_pubseq_clonal_percentage/"
folder_beta_pubseq_percentage="/Users/yingy_adm/Documents/TCRdata/all_results/percentage/report_matching_sequences_percentage_mixcr_umi_no_filter_only_beta_pubseq_percentage/"
counts_matched = read.csv(paste(folder_beta_pubseq_percentage,"matching_sequence_overview.tsv",sep = ""),sep = '\t')
counts_matched$patient<-as.character(counts_matched$patient)

counts_matched <- list_matchNum[['matchNum_no_filter_only_beta_pubseq_percentage']]
counts_matched$logical_CD_status<-as.numeric(as.logical(counts_matched$CD_status))

#logit_fit<-function(df){glm(formula = logical_CD_status ~ count_of_sequences_matched, family=binomial, data = df)}
logit_fit<-function(df){glm(formula = logical_CD_status ~ proportion_clones_matched, family=binomial, data = df)}
foo <- function(counts_matched){
  test_ind<-sample(1:nrow(counts_matched),1)
  train_indices<-sample(setdiff(1:nrow(counts_matched),test_ind),15,replace = T)
  dt<-counts_matched[train_indices,]
  model <- logit_fit(dt)
  pred <- predict(model, counts_matched[test_ind,], type="response") #type="prob"
  return(c(donor=as.character(counts_matched[test_ind,'patient']),
              pred=pred,target=as.integer(as.logical(counts_matched[test_ind,'CD_status'])),
           repertoire_size_clonecounts=counts_matched[test_ind,'repertoire_size_clonecounts']))
}

bootstrap_rep<-1000
bootstrap_result<-matrix(NA,nrow = bootstrap_rep,ncol = 4)
for (i in 1:bootstrap_rep){
  bootstrap_result[i,]<-foo(counts_matched)
}
bootstrap_result<-as.data.frame(bootstrap_result)
colnames(bootstrap_result)<-c("donor","pred","target","repertoire_size_clonecounts")

bootstrap_result[,c("pred","target","repertoire_size")] = apply(bootstrap_result[,c("pred","target","repertoire_size")], 2, function(x) as.numeric(as.character(x)));
bootstrap_result$pred_binomial<-round(bootstrap_result$pred)
bootstrap_result_donor_list<-split(bootstrap_result,bootstrap_result$donor)

calculate_accuracy<-function(pred_df){
  accuracy=sum(pred_df$pred_binomial==pred_df$target)/nrow(pred_df)
  #sensitivity=sum(pred_df$pred_binomial==1 & pred_df$target==1)/sum(pred_df$target==1)
  #specificity=sum(pred_df$pred_binomial==0 & pred_df$target==0)/sum(pred_df$target==0)
  #precision=sum(pred_df$pred_binomial==1 & pred_df$target==1)/sum(pred_df$pred_binomial==1)
  #return(c(accuracy,precision,sensitivity,specificity))
  return (accuracy)
}

accuracy_list_donor<-unlist(lapply(bootstrap_result_donor_list,calculate_accuracy))
counts_matched$accuracy<-accuracy_list_donor[match(counts_matched$patient, names(accuracy_list_donor))]

library(ggplot2)
library(ggrepel)
ggplot(counts_matched,aes(x=repertoire_size,y=accuracy,color=CD_status))+geom_point()+
  geom_label_repel(aes(label = patient),box.padding   = 0.35, point.padding = 0.5)
  #stat_summary(fun.data=mean_cl_normal) + 
  #geom_smooth(method='lm',formula=y~x)
 
#matching_seq_file_list<-list.files(folder)[1:16]
