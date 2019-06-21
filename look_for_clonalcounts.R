#load main functions 
library(stringr)
source('~/Documents/TCRdata/script/main_functions.R') 
options(stringsAsFactors = FALSE)

#read in all "matching_sequence_overview.tsv" in all settings as data frame "matchNum_*"
#settings<-list.dirs("/Users/yingy_adm/Documents/TCRdata/all_results/percentage",recursive = F)[7:12]
settings<-list.dirs("/Users/yingy_adm/Documents/TCRdata/all_result_new/analysis_runs/gut_data/fixed_results/count",recursive = F)[3:6]
for (i in settings){
  filename = paste("matchNum",paste(str_split(basename(i),"_")[[1]][7:length(str_split(basename(i),"_")[[1]])],collapse = '_'),sep = "_")
  x = read.csv(paste(i,"/matching_sequence_overview.tsv",sep = ""),sep = '\t')
  x$chain = ifelse(grepl("merged",filename),"merged","beta")
  x$data = ifelse(grepl("no_filter",filename),"umi_no_filter","umi")
  x$ref = str_split(basename(i),"_")[[1]][length(str_split(basename(i),"_")[[1]])-1]
  x <- subset(x,donor!="CD1426")
  assign(filename,x)
}

#apply function "create_df_clonecount" to produce multiple dataframes names "clonecount_(setting)"
# folder_report_matching_sequence_settings <-"/Users/yingy_adm/Documents/TCRdata/all_results/percentage"
# report_matching_seq_settings <- list.dirs(folder_report_matching_sequence_settings,recursive = F)[7:12]
# for (setting in report_matching_seq_settings){
#   create_df_clonecount(setting)
# }
for (setting in settings){
  create_df_clonecount(setting)
}

#get numbers for Fig1
nrow(unique(clonecount_no_filter_only_beta_pubseq_count[,1:3]))#39
nrow(unique(clonecount_no_filter_merged_pubseq_count[,1:3]))#54
nrow(unique(clonecount_no_filter_merged_scTCR_IRIS_allClonotypes_count[,1:3]))#151
nrow(unique(clonecount_no_filter_only_beta_scTCR_IRIS_allClonotypes_count[,1:3]))#93

#Add columns ("count_of_sequences_matched","count_of_clones_matched") to data frames "matchNum*" with different settings.
#NB: "clonecount_*" should be produced and the suffix(setting) should match with "matchNum_*"
dfnames_match<-ls(pattern = "matchNum_no_filter_")
for (dfname in dfnames_match){
  dfname_clonecount<-gsub("matchNum","clonecount",dfname)
  df_clonecount<-get(dfname_clonecount)
  df_match<-get(dfname)
  df_match$count_of_sequences_matched = NA
  df_match$count_of_clones_matched = NA
  for (i in 1:nrow(df_match)){
    patient = df_match$repertoire_identifier[i]#!!!!!!!!!df_match$patient[i]
    chain = df_match$chain[i]
    df_match$count_of_sequences_matched[i] = ifelse(chain=="merged",sum(df_clonecount$donor==patient),
                                                 sum(df_clonecount$donor==patient & df_clonecount$chain=="B"))
    df_match$count_of_clones_matched[i] = ifelse(chain=="merged",sum(df_clonecount$cloneCount[df_clonecount$donor==patient]),
                                              sum(df_clonecount$cloneCount[df_clonecount$donor==patient & df_clonecount$chain=="B"]))
  }
  assign(dfname,df_match)
}

library(readxl)
library(data.table)
#read in repertoire_size_clonecounts info for all repertoires (df:donor+chain+repertoire_size_clonecounts)
repertoire_size_clonecounts <- read_excel("~/Documents/TCRdata/repertoire_size_clonecounts.xlsx")#alpha & beta
repertoire_size_clonecounts_merged<-repertoire_size_clonecounts %>% setDT %>% 
  .[,.(chain="merged",repertoire_size_clonecounts=sum(repertoire_size_clonecounts)),by=patient]#meged
repertoire_size_clonecounts<<-repertoire_size_clonecounts  %>% .[chain=="beta",,] %>% rbind(repertoire_size_clonecounts_merged)#beta & merged
repertoire_size_clonecounts_merged<-subset(repertoire_size_clonecounts_merged,patient!="CD1426")
repertoire_size_clonecounts<-subset(repertoire_size_clonecounts,patient!="CD1426")
colnames(repertoire_size_clonecounts)[1]=colnames(repertoire_size_clonecounts_merged)[1]="repertoire_identifier"#!!!!!!

#creat list_matchNum
list_matchNum <- lapply(as.list(ls(pattern = "matchNum_no_filter")),FUN = get)
names(list_matchNum)<-ls(pattern = "matchNum_no_filter")

#merge_matchNum_with_repertoire_size_clonecounts for all elements in list_matchNum
list_matchNum<-lapply(list_matchNum,merge_matchNum_with_repertoire_size_clonecounts)

#add 2 cloumns "proportion_" for all elements in list_matchNum
list_matchNum<-lapply(list_matchNum,get_proportion)

#add 2 cloumns "pred_cd_" for all elements in list_matchNum
list_matchNum<-lapply(list_matchNum, rename_CD)
#plot_paper.R Fig3 can be done here

#########################################################TCR papar result stops here############
list_matchNum<-lapply(list_matchNum, loocv_fit)

#calculate "accuracy","precision","sensitivity","specificity" for all elements in list_matchNum
list_evaluation_seq<-lapply(list_matchNum, evaluate_pred_seq)
list_evaluation_seq<-t(matrix(data = unlist(list_evaluation_seq),ncol=6))
list_evaluation_clone<-lapply(list_matchNum, evaluate_pred_clone)
list_evaluation_clone<-t(matrix(data = unlist(list_evaluation_clone),ncol=6))

#combine above as dataframe "df_evaluation"
df_evaluation<-as.data.frame(rbind(list_evaluation_seq,list_evaluation_clone),
                             row.names = c(paste(names(list_matchNum),"seq",sep = "_"),paste(names(list_matchNum),"clone",sep = "_")))
rownames(df_evaluation)<-str_replace_all(rownames(df_evaluation),c("matchNum_"="","no_filter"="NoFilter",
                                                                   "only_"="","scTCR_IRIS_"="","percentage_"=""))
df_evaluation <- cbind(df_evaluation,str_split_fixed(rownames(df_evaluation),"_",4))
colnames(df_evaluation)<-c("accuracy","precision","sensitivity","specificity","data","chain","reference","match_quant")

#plot df_evaluation
ggplot(df_evaluation,aes(reference,accuracy))+ 
  geom_bar(aes(fill = chain), stat = "identity", color = "white", position = position_dodge(0.9) )+ 
  facet_wrap(~match_quant)

#sort all "matchNum_no_filter" object by column "patient" 
for (obj_name in ls(pattern = "matchNum_no_filter_only_beta")){
  obj=get(obj_name)
  assign(obj_name,obj[order(obj$patient),])  
}

#sum "matchNum" and "matchCLONE" for A, B chain as dataframe "AB_merged*" 
for (merged_name in ls(pattern = "matchNum_no_filter_merged")){
  merged = get(merged_name)
  beta_name = gsub("merged","only_beta",merged_name)
  beta = get(beta_name)
  MATCH_NUM_AB = gsub("matchNum_no_filter","AB",merged_name)
  assign(MATCH_NUM_AB,data.frame(patient=merged$patient,
                       CD_status=merged$CD_status,
                       matchNum_A=merged$count_of_sequences_matched-
                         beta$count_of_sequences_matched,
                       matchCLONE_A=merged$count_of_clones_matched-
                         beta$count_of_clones_matched,
                       matchNum_B=beta$count_of_sequences_matched,
                       matchCLONE_B=beta$count_of_clones_matched))
}

# matched alpha chains in control
matched_alpha_in_control<-data.frame()
for (obj_name in ls(pattern = "AB_merged_")){
  obj=get(obj_name)#AB_merged_pubseq_percentage
  setting=gsub("AB_merged_","",obj_name)
  donors = obj$patient[which(obj$CD_status=="False" & obj$matchNum_A>0)]
  count_name=gsub("AB","clonecount_no_filter",obj_name)#"clonecount_no_filter_merged_pubseq_percentage"
  count_df=get(count_name)
  count_df<-count_df[count_df$donor%in%donors,c("donor","sequence","v_gene","j_gene","cloneCount")]
  count_df$setting=setting
  matched_alpha_in_control<-rbind(matched_alpha_in_control,count_df)
}
#look for paired seq for the above matched alpha chains in control
pubseq <- read.csv("~/Documents/TCRdata/pubseq.csv", sep=";",na.strings=c("","NA"))
extract_pubseq <- read.csv2("~/Documents/TCRdata/extract_pubseq.csv")
#if a chain is matched, is the other chain from sc_pub_seq is also matched
matched_alpha_in_control$other_seq<-NA
matched_alpha_in_control$other_seq_patient<-NA
matched_alpha_in_control$other_seq_patient_group<-NA

list_matched_alpha_in_control <- split(matched_alpha_in_control,matched_alpha_in_control$setting)

get_donor_group<-function(vec_donors){
  vec_donors_group<-rep(NA,length(vec_donors))
  for (j in 1:length(vec_donors)){
    vec_donors_group[j]<-c("DQ2_CD","DQ2_HC","DQ8_HC")[c(is.element(vec_donors[j],DQ2_CD),is.element(vec_donors[j],DQ2_HC),is.element(vec_donors[j],DQ8_HC))]
    }
  return(vec_donors_group)
}

get_paired_beta<-function(matched_alpha){
  for (i in 1:dim(matched_alpha)[1]){
  setting<-matched_alpha$setting[i]
  clonecount_df<-get(paste("clonecount_no_filter_only_beta",setting,sep = "_"),envir = .GlobalEnv)#!!!get df from outside
  print(dim(clonecount_df))
  seq=matched_alpha$sequence[i]
  id_pub<-which(pubseq$Chain..TRA..1.==seq |pubseq$Chain..TRA..2.==seq)
  other_seq<-gsub("NA,?","",paste(unique(pubseq$Chain..TRB..1.[id_pub]),collapse = ","))
  matched_alpha$other_seq[i]<-other_seq
  if (length(str_split(other_seq,","))>0){#get patient_id from the corresponding clonecount_df
    id_clonecount_df<-clonecount_df$sequence %in% unlist(str_split(other_seq,","))
    id_clonecount_df<-id_clonecount_df[!is.na(id_clonecount_df)]
    if(sum(id_clonecount_df)>0){
      matched_alpha$other_seq_patient[i]<-paste(clonecount_df$donor[id_clonecount_df],collapse = ",")
      other_seq_donor_group<-get_donor_group(clonecount_df$donor[id_clonecount_df])#"sameHC","diffHC","CD"
      matched_alpha$other_seq_patient_group[i]<-paste(other_seq_donor_group,collapse = ",")
      }
    }
  }
#matched_alpha$paired_match<-str_split_fixed(matched_alpha$other_seq,",")
  return(matched_alpha)
}
lapply(list_matched_alpha_in_control, get_paired_beta)

#look at age
library(readxl)
sample_info_TCR <- read_excel("~/Documents/TCRdata/sample_info_TCR.xlsx", sheet = "Sheet1")
only_beta_pubseq<-list_matchNum[[5]]
only_beta_pubseq$age<-sample_info_TCR$Age
#only_beta_pubseq<-only_beta_pubseq[only_beta_pubseq$CD_status=="True",c('patient','age','proportion_clones_matched')]
cor.test(only_beta_pubseq$proportion_clones_matched[only_beta_pubseq$CD_status=="True"],
         only_beta_pubseq$age[only_beta_pubseq$CD_status=="True"])#p-value = 0.09906

#plot ROC curve for LOOCV
beta<-list_matchNum[[5]]
logit_fit<-function(df){glm(formula = as.logical(CD_status) ~ proportion_clones_matched, family=binomial, data = df)}
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
loocv_result<-matrix(NA,nrow = 16,ncol = 4)
for (i in 1:16){
  loocv_result[i,]<-foo(beta)
}
loocv_result<-as.data.frame(loocv_result)
colnames(loocv_result)<-c("donor","pred","target","repertoire_size_clonecounts")

library(ROCR)
pred<-prediction(as.numeric(loocv_result$pred),loocv_result$target)
roc.perf = performance(pred,measure = "tpr",x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1)
acc.perf = performance(pred, measure = "acc")
plot(acc.perf)

#plot ROC curve for bootstrape LOOCV
pred<-prediction(as.numeric(bootstrap_result$pred),bootstrap_result$target)
predictions <- read.csv("~/Downloads/beta_pubseq_regression/proportion_clones/predictions.csv", sep=";")
pred<-prediction(predictions$CD_True_proba,predictions$CD_true_class)
roc.perf = performance(pred,measure = "tpr",x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1)


