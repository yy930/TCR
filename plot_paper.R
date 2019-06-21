library(ggplot2)
library(dplyr)
#Fig2.
#summary_over_settings_new_data <- read.csv("~/Downloads/all_result_new/summary_over_settings_new_data.csv", sep=";")
#df_evaluation<-summary_over_settings_new_data%>%subset(data=="noFilter")%>%subset(!ref%in%c("pubseqJCI","12kpatent"))
#df_evaluation$chain<-df_evaluation$chain%>%gsub("merged","TRA+TRB",.)%>%gsub("beta","TRB",.)
df_evaluation <- read.csv("~/Documents/TCRdata/results_without_CD1426/summary_over_settings_new_data.csv", sep=";")#!
df_evaluation$chain<-df_evaluation$chain%>%gsub("merged","TRA+TRB",.)%>%gsub("beta","TRB",.)#!

df_evaluation$ref<-df_evaluation$ref%>%gsub("allClonotypes","All gluten-specific TCR",.)%>%gsub("pubseq","Public TCR ",.)
df_evaluation$match_quant<-df_evaluation$match_quant%>%gsub("clonalPercentage","with clonecount",.)%>%gsub("percentage","without clonecount",.)
#plot df_evaluation
pdf("/Users/yingy_adm/Documents/TCRdata/Fig2.pdf",height=10,width=14)
ggplot(data=df_evaluation, aes(x=ref, y=accuracy, fill=chain)) + 
  geom_bar(position = 'dodge', stat='identity') +
  geom_text(aes(label=round(accuracy,2)), position=position_dodge(width=0.9), vjust=-0.25,show.legend = FALSE) +
  facet_wrap(~match_quant) +
  labs(x="", y="balanced accuracy",fill = "repertoire type")+
  theme_bw()
ggplot(data=df_evaluation, aes(x=ref, y=accuracy, fill=chain)) + 
  geom_bar(position = 'dodge', stat='identity') +
  geom_text(aes(label=round(accuracy,2)), position=position_dodge(width=0.9), vjust=-0.25,show.legend = FALSE) +
  facet_wrap(~match_quant) +
  labs(x="", y="balanced accuracy",fill = "repertoire type")+
  theme_bw()
dev.off()  

#Fig2.alt: plot df_evaluation only beta
df_evaluation<-summary_over_settings_new_data%>%subset(data=="noFilter")%>%subset(!ref%in%c("pubseqJCI","12kpatent"))%>%subset(chain=="beta")
colnames(df_evaluation)[1:3]<-c("data","repertoire type","reference")
ggplot(data=df_evaluation, aes(x=reference, y=accuracy)) + 
  geom_bar(position = 'dodge', stat='identity') +
  geom_text(aes(label=accuracy), position=position_dodge(width=0.9), vjust=-0.25,show.legend = FALSE) +
  facet_wrap(~match_quant) +
  theme_classic() 

# #Fig3.
#plot for fitted logistic regression with CV (result from immuneML)
predictions <- read.csv("~/Documents/TCRdata/results_without_CD1426/clonal_percentage/machine_learning_lr_mixcr_umi_no_filter_only_beta_pubseq_clonal_percentage/loocv/predictions.csv")
predictions$proportion_clones_matched<-list_matchNum[['matchNum_no_filter_only_beta_pubseq_count']]$proportion_clones_matched
predictions$repertoire_size_clonecounts<-list_matchNum[['matchNum_no_filter_only_beta_pubseq_count']]$repertoire_size_clonecounts
fit<-glm(formula = as.numeric(as.logical(CD_true_class)) ~ proportion_clones_matched, family=binomial, data = predictions)
predictions$pred_cd = predict(fit, newdata=predictions, type="response")
newdat <- data.frame(proportion_clones_matched=seq(min(predictions$proportion_clones_matched), max(predictions$proportion_clones_matched),len=100))
newdat$CD_true_class = predict(fit, newdata=newdat, type="response")
ggplot(predictions, aes(x=proportion_clones_matched, y=CD_True_proba,color=repertoire_size_clonecounts, shape = CD_true_class)) +
  geom_point(size = 3) +
  scale_color_gradient(low="blue", high="red") +
  geom_line(data=newdat,aes(x=proportion_clones_matched,y=CD_true_class),inherit.aes = FALSE)

#PLOT with CV
a<-list_matchNum[[3]]
a$group<-a$CD_status%>%gsub("True","UCD",.)%>%gsub("False","HC",.)
train_control <- trainControl(method="LOOCV",classProbs = TRUE,savePredictions=T,summaryFunction = twoClassSummary)# define training control
model_clones <- train(CD_status~proportion_clones_matched, data=a, method="glm",family=binomial(),trControl=train_control)
a$pred_cd_proportion_clones<-model_clones$pred$True
ggplot(a, aes(x=proportion_clones_matched, y=pred_cd_proportion_clones, color=repertoire_size_clonecounts, shape = group)) +
  geom_point(size = 2) + scale_color_gradient(low="blue", high="red") +
  labs(x="predicted probability", y="proportion of clones matching public TCR",color = "repertoire size")

#PLOT without CV
model<-glm(formula = as.logical(CD_status) ~ proportion_clones_matched, family=binomial, data = a)
a$pred_cd_proportion_clones<-predict(model, newdata=a, type="response")
newdat <- data.frame(proportion_clones_matched=seq(min(a$proportion_clones_matched), max(a$proportion_clones_matched),len=100))
newdat$CD_status = predict(model, newdata=newdat, type="response")
ggplot(a, aes(x=proportion_clones_matched, y=pred_cd_proportion_clones, color=repertoire_size_clonecounts, shape = group)) +
  geom_point(size = 3) + scale_color_gradient(low="blue", high="red") +
  labs(x="predicted probability", y="proportion of clones matching public TCR",color = "repertoire size")+
  geom_line(data=newdat,aes(x=proportion_clones_matched,y=CD_status),inherit.aes = FALSE)

#Fig4 count_of_sequences_matched in TRA,TRB,TRA+TRB
rm(a)
a<-list_matchNum[['matchNum_no_filter_merged_pubseq_count']][,c(1,2,3,4,6,9:12)]
a$chain<-"alpha"
a$count_of_sequences_matched=a$count_of_sequences_matched - list_matchNum[['matchNum_no_filter_only_beta_pubseq_count']]$count_of_sequences_matched
a$repertoire_size=a$repertoire_size-list_matchNum[['matchNum_no_filter_only_beta_pubseq_count']]$repertoire_size
a$count_of_clones_matched = a$count_of_clones_matched  - list_matchNum[['matchNum_no_filter_only_beta_pubseq_count']]$count_of_clones_matched
a$repertoire_size_clonecounts=a$repertoire_size_clonecounts-list_matchNum[['matchNum_no_filter_only_beta_pubseq_count']]$repertoire_size_clonecounts
a$proportion_sequences_matched = a$count_of_sequences_matched/a$repertoire_size
a$proportion_clones_matched = a$count_of_clones_matched/a$repertoire_size_clonecounts

#a<-rbind(a,list_matchNum[['matchNum_no_filter_only_beta_pubseq_count']][, c(1,2,3,4,6,9:12)], list_matchNum[['matchNum_no_filter_merged_pubseq_count']][, c(1,2,3,4,6,9:12)])
a<-rbind(a[, c(1:5,9:11)],list_matchNum[['matchNum_no_filter_only_beta_pubseq_count']][, c(1:4,6,12:14)], list_matchNum[['matchNum_no_filter_merged_pubseq_count']][, c(1:4,6,12:14)])
a$chain<-a$chain%>%gsub("merged","TRA+TRB",.)%>%gsub("beta","TRB",.)%>%gsub("alpha","TRA",.)
#a$CD_status<-a$CD_status%>%gsub("True","UCD",.)%>%gsub("False","HC",.)
a$CD<-a$CD%>%gsub("True","UCD",.)%>%gsub("False","HC",.)

p<-ggplot(transform(a,chain=factor(chain,levels=c("TRA","TRB","TRA+TRB"))),
          aes(x=CD,y=proportion_sequences_matched))+#CD_status
  geom_jitter(aes(shape=CD,colour=repertoire_size),width = 0.25)+#!CD_status
  facet_wrap(~chain) +
  labs(x="", y="proportion of public sequences",color = "repertoire size") +
  theme_bw() 
p+guides(shape = FALSE)

p<-ggplot(transform(a,chain=factor(chain,levels=c("TRA","TRB","TRA+TRB"))),
          aes(x=CD,y=proportion_clones_matched))+#!CD_status
  geom_jitter(aes(shape=CD,colour=repertoire_size_clonecounts),width = 0.25)+#!CD_status
  facet_wrap(~chain) +
  labs(x="", y="proportion of public sequences",color = "repertoire size") +
  theme_bw() 
p+guides(shape = FALSE)

#Likelihood ratio test
library(lmtest)
model<-glm(formula = as.logical(CD_status) ~ proportion_clones_matched, family=binomial, data = subset(a,chain=="TRB"))
model0<-glm(formula = as.logical(CD_status) ~ 1, family=binomial, data = subset(a,chain=="TRB"))
lrtest(model, model0)

#Fig5 count_of_sequences_matched in TRB with label
library(ggrepel)
p<-ggplot(a[a$chain=="TRB",],aes(x=CD_status,y=proportion_clones_matched))+geom_point(aes(colour=CD_status,shape=CD_status))+
  geom_label_repel(aes(label = patient),box.padding = 0.5, point.padding = 0.5)+
  labs(x="", y="proportion of public sequences",color = "repertoire size") +
  theme_classic()
p+guides(shape = FALSE, color=FALSE)