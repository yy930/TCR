#rbind matching_sequence_overview.tsv for all settings
for (i in settings){
  #filename = paste(str_split(basename(i),"_")[[1]][3:5],collapse = '_')
  filename = paste("matchNum",paste(str_split(basename(i),"_")[[1]][7:length(str_split(basename(i),"_")[[1]])],collapse = '_'),sep = "_")
  #filename = gsub("iris","",filename)
  x = read.csv(paste(i,"/matching_sequence_overview.tsv",sep = ""),sep = '\t')
  # x$chain = str_split(basename(i),"_")[[1]][4]
  # x$data = str_split(basename(i),"_")[[1]][3]
  # x$ref = gsub("iris","",str_split(basename(i),"_")[[1]][5])
  x$chain = ifelse(grepl("merged",filename),"merged","beta")
  x$data = ifelse(grepl("no_filter",filename),"umi_no_filter","umi")
  x$ref = str_split(basename(i),"_")[[1]][length(str_split(basename(i),"_")[[1]])-1]
  assign(filename,x)
}

ExactMatch<-data.frame()
for (i in ls(pattern = "clone*")){
  ExactMatch = rbind(ExactMatch,get(i))
}
for (i in ls(pattern = "umi*")){
  ExactMatch = rbind(ExactMatch,get(i))
}
#ExactMatch$CD_status<-as.logical(ExactMatch$CD_status)
table(ExactMatch$chain,ExactMatch$data,ExactMatch$ref)#check dim

#function to sum clone count by matching v j cdr3
sumclonecount<-function(logical_inds_mixcr_fileB,mixcr_fileB){
  sub_mixB<-mixcr_fileB[logical_inds_mixcr_fileB,]
  sub_mixB$allVHitsWithScore<-str_split_fixed(sub_mixB$allVHitsWithScore,"[*]",2)[,1]
  sub_mixB$allJHitsWithScore<-str_split_fixed(sub_mixB$allJHitsWithScore,"[*]",2)[,1]
  if (length(unique(sub_mixB$allVHitsWithScore))==1 & length(unique(sub_mixB$allJHitsWithScore))==1)
    {return(sum(sub_mixB$cloneCount))}
  else{print ("V or J not matching")}
}

#beta_pubseq_percentage
folder="/Users/yingy_adm/Documents/TCRdata/all_results/percentage/report_matching_sequences_percentage_mixcr_umi_no_filter_only_beta_pubseq_percentage"
#merged_pubseq_percentage
folder="/Users/yingy_adm/Documents/TCRdata/all_results/percentage/report_matching_sequences_percentage_mixcr_umi_no_filter_merged_pubseq_percentage"
#beta_pubseq_clonal_percentage
folder="/Users/yingy_adm/Documents/TCRdata/all_results/clonal_percentage/report_matching_sequences_percentage_mixcr_umi_no_filter_only_beta_pubseq_clonal_percentage"

#print the number of match
file16<-list.files(folder_report_matching_sequence_onesetting,pattern="CD*")
for (icd in file16){
  x = read.csv(paste(folder,"/",icd,sep = ""),sep = '\t')
  x = x[!is.na(x$matching_sequences),]
  x = x[x$matching_sequences!="",]
  print (nrow(x))#print the number of match
}

#print repertore_size_clonecounts/num_clonetypes
mixcr_folder<-"/Users/yingy_adm/Documents/TCRdata/mixcr_umi_no_filter"
donors=c("CD1357","CD1358","CD1364","CD1368","CD1370","CD1386","CD1390","CD1393","CD1408","CD1409","CD1422","CD1426","CD1431","CD1450","CD1451","CD1453")
for (donor in donors){
  cdgroup<-names(CD_list)[unlist(lapply(CD_list, function(x) donor%in%x))]####
  mixcr_fileB_path<-paste(mixcr_folder,"/",cdgroup,"/",donor,"/TRB/",donor,"_clones_TRB.csv",sep = "")
  mixcr_fileB <- read.csv(mixcr_fileB_path,sep = '\t')###
  mixcr_fileA_path<-paste(mixcr_folder,"/",cdgroup,"/",donor,"/TRA/",donor,"_clones_TRA.csv",sep = "")
  mixcr_fileA <- read.csv(mixcr_fileA_path,sep = '\t')###
  #print(c(donor, "alpha", as.character(sum(mixcr_fileA$cloneCount))))####clonecounts
  #print(c(donor, "beta", as.character(sum(mixcr_fileB$cloneCount))))#####clonecounts
  print(c(donor, "alpha", as.character(nrow(mixcr_fileA))))####num_clonetypes
  print(c(donor, "beta", as.character(nrow(mixcr_fileB))))#####num_clonetypes
}
# manually save repertore_size_clonecounts as ("~/Documents/TCRdata/repertore_size_clonecounts.xlsx")

#plot for repertoire_size for A,B,merged
a<-list_matchNum[[2]][,c(1,2,4,6,10)]
a$chain<-"alpha"
a$repertoire_size=a$repertoire_size-list_matchNum[[5]][,4]
a$repertoire_size_clonecounts=a$repertoire_size_clonecounts-list_matchNum[[5]][,10]
a<-rbind(a,list_matchNum[[5]][,c(1,2,4,6,10)],list_matchNum[[2]][,c(1,2,4,6,10)])
ggplot(a,aes(x=CD_status,y=repertoire_size))+geom_jitter(aes(colour=CD_status,shape=CD_status),width = 0.25)+
  facet_wrap(~chain)

#plot for proportion_[clones/sequences]_matched for A,B,merged
a<-list_matchNum[[2]][,c(1,2,3,4,6,9:12)]
#a$chain<-"alpha"
a$count_of_sequences_matched=a$count_of_sequences_matched - list_matchNum[[5]]$count_of_sequences_matched
a$repertoire_size=a$repertoire_size-list_matchNum[[5]]$repertoire_size
a$count_of_clones_matched = a$count_of_clones_matched  - list_matchNum[[5]]$count_of_clones_matched
a$repertoire_size_clonecounts=a$repertoire_size_clonecounts-list_matchNum[[5]]$repertoire_size_clonecounts
a$proportion_sequences_matched = a$count_of_sequences_matched/a$repertoire_size
a$proportion_clones_matched = a$count_of_clones_matched/a$repertoire_size_clonecounts
ggplot(a,aes(x=CD_status,y=proportion_sequences_matched))+geom_jitter(aes(colour=CD_status,shape=CD_status),width = 0.25)+
  facet_wrap(~chain)
ggplot(a,aes(x=CD_status,y=proportion_clones_matched))+geom_jitter(aes(colour=CD_status,shape=CD_status),width = 0.25)+
  facet_wrap(~chain)

#plot for fitted logistic regression
DTB<-list_matchNum[['matchNum_no_filter_only_beta_pubseq_percentage']]
DTB$correct_pred<-as.factor(DTB$CD_status==DTB$pred_cd_proportion_sequences)
DTB$CD_status_numric<-as.numeric(as.logical(DTB$CD_status))
fit<-glm(formula = CD_status ~ proportion_sequences_matched, family=binomial, data = DTB)
DTB$pred_cd = predict(fit, newdata=DTB, type="response")
newdat <- data.frame(proportion_sequences_matched=seq(min(DTB$proportion_sequences_matched), max(DTB$proportion_sequences_matched),len=100))
newdat$CD_status = predict(fit, newdata=newdat, type="response")

ggplot(DTB, aes(x=proportion_sequences_matched, y=pred_cd, color=repertoire_size, shape = CD_status)) +
  geom_point(size = 4) +
  scale_color_gradient(low="blue", high="red") +
  geom_line(data=newdat,aes(x=proportion_sequences_matched,y=CD_status),inherit.aes = FALSE)

#Collapse pubseq
pubseq$TRA...V.gene..1.<-str_split_fixed(pubseq$TRA...V.gene..1.,"[*]",2)[,1]
pubseq$TRA...J.gene..1.<-str_split_fixed(pubseq$TRA...J.gene..1.,"[*]",2)[,1]
pubseq$TRB...V.gene..1.<-str_split_fixed(pubseq$TRB...V.gene..1.,"[*]",2)[,1]
pubseq$TRB...J.gene..1.<-str_split_fixed(pubseq$TRB...J.gene..1.,"[*]",2)[,1]
pubseq$TRA...V.gene..2.<-str_split_fixed(pubseq$TRA...V.gene..2.,"[*]",2)[,1]
pubseq$TRA...J.gene..2.<-str_split_fixed(pubseq$TRA...J.gene..2.,"[*]",2)[,1]
pubseq$TRB...V.gene..2.<-str_split_fixed(pubseq$TRB...V.gene..2.,"[*]",2)[,1]
pubseq$TRB...J.gene..2.<-str_split_fixed(pubseq$TRB...J.gene..2.,"[*]",2)[,1]
A1<-pubseq[pubseq$TRA...V.gene..1.!="",3:6]
A2<-pubseq[pubseq$TRA...V.gene..2.!="",7:10]
B1<-pubseq[pubseq$TRB...V.gene..1.!="",11:14]
B2<-pubseq[pubseq$TRB...V.gene..2.!="",15:18]
colnames(A1)=colnames(A2)=colnames(B1)=colnames(B2)=c("Chain","V.gene","D.gene","J.gene")
A<-rbind(A1,A2)
B<-rbind(B1,B2)
rm(A1,A2,B1,B2)
dim(unique(A))#151
dim(unique(B))#226

#collapse clones in our study across donors
setwd("/Users/yingy_adm/Documents/TCRdata/mixcr_umi_no_filter")
file_a<-list.files(recursive = TRUE,pattern ="TRA")
file_b<-list.files(recursive = TRUE,pattern ="TRB")

subset_mixcr_out<-function(mixcr_out){
  mixcr_out$allVHitsWithScore<-str_split_fixed(mixcr_out$allVHitsWithScore,"[*]",2)[,1]
  mixcr_out$allJHitsWithScore<-str_split_fixed(mixcr_out$allJHitsWithScore,"[*]",2)[,1]
  mixcr_out<-mixcr_out[,c('allVHitsWithScore','allJHitsWithScore','aaSeqCDR3')]
  return(mixcr_out)
}

for (i in file_a){
  filename = paste(unlist(strsplit(i,"/"))[2],"A",sep = "_")
  x = read.csv(i,sep = '\t')
  x<-subset_mixcr_out(x)
  assign(filename,x)
}

for (i in file_b){
  filename = paste(unlist(strsplit(i,"/"))[2],"B",sep = "_")
  x = read.csv(i,sep = '\t')
  x<-subset_mixcr_out(x)
  assign(filename,x)
}
rm(list = ls(pattern = "CD1426_"))

file_a<-as.data.frame(matrix(,0,3))
colnames(file_a)<-colnames(CD1393_A)
for (a in ls(pattern = "_A$")){
  afile=get(a)
  file_a<-rbind(file_a,afile)
  #print(dim(file_a))
}

file_b<-as.data.frame(matrix(,0,3))
colnames(file_b)<-colnames(CD1393_B)
for (a in ls(pattern = "_B$")){
  afile=get(a)
  file_b<-rbind(file_b,afile)
  print(dim(file_b))
}
dim(unique(file_a))#16616
dim(unique(file_b))#26097

rm(list = ls(pattern = "CD1"))

#collapse clones from scTCR_IRIS_allClonotypes(ref)
scTCR_IRIS_allClonotypes <- read.delim("~/Downloads/scTCR_IRIS_allClonotypes.txt")
scTCR_IRIS_allClonotypes$TRA...V.gene..1.<-str_split_fixed(scTCR_IRIS_allClonotypes$TRA...V.gene..1.,"[*]",2)[,1]
scTCR_IRIS_allClonotypes$TRA...J.gene..1.<-str_split_fixed(scTCR_IRIS_allClonotypes$TRA...J.gene..1.,"[*]",2)[,1]
scTCR_IRIS_allClonotypes$TRB...V.gene..1.<-str_split_fixed(scTCR_IRIS_allClonotypes$TRB...V.gene..1.,"[*]",2)[,1]
scTCR_IRIS_allClonotypes$TRB...J.gene..1.<-str_split_fixed(scTCR_IRIS_allClonotypes$TRB...J.gene..1.,"[*]",2)[,1]
scTCR_IRIS_allClonotypes$TRA...V.gene..2.<-str_split_fixed(scTCR_IRIS_allClonotypes$TRA...V.gene..2.,"[*]",2)[,1]
scTCR_IRIS_allClonotypes$TRA...J.gene..2.<-str_split_fixed(scTCR_IRIS_allClonotypes$TRA...J.gene..2.,"[*]",2)[,1]
scTCR_IRIS_allClonotypes$TRB...V.gene..2.<-str_split_fixed(scTCR_IRIS_allClonotypes$TRB...V.gene..2.,"[*]",2)[,1]
scTCR_IRIS_allClonotypes$TRB...J.gene..2.<-str_split_fixed(scTCR_IRIS_allClonotypes$TRB...J.gene..2.,"[*]",2)[,1]
A1<-scTCR_IRIS_allClonotypes[scTCR_IRIS_allClonotypes$TRA...V.gene..1.!="",c(3,4,6)]
A2<-scTCR_IRIS_allClonotypes[scTCR_IRIS_allClonotypes$TRA...V.gene..2.!="",c(7,8,10)]
B1<-scTCR_IRIS_allClonotypes[scTCR_IRIS_allClonotypes$TRB...V.gene..1.!="",c(11,12,14)]
B2<-scTCR_IRIS_allClonotypes[scTCR_IRIS_allClonotypes$TRB...V.gene..2.!="",c(15,16,18)]
colnames(A1)=colnames(A2)=colnames(B1)=colnames(B2)=c("Chain","V.gene","J.gene")
A<-rbind(A1,A2)
B<-rbind(B1,B2)
rm(A1,A2,B1,B2)
dim(unique(A))#2929
dim(unique(B))#2662
