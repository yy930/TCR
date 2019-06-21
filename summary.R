library(stringr)
all_report <- read.table("~/Documents/TCRdata/mixcr_result/all_report.txt", quote="", 
                         comment.char="",fill = T,stringsAsFactors = F)
all_report[is.na(all_report)] <- 0
all_report[all_report==""] <- 0
all_report$V2 <- str_split_fixed(str_split_fixed(all_report$V2,"%",2)[,1],"[(]",2)[,2]
all_report$V2 <- as.numeric(all_report$V2)/100
colnames(all_report) <- c("Number_mapped_reads","Proportion_mapped_reads","Number_reads_mapped_to_top_clone",
                          "Proportion_reads_mapped_to_top_clone")
all_report$Number_mapped_reads <- as.numeric(all_report$Number_mapped_reads)
all_report$Number_reads_mapped_to_top_clone <- as.numeric(all_report$Number_reads_mapped_to_top_clone)

all_report$Number_reads <- all_report$Number_mapped_reads/all_report$Proportion_mapped_reads
Number_reads <- all_report$Number_reads[!is.na(all_report$Number_reads)]
Number_reads<-Number_reads[Number_reads!=Inf]

sum(all_report$Proportion_reads_mapped_to_top_umi==1)
Proportion_reads_mapped_to_top<-all_report$Proportion_reads_mapped_to_top_umi
Proportion_reads_mapped_to_top<-Proportion_reads_mapped_to_top[Proportion_reads_mapped_to_top>0]
Proportion_reads_mapped_to_top<-Proportion_reads_mapped_to_top[Proportion_reads_mapped_to_top<1]

library(data.table)
matched <- read.csv("~/Downloads/manually_assembled_matched_info_overview.csv", comment.char="#")
n_seq_persample<-lapply(flat_list_patientAB, function(x) dim(x)[1])
matched$n_seq<-as.numeric(n_seq_persample[paste0(matched$patient_id,"_clones_TR",matched$chains)])
DT<-as.data.table(matched)
#use sum of all sequences
all_tab<-DT[,.(sum(n_seq),sum(sequences_matched)),CD]
all_tab<-t(all_tab[,.(V2,V1-V2)])
fisher.test(all_tab, alternative = "greater")#p-value = 1.345e-08
#test for each patient
a_tab<-DT[chains=="A",.(sum(n_seq),sum(sequences_matched)),CD]
a_tab<-t(a_tab[,.(V2,V1-V2)])
fisher.test(a_tab, alternative = "greater")#p-value = 0.4767
b_tab<-DT[chains=="B",.(sum(n_seq),sum(sequences_matched)),CD]
b_tab<-t(b_tab[,.(V2,V1-V2)])
fisher.test(b_tab, alternative = "greater")#p-value = 1.021e-09



#Logistic regression
DTB<-DT[chains=="B",]
DTB$pct<-DTB$sequences_matched/DTB$n_seq
DTB$CD<-as.numeric(DTB$CD)-1
fit<-glm(formula = CD ~ pct, family=binomial, data = DTB)
DTB$pred_cd = predict(fit, newdata=DTB, type="response")
#to plot
newdat <- data.frame(pct=seq(min(DTB$pct), max(DTB$pct),len=100))
newdat$CD = predict(fit, newdata=newdat, type="response")
plot(CD~pct, data=DTB, col="red4")
lines(CD~pct, newdat, col="green4", lwd=2)
abline(v=(0.95+7.452)/6288.692)

DTB<-list_matchNum[['matchNum_no_filter_only_beta_pubseq_percentage']]
DTB$correct_pred<-as.factor(DTB$CD_status==DTB$pred_cd_proportion_sequences)
DTB$CD_status_numric<-as.numeric(as.logical(DTB$CD_status))
DTB<-subset(DTB,repertoire_size>600)#!!!!alt
fit<-glm(formula = CD_status ~ proportion_sequences_matched, family=binomial, data = DTB)
DTB$pred_cd = predict(fit, newdata=DTB, type="response")#!!!!alt
DTB$pred_cd_proportion_sequences<-(DTB$pred_cd>0.5)#!!!!alt
evaluate_pred_seq(DTB)#0.9285714 1.0000000 0.8571429 1.0000000
newdat <- data.frame(proportion_sequences_matched=seq(min(DTB$proportion_sequences_matched), max(DTB$proportion_sequences_matched),len=100))
newdat$CD_status = predict(fit, newdata=newdat, type="response")

ggplot(DTB, aes(x=proportion_sequences_matched, y=pred_cd, color=repertoire_size, shape = CD_status)) +
  geom_point(size = 3) +
  scale_color_gradient(low="blue", high="red") +
  #geom_label_repel(aes(label = patient),box.padding = 0.15, point.padding = 1) +
  geom_line(data=newdat,aes(x=proportion_sequences_matched,y=CD_status),color="grey",inherit.aes = FALSE)


#Prevalence and frequency is positively associated for the public TCR clonetypes 
library(stringr)
i=list.dirs("/Users/yingy_adm/Documents/TCRdata/all_results/clonal_percentage",recursive = F)[11]
for (icd in list.files(i,pattern="CD*")){
  filename=paste("matchNum",gsub("_['B'].tsv","",icd),sep = "_")
  x = read.csv(paste(i,"/",icd,sep = ""),sep = '\t')
  x$chain = str_split(basename(i),"_")[[1]][10]
  x$ref = gsub("iris","",str_split(basename(i),"_")[[1]][11])
  x$patient_id = gsub("_['B'].tsv","",icd)
  x = x[!is.na(x$matching_sequences),]
  x = x[x$matching_sequences!="",]
  x$matching_sequences = unlist(lapply(x$matching_sequences,FUN = function(x) str_split_fixed(x,",",2)[1]))
  assign(filename,x)}

seq_ExactMatch<-data.frame()
for (i in ls(pattern = "matchNum_CD")){
  seq_ExactMatch = rbind(seq_ExactMatch,get(i))
}

pubseq <- read.csv("~/Documents/TCRdata/pubseq.csv", sep=";")
extract_pubseq <- read.csv2("~/Documents/TCRdata/extract_pubseq.csv")
DQ2_HC <- c("CD1331", "CD1386", "CD1408", "CD1409", "CD1440", "CD1450")#to drop:"CD1331""CD1440"
DQ8_HC <- c("CD1370", "CD1390", "CD1443", "CD1453", "CD1465")#to drop:CD1443,CD1465
DQ2_CD <- c("CD1357", "CD1358", "CD1364", "CD1368", "CD1393", "CD1422", "CD1426", "CD1451")#to drop:CD1422,CD1451

#if a chain is matched, is the other chain from sc_pub_seq is also matched
seq_ExactMatch$paired_match<-NA
seq_ExactMatch$paired_match_patient<-NA
seq_ExactMatch$other_seq<-NA
seq_ExactMatch$patient_group<-ifelse(seq_ExactMatch$patient_id %in% DQ2_CD,"CD","HC")
seq_ExactMatch$patients_shared_across<-NA

for (i in 1:dim(seq_ExactMatch)[1]){
  seq=seq_ExactMatch$matching_sequences[i]
  ab=seq_ExactMatch$chain[i]
  if (ab=="alpha"|ab=="A"){
    id_pub<-which(pubseq$Chain..TRA..1.==seq |pubseq$Chain..TRA..2.==seq)
    other_seq<-pubseq$Chain..TRB..1.[id_pub]
    seq_ExactMatch$other_seq[i]<-paste(unique(other_seq),collapse = ",")
    if(length(other_seq)>0){seq_ExactMatch$paired_match[i]<-length(other_seq)}#both A and B are found in single cell pubseq
    id_ext_pub<-which(extract_pubseq$Defining_sequence==seq & extract_pubseq$Chain=="TRA")
    seq_ExactMatch$patients_shared_across[i]<-sum(extract_pubseq$Shared_across_items[id_ext_pub])
  }
  else if (ab=="beta"|ab=="B"){
    id_pub<-which(pubseq$Chain..TRB..1.==seq |pubseq$Chain..TRB..2.==seq)
    other_seq<-pubseq$Chain..TRA..1.[id_pub]
    seq_ExactMatch$other_seq[i]<-paste(unique(other_seq),collapse = ",")
    if(length(other_seq)>0){seq_ExactMatch$paired_match[i]<-length(other_seq)}
    id_ext_pub<-which(extract_pubseq$Defining_sequence==seq & extract_pubseq$Chain=="TRB")
    seq_ExactMatch$patients_shared_across[i]<-sum(extract_pubseq$Shared_across_items[id_ext_pub])
  }
  if (length(other_seq)>0){
    seq_ExactMatch$paired_match_patient[i]<-paste(seq_ExactMatch$patient_id[match(unique(other_seq),seq_ExactMatch$matching_sequences)],collapse = ",")
  }
}

tab_Nmatch_Npatient<-as.data.frame(table(seq_ExactMatch$matching_sequences,seq_ExactMatch$patient_id))%>%subset(Freq>0)
tab_Nmatch_Npatient$Npatient<-seq_ExactMatch$patients_shared_across[match(tab_Nmatch_Npatient$Var1,seq_ExactMatch$matching_sequences)]
cor.test(tab_Nmatch_Npatient$Freq, tab_Nmatch_Npatient$Npatient,method = "kendall", alternative = "greater")#p-value = 0.0002131
cor.test(tab_Nmatch_Npatient$Freq, tab_Nmatch_Npatient$Npatient,method = "spearm", alternative = "g")#p-value = 0.0001425
                                                                                                                                 