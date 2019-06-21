#if a chain is matched, is the other chain from sc_pub_seq is also matched

seq_ExactMatch<-clonecount_no_filter_merged_pubseq_percentage
seq_ExactMatch$patient_group<-ifelse(seq_ExactMatch$donor %in% DQ2_CD,"CD","HC")
seq_ExactMatch$patients_shared_across<-NA
seq_ExactMatch$paired_match<-NA
seq_ExactMatch$other_seq<-NA
seq_ExactMatch$other_seq_patient<-NA
seq_ExactMatch$other_seq_patient_group<-NA
for (i in 1:dim(seq_ExactMatch)[1]){
  seq=seq_ExactMatch$matching_sequences[i]
  ab=seq_ExactMatch$chain[i]
  if (ab=="A"){
    id_pub<-which(pubseq$Chain..TRA..1.==seq |pubseq$Chain..TRA..2.==seq)
    other_seq<-gsub("NA,?","",paste(unique(pubseq$Chain..TRB..1.[id_pub]),collapse = ","))
    seq_ExactMatch$other_seq[i]<-other_seq
    if(length(unique(pubseq$Chain..TRB..1.[id_pub])>0)){seq_ExactMatch$paired_match[i]<-length(unique(pubseq$Chain..TRB..1.[id_pub]))}
       #both A and B are found in single cell pubseq
    id_ext_pub<-which(extract_pubseq$Defining_sequence==seq & extract_pubseq$Chain=="TRA")
    seq_ExactMatch$patients_shared_across[i]<-sum(extract_pubseq$Shared_across_items[id_ext_pub])
  }
  else if (ab=="B"){
    id_pub<-which(pubseq$Chain..TRB..1.==seq |pubseq$Chain..TRB..2.==seq)
    other_seq<-gsub("NA,?","",paste(unique(pubseq$Chain..TRA..1.[id_pub]),collapse = ","))
    seq_ExactMatch$other_seq[i]<-other_seq
    if(length(unique(pubseq$Chain..TRA..1.[id_pub])>0)){seq_ExactMatch$paired_match[i]<-length(unique(pubseq$Chain..TRA..1.[id_pub]))}
    id_ext_pub<-which(extract_pubseq$Defining_sequence==seq & extract_pubseq$Chain=="TRB")
    seq_ExactMatch$patients_shared_across[i]<-sum(extract_pubseq$Shared_across_items[id_ext_pub])
  }
  if (length(str_split(other_seq,","))>0){
    id_clonecount_df<-clonecount_no_filter_merged_pubseq_percentage$sequence %in% unlist(str_split(other_seq,","))
    id_clonecount_df<-id_clonecount_df[!is.na(id_clonecount_df)]
    if(sum(id_clonecount_df)>0){
      seq_ExactMatch$other_seq_patient[i]<-paste(clonecount_no_filter_merged_pubseq_percentage$donor[id_clonecount_df],collapse = ",")
      other_seq_donor_group<-get_donor_group(clonecount_no_filter_merged_pubseq_percentage$donor[id_clonecount_df])#"sameHC","diffHC","CD"
      seq_ExactMatch$other_seq_patient_group[i]<-paste(other_seq_donor_group,collapse = ",")
    }
  }
}
seq_ExactMatch%>%subset(chain=="B")%>% nrow #77
seq_ExactMatch%>%subset(chain=="B" & other_seq!="")%>%nrow #60 out of 77 matched beta chains has a paired alpha chain in ref
seq_ExactMatch_b<-seq_ExactMatch%>%subset(chain=="B")
seq_ExactMatch_a<-seq_ExactMatch%>%subset(chain=="A")

#function to get "sameCD":if paired other_seq was found in the same donor
get_sameCD<-function(seq_ExactMatch){
  list_other_seq_patient<-str_split(seq_ExactMatch$other_seq_patient,",")
  seq_ExactMatch$sameCD<-NA
  for (i in 1:nrow(seq_ExactMatch)){
    donor<-seq_ExactMatch$donor[i]
    seq_ExactMatch$sameCD[i] <-donor %in% list_other_seq_patient[[i]]}
  return(seq_ExactMatch)
}

seq_ExactMatch <- get_sameCD(seq_ExactMatch)
seq_ExactMatch_b <- get_sameCD(seq_ExactMatch_b)
seq_ExactMatch_a <- get_sameCD(seq_ExactMatch_a)
table(seq_ExactMatch$patient_group,seq_ExactMatch$sameCD)
table(seq_ExactMatch_b$patient_group,seq_ExactMatch_b$sameCD)
table(seq_ExactMatch_a$patient_group,seq_ExactMatch_a$sameCD)
