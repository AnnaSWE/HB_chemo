#-------------------------------------------------------------------------------------
#
# Anna Wenger
#
#
# Code to calculate the probability of an SNV being caused by a particular SBS signature
# Note: Requires info on which signatures that are present in a sample

# Formula for calculation taken from S1 from
# https://ashpublications.org/blood/article/137/21/2992/475247/Clonal-hematopoiesis-and-therapy-related-myeloid

#-------------------------------------------------------------------------------------
#
# Libraries
#
#-------------------------------------------------------------------------------------

library(ggplot2)
library(gplots)
library(ggpubr)
library(GenomicRanges)
library(Rsamtools)
library(MASS)
library(viridis)
library(ComplexHeatmap)
library(plyr)
library(dplyr)

#-------------------------------------------------------------------------------------
#
# Read in samples
#
#-------------------------------------------------------------------------------------

samples=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/t-NS/Table_S1.txt", 
                   header = T, sep = '\t', stringsAsFactors = F)

#-------------------------------------------------------------------------------------
#
# Read in COSMIC reference signature and (optional) add your novel signature to it
#
#-------------------------------------------------------------------------------------

# Read in extracted novel signature from table S5
sig_hdp=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/t-NS/03_DNA_analyses/mSigHdp/02_Output/RunHdpxParallel_out/extracted.signatures.csv", 
                   header=TRUE, sep=",", stringsAsFactors = F)
rownames(sig_hdp)=paste0(substr(sig_hdp$Trinucleotide,1,1),
                         "[",
                         sig_hdp$Mutation.type,
                         "]",
                         substr(sig_hdp$Trinucleotide,3,3))
sig_hdp=sig_hdp[,c("Mutation.type","hdp.1")] #NOTE! Hardcoded for our novel sig
colnames(sig_hdp)=c("Mutation.type", "novel") #NOTE! Hardcoded for our novel sig


# Read in COSMIC REF signature
cosmic_37=read.table("/lustre/scratch126/casm/team274sb/aw35/targeted_NanoSeq/drivers/mSigHdp/01_Input/COSMIC_v3.4_SBS_GRCh37.txt", 
                     header=TRUE, row.names = 1, sep = '\t', stringsAsFactors = F)

# Put the mutation types into the same order in the cosmic_ref and the novel sig
reorder_idx=match(row.names(sig_hdp), row.names(cosmic_37))
cosmic_37=cosmic_37[reorder_idx, ]

# Combine COSMIC ref with the novel sig
SBS_ref=cbind(cosmic_37, sig_hdp[,"novel"])
colnames(SBS_ref)[ncol(SBS_ref)]="novel"


#-------------------------------------------------------------------------------------
#
# Read in signature contributions for the samples of interest and (if needed) 
# break down into SBS signatures (+ novel signature if applicable)
#
#-------------------------------------------------------------------------------------

#Read in decomposition of the extracted signatures into SBS
sigs_fraction=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/t-NS/03_DNA_analyses/mSigHdp/02_Output/signature_fractionR2_0.1.txt", header=T)

sigs_fraction=sigs_fraction[,2:3] #remove the novel signature (will not be decomposed into SBS signatures)
sigs_fraction=sigs_fraction[!rowSums(sigs_fraction)==0,]
colnames(sigs_fraction)=c("hdp.2", "hdp.3")

#Read in exposures for the samples
sigs=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/t-NS/03_DNA_analyses/mSigHdp/02_Output/RunHdpxParallel_out/inferred.exposures.csv", 
                header=T, sep=",", row.names = 1)
sigs=as.data.frame(t(as.matrix(sigs)))

#Calculate proportion of sigs
sigs_prop=sigs
sigs_prop$total=rowSums(sigs[,1:3])
for (i in 1:3) {
  sigs_prop[,i]=sigs_prop[,i]/sigs_prop$total
}

sigs_prop=sigs_prop[,1:3] #remove "total" column. It's done its job

sigs_prop=sigs_prop[,c(2,3,1)] #put novel sig last (to get correct calculations in the loop below)

SBS_samples=as.data.frame(matrix(0, nrow(sigs_prop),nrow(sigs_fraction))) #ncol=number of SBS signatures that has been extracted
colnames(SBS_samples)=row.names(sigs_fraction)



for(i in 1:nrow(sigs_prop)) {
  sample=row.names(sigs_prop)[i]
  row.names(SBS_samples)[i]=sample
  
  for(j in 1:(ncol(sigs_prop)-1)) { 
     #loop through the extracted signatures EXCEPT the novel signature
    print(colnames(sigs_prop)[j])
    SBS_samples[i,]=sigs_prop[i,j]*sigs_fraction[,j]+SBS_samples[i,]
  } 
}

#add on the proportion of novel sig
SBS_samples$novel=rep(NA, nrow(SBS_samples))
for(i in 1:nrow(SBS_samples)) {
  sample=row.names(SBS_samples)[i]
  
  idx=which(sample==row.names(sigs_prop))
  
  SBS_samples$novel[i]=sigs_prop[idx, "hdp.1"] #hard coded for novel sig

}

rowSums(SBS_samples) #check that rowSums are equal to 1 as a sense check of the calculations above

row.names(SBS_samples)=gsub('_.*', '', row.names(SBS_samples)) #remove tds suffix from sample names

#-------------------------------------------------------------------------------------
#
# Read in drivers of interest and select SNVs only
#
#-------------------------------------------------------------------------------------

annotated_all=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/t-NS/drivers/mutations_annotated_as_driver_candidates_25_03_14.txt",
                           header=TRUE, sep = '\t', stringsAsFactors = F) #table S7

annotated_all$type=ifelse(!annotated_all$impact=="no-SNV", "sub", "indel")
table(annotated_all$type)

#select all the drivers
driver_list_all=annotated_all[annotated_all$driver_candidate=="yes",]

# Only look at the SNVs indicated as drivers for calculating the probability of them arising from the SBS signatures
driver_list=driver_list_all[driver_list_all$type=="sub" & driver_list_all$driver_candidate=="yes",]

# Annotate pyr context of the drivers and use same annotation as COSMIC SBS signatures
genomeFile="/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/genome.fa"
driver_list$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(driver_list$chr, IRanges(as.numeric(driver_list$pos)-1, 
                                                                                       as.numeric(driver_list$pos)+1))))
ntcomp = c(T="A",G="C",C="G",A="T")
driver_list$sub = paste(driver_list$ref,driver_list$mut,sep=">")
driver_list$trinuc_ref_py = driver_list$trinuc_ref
for (j in 1:nrow(driver_list)) {
  if (driver_list$ref[j] %in% c("A","G")) { # Purine base
    driver_list$sub[j] = paste(ntcomp[driver_list$ref[j]],ntcomp[driver_list$mut[j]],sep=">")
    driver_list$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(driver_list$trinuc_ref[j],split="")[[1]])],collapse="")
  }
}

# This is the final annotation to use for the trinuc changes
driver_list$tri_change=paste0(substr(driver_list$trinuc_ref_py,1,1),
                              "[",
                              driver_list$sub,
                              "]",
                              substr(driver_list$trinuc_ref_py,3,3))

#-------------------------------------------------------------------------------------
#
# Loop through each SNV driver and calculate the probability of it being derived from each
# of the SBS signatures that were assigned to it
#
#-------------------------------------------------------------------------------------

#Create columns for the normalised probability for each of the SBS signatures + the novel sig
driver_list$SBS1_prob=rep(NA, nrow(driver_list)) #NOTE! Need to add other signatures manually here and in the "j" for loop below if you want additional signatures
driver_list$SBS5_prob=rep(NA, nrow(driver_list))
driver_list$SBS40a_prob=rep(NA, nrow(driver_list))
driver_list$SBS31_prob=rep(NA, nrow(driver_list))
driver_list$SBS35_prob=rep(NA, nrow(driver_list))
driver_list$novel_prob=rep(NA, nrow(driver_list))

for(i in 1:nrow(driver_list)) {
  print(i)
  
  trinuc=driver_list$tri_change[i]
  
  sample=driver_list$sampleID[i]
  
  #Look up the SBS contributions for the sample
  SBS=as.data.frame(SBS_samples[sample,])
  
  #Loop through each SBS in the sample
  total=0
  SBS_prob=SBS
  SBS_prob[]=0
  for(j in 1:ncol(SBS)) {
    SBS_temp=colnames(SBS)[j]
    SBS_prob[1,SBS_temp]=SBS[1,j]*SBS_ref[trinuc, SBS_temp]
    total=total+SBS_prob[1,SBS_temp]
  }
  
  driver_list$SBS1_prob[i]=SBS_prob[1,"SBS1"]/total
  driver_list$SBS5_prob[i]=SBS_prob[1,"SBS5"]/total
  driver_list$SBS40a_prob[i]=SBS_prob[1,"SBS40a"]/total
  driver_list$SBS31_prob[i]=SBS_prob[1,"SBS31"]/total
  driver_list$SBS35_prob[i]=SBS_prob[1,"SBS35"]/total
  driver_list$novel_prob[i]=SBS_prob[1,"novel"]/total #add more SBS/sigs here if needed
  
}

#Save table with results or plot data as desired
