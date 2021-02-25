library(readstata13)
library(stats)
library(vegan)
library(foreign)
library(lme4)
library(lmerTest)
library(gdata)
library(readxl)
library(Hmisc)
library(ggplot2)
library(survival)
library(survminer)
library(ggfortify)

#set pathway for datafiles
setwd("")

#patient data
patdata<-read.dta13("pat.dta")

#backup of unchanged, raw patient data file
rawpatdata<-patdata

#row names are set equal to ids
row.names(patdata)<-patdata[,1]

#load HLA data for patients
hladata<-read.dta13(".dta")
rownames(hladata)<-hladata[,'ID']
hladata<- hladata[,-1]
colnames(hladata)[1]<- "id"

#make list of all ids of patients with HLA data
hlarows<-rownames(hladata)

#cleanup step
#replace "NA" with actual NA
for(abc in 1:nrow(hladata)){
  for(def in 1:ncol(hladata))
    if(!is.na(hladata[abc, def])){
      if(hladata[abc, def]== "NA"){
        hladata[abc, def]<- NA
      }
    }
}


#generates HLA-A 2 digit alleles from the 4 digit alleles
hladata$id<- rownames(hladata)
hladata$A1<- substr(hladata$A1g, 1,2)
hladata$A2<-substr(hladata$A2g, 1,2)

#set pathway
setwd("")

#load mutations
mutdata<- read.dta13("mutations.dta")
#make string of protein AND mutation name
mutdata$mut_full2<- paste(mutdata$seqtype, mutdata$mut_full, sep="_")


#backup of original mutdata file
rawmutdata<-mutdata

#load file with sequencing test info
testinfo<- read.dta13("test_info.dta")
#backup of original testinfo file
rawtestinfo<- testinfo

#resistance testing sequence length values cleaned up
testinfo$rt_len[is.na(testinfo$rt_len)]<-0
testinfo$pr_len[is.na(testinfo$pr_len)]<-0
testinfo$int_len[is.na(testinfo$int_len)]<-0
testinfo$gp120_len[is.na(testinfo$gp120_len)]<-0
testinfo$gp41_len[is.na(testinfo$gp41_len)]<-0
#total length of sequence from the sum
testinfo$len <- testinfo$rt_len  + testinfo$pr_len + testinfo$int_len + testinfo$gp120_len + testinfo$gp41_len
#only keep sequences of length over 200
testinfo<-testinfo[testinfo$len>200,]

#dummy variable on which rows to keep, 1 is true
testinfo$keep<-1


#sequences ordered by patient id, within each id ordered by sample date
testinfo<- testinfo[order(testinfo$sampledate),]
testinfo<- testinfo[order(testinfo$id),]
rownames(testinfo)<- 1:nrow(testinfo)

#if 2 seqs of same person and same date, keep one with longer length
for(abc in 2:nrow(testinfo)){
  if((testinfo[abc, "sampledate"] == testinfo[abc-1, "sampledate"]) &&
     (testinfo[abc, "id"] == testinfo[abc-1, "id"])){
    
    if(testinfo[abc-1, "len"]< testinfo[abc, "len"]){
      testinfo[abc-1, "keep"]<-0
    }else{
      testinfo[abc, "keep"]<-0
    }
  }
}
testinfo<-testinfo[testinfo$keep==1,]



#list of mutations we are looking at, from stanford
library(readr)
stanmutlist <- read_delim("stanmutlist.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
#make strings: gene plus mutation
resmutlist2<- paste(stanmutlist$gene, stanmutlist$mut, sep="_")
#another one, removing what the specific new amino acid  is
mutlisttrunc<- unique(substr(resmutlist2, 1, nchar(resmutlist2)-1))



#all art drugs
artnamelist<- c("ENF"   , "MVC","ABC" ,"DDL" ,"FTC" , "3TC" , "D4T" , "AZT" , 
                "DDC" ,"TDF" , "TAF" ,"EFV" , "NVP",   "DLV" ,"ETR" ,  "RPV", "DTG" ,"EVG", "RAL" ,
                "APV" , "FPV" ,  "IDV",    "LPV",  "NFV" ,"RTV" ,"SQV","ATV" ,"DRV" , "TPV" ,"COB" ,"RTV" )

#only use patients with hla data
patdata<-patdata[patdata$id %in% hladata$id,]

#only use patients with sequences
patdata<-patdata[patdata$id %in% rawtestinfo$id,]

#variable on ART start date
patdata$artstart<- NA


drugdata<-read.dta13("drug.dta")
#only art drugs, and only drug data for patients with hla data
drugdata<- drugdata[drugdata$id %in% patdata$id,]
drugdata<-drugdata[drugdata$drug %in% artnamelist,]



#keep sequences only for patients indicated above
testinfo<- testinfo[testinfo$id %in% patdata$id,]

#create variable storing sequencing id and date
patdata$testid<-NA
patdata$testdate<-NA


#patdata is no longer one row per patient, but rather one row per sequence (repeated data for each patient)
#hence number of times patient is listed in patdata is how many sequences are in testinfo

mnopq<-0
for(abc in 1:nrow(testinfo)){
  if((abc>1) && (testinfo[abc, "id"]== testinfo[abc-1,"id"])){
   
    mnopq<-mnopq +1
    temprows<- patdata[patdata$id== testinfo[abc, "id"],]
    patdata<- rbind(patdata, temprows[1,])
    
    patdata[nrow(patdata), "testid"]<- testinfo[abc, "header_id"]
    patdata[nrow(patdata), "testdate"]<- testinfo[abc, "sampledate"]
    
  }else{
    patdata[patdata$id== testinfo[abc, "id"], "testid"]<- testinfo[abc, "header_id"]
    patdata[patdata$id== testinfo[abc, "id"], "testdate"]<- testinfo[abc, "sampledate"]
  }
}
patdata<- patdata[order(patdata$testdate),]
patdata<- patdata[order(patdata$id),]
rownames(patdata)<- 1:nrow(patdata)

#cleanup of ethnicity/white variable
patdata$white<-0
patdata$white[patdata$ethnicity == 1]<-1
patdata$white[patdata$ethnicity == 9]<-NA

#artstart is first ART drug start date
for(abc in 1:nrow(patdata)){
  tempmat<- drugdata[drugdata$id == patdata[abc, "id"],]
  if(nrow(tempmat)>=1){
    datelist<- tempmat$drstartd
    patdata[abc, "artstart"]<- min(datelist)
  }
}


#only keep sequences before art start
patdata$keep<-1
for(abc in 1:nrow(patdata)){
  if(!is.na(patdata[abc, 'artstart'])){
    if(patdata[abc, "testdate"]> patdata[abc, "artstart"]){
      patdata[abc, "keep"]<-0
    }
  }
}
patdata<- patdata[patdata$keep ==1,]
rownames(patdata)<- 1:nrow(patdata)


#only keep mutations in the stanford mutation list
mutdata$id<- NA
mutdata<-mutdata[mutdata$mut_full2 %in% resmutlist2,]


#only keep mutations detected from the selected testinfo sequences
for(abc in 1:nrow(mutdata)){

  if(mutdata[abc, "header_id"] %in% testinfo$header_id){
    mutdata[abc,"id"]<- testinfo[testinfo$header_id == mutdata[abc, "header_id"], "id"]  
  }
}

#throwout mutation rows with empty ids
mutdata<-mutdata[!is.na(mutdata$id),]

#date sample taken
mutdata$datelab<- NA


#load dates into mutation dataframe
for(abc in 1:nrow(mutdata)){
  mutdata[abc, "datelab"]<- testinfo[testinfo$header_id == mutdata[abc, "header_id"],"sampledate"]
}

#only mutdata for patients of interest
mutdata<- mutdata[mutdata$id %in% patdata$id,]


#only keep mutations detected before artstart (i.e. ART-naive)
mutdata$pretreat<-0
for(abc in 1:nrow(mutdata)){
  if(is.na(patdata[patdata$id==mutdata[abc, "id"], "artstart"])[1]){
    mutdata[abc, "pretreat"]<-1
  }else{
    if(mutdata[abc, "datelab"]<= patdata[patdata$id==mutdata[abc, "id"], "artstart"][1]){
      mutdata[abc, "pretreat"]<-1
    }
  }
}
mutdata<-mutdata[mutdata$pretreat==1,]


#add dozens of new columns to patdata, one column for each mutation, preset to zero (absent)
for(abc in 1:length(resmutlist2)){
  patdata[, resmutlist2[abc]]<- 0
}


#variables below with name like "First" actually refer to last

#refers to last
mutdata$isfirst<-0
######MUTATIONS FROM LAST PRE-ART SEQUENCE ONLY
for(abc in 1:nrow(mutdata)){
  datelist<- testinfo[testinfo$id == mutdata[abc, "id" ],"sampledate"]
  firstseqdate<- max(datelist)
  
 if(as.Date(mutdata[abc, "datelab"], origin="1970-01-01") == firstseqdate){
    mutdata[abc, "isfirst"]<-1
    patdata[patdata$testid == mutdata[abc, "header_id"], mutdata[abc, "mut_full2"]]<-1
 }
}
mutdatafirstonly<- mutdata[mutdata$isfirst==1,]





#### new variables added to patdata for each mutation, regardless of which exact mutation amino acid it is
for(abc in 1:length(mutlisttrunc)){
  patdata[, mutlisttrunc[abc]]<- 0
  submuts<- resmutlist2
  submuts<- submuts[grepl(mutlisttrunc[abc], submuts)]
  
  for (def in 1:length(submuts)){
    patdata[, mutlisttrunc[abc]]<- (patdata[, mutlisttrunc[abc]] + patdata[, submuts[def]]) 
  }
}



#HLA TYPES
patdata$A1<- NA
patdata$A2<- NA
patdata$B1<- NA
patdata$B2<- NA
patdata$C1<- NA
patdata$C2<- NA


patdata$DRB1<- NA
patdata$DRB2<- NA
patdata$DQA1<- NA
patdata$DQA2<- NA
patdata$DQB1<- NA
patdata$DQB2<- NA
patdata$DPA1<- NA
patdata$DPA2<- NA
patdata$DPB1<- NA
patdata$DPB2<- NA







#load HLA data from hla dataframe
for(abc in 1:nrow(patdata)){
  patdata[abc, "A1"]<- hladata[hladata$id == patdata[abc, "id"], "A1"]
  patdata[abc, "A2"]<- hladata[hladata$id == patdata[abc, "id"], "A2"]
  patdata[abc, "B1"]<- hladata[hladata$id == patdata[abc, "id"], "B1"]
  patdata[abc, "B2"]<- hladata[hladata$id == patdata[abc, "id"], "B2"]
  patdata[abc, "C1"]<- hladata[hladata$id == patdata[abc, "id"], "C1"]
  patdata[abc, "C2"]<- hladata[hladata$id == patdata[abc, "id"], "C2"]
  
  patdata[abc, "DRB1"]<- hladata[hladata$id == patdata[abc, "id"], "DRB1"]
  patdata[abc, "DRB2"]<- hladata[hladata$id == patdata[abc, "id"], "DRB2"]
  patdata[abc, "DQA1"]<- hladata[hladata$id == patdata[abc, "id"], "DQA1"]
  patdata[abc, "DQA2"]<- hladata[hladata$id == patdata[abc, "id"], "DQA2"]
  patdata[abc, "DQB1"]<- hladata[hladata$id == patdata[abc, "id"], "DQB1"]
  patdata[abc, "DQB2"]<- hladata[hladata$id == patdata[abc, "id"], "DQB2"]
  patdata[abc, "DPA1"]<- hladata[hladata$id == patdata[abc, "id"], "DPA1"]
  patdata[abc, "DPA2"]<- hladata[hladata$id == patdata[abc, "id"], "DPA2"]
  patdata[abc, "DPB1"]<- hladata[hladata$id == patdata[abc, "id"], "DPB1"]
  patdata[abc, "DPB2"]<- hladata[hladata$id == patdata[abc, "id"], "DPB2"]
}

patdata$hlaa<- paste(patdata$A1, patdata$A2, sep="__")
patdata$hlab<- paste(patdata$B1, patdata$B2, sep="__")
patdata$hlac<- paste(patdata$C1, patdata$C2, sep="__")

patdata$hlaDRB<- paste(patdata$DRB1, patdata$DRB2, sep="__")
patdata$hlaDQA<- paste(patdata$DQA1, patdata$DQA2, sep="__")
patdata$hlaDQB<- paste(patdata$DQB1, patdata$DQB2, sep="__")
patdata$hlaDPA<- paste(patdata$DPA1, patdata$DPA2, sep="__")
patdata$hlaDPB<- paste(patdata$DPB1, patdata$DPB2, sep="__")







for(tempnum in 1:3)
{
  if(tempnum==1){
    hlaABC<- "A"
  }else if(tempnum==2){
    hlaABC<-"B"
  }else{
    hlaABC<-"C"
  }
  
  # make list of all HLA alleles
  if(hlaABC== "A"){
    hlab1<- unique(patdata$A1)
    hlab2<- unique(patdata$A2)
  }else if(hlaABC == "B"){
    hlab1<- unique(patdata$B1)
    hlab2<- unique(patdata$B2)
  }else if(hlaABC == "C"){
    hlab1<- unique(patdata$C1)
    hlab2<- unique(patdata$C2)
  }
  
  
  
  #make list of all HLA alleles within -A, -B or -C
  hlabtypes<- c(hlab1, hlab2)
  hlabtypes<- unique(hlabtypes)
  hlabtypes<- hlabtypes[hlabtypes!="0"]
  hlabtypes<- hlabtypes[!is.na(hlabtypes)]
  
  
  #selects last ART-naive sequence
  patdata$isMin<-0
  for(abc in 1:nrow(patdata)){
    datelist<- patdata[patdata$id== patdata[abc, "id"], "testdate"]
    if(patdata[abc, "testdate"]== max(datelist)){
      patdata[abc, "isMin"]<-1
    }
  }
  
  #save copy of patdata here, before further steps
  patdataALLseqs<- patdata
  
  #only keep last sequence, one per patient
  patdata<-patdata[patdata$isMin==1,]
  
  
  
  
  
  
  
  ##table to calculate power: all hla alleles X mutations
  powertable <- matrix(, nrow = length(hlabtypes), ncol = length(mutlisttrunc))
  powertable <- as.data.frame(powertable)
  rownames(powertable)<- hlabtypes
  colnames(powertable)<- mutlisttrunc
  
  
  
  
  #count variable keeping track of all combinations
  looptally<-0
  
  #count variable, counting how many with sufficient power
  numpower<-0
  
  for(abc in 1:nrow(powertable)){
    for(def in 1:ncol(powertable)){
      
      looptally<- (looptally+1)
      if(looptally %% 250 == 0){
        Sys.sleep(0.01)
        print(looptally)
      } 
      
      #make strings of hla allele names per the "correct" format
      if(hlaABC== "A"){
        localtable<- table(grepl(hlabtypes[abc],patdata$hlaa), patdata[,mutlisttrunc[def]]>0) 
        
        testname<- paste(colnames(powertable)[def], 
                         paste("HLA-A", rownames(powertable)[abc], sep=""), sep="_" )
      }else if(hlaABC == "B"){
        localtable<- table(grepl(hlabtypes[abc],patdata$hlab), patdata[,mutlisttrunc[def]]>0) 
        
        testname<- paste(colnames(powertable)[def], 
                         paste("HLA-B", rownames(powertable)[abc], sep=""), sep="_" )
      }else if(hlaABC == "C"){
        localtable<- table(grepl(hlabtypes[abc],patdata$hlac), patdata[,mutlisttrunc[def]]>0)
        testname<- paste(colnames(powertable)[def], 
                         paste("HLA-C", rownames(powertable)[abc], sep=""), sep="_" )
      }
      
      
      
      #if all FALSES/absent above, then add dummy column of 0/0 so that next part of code will run
      if(ncol(localtable)==1){
        localtable<- cbind(localtable, c(0,0))
      }
      
      #power calculation
      powervalue <-as.numeric( bpower(p1= (localtable[1,2]/(localtable[1,2] + localtable[1,1])),
                                      p2= (localtable[2,2]/(localtable[2,2] + localtable[2,1])),
                                      odds.ratio = 3,
                                      n1=(localtable[1,2] + localtable[1,1]),
                                      n2=(localtable[2,2] + localtable[2,1]) ,
                                      alpha=0.05)   )
      
      
      
      if(is.nan(powervalue)){
        powertable[abc, def]<-NA
      }else{
        if(powervalue>=0.8){    #power threshold: 0.8
          powertable[abc, def]<-powervalue
          numpower<- (numpower+1)
        }else{
          powertable[abc, def]<-NA
        }
      }
    }
  }
  
  
  #make list of candidate pairs
  candidatesA<-NA
  candidatesB<-NA
  candidatesC<-NA
  
  for(abc in 1:nrow(powertable)){
    for(def in 1:ncol(powertable)){
      if(!is.na(powertable[abc, def])){
        if(hlaABC== "A"){
          candidatesA<- c(candidatesA, paste(colnames(powertable)[def],
                                             paste("HLA-A", rownames(powertable)[abc], sep=""), sep="_" ))
        }else if(hlaABC == "B"){
          candidatesB<- c(candidatesB, paste(colnames(powertable)[def],
                                             paste("HLA-B", rownames(powertable)[abc], sep=""), sep="_" ))
        }else{
          candidatesC<- c(candidatesC, paste(colnames(powertable)[def],
                                             paste("HLA-C", rownames(powertable)[abc], sep=""), sep="_" ))
        }
      }
    }
  }
  
  if(hlaABC== "A"){
    candidatesA<- candidatesA[!is.na(candidatesA)]
  }else if(hlaABC == "B"){
    candidatesB<- candidatesB[!is.na(candidatesB)]
  }else{
    candidatesC<- candidatesC[!is.na(candidatesC)]
  }
  
  
  patdata$queriedhla<- 0    #specific hla allele being analyzed, set to 0/not
  patdata$queriedmut<- 0    #specific mutation being analyzed, set to 0/absent
  
  
  #set pathway
  setwd("")
  #separate file with calculations for time since infection
  seqtimes<- read.csv("datesForSequences.csv")
  
  #ignore NAs
  seqtimes<- seqtimes[!is.na(seqtimes$deltaTime),]
  
  #variable if precise time calculation is available
  patdata$preciseseqtime<-1
  
  #patient always had been ART-naive
  patdata$alwaysartnaive<-1
  patdata$alwaysartnaive[patdata$id %in% drugdata$id]<-0
  
  #days infected with HIV
  patdata$daysinfected<-NA
  seqtimes$sampledate2<-NA
  
  
  #time sample was taken, linked from testinfo to seqtimes
  for(abc in 1:nrow(seqtimes)){
    seqtimes[abc, "sampledate2"]<-rawtestinfo[rawtestinfo$header_id == seqtimes[abc, "header_id"], "sampledate"]
  }
  
  
  #loading times since infection
  for(abc in 1:nrow(patdata)){
    timetemp<- seqtimes[seqtimes$header_id == patdata[abc, "testid"], "deltaTime"]
    patdata[abc, "daysinfected"] <- timetemp
  }
  
  #make any negative/infinite calculated values zero
  patdata$preciseseqtime[patdata$daysinfected<0]<-0
  patdata$daysinfected[patdata$daysinfected <0]<-0
  patdata$preciseseqtime[is.infinite(patdata$daysinfected)]<-0
  patdata$daysinfected[is.infinite(patdata$daysinfected)]<-0
  
  
  
  overridecount<-0
  for(abc in 1:nrow(patdata)){
    #use makeshift calculation for days infected if estimation from file has negative/infinite value
    startdatelist<- c(patdata[abc, "consent_date"], patdata[abc, "hiv_posdate"], 
                      patdata[abc, "hiv_posdocdate"])
    startdatelist<- startdatelist[!is.na(startdatelist)]
    posdate<- min(startdatelist)
    
    timelist<- rawtestinfo[rawtestinfo$header_id== patdata[abc, "testid"], "sampledate"]
    
    if ((timelist-posdate)> patdata$daysinfected[abc]){
      patdata[abc, "daysinfected"] <- (timelist-posdate)
      overridecount<- overridecount+1
    }
    
  }
  
  #convert days to years
  patdata$yearsinfected <- patdata$daysinfected / 365.25
  
  #keep track of which patient sequences have exact date
  patdata$exactdate<-FALSE
  for(abc in 1:nrow(patdata)){
    if (seqtimes[seqtimes$header_id == patdata[abc, "testid"], "dateExact"] != 0){
      patdata[abc, "exactdate"]<-TRUE
    }
  }
  
  
  if(hlaABC=="A"){
    candidatetable<-data.frame(matrix(vector(), length(candidatesA), 6,
                                      dimnames=list(c(), c("candidate", "fisher p", "unadj OR", "adj p"
                                                           ,"hla 95", "hla 95"))),
                               stringsAsFactors=F)
    #column names above refer to: candidate pair names, p value of fisher test of contingency table,
    #unadjusted Odds Ratio of having DRM with given HLA type, OR when adjusted for infection time,
    #and 95% CIs for that OR

    candidatetable<-candidatetable[1:length(candidatesA),]
    
    candidates<- candidatesA
  }else if (hlaABC=="B"){
    candidatetable<-data.frame(matrix(vector(), length(candidatesB), 6,
                                      dimnames=list(c(), c("candidate", "fisher p", "unadj OR", "adj p"
                                                           ,"hla 95", "hla 95"))),
                               stringsAsFactors=F)
    candidatetable<-candidatetable[1:length(candidatesB),]
    
    candidates<-candidatesB
  }else{
    candidatetable<-data.frame(matrix(vector(), length(candidatesC), 6,
                                      dimnames=list(c(), c("candidate", "fisher p", "unadj OR", "adj p"
                                                           ,"hla 95", "hla 95"))),
                               stringsAsFactors=F)
    candidatetable<-candidatetable[1:length(candidatesC),]
    
    candidates <-candidatesC
  }
  
  
  
  
  ######################## calculate log regression HRs for all qualifying pairs
  for (pqr in 1:length(candidates)){
    abc<- which(rownames(powertable)== substr(candidates[pqr],nchar(candidates[pqr])-1 , nchar(candidates[pqr])))
    
    def<- which(colnames(powertable)== substr(candidates[pqr],1 , nchar(candidates[pqr])-8))
    
    
    candidatetable[pqr, 1]<- candidates[pqr]
    
    patdata$queriedhla<- 0
    
    for(mno in 1:nrow(patdata)){
      if(hlaABC== "A"){
        if(grepl(hlabtypes[abc], patdata[mno, "hlaa"])){
          patdata[mno, "queriedhla"]<- 1
        }
      }else if(hlaABC == "B"){
        if(grepl(hlabtypes[abc], patdata[mno, "hlab"])){
          patdata[mno, "queriedhla"]<- 1
        }
      }else{
        if(grepl(hlabtypes[abc], patdata[mno, "hlac"])){
          patdata[mno, "queriedhla"]<- 1
        }
      }
    }
    patdata$queriedmut<- patdata[,mutlisttrunc[def]]
    
    table(patdata$queriedhla, patdata$queriedmut>0)
    ftest<-fisher.test(table(patdata$queriedhla, patdata$queriedmut>0))
    candidatetable[pqr, 2] <- round(ftest$p.value, 3)
    
    patdata$queriedmut <-( patdata$queriedmut>0)
    
    modelmut = glm(queriedmut~ queriedhla+yearsinfected, family=binomial(link='logit'),
                   data=patdata)
    modsum<- summary(modelmut)
    candidatetable[pqr, 3] <- round(exp(modsum$coefficients[2,1]),6)
    candidatetable[pqr, 4] <- round(modsum$coefficients[2,"Pr(>|z|)"  ],15)
    candidatetable[pqr, 5] <- signif(exp(confint(modelmut)[2,1]), 6)
    candidatetable[pqr, 6] <- signif(exp(confint(modelmut)[2,2]), 6)
  }
  
  #save these candidate tables
  if(hlaABC=="A"){
    candidatetableA <- candidatetable
  }else if (hlaABC=="B"){
    candidatetableB<- candidatetable
  }else{
    candidatetableC<- candidatetable
  }
}

#Benjamini Hochberg Adjustment for all candidates 
candidatetable<- rbind(candidatetableA, candidatetableB, candidatetableC)
adjcandidates<- candidatetable[order(candidatetable$adj.p),]
adjcandidates$rank <- c(1:nrow(adjcandidates))
adjcandidates$imq <- (adjcandidates$rank/nrow(adjcandidates)) * 0.2  #False Discovery Rate= 0.2
adjcandidates$keep<-FALSE   #too keep track of which candidates to keep

#start from the bottom of ordered candidate pairs, once first pair where I/MQ is greater
#than adjusted P from log regression is reached, then that pair and all above are kept for
#Ben. Hochberg adjustment
firsttrue<-FALSE
for(abc in nrow(adjcandidates):1){
  if(adjcandidates[abc, "adj.p"] < adjcandidates[abc, "imq"]){
    firsttrue<-TRUE
  }
  if(firsttrue== TRUE){
    adjcandidates[abc, "keep"]<-TRUE
  }
}

keptcandidates<-adjcandidates[adjcandidates$keep==TRUE,]



####
####
####generates output table, each row shows cross-sectional analysis of each candidate pair
###
###
outputtable<-data.frame(matrix(vector(), nrow(keptcandidates), 13, 
                                  dimnames=list(c(), c("candidate", "hla", "hla.p", 
                                                       "years.infected", "years.infected.p", 
                                                       "interaction", "interaction.p",
                                                        "hla.CI.low", "hla.CI.up",
                                                  "years.infected.CI.low", "years.infected.CI.up",
                                                  "interaction.CI.low", "interaction.CI.up"))),  stringsAsFactors=F)

for(abc in 1:nrow(outputtable)){
  candidatepair<-keptcandidates[abc,1]
  hlaquery<- substr(keptcandidates[abc,1], (nchar(keptcandidates[abc,1])-1), nchar(keptcandidates[abc,1]))
  mutquery<-substr(keptcandidates[abc,1], 1, (nchar(keptcandidates[abc,1])-8))
  hlaABC<- substr(keptcandidates[abc,1], (nchar(keptcandidates[abc,1])-2), (nchar(keptcandidates[abc,1])-2))

  #keep track of whether or not years infected with HIV at least a year
  patdata$yearsinfectedbin<- patdata$yearsinfected>1
  
  #whether or not patient's HLA type is of queried type
  patdata$queriedhla<- 0
  for(mno in 1:nrow(patdata)){
    if(hlaABC== "A"){
      if(grepl(hlaquery, patdata[mno, "hlaa"])){
        patdata[mno, "queriedhla"]<- 1
      }
    }else if(hlaABC == "B"){
      if(grepl(hlaquery, patdata[mno, "hlab"])){
        patdata[mno, "queriedhla"]<- 1
      }
    }else{
      if(grepl(hlaquery, patdata[mno, "hlac"])){
        patdata[mno, "queriedhla"]<- 1
      }
    }
  }
  
  #sets queried mutation variable equal to DRM of interest
  patdata$queriedmut<- patdata[,mutquery]
  #variable set to TRUE or FALSE depending on if DRM is there
  patdata$queriedmut <-( patdata$queriedmut>0)
  
  
  modelmut <- glm(queriedmut~ queriedhla+yearsinfected + queriedhla*yearsinfected ,  family=binomial(link='logit'),
                  data=patdata)
  outputtable[abc, "candidate"] <- candidatepair
  outputtable[abc, "hla"]<- as.numeric(round(exp(modelmut$coefficients),3)[2])
  outputtable[abc, "years.infected"]<- as.numeric(round(exp(modelmut$coefficients),3)[3])
  outputtable[abc, "interaction"]<- as.numeric(round(exp(modelmut$coefficients),3)[4])
  modsum<-summary(modelmut)
  outputtable[abc, "hla.p"]<- as.numeric(round(modsum$coefficients[,"Pr(>|z|)"  ],5)[2])
  outputtable[abc, "years.infected.p"]<- as.numeric(round(modsum$coefficients[,"Pr(>|z|)"  ],5)[3])
  outputtable[abc, "interaction.p"]<- as.numeric(round(modsum$coefficients[,"Pr(>|z|)"  ],5)[4])
  outputtable[abc, "hla.CI.low"]<- as.numeric(round(exp(confint(modelmut)),3)[2,1])
  outputtable[abc, "hla.CI.up"]<- as.numeric(round(exp(confint(modelmut)),3)[2,2])
  outputtable[abc, "years.infected.CI.low"]<- as.numeric(round(exp(confint(modelmut)),3)[3,1])
  outputtable[abc, "years.infected.CI.up"]<- as.numeric(round(exp(confint(modelmut)),3)[3,2])
  outputtable[abc, "interaction.CI.low"]<- as.numeric(round(exp(confint(modelmut)),3)[4,1])
  outputtable[abc, "interaction.CI.up"]<- as.numeric(round(exp(confint(modelmut)),3)[4,2])
}




#############################################
###################
#### survival, developing of mutation as outcome
############
#############


for(abc in 1:nrow(keptcandidates)){
  candidatepair<-keptcandidates[abc,1]
  hlaquery<- substr(keptcandidates[abc,1], (nchar(keptcandidates[abc,1])-1), nchar(keptcandidates[abc,1]))
  mutquery<-substr(keptcandidates[abc,1], 1, (nchar(keptcandidates[abc,1])-8))
  hlaABC<- substr(keptcandidates[abc,1], (nchar(keptcandidates[abc,1])-2), (nchar(keptcandidates[abc,1])-2))
  

  #make copy from testinfo
  testinfo5<- testinfo
  
  #only keep sequences of patients of interest
  testinfo5<- testinfo5[testinfo5$id %in% patdata$id,]
  testinfo5$keep<-1
  
  #keep sequences that are artnaive
  for(abc in 1:nrow(testinfo5)){
    if(!is.na(patdata[patdata$id == testinfo5[abc,"id"], "artstart"])[1]){
      if(testinfo5[abc, "sampledate"] >patdata[patdata$id == testinfo5[abc,"id"], "artstart"][1]){
        testinfo5[abc, "keep"]<-0
      } 
    }
  }
  testinfo5<- testinfo5[testinfo5$keep==1,]
  
  #make table for survival analysis, only for patients with at least 2 ART naive sequences
  survivalpatients<- table(testinfo5$id)
  survivalpatients<- as.data.frame(survivalpatients)
  survivalpatients<-survivalpatients[survivalpatients$Freq>1,]
  #keeps track of whether patient has HLA of interest
  survivalpatients$hlab<- NA
  
  
  #see if there is event for each patient (i.e. no DRM -> has DRM)
  for(abc in 1:nrow(survivalpatients)){
    if(hlaABC=="A"){
      if(grepl(hlaquery,patdata[patdata$id== survivalpatients[abc, "Var1"], "hlaa" ])){
        survivalpatients[abc,"hlab"]<-1
      }else{
        survivalpatients[abc,"hlab"]<-0
      }
    }else if(hlaABC=="B"){
      if(grepl(hlaquery,patdata[patdata$id== survivalpatients[abc, "Var1"], "hlab" ])){
        survivalpatients[abc,"hlab"]<-1
      }else{
        survivalpatients[abc,"hlab"]<-0
      }
    }else{
      if(grepl(hlaquery,patdata[patdata$id== survivalpatients[abc, "Var1"], "hlac" ])){
        survivalpatients[abc,"hlab"]<-1
      }else{
        survivalpatients[abc,"hlab"]<-0
      }
    }
  }
  
  
  survivalpatients$keep<-1
  mutdata5<-mutdata
  mutdata5<- mutdata5[mutdata5$seqtype == substr(mutquery, 1,2),]
  mutdata5<- mutdata5[mutdata5$pos == substr(mutquery, 5,nchar(mutquery)),]
  
  for (abc in 1:nrow(survivalpatients)){
    tempmat<- testinfo5[testinfo5$id == survivalpatients[abc, "Var1"],]
    if(tempmat[1, "header_id"] %in% mutdata5$header_id){
      #throw out patients that already initially have DRM
      survivalpatients[abc, "keep"]<-0
    }
  }
  
  survivalpatients<- survivalpatients[survivalpatients$keep==1,]
  
  survivalpatients$failure<-0
  survivalpatients$timeatrisk<-NA
  
  #identifying patients with DRM-developing event, and time at risk
  for (abc in 1:nrow(survivalpatients)){
    tempmat<- testinfo5[testinfo5$id == survivalpatients[abc, "Var1"],]
    tempmat$mut<-0
    
    for(def in 1:nrow(tempmat)){
      if(tempmat[def, "header_id"] %in% mutdata5$header_id){
        tempmat[def, "mut"]<-1
      }
    }
    
    if(sum(tempmat$mut)==0){
      survivalpatients[abc, "timeatrisk"]<- (tempmat[nrow(tempmat), "sampledate"]- tempmat[1, "sampledate"])
    }else{
      datemut<- tempmat[tempmat$mut ==1, "sampledate"]
      
      if(length(datemut)>1){
        datemut<- datemut[1]
      }
      
      survivalpatients[abc, "timeatrisk"]<- datemut - tempmat[1, "sampledate"]
      survivalpatients[abc, "failure"]<-1
    }
  }
  
  
  
  
  res.cox <- coxph(Surv(timeatrisk, failure) ~ hlab, data = survivalpatients)
  print(candidatepair)
  print(res.cox)
  
  
  fit <- survfit(Surv(timeatrisk, failure) ~ hlab, data = survivalpatients)
  print(fit)
  print(ggsurvplot(fit, data = survivalpatients, pval = FALSE, size =0.8, palette= "lancet", font.x= c(15,"bold"),
             font.y= c(15,"bold"), font.main= c(15,"bold"), 
             x.lab= "Time (years)",xscale = "d_y", fun = "cumhaz", break.x.by = (5*365.25), ylim= c(0, 0.3)))
  
  print(exp(confint(res.cox)))
  
  print("mean time at risk")
  print(mean(survivalpatients$timeatrisk))
}
