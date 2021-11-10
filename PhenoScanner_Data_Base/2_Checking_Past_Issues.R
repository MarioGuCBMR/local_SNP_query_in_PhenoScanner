##############
#INTRODUCTION#
##############

#This code is to check what is going on with some SNPs and PhenoScanner.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(stringr)

###########
#Functions#
###########

PS_query <- function(chr_pos){
  
  #Getting the chr_pos df ready now
  
  chr_pos_df_original_copy <- as.data.frame(chr_pos)
  chr_pos_df_original <- as.data.frame(chr_pos)
  
  chr_pos_df_original_copy$rsid <- NA
  chr_pos_df_original$rsid <- NA
  
  #Let's get the data from our PS folder.
  
  path_ <- "C:/Users/zlc436/Downloads/dbSNP.txt/dbSNP"
  
  indexing <- seq(0, 298, by = 10)
  
  #We will need to parse this by partitioning the jobs:
  
  for(i in indexing){
    
    my_files <- c()
    
    for(j in seq(i, i+9)){
      
      print(j)
      
      my_files <- c(my_files, paste(path_, j, sep = "")) 
      
    }
    
    files_list <- as.list(my_files)
    l <- lapply(files_list, fread, sep=",")
    dt <- rbindlist( l )
    
    #Perfect!!
    #Now we check...
    
    dt_ <- dt[which(dt$hg19_coordinates%in%chr_pos),] 
    
    dt_ <- dt_[which(dt_$a1 != "-"),]
    dt_ <- dt_[which(dt_$a2 != "-"),]
    
    allele_vect <- c("A", "G", "C", "T")
    
    dt_ <- dt_[which(dt_$a1%in%allele_vect),]
    dt_ <- dt_[which(dt_$a2%in%allele_vect),]
    
    chr_pos_gotcha <- chr_pos[which(chr_pos%in%dt_$hg19_coordinates)]
    
    chr_pos_gotcha <- chr_pos_gotcha[order(match(chr_pos_gotcha, dt_$hg19_coordinates))]
    
    check <- which(chr_pos_gotcha != dt_$hg19_coordinates)
    
    if(is_empty(check) == FALSE){
      
        print("Loop number...")
        print(i)
        print("dind't work")
        
    }
    
    chr_pos_tmp <- chr_pos_df_original[which(chr_pos_df_original$chr_pos%in%chr_pos_gotcha),]
    
    chr_pos_tmp <- chr_pos_tmp[order(match(chr_pos_tmp$chr_pos, chr_pos_gotcha)),]
    
    chr_pos_tmp$rsid <- dt_$rsid
    
    chr_pos_df_old <- chr_pos_df_original[-which(chr_pos_df_original$chr_pos%in%chr_pos_gotcha),]
    
    chr_pos_df_original_ <- rbind(chr_pos_df_old, chr_pos_tmp)
       
    chr_pos_df_original <- chr_pos_df_original_
    
    my_files <- c()

    

}

other_parser <- function(snp){
  
  checkity <- as.character(unlist(strsplit(as.character(snp), "rs")))
  
  if(checkity[1] != ""){
    
      return(as.character(snp))
  }
  
} 

rsid_parser <- function(snp){
  
  checkity <- as.character(unlist(strsplit(as.character(snp), "rs")))
  
  if(checkity[1] == ""){
    
    if(length(checkity) > 2){
      
      print("There is something fishy here...")
      
      print(as.character(snp))
      
      print("Clean these weird rs SNPs before running EGoS and try again ;)")
      
    } else {
      
      return(as.character(snp))
      
    } 
    
  }
  
} 

rsid_parser_df <- function(df, snp_col, chr_col, pos_col){
  
  col_id <- which(colnames(df) == snp_col)
  
  rsid_df <- df[, ..col_id]
  
  rsid_ <- rsid_df$rsid
  
  rsid_index <- which(str_detect(rsid_, "rs") == TRUE)

  other_index <-which(str_detect(rsid_, "rs") == FALSE)
  
  #####################
  #Getting the cool df#
  #####################
  
  rsid_df <- df[rsid_index,]
  
  ###################################
  #Getting the df of the weirds SNPs#
  ###################################
  
  other_df <- df[other_index,]
  
  if(is_empty(other_df) == FALSE){
    
    print("There are weird SNPs...")
    print("We are going to check these bad boys...")
    
    other_df <- other_df[order(other_df$rsid),]
    
    head(other_df$rsid)
    tail(other_df$rsid)
    
    print("Hopefully you didn't get weird things...")
    print("To avoid really weird things we are going to make some bold moves")
    
    print("###REMOVE WEIRD ALLELES###")
    
    print("We are keeping only A, G, T and C"):
      
    effect_allele <- readline("Gimme you Effect Allele Column")
    other_allele <- readline("Gimme you Other Allele Column")
    
    ea_id <- which(colnames(other_df) == effect_allele)
    nea_id <- which(colnames(other_df) == other_allele)
    
    ea_df <- other_df[, ..ea_id]
    nea_df <- other_df[, ..nea_id]
    
    alleles_df <- cbind(ea_df, nea_df)
    
    colnames(alleles_df) <- c("ea", "nea")
    
    #And now we get the index:
    
    yes_vect <- c("A", "T", "C", "G")
    
    keep_index <- which(alleles_df$ea%in%yes_vect & alleles_df$nea%in%yes_vect)
    
    other_keep_df <- other_df[keep_index,]
    
    print("Let's see if the weirds unusable SNPs are gone...")
    
    head(other_keep_df)
    tail(other_keep_df)
    
    print("Let's remove those with EAF < 0.01 or EAF > 0.99")
    print("if you have MAF, we will only remove those with MAF <0.01")
    
    answer_freq <- readline("Do you wanna filter for MAF? [y/n]")
    
    if(answer_freq == "y"){
    
    freq <- readline("Gimme you freq column")

    freq_id <- which(colnames(other_keep_df) == freq)
    
    freq_df <- other_keep_df[, ..freq_id]
    
    colnames(freq_df) <- c("freq")
    
    summary(freq_df$freq)
    
    freq_index <- which(freq_df$freq > 0.01 & freq_df$freq < 0.99)
    
    other_keep_freq_df <- other_keep_df[freq_index,]
    
    #Perfect:
    
    head(other_keep_freq_df)
    tail(other_keep_freq_df)
    
    } else {
      
      other_keep_freq_df <- other_keep_df
      
    }
    
    print("Wanna check for only autosomical chromosomes?")
    
    answer <- readlines("y/n")
    
    if(answer == "y"){
      
      chr_id <- which(colnames(other_keep_freq_df) == chr_col)
      
      chr_df <- other_keep_freq_df[, ..chr_id]
      
      colnames(chr_df) <- c("chr")
      
      summary(chr_df$chr)
      
      yes_vect_chr <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
      
      chr_index <- which(chr_df$chr%in%yes_vect_chr)
      
      other_keep_chr_freq_df <- other_keep_freq_df[chr_index,]
      
      head(other_keep_chr_freq_df)
      tail(other_keep_chr_freq_df)
      
    } else {
      
      other_keep_chr_freq_df <- other_keep_freq_df
      
    }
    
    print("Now we ARE ready to check get the RSIDs of all those SNPs that lack them:")
    print("Getting ready chr:pos vector in PhenoScanner format...")
    
    chr_id <- which(colnames(other_keep_chr_freq_df) == chr_col)
    chr_df <- other_keep_chr_freq_df[, ..chr_id]
    colnames(chr_df) <- c("chr")
    
    pos_id <- which(colnames(other_keep_chr_freq_df) == pos_col)
    pos_df <- other_keep_chr_freq_df[, ..pos_id]
    colnames(pos_df) <- c("pos")
    
    head(chr_df)
    
    print("Do wee need to add a chr so they have the format of:")
    print("chr1")
    
    question_chr <- readline("y/n")
    
    if(question_chr == "y"){
      
      chr_df$chr <- paste("chr", chr_df$chr, sep = "")
      
    } 
    
    chr_pos <- paste(chr_df$chr, pos_df$pos, sep = ":")
    
    head(chr_pos)
    tail(chr_pos)
    
    print("Now we are going to check the SNPs in PhenoScanner")
    print("And check the RSID")
    
    
  }

  
}

EGoS_both <- function(df, snp_col, chr_col, pos_col){
  
  print("#########################################################")
  print("########STEP 2: Get the RSID and CHR and POS#############")
  print("#########################################################")
  
  print("First we are going to check what info do you have:")
  
  rsid_df <- rsid_parser(df, snp_col)
  
  
  
  
  
}

EGoS <- function(df){
  
  print("############################################")
  print("########STEP 1: Get the columns#############")
  print("############################################")
  
  #Let's check those SNPs that have RSID:
  
  print("EGoS is gonna make things easier for you." )
  print("But we need some help from your part." )
  print("We are going to ask you the SNP, Chr and Pos Column" )
  print("If you do not have one of those columns..." )
  print("Then just press enter" )
  print("These are your columns... ")
  
  print(colnames(df))
  
  empty <- readline("Press ENTER to continue")
  
  RSID <- readline("What is the column with RSIDs?")
  CHR <- readline("What is the column with CHR?")
  POS <- readline("What is the column for position?")
  
  empty <- readline("Press ENTER to continue")
  
  print("Great!! Now we have info on what can we do...")
  print("But EGoS is gonna behave differently depending on what YOU want to do.")
  
  empty <- readline("Press ENTER to continue")
  
  print("If you wanna generate a Polygenic Risk Score you need as many SNPs possible...")
  print("Hence, you are gonna need to do some checkity checks to your data to take...")
  print("ALL THE SNPs POSSIBLE")
  print("To do so we will work with chromosome and positions unless... you don't have any.")
  print("This option is really cool to match dataframs to avoid issues.")
  
  
  empty <- readline("Press ENTER to continue")
  
  print("If you wanna generate a genetic correlation or do mendelian randomization and you are certain of your data is all OK...")
  print("You need to work with RSIDs, unfortunately.")
  print("The reference panels are with RSIDs only, so we have to make do.")

  empty <- readline("Press ENTER to continue")
  
  print("SO, how do you want your data?")
  print("If you choose RSID, it will use the RSID column and will try to find RSIDs missing from the RSID column")
  print("If you choose CHR:POS, it will use the Chromosome and Position Columns and will do some checks with you.")
  print("If you choose CHR:POS, and you don't have Chr and Pos columns... it will try to find them in either build 37 or build 38.")
  print("If you choose BOTH, it will do all of the above depending on the data that you have")
  
  empty <- readline("Press ENTER to continue")
  
  Answer <- readline("Answer: RSID/CHR:POS/BOTH")
  
  if(Answer == "RSID"){
    
    #df_ <- EGoS_RSID(df)
    
    #return(df_)
    
  }
  
  if(Answer == "CHR:POS"){
    
    #df_ <- EGoS_CHR:POS(df)
    
    #return(df_)
    
  }
  
  if(Answer == "BOTH"){
    
    #df_ <- EGoS_BOTH(df, RSID, CHR, POS)
    
    #return(df_)
    
  }
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


##############
#Loading data#
##############

#We are going to load:

#One birth weight trait:

FBW <- fread("N:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/MRanalysis_IRandFI_withBW/Outcomes/Fetal_BW_European_meta.NG2019.txt")

#A Doherty trait.

Moderate <- fread("C:/Users/zlc436/Desktop/Leisure_Project/Doherty_Traits/Doherty-2018-NatureComms-moderate.csv")

#A Klimentidis trait:

Vigorous <- fread("C:/Users/zlc436/Desktop/Leisure_Project/Acc425_Model1_BOLTLMM_500K.txt")

#Smoking trait by Liu et al.

SmkInit <- fread("C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/SmokingInitiation.txt")

#Life time smoking trait by Wotton et al.

LifeTimeSmk <- fread("C:/Users/zlc436/Desktop/SMK_WHRAdjBMI/SmokingCombined_Wootton.txt")

#And Fasting Insulin from Lagou et al 2021.

FI <- fread("N:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/FI/FI_combined_1000G_density (1).txt.gz")

#######################################
#Let's get the data that are not rsids#
#######################################


