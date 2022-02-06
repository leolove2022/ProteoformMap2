
#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
ProteoformMap_3 = function(Peptide_FC_P){
  Peptide_FC_P[,"Fold"]=0 # label all = 0 and add the FC>1.5 as 1 in next
  
  for(n in 1:nrow(Peptide_FC_P)){
    Peptide_FC_P[n,"Protein.Id"]# get protein.Id
    test<-subset(Peptide_FC_P,Protein.Id==Peptide_FC_P[n,"Protein.Id"])
    if((min(test$Plasma_FC)<1)&&(Peptide_FC_P[n,"Plasma_FC"]>1.5)&&(Peptide_FC_P[n,"Plasma_p"]<0.05)){# now small means <1
      #print("OK")
      Peptide_FC_P[n,"selected"]=1
      Peptide_FC_P[n,"Fold"]=1 ## new add to label the fold change for further N,C,M analysis
      #print(Peptide_FC_P[n,1])
      
    }
    # use to label the FC<1 as -1
    if((min(test$Plasma_FC)<1)&&(Peptide_FC_P[n,"Plasma_FC"]<1)){# now small means <1
      
      
      Peptide_FC_P[n,"Fold"]=-2 ## new add to label the fold change for further N,C,M analysis
      
      
    }
    # use to label the FC>1 but P>0.05 as -1
    if((min(test$Plasma_FC)<1)&&(Peptide_FC_P[n,"Plasma_FC"]>1)&&(Peptide_FC_P[n,"Plasma_p"]>0.05)){# now small means <1
      
      
      Peptide_FC_P[n,"Fold"]=-1 ## new add to label the fold change for further N,C,M analysis
      
      
    }
    #print(n)
  }
  # add all contain 1 row
  # the Sequence may change, 
  for(n in 1:nrow(Peptide_FC_P)){
    { 
      test<-subset(Peptide_FC_P,Protein.Id==Peptide_FC_P[n,"Protein.Id"])
      if((max(test$selected)==1)){# not know why have a ".1"
        
        Peptide_FC_P[n,"selected"]=1
        #print(Peptide_FC_P[n,"Protein.Id"])
      }
      
    }
  }
  
  for(n in 1:nrow(Peptide_FC_P)){
    Peptide_FC_P[n,"Protein.Id"]# get protein.Id
    test<-subset(Peptide_FC_P,Protein.Id==Peptide_FC_P[n,"Protein.Id"])
    if((min(test$Plasma_FC)<1)&&(Peptide_FC_P[n,"Plasma_FC"]>1.5)&&(Peptide_FC_P[n,"Plasma_p"]<0.05)){# now small means <1
      #print("OK")
      Peptide_FC_P[n,"selected"]=1
      #Peptide_FC_P[n,"Fold"]=0 ## new add to label the fold change for further N,C,M analysis
      #print(Peptide_FC_P[n,1])
      
    }
    #print(n)
  }
  # add all contain 1 row
  # the Sequence may change, 
  for(n in 1:nrow(Peptide_FC_P)){
    { 
      test<-subset(Peptide_FC_P,Protein.Id==Peptide_FC_P[n,"Protein.Id"])
      if((max(test$selected)==1)){# not know why have a ".1"
        
        Peptide_FC_P[n,"selected"]=1
        #print(Peptide_FC_P[n,"Protein.Id"])
      }
      
    }
  }
  
  # run above to label all contain select =1 proteins, the other peptide is alo =1
  #address_Threshold_peptide=("~/Documents/Experiment/EXP 61 NASH TBG TurboID/Peptide/Peptide_FC_P_complete_1_all_thresthod0.csv")
  #address_Threshold_peptide_2=("~/Documents/Experiment/EXP 61 NASH TBG TurboID/Peptide/Peptide_FC_P_completed_Selected_thresthod0.csv")
  #write.csv(Peptide_FC_P,file=address_Threshold_peptide,row.names=FALSE)# transform ok
  Peptide_FC_P_Selected=subset(Peptide_FC_P,selected==1) ############### This will use in deduplicate
  return(Peptide_FC_P_Selected)
}