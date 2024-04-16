#7 OCTOBER 2021
#FUNCTION TO CALCULATE MIC FOR MULTIPLE STRAINS
#Author: Madelaine Boggan


#use grow96 for growth analysis
#make all data into CSV so outside of R so don't need to rerun everysingle time

data <- AnalyseODdata
# need the means file from AnalyseODdata
#uses maxOD for data, not growth rate

#write.csv(data$means , [specify where to save and file name ])

#ICx is inhibitory concentration, leave at 90

#run everything from 15 to 58
MICcalculator <- function(data, ICx = 90) {
  maxConc <- max(data$RIFconc)
  #makes dataframe workable for the function
  maxODData <- select(data, -lag, -r2, -mumax)
  colnames(maxODData)[2] <- "antiBConc" #renames column from users input to add robustness

  #calculate MIC90 values from data
  MIC90 <- filter(maxODData, maxODData$antiBConc == 0) %>%  #finds all Strain data at antiBConc 0
    tibble(strainName = Strain, MIC90 = maxOD * (1 - ICx/100)) %>% #calculates MIC90
    select(strainName, MIC90) #saves MIC90

  #create empty dataframe for loops to fill with MIC data
  MIC_final <- data.frame(
    Strain = character(),
    MIC = double(),
    stringsAsFactors =  FALSE
  )

  #find all Strains and concentrations where maxOD is below 90%
  for (i in 1:nrow(maxODData)) {
    for (j in 1:nrow(MIC90)) {
      if (maxODData$Strain[i] %in% MIC90$strainName[j] && maxODData$antiBConc[i] == maxConc && maxODData$maxOD[i] > MIC90$MIC90[j]){ #finds all Strains with MIC>512
        MIC_final[i, 1] <- maxODData$Strain[i]
        MIC_final[i, 2] <- NA
      }
      if (maxODData$Strain[i] %in% MIC90$strainName[j] && maxODData$maxOD[i] < MIC90$MIC90[j]) { #appends Strains and antiBConcs where maxOD has dropped by min 90%, DOES NOT find MIC
        MIC_final[i, 1] <- maxODData$Strain[i]
        MIC_final[i, 2] <- maxODData$antiBConc[i]
      }
    }
  }

  #cleanup and output
  MIC_final <- MIC_final %>%
    filter(!is.na(Strain) & Strain != "none") %>% # removes NA and none values
    group_by(Strain) %>%
    summarise(MIC = min(MIC), .groups = "drop") #defines the MIC

return(MIC_final)

}

