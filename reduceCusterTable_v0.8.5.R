### 2020/03/31 Manuel Kuhn | University of Maryland
reduceClusterTable <- function(file              = NA, 
                               specialROIs       = NA,
                               specialROILabels  = NA,
                               amyThresh         = 0.25,
                               wholeBrainThresh  = 0.3,
                               outAlternatives   = TRUE,
                               printSROIsToHead  = TRUE,
                               autoInstallPkgs   = TRUE,
                               macroTable        = TRUE){
  
  ## USAGE:
  # file              REQUIRED: path to clusterTable input file!        
  #                   Reads xlsx, csv, and txt as long as it is output from jsClusterTable
  # specialROIs       list of .nii Files of ROI's. Peaks within ROI will be relabeled (use Binary 0/1 only!)
  # specialROILabels  list of label names corresponding to specialROIs. Can be empty and ROI filenames will be used as labels
  # amyThresh         minimum threshold for amygdala probability
  # wholeBrainThresh  minimum threshold for probabilities for brain labels (except amygdala)
  # outAlternatives   output Alternative labels besided the winner-takes-all label
  # printSROIsToHead  print detected peaks of specialROIs to Table header
  # autoInstallPkgs   automatically install R packages if not available
  # macroTable        print as excel table with macro header for sorting
  
  
  #### USAGE: ####
  #Running the script is as easy as that (in R type for example):
  # source("C:/Users/kuhn/Box Sync/Projects/ACC/reduceCusterTable_v0.8.5.R")
  # reduceClusterTable("C:/Users/kuhn/Box Sync/Projects/ACC/AlexData/clusterTables_2mm/THR_Smooth6mm_dGMMask_PvsWAll_con3_spmT1_mFDR_clusterTable_2mm.txt")
  #If you want to use extra binary ROIs such as AAN_PAG (can be also any other .nii such as BST), 
  #you can specify the path to the respective .nii  images (and even give them shorter labels if you want, see bottom)
  # source("C:/Users/kuhn/Box Sync/Projects/ACC/reduceCusterTable_v0.8.5.R")
  # reduceClusterTable(file = "C:/Users/kuhn/Box Sync/Projects/ACC/AlexData/clusterTables_2mm/THR_Smooth6mm_dGMMask_PvsWAll_con3_spmT1_mFDR_clusterTable_2mm.txt", 
  #                   specialROIs  = c("C:/Users/kuhn/Box Sync/Projects/ACC/AAN_PAG_MNI152_1mm_v1p0_20150630.nii"))
  
  #### WHAT THIS SCRIPT DOES: ####
  # In each hemisphere, in each cluster for each peak:
  # Find highest label (and store as alternative in case of equal probabilities)
  # Look in current peak if same label is also winning label for other peaks w/ smaller T's --> remove those loosers (winner-takes-all)
  # CAVE: THIS MEANS: the one peak with highest T Value for one label is the winner irrespective of potential higher probability but smaller T!)
  # This also means, the entire peak is removed if winning label was already winner for other peak with the higher T. All other labels ignored!
  # --> In other words, within one cluster, the peak with the highest T wins per label 
  # --> i.e. aCG1 prob = 0.36 T = 5.24 will be favored over aCG2 prob = .95, T= 4.85 and acG2 will be deleted
  
  ### When using special masks
  # Any MNI coordinate included in special ROI Mask will survive as the probability is set to 999 (even if H0 is less then for example 0.25)
  # --> Remember: If Amygdala is part of Special ROI mask then HO probability will be overriden with 999!
  
  ############################# PREPARATIONS ############################
  ### LOAD AND IF NEEDED INSTALL REQUIRED PACKAGES
  # Packages Needed
  packages <- c("openxlsx","oro.nifti") #"rstudioapi"
  # Install packages if packages are not existing and autoInstall is flagged TRUE
  if (autoInstallPkgs){
    if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
      install.packages(setdiff(packages, rownames(installed.packages())))  
    }
  }
  # Load packages
  pkgsAvail <-lapply(packages, require, character.only = TRUE)
  # If loading of one of the packages does not work, suggest autoInstall = TRUE
  if (any(pkgsAvail == F)) stop("Missing packages, please use 'autoInstallPkgs = T' to install missing packages")  
  
  ### Checks
  # Check if input file is provided
  if (is.na(file)) {file = readline(prompt="File missing. Provide path to InputFile: ")}
  # If special ROIs are provided and also labels, # of ROIs and labels needs to correspond
  if ((length(specialROIs) != length(specialROILabels)) && (!is.na(specialROILabels))) stop("Lengths of ROIs and Labels are not the same!")
  # Initialize vector to collect possible warnings
  warnText <- character()
  
  ### Read in data
  # find input file 
  # Get filetype from file extension and use corresponding 'read' function
  # Top if no read function for filetyp is implemented here
  fileType <- strsplit(file,"\\.")[[1]][length(strsplit(file,"\\.")[[1]])]
  fileName <- strsplit(file, .Platform$file.sep)
  filePath <- paste0(fileName[[1]][1:(length(fileName[[1]])-1)], sep=.Platform$file.sep, collapse = "")
  fileName <- fileName[[1]][length(fileName[[1]])]
  if (fileType == "xlsx"){
    dat <- read.xlsx(file, startRow = 5, colNames = F )         
  }else if (fileType %in% c("csv", "txt")){
    dat =read.csv(file, header = F, skip = 4, sep = "\t", na.strings = "", stringsAsFactors = F)
  }else{
    stop(paste('Unsupported file type: ', fileType))
  }
  
  ### Get special-ROIs
  # Read in all specialROI NIFTIS and generate one MNI coordinate table per ROI
  if (all(!is.na(specialROIs))){
    sROIMNIs <- list(0)
    for (s in 1:length(specialROIs)){
      currSROI <-  specialROIs[s]
      ### TRY 
      img <- readNIfTI(currSROI)
      #CATCH HERE
      #  readline(prompt="Enter path to where the .nii lives: ")
      if (any(!(img_data(img) %in% c(0,1)))) stop('Use Binary Masks for Special ROIS only!')
      idx <-  which(img_data(img)==1, arr.ind = T)
      sROIMNIs[[s]] <- t(apply(idx, 1, function(x) translateCoordinate(x,nim = img)))
      if (!is.na(specialROILabels)){
        names(sROIMNIs)[s] <- specialROILabels[s]
      }else{
        names(sROIMNIs)[s] <- currSROI
      }
    }  
  }
  #######################################################################
  
  ######################### Correct HO Typos ############################
  dat[which(grepl("Precuneous", dat[,4])),4] <-  gsub("Precuneous","Precuneus", dat[which(grepl("Precuneous", dat[,4])),4])
  #######################################################################
  
  ######################### Start actual data reduction #################
  #Remove White Matter and Cereb Cortex, and ventricle
  origSizeDat <- dim(dat)[1]
  dat <- dat[!grepl('Cerebral White Matter', dat[,4]),]
  dat <- dat[!grepl('Cerebral Cortex', dat[,4]),]
  dat <- dat[!grepl('Ventrical', dat[,4]),]
  dat <- dat[!grepl('Ventricle', dat[,4]),]
  
  ### Store volume in mm3 for all clusters - will be appended later
  volInfo <- data.frame("Cluster" = gsub('Cluster ', "", dat[!is.na(dat[,1]),1]), "mm3" = dat[!is.na(dat[,1]),3] )
  volInfo$mm3 <- gsub('\\(mmm\\)' ,'' ,gsub('Volume ','', volInfo$mm3))
  
  ### Extract list of clusters
  # Identify Cluster by cluster number and generate a list with all Clusters
  clustIdx    <- which(!is.na(dat[,1]))
  if (any(clustIdx == 0)){stop('One cluster does not have surviving peaks')}
  clustList <- list()
  for (g in 1:length(clustIdx)) {
    if (g < length(clustIdx)) {
      clustList[g] <- list(dat[(clustIdx[g]+1):(clustIdx[g+1]-1),])
      names(clustList)[g] <- dat[clustIdx[g],1]
    }else{
      clustList[g] <- list(dat[(clustIdx[g]+1):dim(dat)[1],])
      names(clustList)[g] <- dat[clustIdx[g],1]
    }
  }
  
  ### Create PeakList per Cluster 
  if (exists('allLabelsDF')){rm(allLabelsDF)}
  for (p in 1:length(clustList)){
    ### Find peaks in cluster and create info list
    curDat <- clustList[[p]]
    peakIdx <- which(grepl('Peak', curDat[,2]))
    peakList <- list()
    peakTValsList <- list()
    peakCoordsList <- list()
    for (g in 1:length(peakIdx)){
      if ((diff(peakIdx)[g] != 1) & (g != length(peakIdx))){
        peakList[g] <- list(curDat[(peakIdx[g]+1):(peakIdx[g+1]-1),4])
        names(peakList)[g] <- curDat[peakIdx[g],2]
        peakCoordsList[g] <- gsub("MNI", "", curDat[(peakIdx[g]),3])
        peakTValsList[g] <- as.numeric(unlist(strsplit(curDat[(peakIdx[g]),4],"Value = ")[[1]][2]))
      }else{
        if (g == length(peakIdx)){
          if (peakIdx[g]+1 <= dim(curDat)[1]){
            peakList[g] <- list(curDat[(peakIdx[g]+1):dim(curDat)[1],4])
            names(peakList)[g] <- curDat[peakIdx[g],2]
          }
        }
        peakTValsList[g] <- as.numeric(unlist(strsplit(curDat[(peakIdx[g]),4],"Value = ")[[1]][2]))
        peakCoordsList[g] <- gsub("MNI", "", curDat[(peakIdx[g]),3])
      }
    }
    if (length(peakList) > 0){
      for (l in 1:length(peakList)){
        if (length(peakList[[l]]) > 0){
          curPeak <- peakList[[l]] 
          curSplit <- strsplit(curPeak,' = ')   
          for (y in 1:length(curSplit)) {
            curDF <- data.frame("Cluster" = as.numeric(gsub('Cluster ', "",names(clustList)[p])), "PeakNumber" = l,  
                                "MNI" = peakCoordsList[[l]], 
                                "Probability" = as.numeric(curSplit[[y]][2]), 
                                "Left/Right" =  ifelse(grepl(' -', strsplit(as.character(peakCoordsList[[l]]), ',')[[1]][1]), "L", "R"),
                                "Label" = curSplit[[y]][1], 
                                "T" = peakTValsList[[l]])
            if (exists("allLabelsDF")){
              allLabelsDF <- rbind(allLabelsDF, curDF)
            }else{
              allLabelsDF <- curDF
            }
          }
        }
      }
    }
  }
  ### Add Clusters Volume Info
  allLabelsDF <-merge(allLabelsDF, volInfo, all.x = T, by = "Cluster")
  ### Format MNI Labels
  allLabelsDF$MNIx <- as.numeric(gsub('\\[','',sapply(strsplit(as.character(allLabelsDF$MNI), ','),function(x) x[1])))
  allLabelsDF$MNIy <- as.numeric(sapply(strsplit(as.character(allLabelsDF$MNI), ','),function(x) x[2]))
  allLabelsDF$MNIz <- as.numeric(gsub('\\]','',sapply(strsplit(as.character(allLabelsDF$MNI), ','),function(x) x[3])))
  ### Make label names character!
  allLabelsDF$Label <- as.character(allLabelsDF$Label)
  
  
  ##### HERE YOU COULD PRINT OUT THE FULL TABLE
  ##### THIS IS ONE OF THE REASONS WHY IT'S DONE IN A TWO STEP PROCESS
  
  ################ CHECK and DEAL WITH SPECIAL ROIs  ####################
  if (all(!is.na(specialROIs))){
    tableCoords <-  paste(allLabelsDF$MNIx, allLabelsDF$MNIy, allLabelsDF$MNIz)
    storeSMNICoords <- list()
    for (sc in 1:length(specialROIs)){
      smniCoords <-  paste(sROIMNIs[[sc]][,1],sROIMNIs[[sc]][,2], sROIMNIs[[sc]][,3])
      smniIdx <- tableCoords %in% smniCoords
      # StoreAway all sROI coords and later print in Header!
      # Rename all labels for detected peaks within ROI to name of ROI
      if (any(smniIdx)){
        storeSMNICoords[[sc]] <- tableCoords[smniIdx]
        allLabelsDF$Label[smniIdx] <- names(sROIMNIs[sc])
        # Set probabilities to 999? As indicator that this is in ROI!
        allLabelsDF$Probability[smniIdx] <- 999
      }else{
        storeSMNICoords[sc] <- paste("No peaks detected within",names(sROIMNIs[sc]))
      }
      names(storeSMNICoords)[sc] <- names(sROIMNIs[sc])
    }
  }
  
  #######################################################################
  ### Remove based on thresholds
  # Remove amydgala < amyThresh
  # If amy is part of special ROI these peaks will survive as the probability is set to 999
  allLabelsDF <- allLabelsDF[!(grepl('mygdala', allLabelsDF$Label) & (allLabelsDF$Probability < amyThresh)),]
  # Remove rest of WB < wholeBrainThresh
  allLabelsDF <- allLabelsDF[!(!(grepl('mygdala', allLabelsDF$Label)) & (allLabelsDF$Probability < wholeBrainThresh)),]
  #######################################################################
  
  
  #################### FINALLY RUN SUPERDOOPER LOGIC!####################
  # seperate for R and L and do per cluster
  hemis <- c("R", "L")
  storeAmyCluster <- list()
  alternaLabels <-list()
  allLabelsDF$Label <- as.character(allLabelsDF$Label)
  redAllLabelsDF <-  allLabelsDF
  redAllLabelsDF <- redAllLabelsDF[with(redAllLabelsDF, order(Left.Right,Cluster, -T, -Probability)),]
  alternaMaxPeaks <- data.frame(0)
  amyCount <- 0
  rm(alternaMaxPeaks)
  for (h in hemis){
    hemiallLabelsDF <- redAllLabelsDF[redAllLabelsDF$Left.Right == h,]
    clusts <- unique(hemiallLabelsDF$Cluster)
    for (c in 1:length(clusts)){
      curClust <- hemiallLabelsDF[hemiallLabelsDF$Cluster == clusts[c],]
      
      # Store away AmyClusters
      amyClustIdx <- grepl('mygdala', curClust$Label)
      if (any(amyClustIdx)){
        amyCount <-  amyCount+1
        storeAmyCluster[[amyCount]]  <- curClust[amyClustIdx,]
      }
      
      origCurClust <- curClust
      peakNames <- unique(curClust$PeakNumber)
      curClustAllLabels <- unique(curClust$Label)
      for (p in 1:length(peakNames)){
        curPeak <- curClust[curClust$PeakNumber == peakNames[p],]
        #if (( peakNames[p] == 1) & clusts[c] == 19) stop()
        ### ADD REMOVE entire peak if firstPeak is removed
        if ((dim(curPeak[curPeak$Label != "removed",])[1] > 0) & (curPeak$Label[1] == "removed")){
          warnText[length(warnText)+1] <- paste0("Maybe a bad decision to exclude a peak. Please check Cluster=", clusts[c], 
                  ", Peak=", peakNames[p], ", x=", curPeak$MNIx[1],", y=",curPeak$MNIy[1] ,", z=",curPeak$MNIz[1])
          warning(warnText[length(warnText)], immediate. = T)
          curClust$Label[curClust$PeakNumber == peakNames[p]] <- "removed"
        }else{
          curPeak <- curPeak[curPeak$Label != "removed",]
          if (dim(curPeak)[1] != 0){
            highestPeakLabel <-  curPeak$Label[curPeak$Probability == max(curPeak$Probability)]
            if (length(highestPeakLabel) > 1){
              if (exists("alternaMaxPeaks")){
                alternaMaxPeaks <- rbind(alternaMaxPeaks, curPeak[curPeak$Label %in% highestPeakLabel[2:length(highestPeakLabel)],])
              }else{
                alternaMaxPeaks <- curPeak[curPeak$Label %in% highestPeakLabel[2:length(highestPeakLabel)],]
              }
              curClust$Label[(curClust$Label %in% highestPeakLabel[2:length(highestPeakLabel)]) & (curClust$Probability == max(curPeak$Probability))] <- "removed"
              highestPeakLabel <- highestPeakLabel[1]
            } 
            alternaLabels[[length(alternaLabels)+1]] <- origCurClust[(((origCurClust$PeakNumber == peakNames[p])) & (!(origCurClust$Label == highestPeakLabel))), ]
            # Remove this highest label from other peaks
            curClust$Label[((!(curClust$PeakNumber == peakNames[p])) & (curClust$Label == highestPeakLabel)) ] <- "removed"
            # Remove alternalables from current peak
            curClust$Label[(((curClust$PeakNumber == peakNames[p])) & (!(curClust$Label == highestPeakLabel)))] <- "removed"
            ### Add remove alternalabel
          }
        }
      }
      redAllLabelsDF[(redAllLabelsDF$Cluster == clusts[c]) & (redAllLabelsDF$Left.Right == h),] <- curClust
    }
  }
  
  redAllLabelsDF <- redAllLabelsDF[!(redAllLabelsDF$Label == "removed"),]
  
  ################## ADD ALTERNATIVES #####################
  if (outAlternatives){
    # Alternative maximum peak labels
    if (exists("alternaMaxPeaks")){
      redAllLabelsDF$AlternaMaxLabel <- 
        if (any(as.numeric(names(table(table(alternaMaxPeaks$Cluster, alternaMaxPeaks$PeakNumber)))) > 1)){
          warning('More than one AlternaMaxPeak. Implement exception',immediate. = T )}
      for (a in 1:dim(alternaMaxPeaks)[1])
        redAllLabelsDF$AlternaMaxLabel[(redAllLabelsDF$Cluster == alternaMaxPeaks$Cluster[a]) & (redAllLabelsDF$PeakNumber == alternaMaxPeaks$PeakNumber[a])] <- 
          alternaMaxPeaks$Label[a]
    }
    
    # Alternative labels for peaks
    sizeAlternaLabel <- numeric()
    for (a2 in 1:length(alternaLabels)[1]){
      sizeAlternaLabel[a2] <- dim(alternaLabels[[a2]])[1]
    }
    if (max(sizeAlternaLabel)>0){
      alternaMat <- data.frame(matrix(data=NA,nrow=dim(redAllLabelsDF)[1],ncol=max(sizeAlternaLabel)))
      names(alternaMat) <- paste("AlternaLabel_", 1:max(sizeAlternaLabel), sep = "")
      for (a2 in 1:length(alternaLabels)[1]){
        curAlternaLabel <- alternaLabels[[a2]]
        if (dim(curAlternaLabel)[1] > 0){
          for (a2a in 1:dim(curAlternaLabel)[1]){
            alClu <- curAlternaLabel$Cluster[a2a]
            alPea <- curAlternaLabel$PeakNumber[a2a]
            alLab <- curAlternaLabel$Label[a2a]
            alProb <- curAlternaLabel$Probability[a2a]
            alternaMat[which((redAllLabelsDF$Cluster == alClu) & (redAllLabelsDF$PeakNumber == alPea)),a2a] <- paste0(alLab, " (", alProb, ")")
          }
        }
      }
      redAllLabelsDF <- cbind(redAllLabelsDF, alternaMat)  
    }
  }
  
  
  #######################################################################
  #################### HEADER INFO GENERATION ###########################
  #######################################################################
  header <- paste('Reduced Cluster Table from: ', fileName)
  header <- rbind(header,paste("generated:", Sys.Date()))
  header <- rbind(header,paste("Reduced table from", origSizeDat, "to", dim(redAllLabelsDF)[1],"lines"))
  header <- rbind(header,paste("################################################ DON'T FORGERT!#######################################################"))
  header <- rbind(header,paste("THESE ARE  PEAK VOXEL COORDINATES ONLY. CHECK VISUALLY FOR non-PEAK ACTIVATION IN YOUR ROIS!"))
  header <- rbind(header,paste("jsClusterTable only reports positive peaks. to obtain the negative peaks, take care to run this code on the appropr. sign-flipped .nii file"))
  header <- rbind(header,paste("BST is NOT included in Harvard Oxford labels"))
  header <- rbind(header,paste("################################################## WARNINGS ##########################################################"))
  header <- rbind(header,paste("No Warranty for Correctness - Check your data manually!"))
  header <- rbind(header,paste("######################################################################################################################"))
  header <- rbind(header,paste("This code does not incorporate a minimum cluster extent."))
  header <- rbind(header,paste("In typical cases, the data analyst should strongly consider omitting clusters that are smaller than several native EPI voxels."))
  header <- rbind(header,paste("In some cases (e.g. spatially unsmoothed data; minimum conjunctions), special considerations may apply "))
  header <- rbind(header,paste("######################################################################################################################"))
  header <- rbind(header,paste("This code omits peaks that have a high probability of lying outside of GM."))
  header <- rbind(header,paste("It reports the most extreme peak for each combination of Cluster, Hemisphere, and HO Regional label. "))
  header <- rbind(header,paste("In cases where the extreme peak lies at the intersection of 2 Regions (i.e. the probabilities are similar), both are reported. "))
  header <- rbind(header,paste("######################################################################################################################"))
  header <- rbind(header,paste("For additional details, please see the header to the code."))
  header <- rbind(header,paste("######################################################################################################################"))
  header <- rbind(header,paste("If Peak Coordinates are within a special ROI, prob is set to 999. This will override any other (more probable) label"))
  header <- rbind(header,paste("######################################################################################################################"))
  
  header <- rbind(header,paste("INFOS:"))
  header <- rbind(header,paste("Amygdala was thresholded by:",amyThresh))
  header <- rbind(header,paste("Rest of WholeBrain was thresholded by:", wholeBrainThresh))
  
  fHeader <- redAllLabelsDF[1,]
  fHeader[1,] <- NA
  for (t in length(header):1){
    fHeader <- rbind(NA, fHeader )
    fHeader[1,1] <- header[t]
  }
  
  ### Amy Header
  if (length(storeAmyCluster) > 0){
    for(ac in 1:length(storeAmyCluster)){
      tmpsStoreAmyClust <- storeAmyCluster[[ac]]
      if (length(tmpsStoreAmyClust) > 0){
        if (ac == 1){
        amyHead <- redAllLabelsDF[1,]
        amyHead[1,] <- NA
        amyHead[1,1] <- paste("AMYGDALA PEAKS:")
        }
        tmp <- redAllLabelsDF[1:dim(tmpsStoreAmyClust)[1],]
        tmp[1:dim(tmpsStoreAmyClust)[1],] <- NA
        tmp[1:dim(tmpsStoreAmyClust)[1],1:dim(tmpsStoreAmyClust)[2]] <- tmpsStoreAmyClust
        amyHead <- rbind(amyHead,tmp)
      }
    }
  }else{
    amyHead <- redAllLabelsDF[1:2,]
    amyHead[1:2,] <- NA
    amyHead[1,1] <- paste("AMYGDALA HARVARD-OXFORD PEAKS:")
    amyHead[2,1] <- paste("No peaks detected within Harvard-Oxford Amygdala mask thresholded at:", amyThresh)
  }
  amyHead[dim(amyHead)[1]+1,] <- NA
  
  ### STORED SPECIAL ROI Header
  if (exists("storeSMNICoords")){
    roiHead <- redAllLabelsDF[1,]
    roiHead[1,] <- NA
    for (l in 1:length(storeSMNICoords)){
      curROI <- redAllLabelsDF[1,]
      curROI[1,] <- NA
      curROI[1,1] <- paste("PEAKS FOR SPECIAL ROI:",names(storeSMNICoords)[l])
      for (m in 1:length(storeSMNICoords[[l]])){
        cMNI <- redAllLabelsDF[1,]
        cMNI[1,]  <- NA
        cMNI[1,1] <- storeSMNICoords[[l]][m]
        curROI <- rbind(curROI, cMNI)
      }
      roiHead <-  rbind(roiHead, curROI)
    }
    tmp <- redAllLabelsDF[1,]
    tmp[1,] <- NA
    roiHead <- rbind(roiHead,tmp)
  }
  
  if (length(warnText) > 0){
    warnHead <- redAllLabelsDF[1,]
    warnHead[1,] <- NA
    warnHead[1,1] <- paste("Collected warnings while running:")
    for (w in 1:length(warnText)){
      tmp <- redAllLabelsDF[1,]
      tmp[1,] <- NA
      tmp[1,1] <-  warnText[w]
      warnHead <- rbind(warnHead,tmp)
    }
    tmp <- redAllLabelsDF[1,]
    tmp[1,] <- NA
    warnHead <- rbind(warnHead,tmp)
  }
  
  #### Bind all headers to dataframe
  fHeader <- rbind(fHeader, amyHead)
  if (exists("storeSMNICoords")){
    fHeader <- rbind(fHeader, roiHead)
  }
  #redAllLabelsDF <- rbind(fHeader,redAllLabelsDF)
  if (length(warnText) > 0){
    fHeader <- rbind(fHeader, warnHead)
  }
  
  ### Appand actual table
  redAllLabelsDF <- rbind(fHeader,redAllLabelsDF)
  
  ### Write file
  outName <- paste0(filePath, "reduced_", strsplit(fileName, "\\.")[[1]][1],".xlsx")
  write.xlsx(redAllLabelsDF,file = outName,asTable = macroTable)
}
