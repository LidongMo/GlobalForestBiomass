#Biomass Calculation Tropical Tree Plots#

library(BIOMASS)
library(stringr)
library(ggplot2)
library(fields)
library(plyr)
library(gtools)
library(parallel)
library(dplyr)
library(RANN)

##########################################################################################
# STEP 1 
##########################################################################################

WoodDensityData = read.csv("/Volumes/CrowtherLabRAID/Lidong_Mo/BiomassEstimation/original_wd_of_all_spsecies_GFBI.csv") #the file path should be changed if you want to apply the code on other machine#
# this the normal approach for reltively small file size csv loading#
#for the loading of the extremely large csv data frame, the most effictive approach is the function
#from the package "data.table" and the function is "fread"#
#l
#this is the data frame of the raw binomials and the TNRS corrected genus and species name#
GFBIRawAndCorrectedName = read.csv("GFBI_raw_and_tnrs_corrected_binomial_df.csv")[,-1]
# read the ecoregion equation data frame
EcoregionEquation = read.csv("pesudo_biomass_based_ecoregion_equations.csv")

# read the original GFBI data
Original_GFBI = fread("Global_treelist_Biome_PLT_TPH.csv")
# remove the rows with NAs, specifical for the plots which located out side of the WWF_Biome mask
NAomit_GFBI = na.omit(Original_GFBI)
# remove the rows which has the WWF_Biome values 99 and 98
CleanedGFBI = subset(NAomit_GFBI,NAomit_GFBI$WWF_Biome!=98&NAomit_GFBI$WWF_Biome!=99)
# kick out those individuals whose diameter is less than 5
FiveDeleteGFBI = subset(CleanedGFBI,CleanedGFBI$DBH>=5)
# split the big data frame into samller pieces
SplitList = split(FiveDeleteGFBI, (as.numeric(rownames(FiveDeleteGFBI))-1) %/% 30000)
# for (i in names(SplitList))
# {
# 	write.csv(SplitList[[as.character(i)]], file = paste("GFBI_Data_",i,"_piece.csv", sep = ""))
# }
# in order to inprove the calcualtion speed, we split the big data frame into smalle pieces and name them
lapply(names(SplitList), function(x){write.csv(SplitList[[as.character(x)]], file = paste("Split_After_Cleaning/Split_GFBI_Data_",x,"_piece.csv", sep = ""))})

# get the file neme list in the split data frame folder
FileList                             <- list.files(path="Split_After_Cleaning",pattern="*.csv")
for (fl in FileList[134:945])
# for (fl in FileList)
{
	IndividualDataDBH                <- fread(paste("Split_After_Cleaning/",fl,sep=""))[,-1]
	# ListElement                      <- which(FileList==fl)
	BiomassCalculation               <- function(individual)
    {
    	#got the ecoregion code for each tree#
	    EcoregionCode                <- IndividualDataDBH[individual,]$WWF_Biome
	    #there is no allometic equations trained for 14 Mangroves and 3 Tropical and Subtropical Coniferous Forests#
	    if (EcoregionCode!=3&EcoregionCode!=14)
	    {
		    #get the equation information for each specific ecoregion code#
		    EquationInformation      <- EcoregionEquation[EcoregionEquation$ID_WWF==EcoregionCode,]
		    #get the individual information for each individual in the  IndividualDataDBH data frame#
		    IndividualInformation    <- IndividualDataDBH[individual,]
		    # get the dbh value for this individual
		    dbh                      <- IndividualInformation$DBH
		    # get the eqaution for the ecoregion where the individual is#
		    EquationFormula          <- tolower(as.character(EquationInformation$Eqaution))
            #Find the "=" symbol in the equation string,the charater after "=" will be the start of the equation#
            StartPosition            <- regexpr("=",EquationFormula)[1]+1
            #Find the end point in the equation string#
            EndPostion               <- nchar(EquationFormula)
            #Get the final equation for following calculations#
            EqautionForCalculation   <- str_sub(EquationFormula, start = StartPosition, end = EndPostion)
            #calculate out the individual biomass#
            IndividualBiomass        <- exp(eval(parse(text=EqautionForCalculation)))
            OutputDataFrame          <- data.frame(select(IndividualInformation,-Year),BiomassEquationBased=IndividualBiomass,BiomassPackageBased=NA)
            #for the tropical trees, we also apply the calculation approach from BIOMASS package#
            if(EcoregionCode!=1&EcoregionCode!=2&EcoregionCode!=7)
            {
        	    # for the individuals which come from the non-tropical region, we do not apply the computeAGB function to get the BIOMASS package based biomass, then set it to NA@
        	    OutputDataFrame$BiomassPackageBased <-NA
            }else
            {
        	    # for the individuals come from the tropical region, we will apply the computeAGB function from BIOMASS
        	    SpeciesBinomial      <- as.character(IndividualInformation$SPCD)
        	    # get the wood density row of the specific species
        	    WoodDensityRow       <- WoodDensityData[WoodDensityData$GFBIRawName==SpeciesBinomial,]
        	    # compute the biomass with the computeAGB equation from BIOMASS
        	    OutputDataFrame$BiomassPackageBased <-computeAGB(dbh,WoodDensityRow$meanWD,coord=cbind(IndividualInformation$LON,IndividualInformation$LAT))*1000
            }
        }else
        # (EcoregionCode==3|14)
        {
    	    # for the individuals come from the two ecoregions which do not hava any equations to get the biomass calculation, we set the biomass value to NA
    	    IndividualInformation    <- IndividualDataDBH[individual,]
    	    OutputDataFrame          <- data.frame(select(IndividualInformation,-Year),BiomassEquationBased=NA,BiomassPackageBased=NA)
        }

        print(paste("---the ",individual,"th line of ",nrow(IndividualDataDBH),",has been calculated---",sep=""))

        return(OutputDataFrame)
    }
    # apply parallel calculation of all the individuals in the GFBI database
    system.time(OutputLists          <- mclapply(1:nrow(IndividualDataDBH), BiomassCalculation, mc.cores = 33, mc.preschedule=FALSE))
    OutputDF                         <- do.call(rbind, OutputLists)
    write.csv(OutputDF,paste("Split_After_Cleaning_Biomass/individual_biomass_for_",fl,sep=""))
    print(paste("-----the ",fl,"has been calculated-----",sep=""))
}



##########################################################################################
# STEP 2
##########################################################################################

# transfer the data from the big computer to my local drive 
# scp  standardlogin@uwis-cx-dock-11-040.ethz.ch:/Users/Shared/Lidong_Mo/Split_Data_GFBI/individual_biomass_of_GFBI_data_for_ 6 _split.csv /Users/LeonidMoore/Desktop/BIOMASS/individual_biomass_of_GFBI_data_for_ 6 _split.csv

# this the work dirctory setting in big computer
setwd("/Volumes/CrowtherLabRAID/Lidong_Mo/BiomassEstimation/")
FileList           <- list.files(path="Split_After_Cleaning_Biomass",pattern="*.csv")
# use lapply to combine all the data into a singe object in R
ReunionDataList                 <- lapply(paste("Split_After_Cleaning_Biomass/",FileList,sep=""),fread)
ReunionData                     <- do.call(rbind,ReunionDataList)
# write the data cleaned biomass calculation result to the local folder
write.csv(ReunionData[,-1:-2],"/Volumes/CrowtherLabRAID/Lidong_Mo/BiomassEstimation/GFBI_Biomass_Calcualtion_Result_Data_Frame_0522.csv")

##########################################################################################
# STEP 3
##########################################################################################

# merge the data frame with TPH with the one just has the biomass value
# Load the table with TPH
DataBiomePLTTPH               <- fread("Global_treelist_Biome_PLT_TPH.csv")[,-1]
# load the biomass data
BiomassBiome                  <- fread("GFBI_Biomass_Calcualtion_Result_Data_Frame_0522.csv")[,-1]

# add new column to the data frame
DataBiomePLTTPH$New           <- paste(DataBiomePLTTPH$LAT,DataBiomePLTTPH$LON,DataBiomePLTTPH$SPCD,DataBiomePLTTPH$DBH,DataBiomePLTTPH$FID,sep="_")
BiomassBiome$New              <- paste(BiomassBiome$LAT,BiomassBiome$LON,BiomassBiome$SPCD,BiomassBiome$DBH,BiomassBiome$FID,sep="_")
# merge the two data frames by the "new" column
MergeData                     <- merge(DataBiomePLTTPH[,c("Year","TPH","New")],BiomassBiome,by="New")
# delete the column "New"
BiomassBiomeTPH               <- MergeData[,-"New"]

write.csv(BiomassBiomeTPH,"GFBI_Biomass_Calcualtion_Result_with_TPH_Data_Frame_20181008.csv")

# load the bimass data frame for plots only have 1 year record
SingleYearTreeLevelData        <- fread("GFBI_Biomass_Calcualtion_Result_with_TPH_Data_Frame_20181008.csv")
# do not do the NA kick out process
ReunionData                    <- SingleYearTreeLevelData
# kick out the data in biome 3 and 14
SubDataFrame                   <- ReunionData[ReunionData$WWF_Biome!=3&ReunionData$WWF_Biome!=14,]
# subset the individual biomass for the tropcial region
TropicalDataFrame              <- SubDataFrame[SubDataFrame$WWF_Biome==1|SubDataFrame$WWF_Biome==2|SubDataFrame$WWF_Biome==7,]
TropicalDataFrame              <- TropicalDataFrame[,c("FID","PLT","LAT","LON","Year","SPCD","DBH","TPH","WWF_Biome","BiomassPackageBased")]
TropicalDataFrame              <- na.omit(TropicalDataFrame)
names(TropicalDataFrame)       <- c("FID","PLT","LAT","LON","Year","SPCD","DBH","TPH","WWF_Biome","Biomass")
write.csv(TropicalDataFrame,"Tropical_Individual_Biomass_Data_Frame_20181012.csv") 

# subset the data frame for the non tropical region and allocate the name Biomass to the Biomasspackagebased column
NonTropicalDataFrame           <- SubDataFrame[SubDataFrame$WWF_Biome!=1&SubDataFrame$WWF_Biome!=2&SubDataFrame$WWF_Biome!=7,]
# select the column of biomass for the nontropical region individuals
NonTropicalDataFrame           <- NonTropicalDataFrame[,c("FID","PLT","LAT","LON","Year","SPCD","DBH","TPH","WWF_Biome","BiomassEquationBased")]
NonTropicalDataFrame           <- na.omit(NonTropicalDataFrame)
# change the name for the biomass column
names(NonTropicalDataFrame)    <- c("FID","PLT","LAT","LON","Year","SPCD","DBH","TPH","WWF_Biome","Biomass")
write.csv(NonTropicalDataFrame,"Non_Tropical_Individual_Biomass_Data_Frame_20181012.csv") 
# rbind the tropical and notropical region individual biomass data 
FinalIndividualBiomass         <- rbind(TropicalDataFrame,NonTropicalDataFrame)
# write the tropical and nontropical individual biomass data to the local directory
write.csv(FinalIndividualBiomass,"GFBI_Biomass_Calcualtion_Result_Reunion_Data_Frame_20181012.csv")



##########################################################################################
# STEP 4
##########################################################################################


library(data.table)
library(stringr)
library(parallel)
library(dplyr)
library(BIOMASS)
library(reshape2)
# because there are some duplicates in the PLT, which is not good for our calculation, We have add our own PLT to the data
# check the unique plots by the lat and lon
# genewrete a new column by the combination of lat and lon
FinalIndividualBiomass         <- fread("GFBI_Biomass_Calcualtion_Result_Reunion_Data_Frame_20181012.csv")[,-1]
FinalIndividualBiomass$new      <- paste(FinalIndividualBiomass$LAT,"_",FinalIndividualBiomass$LON,"_",FinalIndividualBiomass$PLT,sep="")
# find the out the unique lat and lon combination
UniqueCoordinates               <- unique(FinalIndividualBiomass$new)
# ADD NEW PLT INTO EACH COORDINATES PLT
NewPLTAllocating                <- function(uc)
# for (uc in UniqueCoordinates)
{
    # get the order ot the uc in the uniquecorrdinates
    UCOrder                     <- which(UniqueCoordinates==uc)
    # generate new plot id by the order of the uc and allocate this to the data frame
    NewPlotID                   <- paste("PLT_",str_pad(UCOrder, 8, pad = "0"),sep="")
    # subset the data in that coordinates plot
    PerCoordinate               <- FinalIndividualBiomass[FinalIndividualBiomass$new==uc,]
    PerCoordinateNew            <- data.frame(PerCoordinate,NewPLT=NewPlotID)
    PerCoordinateOutput         <- PerCoordinateNew[,c("FID","PLT","LAT","LON","Year","NewPLT","SPCD","DBH","TPH","WWF_Biome","Biomass")] 
    # print(UCOrder)
    # return(PerCoordinateOutput)
    write.csv(PerCoordinateOutput,paste("PerPLT/PLT_",NewPlotID,".csv",sep=""))
    print(UCOrder)
}

system.time(OutputList <- mclapply(UniqueCoordinates[58700:1188771],NewPLTAllocating,mc.cores=28,mc.preschedule=F))
# OutputDataFrame                   <- rbindlist(OutputList)
# write.csv(OutputDataFrame,"GFBI_Biomass_Tree_Level_Biomass_New_PLT_Data_Frame_20181022.csv")
# the noted codes are for the classic calculation

# rbind all the data frame in the PerPLT folder
# list all the data frame in the perplt folder
DataFrameList                    <- list.files(path="PerPLT",pattern=".csv")
# load all the data frame into the memory
ReunionDataList                 <- lapply(paste("PerPLT/",DataFrameList ,sep=""),read.csv)
ReunionData                     <- rbindlist(ReunionDataList)
# write the data cleaned biomass calculation result to the local folder
fwrite(ReunionData[,-1],"GFBI_Biomass_Tree_Level_Biomass_New_PLT_Data_Frame_20181024.csv")

##########################################################################################
# STEP 5
##########################################################################################


# HPC_Duplicates_PLT_LAT_LON.r
library(data.table)
library(parallel)
library(dplyr)
# setwd("/Volumes/CrowtherLabRAID/Lidong_Mo/BiomassEstimation/")
BiomassBiomeNewPLT                 <- fread("GFBI_Biomass_Tree_Level_Biomass_New_PLT_Data_Frame_20181024.csv")[,-1]
# add a new column which is the combination of lat and long
BiomassBiomeNewPLT$LATLON          <- paste(BiomassBiomeNewPLT$LAT,"_",BiomassBiomeNewPLT$LON,sep="")
# get the unique LATLON 
UniqueLATLON                       <- unique(BiomassBiomeNewPLT$LATLON)
SingleLatLONPLTFinding <- function(ll)
# for (ll in UniqueLATLON)
{
    # subse the data frame by the latlon combination
    SubsetDF                       <- BiomassBiomeNewPLT[BiomassBiomeNewPLT$LATLON==ll,]
    # check if there are more than 1 PLT, then just kick it out
    if (length(unique(SubsetDF$NewPLT))==1)
    {
        print(paste("--- the Plot ",unique(SubsetDF$NewPLT)," has been returned---",sep=""))
        return(SubsetDF)
    }
}
system.time(OutputList  <- mclapply(UniqueLATLON[58000:611627],SingleLatLONPLTFinding,mc.cores = 20, mc.preschedule=FALSE))
OutputDF                        <- rbindlist(OutputList)
write.csv(OutputDF,"GFBI_Plots_with_Single_PLT_LAT_LON_Records_20181024.csv")


##########################################################################################
# STEP 6 get the plots data with the last year
##########################################################################################

# for some time series plot, just keep the last year data there
library(data.table)
library(parallel)
library(dplyr)

# set the working directory by cmd 
# set the R version
# module load new gcc/4.8.2 r/3.5.1
# read the data with TPH
BiomassBiomeTPH                 <- fread("GFBI_Plots_with_Single_PLT_LAT_LON_Records_20181024.csv")[,-1]
# get the plot names
PlotNames                       <- as.vector(unique(BiomassBiomeTPH$NewPLT))
# in order to get the plot with time series, we just to check the plot with multiple years records
LastYearPlotFind            <- function(pn)
{
    # subset the dataframe for each plot
    PerplotDF                   <- BiomassBiomeTPH[BiomassBiomeTPH$NewPLT==pn,]
    # check how many years in the records of this plot
    # if the more than one year in the data frame, then return the last year of this plot,
    # get the years in the data frame
    UniqueYears                 <- unique(PerplotDF$Year)
    # get the max year in the plot
    SelectedYear                <- max(UniqueYears)
    # subset the per plot data frame by the max year
    LastYearRecords             <- PerplotDF[PerplotDF$Year==SelectedYear,]
    print(paste("---The ",pn," plot's last year data been returned! ---",sep=""))
    return(LastYearRecords)
}
system.time(OutputList  <- mclapply(PlotNames,LastYearPlotFind,mc.cores = 28, mc.preschedule=FALSE))
OutputDF                        <- rbindlist(OutputList)
write.csv(OutputDF,"GFBI_Plots_with_Single_Year_Records_single_lat_lon_20181024.csv")


##########################################################################################
# STEP 7 
##########################################################################################
# load the data frame negerated from the STEP_1
library(data.table)
library(parallel)
library(dplyr)
library(e1071)

# load the data frame
singleYearPlotData = fread("GFBI_Plots_with_Single_Year_Records_single_lat_lon_20181024.csv")[,-1]
# get the unique plots names 
plotNames = as.vector(unique(BiomassBiomeTPH$NewPLT))
# load the wood density data table
woodDensityData = read.csv("original_wd_of_all_spsecies_GFBI.csv") 


# work out a function to generate the metadata at plot level
plotLevelOperationFunc = function(pn=PlotName)
{
    
    # get the subset data for each plot by plot name
    perPlotDataFrame = singleYearPlotData %>% dplyr::filter(NewPLT == pn)
    # get the total tree number in per plot
    treeNumber = nrow(perPlotDataFrame)
    # plot size by hectare
    plotArea = 1/(perPlotDataFrame$TPH[1])
    # total biomass
    totalBiomass = sum(perPlotDataFrame$Biomass)
    # biomass density kg/ha
    biomassDensity = totalBiomass/plotArea
    # mean biomass per tree
    meanTreeBiomass = totalBiomass/ treeNumber
    # sd of biomass per tree
    sdBiomass = sd(perPlotDataFrame$Biomass)
    # sd of DBH per tree
    sdDBH = sd(perPlotDataFrame$DBH)
    # mean DBH in each plot
    meanDBH = mean(perPlotDataFrame$DBH)
    woodDenstiyAllocateFUnc = function(x)
    {
        # for the individuals come from the tropical region, we will apply the computeAGB function from BIOMASS
        speciesBinomial = as.character(x[5])
        # get the wood density row of the specific species
        woodDensityRow = woodDensityData[woodDensityData$GFBIRawName==speciesBinomial,]
        woodDensity = mean(woodDensityRow$meanWD)
        return(woodDensity)
    }
    # get the wood density data for each tree in the plot
    woodDensityVector = apply (perPlotDataFrame,1,woodDenstiyAllocateFUnc)
    # get the mean of wood density for each plot
    meanWoodDensity = mean(woodDensityVector)
    # get the sd of the wood density in each plot
    sdWoodDensity = sd(woodDensityVector)
    # get the lat and lon for that tph in that year
    Lati = unique(perPlotDataFrame$LAT)
    Long = unique(perPlotDataFrame$LON)
    yr = unique(perPlotDataFrame$Year)
    # return the information row out
    outputRow = data.frame(PLT=pn,
                          LAT=Lati,
                          YEAR=yr,
                          LON=Long,
                          treeNumber = treeNumber,
                          plotArea = plotArea,
                          totalBiomass = totalBiomass,
                          biomassDensity = biomassDensity,
                          meanTreeBiomass = meanTreeBiomass,
                          sdBiomass = sdBiomass,
                          sdDBH = sdDBH,
                          meanDBH = meanDBH,
                          meanWoodDensity = meanWoodDensity,
                          sdWoodDensity = sdWoodDensity)
    print(paste("--- the metadata for plot ",pn," has been calculated ---"))
    return(outputRow)   
}
# parallel running of this function
system.time(OutputList  <- mclapply(plotNames,plotLevelOperationFunc,mc.cores = 24, mc.preschedule=FALSE))
OutputDF                        <- rbindlist(OutputList)
write.csv(OutputDF,"GFBI_Plot_Level_Metadata_with_Single_Year_Iris_20190306.csv")

##########################################################################################
# STEP 8 
##########################################################################################

library(data.table)
library(parallel)
library(dplyr)
library(e1071)

# load the data frame
singleYearPlotData = fread("GFBI_Plots_with_Single_Year_Records_single_lat_lon_20181024.csv")[,-1]
# get the unique plots names 
plotNames = as.vector(unique(singleYearPlotData$NewPLT))
# load the wood density data table
woodDensityData = read.csv("original_wd_of_all_spsecies_GFBI.csv") 

# work out a function to generate the metadata at plot level
plotLevelRandomFunc = function(pn)
{
    
    # get the subset data for each plot by plot name
    perPlotDataFrame = singleYearPlotData %>% filter(NewPLT == pn)
    rawMeanCalc = plotLevelOperationFunc(perPlotDataFrame)
    processTable = data.frame()
    for (seed in c(1:100))
    {
    	# get the total tree number in per plot
    	set.seed(seed)
        perRandomTable = perPlotDataFrame %>%sample_n(nrow(perPlotDataFrame), replace = T) %>% plotLevelOperationFunc()
        processTable = rbind(processTable,perRandomTable)
    }
    # calculate the sd 
    outputTable = data.frame(PLT = unique(perPlotDataFrame$NewPLT)[1],
    	                     Lat = unique(perPlotDataFrame$LAT)[1],
                             Lon = unique(perPlotDataFrame$LON)[1],
                             MeanTreeBiomass = rawMeanCalc$meanTreeBiomass,
                             SDTreeBiomass = sd(processTable$meanTreeBiomass),
                             SkewnessBiomass = skewness(processTable$meanTreeBiomass),
                             MeanTreeBiomassDensity = rawMeanCalc$biomassDensity,
                             SDTreeBiomassDensity = sd(processTable$biomassDensity),
                             SkewnessBiomassDensity = skewness(processTable$biomassDensity))
    return(outputTable)  
}

plotLevelOperationFunc = function(inputTable)
{
    
    # get the subset data for each plot by plot name
    perPlotDataFrame = inputTable
    # get the total tree number in per plot
    treeNumber = nrow(perPlotDataFrame)
    # plot size by hectare
    plotArea = 1/(perPlotDataFrame$TPH[1])
    # total biomass
    totalBiomass = sum(perPlotDataFrame$Biomass)
    # biomass density kg/ha
    biomassDensity = totalBiomass/plotArea
    # mean biomass per tree
    meanTreeBiomass = totalBiomass/ treeNumber
    # sd of biomass per tree
    sdBiomass = sd(perPlotDataFrame$Biomass)
    # sd of DBH per tree
    # get the lat and lon for that tph in that year
    Lati = unique(perPlotDataFrame$LAT)
    Long = unique(perPlotDataFrame$LON)
    yr = unique(perPlotDataFrame$Year)
    # return the information row out
    outputRow = data.frame(PLT=pn,
                          LAT=Lati,
                          YEAR=yr,
                          LON=Long,
                          treeNumber = treeNumber,
                          plotArea = plotArea,
                          biomassDensity = biomassDensity,
                          meanTreeBiomass = meanTreeBiomass)
    print(paste("--- the metadata for plot ",unique(inputTable$NewPLT)," has been calculated ---"))
    return(outputRow)   
}
# parallel running of this function
system.time(OutputList  <- mclapply(plotNames,plotLevelRandomFunc,mc.cores = 28, mc.preschedule=FALSE))
OutputDF = rbindlist(OutputList)
# write.csv(OutputDF,"GFBI_Plot_Level_with_Single_Year_Bootstrapping_subsampled_Mean_SD_Iris_20211207.csv")
write.csv(OutputDF,"GFBI_Plot_Level_with_Single_Year_Bootstrapping_subsampled_Mean_SD_Iris_20211222.csv")