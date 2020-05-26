
rm(list=ls())
try(dev.off, silent=TRUE)

library(rworldmap)
library(raster)
library(sp)
library(rgdal)

dir.create("C:/Phytophthora_SDM")
setwd("C:/Phytophthora_SDM")

#setwd("C:/Users/Dan Chapman/Box/research/phytothreats/niche models")

# load modeling functions...
source("C:/Users/Dan Chapman/Box/research/phytothreats/niche models/model code/biomod_functions_Phytophthora.R")

allSpp = c("Phytophthora cactorum",
           "Phytophthora cambivora",
           "Phytophthora cinnamomi",
           "Phytophthora cryptogea",
           "Phytophthora gonapodyides",
           "Phytophthora lacustris",
           "Phytophthora plurivora",
           "Phytophthora ramorum",
           "Phytophthora x alni")

for(spp in allSpp[-(1:6)]){
  
  source("C:/Users/Dan Chapman/Box/research/phytothreats/niche models/model code/biomod_functions_Phytophthora.R")
  
  # set the species to model...
  message("\n********************************************************************************\n", 
          spp, "\n")
  sppLab = gsub(" ", ".", spp, fixed=TRUE) # spp label
  modelName = paste0(sppLab, ".v2") # label model version for trying different options
  imageName = paste(modelName, "/", modelName, ".RData", sep="") # RData file to save results
  #load(imageName)
  
  
  #################################################################################
  # read predictor rasters...
  
  rP = stack("C:/Users/Dan Chapman/Box/research/phytothreats/niche models/rasters/predictors.grd")
  names(rP)
  rP = subset(rP, c("bio10", "bio6", "bio15", "moisture", "forest", "urban_hort", "agric"))
  names(rP)
  
  # check distributions and log transform skewed ones
  #hist(rP)
  #rP$forest = log1p(rP$forest)
  #rP$urban_hort = log1p(rP$urban_hort)
  #rP$agric = log1p(rP$agric)
  #hist(rP)
  
  #########################################################################################
  # read in species records aligned to model grid...
  
  trainingOcc = read.csv(paste0("C:/Users/Dan Chapman/Box/research/phytothreats/niche models/rasters/train_", gsub(" ", "_", spp), ".csv")) # file name
  validationOcc = read.csv(paste0("C:/Users/Dan Chapman/Box/research/phytothreats/niche models/rasters/validation_", gsub(" ", "_", spp), ".csv")) # file name
  
  # convert to Spatial
  trainingOcc = SpatialPoints(trainingOcc, proj4string=rP@crs)
  validationOcc = SpatialPoints(validationOcc, proj4string=rP@crs)
  
  # extract predictor values at occurrences
  x_tr = data.frame(extract(rP, trainingOcc))
  x_val = data.frame(extract(rP, validationOcc))
  
  # remove any occurrences with no environmental covariate data
  trainingOcc = trainingOcc[!is.na(rowSums(x_tr))]
  validationOcc = validationOcc[!is.na(rowSums(x_val))]
  
  message(length(trainingOcc), " training occurrences")
  message(length(validationOcc), " validation occurrences")
  
  
  ############################################
  # read accessible region rasters (100 km radius)...
  
  rA = raster(paste0("C:/Users/Dan Chapman/Box/research/phytothreats/niche models/rasters/acc_", gsub(" ", "_", spp), "_100.tif"))
  
  #################################################################
  # Define unsuitable region
  
  x_tr = x_tr[!is.na(rowSums(x_tr)),]
  
  # compare values in occurrences with global values...
  #par(mfrow=c(3,3))
  #for(i in 1:nlayers(rP)){
  #  hist(rP[[i]], br=100, freq=FALSE, main=paste(spp, names(rP)[i]))
  #  lines(density(x[,i]), col="red", lwd=2)
  #}
  
  # compare values in occurrences with accessible values...
  #par(mfrow=c(3,3))
  #for(i in 1:nlayers(rP)){
  #  vals = getValues(rP[[i]])[ rA[] == 1 ]
  #  hist(vals, br=100, freq=FALSE, main=paste(spp, names(rP)[i]))
  #  lines(density(x[,i]), col="red", lwd=2)
  #}
  
  # back-transform any transformed
  #x$forest = exp(x$forest)-1
  #...
  
  # view quantiles...
  print(t(apply(x_tr, 2, quantile, c(0, 0.005, 0.01, 0.5, 0.99, 0.995, 1) )))
  print(t(apply(x_val, 2, quantile, c(0, 0.005, 0.01, 0.5, 0.99, 0.995, 1) )))
  
  print(table(x_val$bio10 > max(x_tr$bio10)))
  print(table(x_val$bio10 < min(x_tr$bio10)))
  print(table(x_val$moisture < min(x_tr$moisture)))
  print(table(x_val$bio15 > max(x_tr$bio15)))
  
  # create unsuitable raster
  rU = (rP$bio10 > max(x_tr$bio10)) + # summer too hot
    (rP$bio10 < min(x_tr$bio10)) + # summer too cold?
    (rP$moisture < min(x_tr$moisture)) + # too dry?
    #(rP$bio14 < min(x_tr$bio14)) + # too much drought stress?
    (rP$bio15 > max(x_tr$bio15)) # too seasonal precipitation?
  rU[rU[]>0] = 1
  try(dev.off, silent=TRUE)
  #plot(rU)
  
  table(extract(rU, trainingOcc))
  message(round(100*mean(extract(rU, trainingOcc)==1),1), "% of training records in unsuitable") # % in unsuitable
  
  table(extract(rU, validationOcc))
  message(round(100*mean(extract(rU, validationOcc)==1),1), "% of validation records in unsuitable") # % in unsuitable
  
  
  #####################################################
  # read recording effort
  rE = raster(paste0("C:/Users/Dan Chapman/Box/research/phytothreats/niche models/rasters/tg_", gsub(" ", "_", spp), ".tif"))
  
  numAccessible = sum(getValues(rE)[which(getValues(rA)==1)] > 0) # number of grid cells with records in accessible region
  message(numAccessible, " potential grid cells in accessible area")
  #plot(log10(rE))
  
  ##########################################################
  # start modelling...
  
  # make background samples...
  mySample = sample_background(
    occ=trainingOcc,
    trueAbsences=NULL, # no true absences to use
    rAccessible=rA, 
    rUnsuitable=rU,
    recEffort=rE, 
    rEnv = rP, # predictor layers
    nb.unsuitable=5000, # take 5000 unsuitable background samples (could increase)
    nb.accessible=5000, #5*length(trainingOcc), # same number of accessible samples as records
    nb.replicates=2) # change to 10 replicate samples for full model
  str(mySample)
  
  #save.image(imageName)
  #load(imageName)
  
  # set where to store temporary raster files...
  oldTempDir = rasterOptions()$tmpdir
  myTempDir = paste("c:/temp/", modelName, sep="")
  assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
  
  # model fitting...
  myBIOMOD = biomod_fitting(modelName=modelName, 
                            response=mySample, 
                            predictors=rP, 
                            EMreject="outlier",# reject outlying bad models from the ensemble
                            path_to_maxent="C:/" ) # directory where maxent files are saved
  myBIOMOD$evaluation_table # model evaluation
  myBIOMOD$importance_table # variable importance
  myBIOMOD$EM_algorithms # algorithms used for ensemble model
  save.image(imageName)
  #load(imageName)
  
  # save background samples to file...
  writeOGR(cbind(mySample$resp.var, mySample$PA.table), 
           dsn=modelName, layer=paste0("background_", modelName),
           driver="ESRI Shapefile",  overwrite_layer = TRUE)
  
  # save accessible and unsuitable rasters to file...
  writeRaster(rA, file=paste0(modelName, "/", modelName, "_accessible.tif"),
              format="GTiff", overwrite=TRUE)
  writeRaster(rU, file=paste0(modelName, "/", modelName, "_unsuitable.tif"),
              format="GTiff", overwrite=TRUE)
  #writeRaster(rU==1 & rA==0, file=paste0(modelName, "/", modelName, "_unsuitable_and_inaccessible.tif"),
  #            format="GTiff", overwrite=TRUE)
  
  # write evaluation stats to file...
  eval = data.frame(species=sppLab, algorithm=row.names(myBIOMOD$evaluation_table),
                    selected = ifelse(rownames(myBIOMOD$evaluation_table) %in% myBIOMOD$EM_algorithms, "yes", "no"),
                    myBIOMOD$evaluation_table, 
                    rbind(t(myBIOMOD$importance_table),
                          Ensemble=apply(myBIOMOD$importance_table[,myBIOMOD$EM_algorithm], 1, 
                                         weighted.mean, w=myBIOMOD$evaluation_table[myBIOMOD$EM_algorithm,"TSS"])))
  print(eval) 
  write.csv(eval, paste0(modelName, "/", modelName, "_evaluation_table.csv"), row.names=FALSE)
  
  # make response function plot...
  rp_df = plot_responses(M=myBIOMOD, modelName=modelName)
  save.image(imageName)
  
  # global projection in current climate...
  myProj_current = project_models(M=myBIOMOD, occPts=trainingOcc,
                                  modelName=modelName, response=mySample,
                                  predictors=rP, projectionName="current", 
                                  myCutOff=NULL, myBreaks=seq(0,1,0.01), 
                                  EMalgorithms=myBIOMOD$EM_algorithms)
  
  # plot nicer projection maps...
  rProj = raster(paste0(modelName, "/", modelName, 
                        "_ensemble_projection_current.tif"), band=2)
  rCV = raster(paste0(modelName, "/", modelName, 
                      "_ensemble_projection_current.tif"), band=4)
  plotGlobal(R=rProj, RCV=rCV, binary=FALSE, outputRes = 4*res(rP), 
             occPts=rbind(trainingOcc, validationOcc), myCutOff=myProj_current$cut_off,
             plotFile=paste(modelName, "/", modelName, "_projection_current_PRA", sep=""))
  plotEurope(R=rProj, binary=FALSE, occPts=rbind(trainingOcc, validationOcc), 
             myCutOff=myProj_current$cut_off,
             plotFile=paste(modelName, "/", modelName, "_projection_current_PRA_Europe", sep=""))
  plotUK(rProj, binary=FALSE, occPts=validationOcc, 
         myCutOff=myProj_current$cut_off,
         plotFile=paste(modelName, "/", modelName, "_projection_current_PRA_UK", sep=""))
  save.image(imageName)
  
  # map limiting factors...
  rLF = limiting_factor(M=myBIOMOD, occPts=rbind(trainingOcc, validationOcc),
                        modelName=modelName, 
                        predictors=rP)
  save.image(imageName)
  
  # make plots of the species distribution and the background regions...
  plot_distribution(modelName=modelName,  tr_occ=trainingOcc, val_occ=validationOcc,
                    closeup=TRUE, effort=rE)
  
  plot_background(modelName=modelName)
  
  # reset the temporary raster files and delete any created in the function
  assign("rasterOptions", rasterOptions(tmpdir=oldTempDir), envir = .GlobalEnv)
  unlink(myTempDir, recursive=TRUE)
  
}

# Done!

