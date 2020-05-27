

thin_records = function(occs, templateGrid){
  require(raster)
  require(sp)
  
  # set where to store temporary raster files...
  oldTempDir = rasterOptions()$tmpdir
  myTempDir = paste(oldTempDir, "thin_records", sep="/")
  assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
  
  # convert records to presence on the template grid...
  rownames(occs) = 1:nrow(occs)
  R = rasterize(occs[,1:2], templateGrid, field=1, max)
  xy = coordinates(R)[!is.na(getValues(R)),] # presence grid cell coords
  
  # reset the temporary raster files and delete any created in the function
  assign("rasterOptions", rasterOptions(tmpdir=oldTempDir), envir = .GlobalEnv)
  unlink(myTempDir, recursive=TRUE)  
  
  # return SpatialPoints of the presence grid cell coords...
  SpatialPoints(xy, proj4string=CRS(projection(templateGrid)))
  
}

make_native_accessible = function(occ, bufferDist, rasterRes=0.25){
  
  require(raster)
  require(rworldmap)
  require(sp)
  require(rgeos)
  require(adehabitatHR)
  
  # MCP around native points
  nativePolygon = adehabitatHR::mcp(occ, percent = 100) # minimum convex polygon
  
  # find how much to extend the extent of the nativePolygon
  # maximum cell diagonal distance...
  cellD = max(apply(coordinates(occ), 1, function(xy){
    spDists(rbind(xy, xy+rasterRes), longlat=TRUE)[1,2]
  }))
  
  # template raster larger than polygon...
  myExt = round(extent(nativePolygon) + ceiling(rasterRes*4*bufferDist/cellD))
  rTemp = raster(ext=myExt, res=rasterRes)
  
  # rasterize MCP...
  rMCP = rasterize(nativePolygon, rTemp, field=1, background=NA)
  
  # find edge cells of MCP raster
  rEdge = boundaries(rMCP, type="inner", asNA=TRUE)
  
  # calculate minimum distance from every point in rtemp to edge of MCP...
  allXY = coordinates(rTemp)
  edgeXY = allXY[which(rEdge[] == 1),]
  d = apply(allXY, 1, function(xy){
    min(spDistsN1(pts=edgeXY, pt=xy, longlat=TRUE)) })
  
  # rasterise, and set interior to 0 distance...
  rDist = setValues(rMCP, d)
  rDist[!is.na(rMCP[])] = 0
  
  # convert to the fixed distance buffer...
  rBuffer = rDist <= bufferDist
  plot(rBuffer)
  plot(nativePolygon, add=TRUE)
  points(occ)
  
  # extend to the size of the world...
  rBuffer = extend(rBuffer, raster(res=rasterRes))
  rBuffer[is.na(rBuffer)] = 0 # replace NA with 0
  
  # return...
  rBuffer
  
}

make_invaded_accessible = function(occ=invadedOcc, bufferDist=40, rLand){
  
  require(raster)
  require(rworldmap)
  require(sp)
  require(rgeos)
  
  landXY = coordinates(rLand)[!is.na(rLand[]),]
  d = spDists(x=landXY, y=coordinates(occ), longlat=TRUE)
  minD = apply(d, 1, min)
  
  # rasterize distances...
  rD = setValues(rLand, NA)
  idx = cellFromXY(rD, landXY)
  rD[idx] = minD
  rBuffer = rD <= bufferDist
  zoom(rBuffer, extent(occ)*1.1, new=FALSE)
  plot(countriesCoarseLessIslands, add=TRUE)
  points(occ)
  
  rBuffer
}

sample_background = function(occ, trueAbsences=NULL,
                             rAccessible, recEffort=NULL, rUnsuitable=NULL,
                             rEnv,
                             nb.unsuitable=1000, 
                             nb.accessible=length(occ), nb.replicates=2){
  
  require(sp)
  require(rworldmap)
  require(raster)
  require(nnet)
  
  # set where to store temporary raster files...
  oldTempDir = rasterOptions()$tmpdir
  myTempDir = paste(oldTempDir, "sample_background", sep="/")
  assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
  
  # compile all the data for sampling...
  X = data.frame(coordinates(rAccessible), accessible=getValues(rAccessible))
  X$recEffort = if(!is.null(recEffort)) { getValues(recEffort) } else 1
  X$unsuitable = if(!is.null(rUnsuitable)) { getValues(rUnsuitable) } else 0
  X = cbind(X, getValues(rEnv)) # add covariates for masking
  X = X[!is.na(rowSums(X)),] # remove NAs
  
  # sample each replicate...
  backgroundSamples = lapply(1:nb.replicates, function(r){
    # sample from accessible area...
    accIdx = which(X$accessible==1)
    idx = sample(accIdx, size=nb.accessible, replace=TRUE, prob=X$recEffort[accIdx]) # sample with replacement!!!
    bType = rep("accessible", length(idx))
    # add samples from inaccessible but unsuitable area...
    if(!is.null(rUnsuitable) & nb.unsuitable>0){
      unsIdx = which(X$accessible==0 & X$unsuitable==1)
      idx = c(idx, sample(unsIdx, size=nb.unsuitable, replace=FALSE)) 
      bType = c(bType, rep("unsuitable", nb.unsuitable))
    }
    data.frame(rep=r, idx=idx, bType=bType) # rep(r,length(idx))
  })
  backgroundSamples = do.call(rbind, backgroundSamples)
  
  # convert to BIOMOD pseudo-absence table
  #PA.table = as.matrix(table(backgroundSamples$idx, backgroundSamples$rep)==1) # pseudo-absences per replicate
  PA.table = class.ind(backgroundSamples$rep) == 1
  PA.table = rbind(matrix(TRUE, ncol=nb.replicates, nrow=length(occ)), PA.table) # add presences to the top
  row.names(PA.table) = 1:nrow(PA.table)
  
  # make BIOMOD response as a SPDF...
  #pa.XY = X[sort(unique(backgroundSamples$idx)), 1:2] # grid cell coordinates of pseudo-absences
  pa.XY = X[backgroundSamples$idx, 1:2] # grid cell coordinates of pseudo-absences
  all.XY = rbind(coordinates(occ), pa.XY)
 # all.df = data.frame(occ=c(rep(1, length(occ)), rep(0, nrow(pa.XY))),
#                      ptType = c(rep("presence", length(occ)), 
#                                 ifelse(table(backgroundSamples$idx, backgroundSamples$bType)[,1]>0, "accessible", "unsuitable")) )
  all.df = data.frame(occ=c(rep(1, length(occ)), rep(0, nrow(pa.XY))),
                      ptType = c(rep("presence", length(occ)), backgroundSamples$bType) )
  resp.var = SpatialPointsDataFrame(coords=all.XY, data=all.df, proj4string=occ@proj4string)
  
  # reset the temporary raster files and delete any created in the function
  assign("rasterOptions", rasterOptions(tmpdir=oldTempDir), envir = .GlobalEnv)
  unlink(myTempDir, recursive=TRUE)  
  
  # return the sample
  list(resp.var=resp.var, PA.table=PA.table)
}

biomod_fitting = function(modelName, response, predictors, 
                          factorPredictors=NULL,
                          path_to_maxent, EMreject){
  # modelName = string for the name of the model.
  # occPts = SpatialPoints for the species presences
  # backPts = SpatialPoints for the background (aka pseudo-absences)
  # predictors = RasterStack for all the environmental predictors
  # fittingWts = case weights for model fitting
  
  require(biomod2)
  require(sp)
  require(raster)
  require(PresenceAbsence)
  
  # set where to store temporary raster files...
  oldTempDir = rasterOptions()$tmpdir
  myTempDir = paste(oldTempDir, "biomod_fitting", sep="/")
  assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
  
  # delete any old versions of the directory to save to...
  if(file.exists(modelName)) unlink(modelName, recursive=TRUE, force = TRUE)
  
  # set Biomod options to find MAXENT...
  BIOMOD_options = BIOMOD_ModelingOptions()
  BIOMOD_options@MAXENT.Phillips$path_to_maxent.jar = path_to_maxent # set to where maxent.jar file is saved
  BIOMOD_options@GAM$k = 4 # set GAM to 4 d.f.
  BIOMOD_options@GLM$test = "none" # no stepwise selection for GLM
  
  # make data.frame of explanatory variables...
  xVar = data.frame(extract(predictors, response$resp.var))
  if(!is.null(factorPredictors))
    for(v in factorPredictors) xVar[,v] = as.factor(xVar[,v])
  
  # select only the occurrences (1st column)
  response$resp.var = response$resp.var[,1]
  
  # format the data for input to the models
  BIOMOD_data = BIOMOD_FormatingData(
    resp.var=response$resp.var,
    expl.var = xVar,
    resp.name = modelName,
    PA.strategy = 'user.defined',
    PA.table = response$PA.table)
  
  # fit the models on each pseudo-absence sample...
  message("Running BIOMOD model fitting...")
  BIOMOD_model = BIOMOD_Modeling(
    data=BIOMOD_data, # data
    models = c('GLM','GAM','ANN', 'GBM', 'MARS', #'FDA',
               'RF','MAXENT.Phillips'), ## list of models to run
    models.options = BIOMOD_options, # Biomod options
    NbRunEval=1, # number of evaluation runs per psuedo-absence set (with cross-validation)
    DataSplit=80, # 80, # % of data used to calibrate ()
    VarImport=3, # number of permutations for estimating variable importance
    models.eval.meth = c('KAPPA', 'TSS', 'ROC'), # how to evaluate the models
    SaveObj = TRUE, # if TRUE, saves models to the hard drive
    rescal.all.models = TRUE, # if TRUE models are re-scaled with a GLM
    do.full.models=FALSE # if TRUE models are fitted to all the data (not desirable)
  )
  message("Done BIOMOD model fitting!")
  print(BIOMOD_model)
  
  # get evaluation statistics (AUC, KAPPA, TSS), averaged for each algorithm
  BIOMOD_eval = get_evaluations(BIOMOD_model)
  # average evaluation statistics on testing data across all model runs, for each algorithm 
  evaluation_table = apply(BIOMOD_eval[,"Testing.data",,,], c(2,1), mean)
  
  # for ensemble prediction, discard the worst performing algorithms
  if(is.numeric(EMreject)) # reject 'EMreject' worst algorithms...
    idx = rank(evaluation_table[,"ROC"], ties.method="random") > EMreject
  # citation: Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
  if(EMreject == "outlier"){# reject low outliers using MAD...
    AUC = evaluation_table[,"ROC"] # AUC of each algorithm
    MAD = median(abs(AUC - median(AUC))) # median absolute deviation (MAD)
    modifiedZ = 0.6745*(AUC-median(AUC))/MAD # modified z-scores
    idx = which(modifiedZ >= -1) # -3.5 is recommended for outliers, but -1 is stricter
  }
  if(EMreject == "none") # use everything
    idx = 1:nrow(evaluation_table)
  
  myAlgorithms = rownames(evaluation_table)[idx]
  myModels = get_built_models(BIOMOD_model)
  myModelsToUse = grep(paste(myAlgorithms, collapse="|"), myModels, value=T)
  
  # ensemble modelling
  # Average proportional to AUC
  message("Producing ensemble model prediction...")
  BIOMOD_EM = BIOMOD_EnsembleModeling(
    modeling.output = BIOMOD_model, # modelling output object
    chosen.models = myModelsToUse, # select from which models
    eval.metric = 'ROC', # evaluation metric
    eval.metric.quality.threshold = NULL, # minimum quality to accept
    prob.mean = FALSE, # Estimate the mean probabilities across predictions
    prob.cv = FALSE, # Estimate the coefficient of variation across predictions
    prob.ci = FALSE, # Estimate the confidence interval around the prob.mean
    prob.median = FALSE, # Estimate the median of probabilities
    committee.averaging = FALSE, # Estimate the committee averaging across predictions
    prob.mean.weight = TRUE, # Estimate the weighted sum of probabilities
    prob.mean.weight.decay = 'proportional', # weights are proportional to the evaluation scores
    em.by="PA_dataset+repet" # ensemble models are evaluated on the same part of the data as the individual models they are made with
  )
  message("Done producing ensemble model prediction!")
  
  # add evaluation stats for the ensemble model...
  BIOMOD_eval2 = get_evaluations(BIOMOD_EM)
  ensemble_eval = colMeans(do.call(rbind, lapply(BIOMOD_eval2, function(x){ x[,"Testing.data"]})))
  evaluation_table = rbind(evaluation_table, ensemble_eval)
  rownames(evaluation_table)[nrow(evaluation_table)] = "Ensemble" 
  
  # get average Variable Importance per variables and algorithm
  VI = get_variables_importance(BIOMOD_model)
  importance_table = apply(VI, c(2,1), mean) # average of each model across PA replicates and cross-validation splits
  # normalise to sum to 1 for each algorithm...
  importance_table = apply(importance_table, 1, function(x) { x / sum(x) })
  
  message("BIOMOD workflow completed successfully!")
  
  # return outputs...
  res = list("model_fitting" = BIOMOD_model, 
             "model_options" = BIOMOD_options,
             "ensemble_model" = BIOMOD_EM,
             "evaluation_table" = evaluation_table, 
             "importance_table" = importance_table,
             "EM_algorithms" = myAlgorithms,
             "modified_Z" = if(EMreject == "outlier") { modifiedZ } else NA)
  
  # save outputs to the directory
  save(res, file=paste(modelName, "/", modelName, "_BIOMOD_model.RData", sep=""))
  
  # reset the temporary raster files and delete any created in the function
  assign("rasterOptions", rasterOptions(tmpdir=oldTempDir), envir = .GlobalEnv)
  unlink(myTempDir, recursive=TRUE)  
  
  # return results...
  return(res)  
}

plot_responses = function(M, modelName){
  
  require(ggplot2)
  require(RColorBrewer)
  
  # set where to store temporary raster files...
  oldTempDir = rasterOptions()$tmpdir
  myTempDir = paste(oldTempDir, "biomod_fitting", sep="/")
  assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
  
  # load all the individuals models
  message("Loading the models...")
  myMod = BIOMOD_LoadModels(M$model_fitting) 
  for(i in 1:length(myMod)) # assign models to the global environment
    assign(myMod[i],  get(myMod[i]), envir=.GlobalEnv) # load models to the global environment
  
  # extract term names
  myTerms = get_formal_data(M$model_fitting, 'expl.var.names')
  
  # extract AUC for algorithms in ensemble
  EMalgorithms = M$EM_algorithms
  myW = M$evaluation_table[EMalgorithms, "ROC"] # weightings
  # myAUC = myAUC[names(myAUC)!="Ensemble"] # exclude the ensemble model
  
  # reject worst algorithms values
  #myAUC = myAUC[rank(myAUC, ties.method="random") > EMreject]
  #myAlgorithms = names(myAUC)
  myModSelection = grep(paste(EMalgorithms, collapse="|"), myMod, value=T)
  
  # sort terms in order of overall importance for the selected algorithms
  meanImportance = apply(M$importance_table[,EMalgorithms], 1, 
                         weighted.mean, w=myW)
  myTerms = myTerms[order(meanImportance, decreasing=T)]
  meanImportance = meanImportance[myTerms]
  
  # create response plot object
  plotDat = lapply(myTerms, function(term){
    message("Creating plotting data for ", term)
    myRespPlot = response.plot2(models=myModSelection,
                                Data = get_formal_data(M$model_fitting, 'expl.var'), 
                                show.variables= term,
                                do.bivariate = FALSE,
                                fixed.var.metric = 'median', # fixes other variables at their median
                                plot=FALSE,
                                data_species = get_formal_data(M$model_fitting,'resp.var'))
    myRespPlot = myRespPlot[[1]]
    xv = myRespPlot[,1] # value of the predictor
    myRespPlot = myRespPlot[,-1] # remove predictor
    
    # mean responses per algorithm...
    myRespPlot = cbind(myRespPlot, myRespPlot) # make sure more than 1 representation (needed for below)
    algorithms = unlist(lapply(strsplit(names(myRespPlot), "_"), tail, 1))
    yv = data.frame(do.call(cbind, tapply(1:ncol(myRespPlot), algorithms, function(i){
      rowMeans(myRespPlot[,i]) })))
    res = cbind(xv, yv)
    colnames(res)[1] = "value"
    res$variable = term
    return(res)
  })
  
  ensemblePlotDat = data.frame(do.call(rbind, lapply(plotDat, function(x){
    xv = as.numeric(x[,1])
    yv = x[,2:(ncol(x)-1)]  
    YV = as.numeric(apply(yv, 1, weighted.mean, w=myW[colnames(yv)]))
    #print(length(xv))
    #print(length(YV))
    data.frame(variable=rep(x$variable[1],length(xv)), x=xv, y=YV)
  })))
  
  # convert plotDat to data.frame for ggplot
  plotDat = do.call(rbind, plotDat) # combine into 1 data.frame
  plotDat$variable = factor(plotDat$variable, levels=unique(plotDat$variable))
  plotDat2 = melt(plotDat, id.vars=c("variable", "value"))
  names(plotDat2) = c("variable", "xvalue", "Model", "yvalue")
  
  plotdim = if(length(myTerms)<= 3) c(1,1) else n2mfrow(length(myTerms)) # number of rows and cols in facet plot
  
  g = ggplot(data=plotDat2, aes(x=xvalue, y=yvalue, colour=Model)) +
    facet_wrap(~variable, scales="free_x",  ncol=plotdim[2], nrow = plotdim[1]) +
    geom_line() +
    geom_line(data=ensemblePlotDat, aes(x=x, y=y), inherit.aes=FALSE, size=1) +
    ylim(0,1) +
    ylab("Partial suitability") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          strip.text = element_text(size=12),
          axis.title.y = element_text(size=12))
  pdf(paste(modelName, "/", modelName, "_response_plots.pdf", sep=""),
      width=7, height=plotdim[1]*(7/2.5))
  print(g)
  dev.off()
  
  rm(list=myMod, envir=.GlobalEnv) # remove models from the environment
  
  # reset the temporary raster files and delete any created in the function
  assign("rasterOptions", rasterOptions(tmpdir=oldTempDir), envir = .GlobalEnv)
  unlink(myTempDir, recursive=TRUE)  
  
  #return(list(algDat=plotDat, ensembleDat=ensemblePlotDat))
  return(ensemblePlotDat)
  
}

project_models = function(M, occPts,
                          modelName, response,
                          predictors=WC, projectionName="current", 
                          myCutOff=NULL, myBreaks=NULL, EMalgorithms){
  
  require(RColorBrewer)
  require(biomod2)
  require(raster)
  require(rworldmap)
  require(PresenceAbsence)
  data("countriesCoarse")
  data("countriesLow")
  rasterOptions(progress="")   
  
  # set where to store temporary raster files...
  oldTempDir = rasterOptions()$tmpdir
  myTempDir = paste(oldTempDir, "project_models", sep="/")
  assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
  
  # discard the worst performing algorithms...
  myModels = get_built_models(M$model_fitting)
  myModelsToUse = grep(paste(EMalgorithms, collapse="|"), 
                       myModels, value=T)
  myAUC = M$evaluation_table[EMalgorithms, "ROC"]
  
  # project each individual retained model...
  BIOMOD_proj = BIOMOD_Projection(
    modeling.output = M$model_fitting,
    new.env = predictors, # rasterStack of predictors
    proj.name = projectionName, # name of the projection
    selected.models = myModelsToUse,
    binary.meth = NULL, # 'ROC', NULL = don't transform to binary
    compress = 'xz',
    build.clamping.mask = FALSE
  )
  
  # re-do ensemble assembly...
  BIOMOD_EM = BIOMOD_EnsembleModeling(
    modeling.output = M$model_fitting, # modelling output object
    chosen.models = myModelsToUse, # select from which models
    eval.metric = 'ROC', # evaluation metric
    eval.metric.quality.threshold = NULL, # minimum quality to accept
    prob.mean = FALSE, # Estimate the mean probabilities across predictions
    prob.cv = TRUE, # Estimate the coefficient of variation across predictions
    prob.ci = FALSE, # Estimate the confidence interval around the prob.mean
    prob.median = FALSE, # Estimate the median of probabilities
    committee.averaging = FALSE, # Estimate the committee averaging across predictions
    prob.mean.weight = TRUE, # Estimate the weighted sum of probabilities
    prob.mean.weight.decay = 'proportional', # weights are proportional to the evaluation scores
    em.by="PA_dataset+repet" # ensemble models are evaluated on the same part of the data as the individual models they are made with
  )
  
  # produce an ensemble projection across all models...
  BIOMOD_ensemble_proj = BIOMOD_EnsembleForecasting(
    EM.output=BIOMOD_EM, projection.output = BIOMOD_proj)
  
  # extract the projection as a RasterStack on each PA dataset
  rProj = BIOMOD_ensemble_proj@proj@val
  
  # average of the ensemble projection (converted to "probabilities")
  rEnsembleProj = mean(rProj[[grep("EMwmeanByROC", names(rProj))]]) / 1000
  
  # coefficient of variation
  rEnsembleCV = mean(rProj[[grep("EMcvByROC", names(rProj))]]) / 1000
  
  # 95% CI
  #rEnsembleCI = mean(rProj[[grep("EMciInfByROC", names(rProj))]]) / 1000
  
  # create a 'clamping mask'...
  message("Creating clamping mask...")
  rClamp = setValues(predictors[[1]], 1)
  myLimits = apply(get_formal_data(M$model_fitting, "expl.var"), 2, range)
  for(i in 1:nlayers(predictors)){
    message("  Clamping ", names(predictors)[i], "...")
    x = getValues(predictors[[i]])
    xx = as.numeric(x>=myLimits[1,i] & x<=myLimits[2,i])
    R = setValues(rClamp, xx)
    rClamp = rClamp * R
    #plot(rClamp, main=i, col=c("red", "grey95"))	
  }
  rClamp[rClamp[] == 0] = NA
  
  # mask the ensemble...
  rEnsembleProj_clamped = mask(rEnsembleProj, rClamp)
  
  message("Plotting projection maps...")
  
  # make colour ramp for continuous projection, if not supplied...
  if(is.null(myBreaks)){
    myBreaks = unique(quantile(rEnsembleProj_clamped, seq(0,1,0.01), na.rm=T))
    #myBreaks[length(myBreaks)] = myBreaks[length(myBreaks)] + 
    #	abs(myBreakls[length(myBreaks)]*0.001)
    myBreaks[1] = 0
    myBreaks[length(myBreaks)] = 1
  }
  
  # find optimal threshold for binary projection if not supplied...
  # cut-off done with min ROC distance
  if(is.null(myCutOff)){
    
    X = data.frame(type=response$resp.var$ptType,
                   obs=response$resp.var$occ,
                   predAll=extract(rEnsembleProj_clamped, response$resp.var))
    th = suppressWarnings(optimal.thresholds(X, opt.methods=1:12, req.sens=0.9, req.spec=0.9))
    write.csv(th, paste0(modelName, "/", modelName, "_threshold_table.csv"))
    message("Potential thresholds for binary conversion...")
    print(th)
    myCutOff =  th[9,2] # minRocDist seems OK?
    message("Chosen threshold value = ", round(myCutOff, 4))
  } else {
    message("User-defined threshold value = ", round(myCutOff, 4))
  }
  
  # plot global projection...
  plotFile = paste(modelName, "/", modelName, "_projection_", 
                   projectionName, sep="")
  map = spTransform(countriesCoarse, rEnsembleProj_clamped@crs)
  #pdf(paste(plotFile, ".png", sep=""), h=4)
  png(paste(plotFile, ".png", sep=""), width=7, height=4, units = 'in', res = 600)
  polys = list("sp.lines", as(map, "SpatialLines"), lwd=0.5)
  points = list("sp.points", occPts, pch=24, cex=0.5, col="black", fill="white", lwd=0.5)
  print(spplot(rEnsembleProj_clamped, 
               maxpixels=ncell(rEnsembleProj_clamped),
               scales = list(draw=FALSE),
               sp.layout=list(polys),
               at=myBreaks,
               colorkey=list(space="bottom") ))
  dev.off()
  
  # plot a binary projection...
  #pdf(paste(plotFile, "_binary.pdf", sep=""), h=4)
  png(paste(plotFile, "_binary.png", sep=""), width=7, height=4, units = 'in', res = 600)
  print(spplot(rEnsembleProj_clamped, maxpixels=ncell(rEnsembleProj_clamped),
               scales = list(draw=FALSE),
               sp.layout=list(polys, points),
               at=c(0,myCutOff,1), col.regions=c("grey95", "tomato"),
               colorkey=list(space="bottom") ))
  dev.off()
  
  # plot Europe projection...
  #pdf(paste(plotFile, "_Europe.pdf", sep=""))
  png(paste(plotFile, "_Europe.png", sep=""), width=7, height=7, units = 'in', res = 600)
  polys = list("sp.lines", as(map, "SpatialLines"), lwd=0.5)
  print(spplot(rEnsembleProj_clamped, maxpixels=ncell(rEnsembleProj_clamped),
               scales = list(draw=FALSE),
               sp.layout=list(polys),
               at=myBreaks,
               colorkey=list(space="bottom"),
               xlim=c(-1635172, 3435823), ylim=c(4064677, 8146210)))
  dev.off()
  
  #pdf(plotFile, "_Europe_binary.pdf", sep=""))
  png(paste(plotFile, "_Europe_binary.png", sep=""), width=7, height=7, units = 'in', res = 600)
  
  print(spplot(rEnsembleProj_clamped, maxpixels=ncell(rEnsembleProj_clamped),
               scales = list(draw=FALSE),
               sp.layout=list(polys, points),
               at=c(0,myCutOff,1), col.regions=c("grey95", "tomato"),
               colorkey=list(space="bottom"),
               xlim=c(-1635172, 3435823), ylim=c(4064677, 8146210) ))
  dev.off()
  
  # write to file...
  rOut = stack(rEnsembleProj, rEnsembleProj_clamped, rClamp, rEnsembleCV)
  names(rOut) = c("ensembleProj", "ensembleProj_masked", "clamp", "ensemble_CV")
  outFile = paste(modelName, "/", modelName, "_ensemble_projection_", projectionName, sep="")
  writeRaster(x=rOut, filename=outFile, 
              overwrite=TRUE, format="GTiff", bylayer=FALSE)
  
  # assess accuracy in accessible area...
  aucDat = data.frame(ptType=response$resp.var$ptType,
                      occ=response$resp.var$occ, 
                      pv=extract(rEnsembleProj, response$resp.var))
  aucDat = aucDat[aucDat$ptType!="unsuitable",] # remove unsuitable
  evalAccessible = presence.absence.accuracy(aucDat, threshold = myCutOff)
  
  # reset the temporary raster files and delete any created in the function
  assign("rasterOptions", rasterOptions(tmpdir=oldTempDir), envir = .GlobalEnv)
  unlink(myTempDir, recursive=TRUE)  
  
  message("All done!")
  
  # return the projection objects...
  return(list(projection = BIOMOD_proj, 
              ensemble_projection = BIOMOD_ensemble_proj,
              cut_off = myCutOff, 
              breaks=myBreaks, evalAccessible=evalAccessible))
}

plotGlobal = function(R, RCV, plotFile="test", binary=FALSE, myCutOff=0.5,
                      outputRes, occPts=occGridCells){
  # plot global projection...
  
  require(raster)
  require(rworldmap)
  require(grid)
  require(gridExtra)
  data(countriesCoarse)
  
  map = spTransform(countriesLow, R@crs)
  
  R2 = aggregate(R, outputRes/res(R), max, na.rm=TRUE)
  RSD = aggregate(RCV*R, outputRes/res(RCV), mean, na.rm=TRUE)
  
  # convert to binary if needed...
  if(binary) R2 = R2 > myCutOff
  
  myBreaks = if(!binary){ seq(0,1,0.01) } else c(0,myCutOff,1)
  myColours = if(!binary){
    #bpy.colors(length(myBreaks)-1)
    #colorRampPalette(c("grey75", "magenta", "green"))(length(myBreaks)-1)
    col1 = colorRampPalette(rev(brewer.pal("Blues", n=9)))(myCutOff*100)
    col2 = colorRampPalette(brewer.pal("OrRd", n=9)[-(1:2)])(100 - myCutOff*100)
    c(col1, col2)
  } else c("grey80", "red")
  
  
  polys = list("sp.lines", as(map, "SpatialLines"), lwd=0.5, col="grey50")
  points = list("sp.points", occPts, pch=24, cex=0.4, col="black", fill="white", lwd=0.5)
  plot1 = spplot(R2, maxpixels=ncell(R2),
                 scales = list(draw=FALSE),
                 sp.layout=list(polys, points),
                 at=myBreaks, col.regions=myColours )
  plot2 = spplot(RSD, maxpixels=ncell(RSD),
                 scales = list(draw=FALSE),
                 sp.layout=list(polys),
                 col.regions=colorRampPalette(brewer.pal("OrRd", n=9))(100) )
  
  pdf(paste(plotFile, ".pdf", sep=""), width=7, height=6)
  #png(paste(plotFile, ".png", sep=""), width=7, height=6, units = 'in', res = 600)
  grid.arrange(plot1, plot2, ncol=1, nrow=2)
  grid.text("Suitable", x=unit(0.925, "npc"), y=unit(0.985, "npc"), 
            gp=gpar(font=2, col=myColours[length(myColours)]))
  grid.text("Unsuitable", x=unit(0.925, "npc"), y=unit(0.515, "npc"),
            gp=gpar(font=2, col=myColours[1]))
  grid.text("(a) Projected suitability", x=unit(0.05, "npc"), y=unit(0.985, "npc"), 
            gp=gpar(font=1), just = "left")
  grid.text("(b) Standard deviation in projected suitability", x=unit(0.05, "npc"), y=unit(0.5, "npc"), 
            gp=gpar(font=1), just = "left")	
  dev.off()
  
}

plotEurope = function(R, plotFile="test1", binary=FALSE, myCutOff=0.5,
                      occPts=rbind(trainingOcc,validationOcc)){
  # plot global projection...
  
  require(raster)
  require(rworldmap)
  require(lattice)
  require(grid)
  data(countriesLow)
  
  map = spTransform(countriesLow, R@crs)
  
  extEU = c(-1635172, 3435823, 4064677, 8146210)
  
  R2 = crop(R, extEU)		
  #w = focalWeight(x=R, d=0.1, type="Gauss")
  #R3 = focal(R2, w, na.rm=TRUE)
  #R3 = mask(R3, R2)
  R3 = R2
  
  # convert to binary if needed...
  if(binary) R3 = R3 > myCutOff
  
  myBreaks = if(!binary){ seq(0,1,0.01) } else c(0,myCutOff,1)
  myColours = if(!binary){
    #bpy.colors(length(myBreaks)-1)
    #colorRampPalette(c("grey75", "magenta", "green"))(length(myBreaks)-1)
    col1 = colorRampPalette(rev(brewer.pal("Blues", n=9)))(myCutOff*100)
    col2 = colorRampPalette(brewer.pal("OrRd", n=9)[-(1:2)])(100 - myCutOff*100)
    c(col1, col2)
  } else c("grey80", "red")
  
  pdf(paste(plotFile, ".pdf", sep=""), width=6, height=4.5)
  #png(paste(plotFile, ".png", sep=""), width=6, height=4.5, units = 'in', res = 600)
  
  polys = list("sp.lines", as(map, "SpatialLines"), lwd=0.5, col="grey50")
  points = list("sp.points", occPts, pch=24, cex=0.4, col="black", fill="white", lwd=0.5)
  
  print(spplot(R3, maxpixels=ncell(R3),
               scales = list(draw=FALSE),
               sp.layout=list(polys, points),
               at=myBreaks, col.regions=myColours ), position = c(0,0,0.95,1))
  grid.text("Suitable", x=unit(0.975, "npc"), y=unit(0.75, "npc"), rot=-90, 
            gp=gpar(font=2, col=myColours[length(myColours)]))
  grid.text("Unsuitable", x=unit(0.975, "npc"), y=unit(0.25, "npc"), rot=-90, 
            gp=gpar(font=2, col=myColours[1]))
  
  #  	print(spplot(R3, maxpixels=ncell(R2),
  #         	scales = list(draw=FALSE),
  #         	sp.layout=list(polys, points),
  #         	at=myBreaks, col.regions=myColours,
  #          colorkey = list(labels = list(at=c(0.05, 0.05+(myCutOff-0.05)/2, myCutOff, myCutOff+(0.95-myCutOff)/2, 0.95 ), # (0.05,0.95,l=5), 
  #                                       labels = c("Highly\nunsuitable",
  #                                                  "Unsuitable",
  #                                                  "Marginal",
  #                                                  "Suitable", "Highly\nsuitable"),        
  #                                       width = 1, cex = 1), space="right"))) 
  dev.off()
  
}

plotUK = function(R, plotFile="test1", binary=FALSE, myCutOff=0.5,
                  occPts=validationOcc){
  # plot global projection...
  
  require(raster)
  require(rworldmap)
  require(lattice)
  require(grid)
  data(countriesLow)
  
  countries = c("GBR", "IRL", "FRA", "BEL", "NLD")
  map = spTransform(countriesLow[countriesLow$ISO3 %in% countries,], R@crs)
  
  extUK = c(-903715.3, 308154.3, 5738330, 7056504)
  
  R2 = crop(R, extUK)		
  #w = focalWeight(x=R, d=0.1, type="Gauss")
  #R3 = focal(R2, w, na.rm=TRUE)
  #R3 = mask(R3, R2)
  R3 = R2
  
  # convert to binary if needed...
  if(binary) R3 = R3 > myCutOff
  
  myBreaks = if(!binary){ seq(0,1,0.01) } else c(0,myCutOff,1)
  myColours = if(!binary){
    #bpy.colors(length(myBreaks)-1)
    #colorRampPalette(c("grey75", "magenta", "green"))(length(myBreaks)-1)
    col1 = colorRampPalette(rev(brewer.pal("Blues", n=9)))(myCutOff*100)
    col2 = colorRampPalette(brewer.pal("OrRd", n=9)[-(1:2)])(100 - myCutOff*100)
    c(col1, col2)
  } else c("grey80", "red")
  
  pdf(paste(plotFile, ".pdf", sep=""), width=5, height=4.5)
  #png(paste(plotFile, ".png", sep=""), width=6, height=4.5, units = 'in', res = 600)
  
  polys = list("sp.lines", as(map, "SpatialLines"), lwd=0.5, col="grey50")
  points = list("sp.points", occPts, pch=24, cex=0.4, col="black", fill="white", lwd=0.5)
  
  print(spplot(R3, maxpixels=ncell(R3),
               scales = list(draw=FALSE),
               sp.layout=list(polys, points),
               at=myBreaks, col.regions=myColours ), position = c(0,0,0.95,1))
  grid.text("Suitable", x=unit(0.975, "npc"), y=unit(0.75, "npc"), rot=-90, 
            gp=gpar(font=2, col=myColours[length(myColours)]))
  grid.text("Unsuitable", x=unit(0.975, "npc"), y=unit(0.25, "npc"), rot=-90, 
            gp=gpar(font=2, col=myColours[1]))
  
  dev.off()
  
}


limiting_factor = function(M=myBIOMOD, occPts,
                           modelName, 
                           predictors,  
                           myProj = raster(paste(modelName, "/", modelName, 
                                                 "_ensemble_projection_current.tif", sep=""), band=1),
                           plotFile = paste0(modelName, "/", modelName, "_limiting_factor")){
  require(RColorBrewer)
  require(biomod2)
  require(raster)
  require(rworldmap)
  data("countriesCoarse")
  data("countriesLow")
  rasterOptions(progress="") 
  
  # set where to store temporary raster files...
  oldTempDir = rasterOptions()$tmpdir
  myTempDir = paste(oldTempDir, "project_models", sep="/")
  assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
  on.exit(unlink(myTempDir, recursive=TRUE))
  
  # occurrence values of environmental variables...
  x = extract(predictors, occPts)
  occVals = apply(x, 2, median, na.rm=TRUE)
  
  # discard the worst performing algorithms...
  EMalgorithms=myBIOMOD$EM_algorithms
  myModels = get_built_models(M$model_fitting)
  myModelsToUse = grep(paste(EMalgorithms, collapse="|"), myModels, value=T)
  myAUC = M$evaluation_table[EMalgorithms, "ROC"]
  
  # project for Europe with each variable fixed to its median occurrence value
  partialProj = lapply(names(predictors), function(V){ # for each variable... 
    
    myTempDir = paste(oldTempDir, "project_models", sep="/")
    assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
    
    extEU = c(-1635172, 3435823, 4064677, 8146210)
    
    message(V)
    predictors2 = if(any(grepl("MAXENT.Tsuruoka", myModelsToUse))){ 
      predictors 
    } else crop(predictors, extEU) # crop to Europe (unless using MEMLR)
    predictors2[[V]] = setValues(predictors2[[V]], occVals[V])
    
    predictors2 = stack(predictors2)
    
    # project each individual retained model...
    BIOMOD_proj = BIOMOD_Projection(
      modeling.output = M$model_fitting,
      new.env = predictors2, # rasterStack of predictors
      proj.name = V, # name of the projection
      selected.models = myModelsToUse,
      binary.meth = NULL, # 'ROC', NULL = don't transform to binary
      compress = 'xz',
      build.clamping.mask = FALSE
    )
    
    # ensemble assembly...
    BIOMOD_EM = BIOMOD_EnsembleModeling(
      modeling.output = M$model_fitting, # modelling output object
      chosen.models = myModelsToUse, # select from which models
      eval.metric = 'ROC', # evaluation metric
      eval.metric.quality.threshold = NULL, # minimum quality to accept
      prob.mean = FALSE, # Estimate the mean probabilities across predictions
      prob.cv = FALSE, # Estimate the coefficient of variation across predictions
      prob.ci = FALSE, # Estimate the confidence interval around the prob.mean
      prob.median = FALSE, # Estimate the median of probabilities
      committee.averaging = FALSE, # Estimate the committee averaging across predictions
      prob.mean.weight = TRUE, # Estimate the weighted sum of probabilities
      prob.mean.weight.decay = 'proportional', # weights are proportional to the evaluation scores
      em.by="PA_dataset+repet" # ensemble models are evaluated on the same part of the data as the individual models they are made with
    )
    
    # produce an ensemble projection across all models...
    BIOMOD_ensemble_proj = BIOMOD_EnsembleForecasting(
      EM.output = BIOMOD_EM,
      projection.output = BIOMOD_proj)
    
    # extract the projection as a RasterStack on each PA dataset
    rProj = BIOMOD_ensemble_proj@proj@val
    if(any(grepl("MAXENT.Tsuruoka", myModelsToUse)))
      rProj = crop(rProj, extEU) # crop to Europe if needed
    
    # average of the ensemble projection (converted to "probabilities")
    mean(rProj) / 1000
  })
  
  partialProj = stack(partialProj)
  names(partialProj) = names(predictors)
  plot(partialProj)
  
  # calculate change in suitability...#
  extEU = c(-1635172, 3435823, 4064677, 8146210)
  baseProj = crop(myProj, extEU)
  suitabilityChange = partialProj - baseProj
  names(suitabilityChange) = names(predictors)
  plot(suitabilityChange)
  
  # find which is the limiting factor (highest increase in suitability)...
  x = getValues(suitabilityChange)
  bestVar = apply(x, 1, which.max)
  bestVar[which(sapply(bestVar, length) == 0)] = NA
  limitingFactor = setValues(baseProj, as.numeric(unlist(bestVar)))
  
  writeRaster(limitingFactor, file=paste0(modelName, "/", modelName, "_limiting.tif"),
              format="GTiff", overwrite=TRUE)
  
  pdf(paste(plotFile, ".pdf", sep=""), width=6, height=4.5)
  #png(paste(plotFile, ".png", sep=""), width=6, height=4.5, units = 'in', res = 600)
  
  map = spTransform(countriesLow, predictors@crs)
  polys = list("sp.lines", as(map, "SpatialLines"), lwd=0.5, col="grey50")
  
  #myColours = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  myColours = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                brewer.pal("Set1", n=8), brewer.pal("Dark2", n=8))
  
  print(spplot(limitingFactor, maxpixels=ncell(limitingFactor),
               sp.layout=list(polys),
               col.regions=myColours[1:length(names(predictors))],
               at=0.5+0:length(names(predictors)),
               colorkey = list(
                 labels=list(
                   at = 1:length(names(predictors)),
                   labels = names(predictors)
                 )
               )
  ))
  
  dev.off()
  
  # reset the temporary raster files and delete any created in the function
  assign("rasterOptions", rasterOptions(tmpdir=oldTempDir), envir = .GlobalEnv)
  unlink(myTempDir, recursive=TRUE) 
  
  # delete the projections...
  sapply(paste0(modelName, "/proj_", names(predictors)), unlink, recursive=TRUE, force = TRUE)
  
  return(limitingFactor)
  
}

plot_distribution = function(modelName, tr_occ=trainingOcc, val_occ=validationOcc,
                             closeup=TRUE, effort=rE){
  
  map = spTransform(countriesCoarseLessIslands, tr_occ@proj4string)
  
  # plot the distribution...
  png(filename=paste0(modelName, "/", sppLab, ".distribution.png"),
      units = 'in', width = 5, height = 4.75, 
      res = 400)
  par(mfrow=c(2,1), mar=c(0.6,0.6,1.6,0.6))
  if(!closeup)
    plot(map, col="grey90", border="grey75", ylim=c(-7210547, 9182709))
  if(closeup){
    e = 1.2*extent(rbind(tr_occ,val_occ))
    plot(map, col="grey90", border="grey75", 
         ylim=c(e@ymin, e@ymax), xlim=c(e@xmin, e@xmax))
  }  
  title(main="(a) Species distribution used in modelling",
        cex.main=1, font.main=1, adj=0)
  
  points(tr_occ, pch=21, bg="tomato", col="red4", cex=0.7)
  points(val_occ, pch=21, bg="white", col="blue4", cex=0.7)
  legend("topright", pch=21, bg="white",
         pt.bg=c("tomato","white"),
         col=c("red4","blue4"),
         legend=c("Training", "Validation"), cex=0.7)
  
  plot(map, border=NA, ylim=c(-7210547, 9182709))
  effort = log10(aggregate(effort, 20, na.rm=TRUE, fun=sum))
  plot(effort, axes=FALSE, bty="n", add=TRUE, 
       legend=FALSE, col=colorRampPalette(brewer.pal("OrRd", n=9)[-(1:2)])(100))
  plot(effort, axes=FALSE, bty="n", add=TRUE, 
       legend.only=TRUE, col=colorRampPalette(brewer.pal("OrRd", n=9)[-(1:2)])(100), 
       smallplot= c(0.125, 0.15, 0.2, 0.5))
  plot(map, border="grey75", add=TRUE)
  title(main="(b) Target group record density (log10-scaled)",
        cex.main=1, font.main=1, adj=0)
  dev.off()
}

plot_background = function(modelName){

  acc = raster(paste0(modelName, "/", modelName, "_accessible.tif"))
  uns = raster(paste0(modelName, "/", modelName, "_unsuitable.tif"))
  
  map = spTransform(countriesCoarseLessIslands, acc@crs)
  
  # plot the backgrounds...
  png(filename=paste0(modelName, "/", sppLab, ".background.png"),
      units = 'in', width = 5, height = 3, 
      res = 400)
  par(mar=c(0.6,0.6,0.6,0.6))
  plot(map, border=NA, ylim=c(-7210547, 9182709))
  plot(uns+acc*2, col=c(NA,terrain.colors(3)[2:1]), 
       zlim=c(0,2), legend=FALSE, axes=FALSE, bty="n", add=TRUE)
  plot(map, border="grey50", add=TRUE)
  legend("topleft", bg="white",
         fill=terrain.colors(3)[1:2],
         legend=c("Accessible", "Unsuitable"), cex=0.7)
  dev.off()
} 

