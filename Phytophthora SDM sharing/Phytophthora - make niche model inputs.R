library(raster)
library(data.table)
library(rworldmap)
library(sf)
library(readr)

rasterOptions(progress="text")

rm(list=ls())

setwd("C:/Users/Dan Chapman/Box/research/phytothreats/niche models")

dir.create("rasters")

# project worldclim layers...
R = stack(paste0("C:/Users/Dan Chapman/Downloads/wc2.0_bio_2.5m_", c("05","06","10","14","15"), ".tif"))
R = projectRaster(R, raster("URBAN_HORT_areakm2.tif"))


# add moisture index
moisture = raster("C:/Users/Dan Chapman/Downloads/global-ai_et0/ai_et0/ai_et0.tif")
moisture = aggregate(moisture, 4, fun=mean, na.rm=TRUE)/10000 # aggregate to ~4x4 km
moisture = projectRaster(moisture, R)
compareRaster(moisture, R)
names(moisture) = "moisture"
R = addLayer(R, moisture)
names(R)

# add land cover layers...
R = stack(R, stack(c("FOREST_areakm2.tif", "URBAN_HORT_areakm2.tif", "AGRIC_areakm2.tif")))
names(R) = c("bio5", "bio6", "bio10", "bio14", "bio15", "moisture", "forest", "urban_hort", "agric")

# apply a land mask...
rLand = raster("land_mask.tif")
#plot(rLand)
#R = mask(R, rLand)
writeRaster(R, "rasters/predictors", overwrite=TRUE)

############################################################################################



R = stack("rasters/predictors.grd")


# read in records data...
recs = fread("global_Phytophthora_records.csv")
table(recs$ISO3, useNA="always")
table(recs$species)
table(grep("alni", recs$species, value=TRUE))
sp = "Phytophthora gonapodyides"

# calculate size of training data (non-UK grid cells with presences) for each species...
nCells = sapply(sort(unique(recs$species)), function(sp){
  message(sp)
  xy = as.matrix(recs[species==sp & ISO3!="GBR", c("longitude","latitude"), with=FALSE])
  if(nrow(xy)>0){
    pts = SpatialPoints(xy, proj4string = CRS("+init=epsg:4326"))
    pts = spTransform(pts, R@crs)
    rSp = rasterize(pts, R, fun="count", background=0) > 0
    rSp = mask(rSp, rLand)
    n = sum(rSp[]==1, na.rm=TRUE)
  } else { n = 0 }
  message("   ", n, " non-UK presence grid cells")
  n
})
cbind(sort(nCells))
write.table(cbind(nCells), "clipboard", sep="\t", col.names=FALSE)

# select species for modelling based on number of training grid cells...
sppList = sort(names(nCells)[nCells > 50])

# plot the records for those species...
pdf("distributions.pdf")
par(mfrow=c(3,3), mar=c(0.3, 0.3, 1.6, 0.3))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
for(sp in sppList){
  plot(countriesCoarseLessIslands, col="grey80", border="grey",
       ylim=c(-50,90), main=sp)
  points(recs$longitude, recs$latitude, pch=16, cex=0.5, col=cbPalette[2])
  idx = recs$species == sp
  points(recs$longitude[idx], recs$latitude[idx], pch=16, cex=0.5)
}
try(dev.off(), silent=TRUE)

# extract training data grid cells...
sapply(sppList, function(sp){
  message(sp)
  xy = as.matrix(recs[species==sp & ISO3!="GBR", c("longitude","latitude"), with=FALSE])
  pts = SpatialPoints(xy, proj4string = CRS("+init=epsg:4326"))
  pts = spTransform(pts, R@crs)
  rSp = rasterize(pts, R, fun="count", background=0) > 0
  rSp = mask(rSp, rLand)
  xyGrid = data.frame(coordinates(rSp)[which(rSp[]==1),])
  write_csv(xyGrid, paste0("rasters/train_", gsub(" ", "_", sp), ".csv"))
  #fn = paste0("rasters/train_", gsub(" ", "_", sp), ".tif")
  #writeRaster(rSp, fn, format="GTiff", overwrite=TRUE)
  NULL
})

# extract validation data grid cells...
sapply(sppList, function(sp){
  message(sp)
  xy = as.matrix(recs[species==sp & ISO3=="GBR", c("longitude","latitude"), with=FALSE])
  pts = SpatialPoints(xy, proj4string = CRS("+init=epsg:4326"))
  pts = spTransform(pts, R@crs)
  rSp = rasterize(pts, R, fun="count", background=0) > 0
  rSp = mask(rSp, rLand)
  xyGrid = coordinates(rSp)[which(rSp[]==1),]
  xyGrid = data.frame(matrix(xyGrid, ncol=2))
  names(xyGrid) = c("x", "y")
  write_csv(xyGrid, paste0("rasters/validation_", gsub(" ", "_", sp), ".csv"))
#  fn = paste0("rasters/validation_", gsub(" ", "_", sp), ".tif")
#  writeRaster(rSp, fn, format="GTiff", overwrite=TRUE)
  NULL
})

# make target group density rasters...
try(dev.off(), silent=TRUE)
sapply(sppList, function(sp){
  message(sp)
  # all points of other species...
  xy = as.matrix(recs[species!=sp, c("longitude","latitude"), with=FALSE])
  pts = SpatialPoints(xy, proj4string = CRS("+init=epsg:4326"))
  pts = spTransform(pts, R@crs)
  rPts = rasterize(pts, R, fun="count", background=0)
  rPts = mask(rPts, rLand)
  plot(aggregate(rPts,4)>0, main=paste(sp, "target group"))
  #plot(spTransform(countriesCoarseLessIslands, R@crs), add=TRUE)
  fn = paste0("rasters/tg_", gsub(" ", "_", sp), ".tif")
  writeRaster(rPts, fn, format="GTiff", overwrite=TRUE)
  NULL
})

# make accessible rasters...
try(dev.off(), silent=TRUE)
for(A in c(100,300)) sapply(sppList, function(sp, aw=A){ # aw = accessible buffer width (km)
  message(sp, " ", aw)
  # non-UK points of the species...
  xy = as.matrix(recs[species==sp & ISO3!="GBR", c("longitude","latitude"), with=FALSE])
  pts = SpatialPoints(xy, proj4string = CRS("+init=epsg:4326"))
  pts = spTransform(pts, R@crs)
  pts = as(pts, "sf")
  pts_buf = st_buffer(pts, dist=aw*1000, nQuadSegs=30)
  pts_buf = st_union(pts_buf)
  r_buf = rasterize(as(pts_buf, "Spatial"), R, background=0)
  r_buf = mask(r_buf, rLand)
  plot(r_buf, main=paste("Accessible for", sp, "with", aw, "km buffer"))
  #plot(spTransform(countriesCoarseLessIslands, rTemplate@crs), add=TRUE)
  fn = paste0("rasters/acc_", gsub(" ", "_", sp), "_", aw, ".tif")
  writeRaster(r_buf, fn, format="GTiff", overwrite=TRUE)
  NULL
})

# get values of predictors at occurrence training grid cells...
try(dev.off(), silent=TRUE)
r = getValues(R)
pdf("niche_plots.pdf", w=9)
xVals = data.frame(rbindlist(lapply(sppList, function(sp){
  message(sp)

  xyTrain = fread(paste0("rasters/train_", gsub(" ", "_", sp), ".csv"))
  xyValid = fread(paste0("rasters/validation_", gsub(" ", "_", sp), ".csv"))
  
  if(nrow(xyValid) > 5){
    xTrain = data.frame(r[cellFromXY(R, xyTrain),]) # data.frame(extract(R, xyTrain))
    xValid = data.frame(r[cellFromXY(R, xyValid),]) # data.frame(extract(R, xyValid))
  
    aw=300
    rAcc = raster(paste0("rasters/acc_", gsub(" ", "_", sp), "_", aw, ".tif"))
    xAcc = r[which(getValues(rAcc)==1),]
  
    par(mfrow=c(3,3))
    for(i in 1:nlayers(R)){
  
      dTrain = density(xTrain[,i], from=min(xAcc[,i]), to=max(xAcc[,i]), adjust=1)
      if(nrow(xValid)>1) dValid = density(xValid[,i], from=min(xAcc[,i]), to=max(xAcc[,i]), adjust=1)
      dAcc = density(xAcc[,i], from=min(xAcc[,i]), to=max(xAcc[,i]), adjust=1)
    
      hist(xAcc[,i], br=100, col="grey", border=NA, 
             xlab=names(R)[i], 
             freq=FALSE, ylim=c(0, max(c(dAcc$y, dTrain$y, dValid$y))),
             main=sp)
      lines(dAcc)
      lines(dTrain, col="tomato")
      lines(dValid, col="blue")
    }
  } else {
    
    xTrain = data.frame(r[cellFromXY(R, xyTrain),]) # data.frame(extract(R, xyTrain))
    xValid = if(nrow(xyValid)>1){ 
        data.frame(r[cellFromXY(R, xyValid),]) 
      } else data.frame(rbind(r[cellFromXY(R, xyValid),]))
    
    aw=300
    rAcc = raster(paste0("rasters/acc_", gsub(" ", "_", sp), "_", aw, ".tif"))
    xAcc = r[which(getValues(rAcc)==1),]
    
    par(mfrow=c(3,3))
    for(i in 1:nlayers(R)){
      
      dTrain = density(xTrain[,i], from=min(xAcc[,i]), to=max(xAcc[,i]), adjust=1)
      dAcc = density(xAcc[,i], from=min(xAcc[,i]), to=max(xAcc[,i]), adjust=1)
      
      hist(xAcc[,i], br=100, col="grey", border=NA, 
           xlab=names(R)[i], 
           freq=FALSE, ylim=c(0, max(c(dAcc$y, dTrain$y))),
           main=sp)
      lines(dAcc)
      lines(dTrain, col="tomato")
      abline(v=xValid[,i], col="blue")
    }
    
    
  }
  
  data.frame(species=sp, xTrain)
})))
try(dev.off(), silent=TRUE)
rm(r)

