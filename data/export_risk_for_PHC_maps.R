rm(list=ls())

library(readr)
library(dplyr)
library(reshape2)
library(htmlTable)
library(countrycode)
library(tidyverse)


####################### generate predictions of export risk for PHC maps ######################
load("trade_model_allspp_zinfl.rData")
load("arrivals_unscaled.rData")
trade_model <- trade_model_allspp_zinfl


# prepare the data needed to unscale the variables for plotting
trade_scale <- 2 * sd(log(1 + arrivals_unscaled$trade), na.rm = TRUE)
trade_center <- mean(log(1 + arrivals_unscaled$trade), na.rm = TRUE)

climate_scale <- 2 * sd(log(arrivals_unscaled$climate_similarity), na.rm = TRUE)
climate_center <- mean(log(arrivals_unscaled$climate_similarity), na.rm = TRUE)

EPPO_reporting_service_scale <- 2 * sd(sqrt(arrivals_unscaled$EPPO_reporting_service), na.rm = TRUE)
EPPO_reporting_service_center <- mean(sqrt(arrivals_unscaled$EPPO_reporting_service), na.rm = TRUE)

impproactive_scale <- 2 * sd(arrivals_unscaled$impproactive, na.rm = TRUE)
impproactive_center <- mean(arrivals_unscaled$impproactive, na.rm = TRUE)

expproactive_scale <- 2 * sd(arrivals_unscaled$expproactive, na.rm = TRUE)
expproactive_center <- mean(arrivals_unscaled$expproactive, na.rm = TRUE)

temp_range_scale <-  2 * sd(arrivals_unscaled$temp_range, na.rm = TRUE)
temp_range_center <- mean(arrivals_unscaled$temp_range, na.rm = TRUE)

oospore_wall_index_scale <-  2 * sd(log(arrivals_unscaled$oospore_wall_index), na.rm = TRUE)
oospore_wall_index_center <- mean(log(arrivals_unscaled$oospore_wall_index), na.rm = TRUE)


# export risk to UK
# fix all other variables at their value for the UK
tradeMatrix <- read_rds("live_tradeMatrix.rds")
sourceMatrix <- read_rds("phytophthora_source_matrix.rds")
arrivalMatrix <- read_rds("phytophthora_arrival_matrix.rds")
# Phytophthora distributions
presenceMatrix <- sourceMatrix + arrivalMatrix
table(presenceMatrix, useNA = "always") # should be present or absent or NA (e.g. no countries should be both source and recipient == 2)
load("MD_mean.rData")
climate_dist <- MD_mean
climate_sim <- 1/climate_dist # inverse for similarity
biosecurityMatrix <- read_rds("biosecurityMatrix.rds")
expbiosecurityMatrix <- read_rds("expbiosecurityMatrix.rds")
impbiosecurityMatrix <- read_rds("impbiosecurityMatrix.rds")
phytosanitary <- read.csv("biosecurity_EPPORS.csv", stringsAsFactors = FALSE) %>% rename(exporter_iso3 = importer_iso3)
rownames(phytosanitary) <- phytosanitary$exporter_iso3


# composite export risk to UK based on all risk factors (excluding species-level variance - traits, phylogeny)
# scale all variables by original mean and sd
composite_newdata <- full_join(full_join(((log(1 + tradeMatrix[,"GBR"]) - trade_center) / trade_scale) %>% 
                                           as.data.frame() %>% rownames_to_column(var = "exporter_iso3") %>% rename(trade = 2),
                                         ((log(1 + climate_sim[,"GBR"]) - climate_center) / climate_scale) %>% 
                                           as.data.frame() %>% rownames_to_column(var = "exporter_iso3") %>% rename(climate_similarity = 2)),
                               ((expbiosecurityMatrix[,"GBR"] - expproactive_center) / expproactive_scale) %>% 
                                 as.data.frame() %>% rownames_to_column(var = "exporter_iso3") %>% rename(expproactive = 2)
) %>% mutate(impproactive = (impbiosecurityMatrix["GBR",] - impproactive_center) / impproactive_scale,
             importer_iso3 = "GBR",
             EPPO_reporting_service = (sqrt(phytosanitary["GBR","EPPO_reporting_service"]) - EPPO_reporting_service_center) / EPPO_reporting_service_scale,
             temp_range = mean(trade_model$data$temp_range),
             oospore_wall_index = mean(trade_model$data$oospore_wall_index),
             chlamydospores = 0)


# use this data to predict export risk posed to UK by country

export_risk_composite <- as.data.frame(fitted(trade_model, 
                                              newdata = composite_newdata, 
                                              sample_new_levels = "gaussian",
                                              allow_new_levels = TRUE
)) %>% 
  rename(composite_export_risk = Estimate, compositeQ2.5 = Q2.5, compositeQ97.5 = Q97.5) %>% 
  select(1,3,4) %>% mutate(importer_iso3 = composite_newdata$importer_iso3, 
                           exporter_iso3 = composite_newdata$exporter_iso3)

export_risk_composite <- 
  full_join(export_risk_composite,
            full_join(tradeMatrix[,"GBR"]  %>% as.data.frame() %>% rownames_to_column(var = "exporter_iso3") %>% rename(trade_raw = 2),
                      full_join(climate_sim[,"GBR"]  %>% as.data.frame() %>% rownames_to_column(var = "exporter_iso3") %>% rename(climate_raw = 2),
                                full_join(expbiosecurityMatrix[,"GBR"] %>% as.data.frame() %>% rownames_to_column(var = "exporter_iso3") %>% rename(expproactive_raw = 2),
                                          full_join(phytosanitary,
                                                    rowSums(presenceMatrix, na.rm = TRUE) %>% 
                                                      as.data.frame() %>% 
                                                      rownames_to_column(var = "exporter_iso3") %>% 
                                                      rename(Phytophthora_richness = 2))))))


# marginal risk: trade
trade_newdata <-  data.frame(exporter_iso3  = rownames(tradeMatrix),
                                                  importer_iso3 = "GBR",
                                                  trade = (log(1 + tradeMatrix[,"GBR"]) - trade_center) / trade_scale,
                                                  climate_similarity = mean(trade_model$data$climate_similarity),
                                                  impproactive = unique(trade_model$data$impproactive[trade_model$data$importer_iso3 == "GBR"]),
                                                  expproactive = mean(trade_model$data$expproactive),
                                                  EPPO_reporting_service = unique(trade_model$data$EPPO_reporting_service[trade_model$data$importer_iso3 == "GBR"]),
                                                  temp_range = mean(trade_model$data$temp_range),
                                                  oospore_wall_index = mean(trade_model$data$oospore_wall_index),
                                                  chlamydospores = 0, #reference level
                                                  phylo = NA,
                                                  spp = NA)

# predict export risk based on the volume of live plant trade from each country to the UK
trade_risk <- as.data.frame(fitted(trade_model,
                                          newdata = trade_newdata,
                                          allow_new_levels = TRUE,
                                          sample_new_levels = "gaussian",
                                          robust = TRUE)) %>%
  rename(trade_risk = Estimate, tradeQ2.5 = Q2.5, tradeQ97.5 = Q97.5) %>%
  select(1,3,4) %>% mutate(exporter_iso3 = trade_newdata$exporter_iso3)

# marginal risk: climate matching
climate_newdata <-  data.frame(exporter_iso3  = rownames(climate_sim),
                             importer_iso3 = "GBR",
                             trade = mean(trade_model$data$trade),
                             climate_similarity = (log(climate_sim[,"GBR"]) - climate_center) / climate_scale,
                             impproactive = unique(trade_model$data$impproactive[trade_model$data$importer_iso3 == "GBR"]),
                             expproactive = mean(trade_model$data$expproactive),
                             EPPO_reporting_service = unique(trade_model$data$EPPO_reporting_service[trade_model$data$importer_iso3 == "GBR"]),
                             temp_range = mean(trade_model$data$temp_range),
                             oospore_wall_index = mean(trade_model$data$oospore_wall_index),
                             chlamydospores = 0, #reference level
                             phylo = NA,
                             spp = NA)

# predict export risk based on the volume of live plant trade from each country to the UK
climate_risk <- as.data.frame(fitted(trade_model,
                                   newdata =climate_newdata,
                                   allow_new_levels = TRUE,
                                   sample_new_levels = "gaussian",
                                   robust = TRUE)) %>%
  rename(climate_risk = Estimate, climateQ2.5 = Q2.5, climateQ97.5 = Q97.5) %>%
  select(1,3,4) %>% mutate(exporter_iso3 = climate_newdata$exporter_iso3)


# marginal risk: export biosecurity
expproactive_newdata <-  data.frame(exporter_iso3  = rownames(expbiosecurityMatrix),
                             importer_iso3 = "GBR",
                             trade = mean(trade_model$data$trade),
                             climate_similarity = mean(trade_model$data$climate_similarity),
                             impproactive = unique(trade_model$data$impproactive[trade_model$data$importer_iso3 == "GBR"]),
                             expproactive = (expbiosecurityMatrix[,"GBR"] - expproactive_center) / expproactive_scale,
                             EPPO_reporting_service = unique(trade_model$data$EPPO_reporting_service[trade_model$data$importer_iso3 == "GBR"]),
                             temp_range = mean(trade_model$data$temp_range),
                             oospore_wall_index = mean(trade_model$data$oospore_wall_index),
                             chlamydospores = 0, #reference level
                             phylo = NA,
                             spp = NA)

# predict export risk based on the volume of live plant trade from each country to the UK
expproactive_risk <- as.data.frame(fitted(trade_model,
                                   newdata = expproactive_newdata,
                                   allow_new_levels = TRUE,
                                   sample_new_levels = "gaussian",
                                   robust = TRUE)) %>%
  rename(expproactive_risk = Estimate, expproactiveQ2.5 = Q2.5, expproactiveQ97.5 = Q97.5) %>%
  select(1,3,4) %>% mutate(exporter_iso3 = expproactive_newdata$exporter_iso3)



presenceMatrix <- arrivalMatrix + sourceMatrix
presenceDF <- melt(presenceMatrix) %>% rename(exporter_iso3 = Var1, species = Var2, presence = value)
presenceDF <- presenceDF[complete.cases(presenceDF),]
presenceDF <- presenceDF[presenceDF$presence == 1,]

UK_species <- read.csv("UK_first_records.csv", stringsAsFactors = FALSE)

UK_risk_register <- 
  c("Phytophthora acerina",
    "Phytophthora alni",
    "Phytophthora austrocedri",
    "Phytophthora bishii",
    "Phytophthora chrysanthemi",
    "Phytophthora foliorum",
    "Phytophthora fragariae",
    "Phytophthora fragariaefolia",
    "Phytophthora infestans",
    "Phytophthora kernoviae",
    "Phytophthora lateralis",
    "Phytophthora pinifolia",
    "Phytophthora pluvialis",
    "Phytophthora polonica",
    "Phytophthora pseudosyringae",
    "Phytophthora ramorum",
    "Phytophthora rubi",
    "Phytophthora siskiyouensis")

traits <- read.csv("2018-01-29_Phytophthora_trait_database_merged.csv", 
                   stringsAsFactors = FALSE,
                   na.strings = c("NA", "#N/A"))


spp_risk_newdata <- left_join(presenceDF %>% select(-presence), traits[,c("species_name", intersect(colnames(traits), colnames(trade_model$data)))] %>% 
                           transmute(species = species_name,
                                     oospore_wall_index = (log(oospore_wall_index) - oospore_wall_index_center) / oospore_wall_index_scale,
                                     temp_range = (temp_range - temp_range_center) / temp_range_scale,
                                     chlamydospores = ifelse(grepl("^N", chlamydospores), yes = 0, no = 1) ))




spp_risk_newdata <- left_join(spp_risk_newdata, 
                              composite_newdata %>% select(-c(temp_range, chlamydospores, oospore_wall_index)))

# if a species is already present in UK, then replace arrival risk estimate with NA

spp_data <- as.data.frame(fitted(trade_model,
                                          newdata = spp_risk_newdata,
                                          allow_new_levels = TRUE,
                                          sample_new_levels = "uncertainty",
                                          robust = TRUE)) %>%
  rename(spp_risk = Estimate, spp_riskQ2.5 = Q2.5, spp_riskQ97.5 = Q97.5) %>%
  select(1,3,4) %>% mutate(exporter_iso3 = spp_risk_newdata$exporter_iso3,
                           species = spp_risk_newdata$species,
                           inUK = ifelse(spp_risk_newdata$species %in% UK_species$species, yes = "yes", no = "no"),
                           risk_register = ifelse(spp_risk_newdata$species %in% UK_risk_register, yes = "yes", no = "no"))

spp_data$format_risk <- paste0(sprintf("%.3f", spp_data$spp_risk), " (", sprintf("%.3f", spp_data$spp_riskQ2.5), ", ", sprintf("%.3f", spp_data$spp_riskQ97.5), ")")

spp_data$format_risk[spp_data$inUK == "yes"] <- NA
spp_data$format_risk[grepl("^N", spp_data$format_risk)] <- NA




species_list_by_country <- split(spp_data, f = spp_data$exporter_iso3)
for (i in 1:length(species_list_by_country)){
  species_list_by_country[[i]] <- 
    species_list_by_country[[i]][order(species_list_by_country[[i]]$spp_risk, decreasing = TRUE),]
}





# add pop up tables with species info for each country
export_risk_composite$species_table <- NA
for (i in export_risk_composite$exporter_iso3){
  if(!is.null(species_list_by_country[[i]])){
    export_risk_composite$species_table[export_risk_composite$exporter_iso3 == i] <- 
      gsub('<table', '<table width=400', gsub('<td', '<td nowrap="nowrap"; ', htmlTable(species_list_by_country[[i]][,5:8],
                                                     header = c("species", "in UK", "UK risk register", "arrival risk (lower, upper)"),
                                                     caption = paste0(countrycode(sourcevar = i,
                                                                                  origin = "iso3c",
                                                                                  destination = "country.name"), " - species reported"),
                                                     rnames = FALSE)))
  }
  if(is.null(species_list_by_country[[i]])) {
    export_risk_composite$species_table[export_risk_composite$exporter_iso3 == i] <- NA
  }
}

save(export_risk_composite, file = "C:/Users/loubar/OneDrive - NERC/PHC_tools/data/export_risk_composite.rData")









