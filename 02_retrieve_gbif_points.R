# Send names to GBIF and download points
# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
library(rgbif)

wcvp_climbers <- readRDS("wcvp_nt_climbers_final.Rdata")

# check reference table for GBIF names
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
reference_table$gbif_name <- fix.names.taxize(reference_table$gbif_name)
subset_reference_table <- subset(reference_table, reference_table$wcvp_name %in% wcvp_climbers$taxon_name)

names_to_retrieve <- subset_reference_table$gbif_name

# Now we send a request to GBIF to download the points for this list of species 
user <- "" # username
pwd <- "" # password
email <- "@gmail.com" # email

rgbif::occ_download(rgbif::pred_in("scientificName", names_to_retrieve),
                    pred_in("basisOfRecord", 'PRESERVED_SPECIMEN'),
                    pred("hasCoordinate", TRUE),
                    format = "SIMPLE_CSV", user=user,pwd=pwd,email=email) # Sending request to GBIF
