###### combine flow data into a single file

###### perform downsampling if necessary
#
#### start with MLN

lab_files <- list.files("FACS_data/Lab_Export_CSVs/MLN", recursive = T, pattern = ".csv", full.names = T)
wild1_files <- list.files("FACS_data/Wild_w1_Export_CSVs/MLN", recursive = T, pattern = ".csv", full.names = T)
wild2_files <- list.files("FACS_data/Wild_w2_Export_CSVs/MLN", recursive = T, pattern = ".csv", full.names = T)

MLN_files <- c(lab_files, wild1_files, wild2_files)

MLN_df <- data.frame(id = NA, read.table(MLN_files[1], T, ',')[1,])
for (i in 1:length(MLN_files)){
  curr <- read.table(MLN_files[i], T, sep=',')
  if(grepl("Wild", MLN_files[i])){
    id <- gsub(" MLN.*", "", basename(MLN_files[i]))
    id <- gsub(".* .* ", "", id)
  } else {
    id <- gsub(".*MLN ", "", basename(MLN_files[i]))
    id <- gsub(" .*", "", id)
  }
  MLN_df <- rbind(MLN_df, data.frame(id=rep(id,nrow(curr)), curr))
}
MLN_df <- MLN_df[-1,]

metadata <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')

MLN_full <- merge(data.frame(id=unique(MLN_df$id)), metadata, all.x = T, by.x = "id", by.y = "mouse_id")
MLN_full$id <- factor(MLN_full$id, levels = MLN_full$id)

write.table(MLN_df, "inputs/MLN_df.txt", sep='\t', row.names=F, col.names=F)

##
#

## start with MLN myeloid
lab_files <- list.files("FACS_data/Lab_Export_CSVs/MLN Myeloid", recursive = T, pattern = ".csv", full.names = T)
wild1_files <- list.files("FACS_data/Wild_w1_Export_CSVs/MLN Myeloid", recursive = T, pattern = ".csv", full.names = T)
wild2_files <- list.files("FACS_data/Wild_w2_Export_CSVs/MLN Myeloid", recursive = T, pattern = ".csv", full.names = T)

MLN_files <- c(lab_files, wild1_files, wild2_files)

MLN_df <- data.frame(id = NA, read.table(MLN_files[1], T, ',')[1,])
for (i in 1:length(MLN_files)){
  curr <- read.table(MLN_files[i], T, sep=',')
  if(grepl("Wild", MLN_files[i])){
    id <- gsub(" MLN.*", "", basename(MLN_files[i]))
    id <- gsub(".* .* ", "", id)
  } else {
    id <- gsub(".*MLN ", "", basename(MLN_files[i]))
    id <- gsub(" .*", "", id)
  }
  MLN_df <- rbind(MLN_df, data.frame(id=rep(id,nrow(curr)), curr))
}
MLN_df <- MLN_df[-1,]

metadata <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')

MLN_full <- merge(data.frame(id=unique(MLN_df$id)), metadata, all.x = T, by.x = "id", by.y = "mouse_id")
MLN_full$id <- factor(MLN_full$id, levels = MLN_full$id)

write.table(MLN_df, "inputs/MLN_myeloid_df.txt", sep='\t', row.names=F, col.names=F)


##
#### start with Blood

lab_files <- list.files("FACS_data/Lab_Export_CSVs/Blood", recursive = T, pattern = ".csv", full.names = T)
wild1_files <- list.files("FACS_data/Wild_w1_Export_CSVs/Blood", recursive = T, pattern = ".csv", full.names = T)
wild2_files <- list.files("FACS_data/Wild_w2_Export_CSVs/Blood", recursive = T, pattern = ".csv", full.names = T)

Blood_files <- c(lab_files, wild1_files, wild2_files)

Blood_df <- data.frame(id = NA, read.table(Blood_files[1], T, ',')[1,])
for (i in 1:length(Blood_files)){
  curr <- read.table(Blood_files[i], T, sep=',')
  if(grepl("Wild", Blood_files[i])){
    id <- gsub(" PBMC.*", "", basename(Blood_files[i]))
    id <- sub(".* ", "", id)
  } else {
    id <- gsub(".* PBMC ", "", basename(Blood_files[i]))
    id <- gsub(".* PBMC", "", id)
    id <- sub(" .*", "", id)
  }
  Blood_df <- rbind(Blood_df, data.frame(id=rep(id,nrow(curr)), curr))
}
Blood_df <- Blood_df[-1,]

metadata <- read.table("mice_metadata.11.19_mouse_id.txt", T, '\t')

Blood_full <- merge(data.frame(id=unique(Blood_df$id)), metadata, all.x = T, by.x = "id", by.y = "mouse_id")
Blood_full$id <- factor(Blood_full$id, levels = Blood_full$id)

write.table(Blood_df, "inputs/Blood_df.txt", sep='\t', row.names=F, col.names=F)
