library(dplyr)
library(glue)

setwd("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/")

Helmsley_sample <- read.table("Helmsley_sample.txt", sep = "\t", header = F)

# pick_df <- lapply(Helmsley_sample[,1], function(idx){
#   pickx <- read.table(glue("02_clean/{idx}_blank_pick.txt"), header = T, sep = "\t")
# }) %>% Reduce("rbind", .) %>% as.data.frame()
# rownames(pick_df) <- Helmsley_sample[,1]
# write.table(pick_df, file = "Helmsley_dsb.txt", col.names = T, row.names = T, quote = F, sep = "\t")

pick_manual <- read.table("thresh_dsb_pick_manual.txt", header = F, sep = "\t")
pick_df <- cbind(Helmsley_sample,
                 pick_manual[match(Helmsley_sample[,1], pick_manual[,1]), 2:3])
write.table(pick_df, file = "Helmsley_pre_dsb.txt", col.names = F, row.names = F, quote = F, sep = "\t")

