library(SAVER)


data.path <- "E:/����ֽ⽵ά/SCEnddemo/traj/samp/cellmix4.csv"

raw.data <- read.csv(data.path,header=T, row.names=1, sep=",")




#data_raw = data_s4@assays[[".->data"]]@listData[["counts"]]
res_saver <- SAVER::saver(raw.data, estimates.only = TRUE , ncores = 12)
# write.table(res_saver, "cellcircle/3Line-qPCR_SAVER.csv",sep=",")
saveRDS(res_saver, file = "E:/����ֽ⽵ά/SCEnddemo/traj/SAVER_samp/cellmix4.rds")


