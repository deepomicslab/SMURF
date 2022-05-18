# Mo Huang, mohuang@wharton.upenn.edu
# Perform Melanoma Drop-seq, FISH analysis

library(SAVER)

dropseq <- readRDS("melanoma_dropseq.rds")

fish <- read.table("fishSubset.txt", header = TRUE, row.names = 1)

fish.filt <- fish[fish[, "GAPDH"] < quantile(fish[, "GAPDH"], 0.9) & 
                    fish[, "GAPDH"] > quantile(fish[, "GAPDH"], 0.1), ]
fish.norm <- 
  sweep(fish.filt, 1, fish.filt[, "GAPDH"]/mean(fish.filt[, "GAPDH"]), "/")

n.genes <- nrow(dropseq)
n.cells <- ncol(dropseq)
gene.names <- rownames(dropseq)
cell.names <- colnames(dropseq)

genes <- which(gene.names %in% colnames(fish.norm))

sf <- colSums(dropseq)/mean(colSums(dropseq))

normalize.gapdh <- function(x) {
  x.filt <- x[, x["GAPDH", ] < quantile(x["GAPDH", ], 0.9) &
                x["GAPDH", ] > quantile(x["GAPDH", ], 0.1)]
  x.norm <- sweep(x.filt, 2, x.filt["GAPDH", ]/mean(x.filt["GAPDH", ]), "/")
}

dropseq.filt <- dropseq[genes, dropseq["GAPDH", ] < quantile(dropseq["GAPDH", ], 0.9) &
                          dropseq["GAPDH", ] > quantile(dropseq["GAPDH", ], 0.1)]


dropseq.norm <- normalize.gapdh(dropseq[genes, ])

saver <- readRDS("melanoma_dropseq_saver.rds")

saver.filt <- saver$estimate[, saver$estimate["GAPDH", ] < 
                               quantile(saver$estimate["GAPDH", ], 0.9) &
                               saver$estimate["GAPDH", ] >
                               quantile(saver$estimate["GAPDH", ], 0.1)]

magic <- readRDS("melanoma_dropseq_magic.rds")
magic.filt <- magic[, magic["GAPDH", ] < quantile(magic["GAPDH", ], 0.9) &
                      magic["GAPDH", ] > quantile(magic["GAPDH", ], 0.1)]
magic.norm <- normalize.gapdh(magic)

scimpute <- readRDS("SAVER-data/melanoma_dropseq_scimpute.rds")
scimpute.filt <- scimpute[, scimpute["GAPDH", ] < quantile(scimpute["GAPDH", ], 0.9) &
                      scimpute["GAPDH", ] > quantile(scimpute["GAPDH", ], 0.1)]
scimpute.norm <- normalize.gapdh(scimpute)

knnsmooth <- read.csv("melanoma_knnsmooth.csv", header = TRUE, row.names = 1)
knnsmooth.filt <- knnsmooth[knnsmooth[, "GAPDH"] < quantile(knnsmooth[, "GAPDH"], 0.9) & 
                    knnsmooth[, "GAPDH"] > quantile(knnsmooth[, "GAPDH"], 0.1), ]
knnsmooth.norm <- 
  sweep(knnsmooth.filt, 1, knnsmooth.filt[, "GAPDH"]/mean(knnsmooth.filt[, "GAPDH"]), "/")


smurf <- read.csv("melanoma_SMURF_CV.csv", header = TRUE, row.names = 1)
smurf.filt <- smurf[smurf[, "GAPDH"] < quantile(smurf[, "GAPDH"], 0.9) & 
                              smurf[, "GAPDH"] > quantile(smurf[, "GAPDH"], 0.1), ]
smurf.norm <- 
  sweep(smurf.filt, 1, smurf.filt[, "GAPDH"]/mean(smurf.filt[, "GAPDH"]), "/")



mult.factor <- sapply(1:16, function(x)
  mean(fish.filt[, gene.names[genes[x]]], na.rm = TRUE)/mean(dropseq.filt[x, ]))
names(mult.factor) <- gene.names[genes]

mult.factor2 <- sapply(1:16, function(x)
  mean(fish.filt[, gene.names[genes[x]]], na.rm = TRUE)/mean(saver.filt[x, ]))
names(mult.factor2) <- gene.names[genes]

mult.factor3 <- sapply(1:16, function(x)
  mean(fish.filt[, gene.names[genes[x]]], na.rm = TRUE)/mean(magic.filt[x, ]))
names(mult.factor3) <- gene.names[genes]

mult.factor4 <- sapply(1:16, function(x)
  mean(fish.filt[, gene.names[genes[x]]], na.rm = TRUE)/mean(scimpute.filt[x, ]))
names(mult.factor4) <- gene.names[genes]

mult.factor5 <- sapply(1:16, function(x)
  mean(fish.filt[, gene.names[genes[x]]], na.rm = TRUE)/mean(knnsmooth.filt[, gene.names[genes[x]]]))
names(mult.factor5) <- gene.names[genes]

mult.factor6 <- sapply(1:16, function(x)
  mean(fish.filt[, gene.names[genes[x]]], na.rm = TRUE)/mean(smurf.filt[, gene.names[genes[x]]]))
names(mult.factor6) <- gene.names[genes]
set.seed(5)

lambda.samp <- sample.saver(saver)

lambda.samp.filt <- 
  lambda.samp[, lambda.samp["GAPDH", ] < quantile(lambda.samp["GAPDH", ], 0.9) &
                lambda.samp["GAPDH", ] > quantile(lambda.samp["GAPDH", ], 0.1)]
lambda.samp.norm <- 
  sweep(lambda.samp.filt, 2,
        lambda.samp.filt["GAPDH", ]/mean(lambda.samp.filt["GAPDH", ]), "/")


###############################################################################
## Gini analysis

library(reldist)
fish.gini <- apply(fish.norm[, gene.names[genes]], 2, function(x) 
  gini(x[complete.cases(x)]))
knnsmooth.gini <- apply(knnsmooth.norm[, gene.names[genes]], 2, function(x) 
  gini(x[complete.cases(x)]))
smurf.gini <- apply(smurf.norm[, gene.names[genes]], 2, function(x) 
  gini(x[complete.cases(x)]))
dropseq.gini <- apply(dropseq.norm, 1, gini)
saver.lambda.gini <- apply(lambda.samp.norm, 1, gini)
names(saver.lambda.gini) <- gene.names[genes]
magic.gini <- apply(magic.norm, 1, gini)



# pdf("plots/fig1b_gini1.pdf", 7.5, 4.25)
# 
# par(mfrow = c(1, 2), cex.main = 1.5, mar = c(4, 4, 1, 2) + 0.1, oma = c(3, 2, 3, 0), 
#     mgp = c(3.5, 1, 0),
#     cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
# 
# plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), pch = 21, xlab = "",
#      ylab = "", bg = "#bdd7e7", cex = 1.5)
# abline(0, 1, col = "gray", lty = 2)
# points(fish.gini[-6], dropseq.gini[-6], pch = 21, bg = "#bdd7e7", cex = 3)
# par(las = 0)
# mtext("Drop-seq Gini", side = 2, line = 3.5, cex = 1.5)
# mtext("FISH Gini", side = 1, line = 3, cex = 1.5)
# points(fish.gini["LMNA"], dropseq.gini["LMNA"], pch = 21, cex = 3,
#        bg = "#a63603")
# points(fish.gini["CCNA2"], dropseq.gini["CCNA2"], pch = 21, cex = 3,
#        bg = "#31a354")
# text(0.42, 0.82, "LMNA", cex = 1.2)
# text(0.72, 0.9, "CCNA2", cex = 1.2)
# text(.7, .1, paste("RMSE = ", round(sqrt(mean((fish.gini[-6]-dropseq.gini[-6])^2)), 2)),
#      cex = 1.2)
# text(0.78, 0.22, paste("r = ", round(cor(fish.gini[-6], dropseq.gini[-6]), 2)),
#      cex = 1.2)
# 
# 
# par(las = 1)
# plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), pch = 21, xlab = "",
#      ylab = "", bg = "#9ecae1", cex = 1.5)
# abline(0, 1, col = "gray", lty = 2)
# points(fish.gini[-6], saver.lambda.gini[-6], pch = 21, bg = "#bdd7e7", cex = 3)
# par(las = 0)
# mtext("SAVER Gini", side = 2, line = 3.5, cex = 1.5)
# mtext("FISH Gini", side = 1, line = 3, cex = 1.5)
# points(fish.gini["LMNA"], saver.lambda.gini["LMNA"], pch = 21, cex = 3,
#        bg = "#a63603")
# points(fish.gini["CCNA2"], saver.lambda.gini["CCNA2"], pch = 21, cex = 3,
#        bg = "#31a354")
# text(0.13, 0.40, "LMNA", cex = 1.2)
# text(0.76, 0.48, "CCNA2", cex = 1.2)
# text(.7, .1, paste("RMSE = ", round(sqrt(mean((fish.gini[-6]-saver.lambda.gini[-6])^2)), 2)),
#      cex = 1.2)
# text(0.8, 0.22, paste("r = ", round(cor(fish.gini[-6], saver.lambda.gini[-6]), 2)),
#      cex = 1.2)
# 
# dev.off()
# 
# ### supplementary gini
# pdf("plots/knnsmooth_gini.pdf", 4.5, 4.5)
# par(mfrow = c(1, 1), cex.main = 1.5, mar = c(4, 4, 1, 2) + 0.1, oma = c(3, 2, 3, 0), 
#     mgp = c(3.5, 1, 0),
#     cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
# 
# plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), pch = 21, xlab = "",
#      ylab = "", bg = "#bdd7e7", cex = 1.5)
# abline(0, 1, col = "gray", lty = 2)
# points(fish.gini[-6], knnsmooth.gini[-6], pch = 21, bg = "#bdd7e7", cex = 3)
# par(las = 0)
# mtext("KNN-Smoothing Gini", side = 2, line = 3.5, cex = 1.5)
# mtext("FISH Gini", side = 1, line = 3, cex = 1.5)
# points(fish.gini["LMNA"], knnsmooth.gini["LMNA"], pch = 21, cex = 3,
#        bg = "#a63603")
# points(fish.gini["CCNA2"], knnsmooth.gini["CCNA2"], pch = 21, cex = 3,
#        bg = "#31a354")
# text(0.1, 0.18, "LMNA", cex = 1.2)
# text(0.7, 0.07, "CCNA2", cex = 1.2)
# text(.7, .8, paste("RMSE = ", round(sqrt(mean((fish.gini[-6]-knnsmooth.gini[-6])^2)), 2)),
#      cex = 1.2)
# text(0.78, 0.92, paste("r = ", round(cor(fish.gini[-6], knnsmooth.gini[-6]), 2)),
#      cex = 1.2)
# dev.off()
# 
# pdf("plots/smurf_gini.pdf", 4.5, 4.5)
# par(mfrow = c(1, 1), cex.main = 1.5, mar = c(4, 4, 1, 2) + 0.1, oma = c(3, 2, 3, 0), 
#     mgp = c(3.5, 1, 0),
#     cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
# 
# plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), pch = 21, xlab = "",
#      ylab = "", bg = "#bdd7e7", cex = 1.5)
# abline(0, 1, col = "gray", lty = 2)
# points(fish.gini[-6], smurf.gini[-6], pch = 21, bg = "#bdd7e7", cex = 3)
# par(las = 0)
# mtext("SCend Gini", side = 2, line = 3.5, cex = 1.5)
# mtext("FISH Gini", side = 1, line = 3, cex = 1.5)
# points(fish.gini["LMNA"], smurf.gini["LMNA"], pch = 21, cex = 3,
#        bg = "#a63603")
# points(fish.gini["CCNA2"], smurf.gini["CCNA2"], pch = 21, cex = 3,
#        bg = "#31a354")
# text(0.1, 0.18, "LMNA", cex = 1.2)
# text(0.7, 0.07, "CCNA2", cex = 1.2)
# text(.7, .8, paste("RMSE = ", round(sqrt(mean((fish.gini[-6]-smurf.gini[-6])^2)), 2)),
#      cex = 1.2)
# text(0.78, 0.92, paste("r = ", round(cor(fish.gini[-6], smurf.gini[-6]), 2)),
#      cex = 1.2)
# dev.off()
# 
# 

###############################################################################
## LMNA and CCNA2 density

set.seed(1)

Y.samp <- t(sapply(1:16, function(i)
  rnbinom(n.cells, mu = saver$estimate[i, ]*mult.factor2[i],
          size = saver$estimate[i, ]^2/saver$se[i, ]^2)))

rownames(Y.samp) <- gene.names[genes]

Y.samp.filt <- Y.samp[, Y.samp["GAPDH", ] < quantile(Y.samp["GAPDH", ], 0.9) &
                        Y.samp["GAPDH", ] > quantile(Y.samp["GAPDH", ], 0.1)]

Y.samp.norm <- sweep(Y.samp.filt, 2, Y.samp.filt["GAPDH", ]/
                       mean(Y.samp.filt["GAPDH", ]), "/")


pdf("plots/fig1c_dens3.pdf", 4.75, 6.5)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(1, 2, 1.5, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.04), 
     xlim = c(-10, 500), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 500, by = 100))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("LMNA", side = 3, line = -2, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "LMNA"], na.rm = TRUE)$bw
lines(density(fish.norm[, "LMNA"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["LMNA", ]*mult.factor["LMNA"], bw = dens.bw), lwd = 3)
lines(density(Y.samp.norm["LMNA", ], bw = dens.bw), lwd = 3, lty = 2)
legend(250, 0.03, c("RNA FISH", "Drop-seq", "SAVER"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)


dens.bw <- density(fish.norm[, "LMNA"], na.rm = TRUE)$bw
de_fish = density(fish.norm[, "LMNA"], na.rm = TRUE)
de_dropseq = density(dropseq.norm["LMNA", ]*mult.factor["LMNA"], bw = dens.bw)
de_saver = density(Y.samp.norm["LMNA", ], bw = dens.bw)
de_magic =  density(magic.norm["LMNA", ]*mult.factor3["LMNA"], bw = dens.bw)
de_scimpute = density(scimpute.norm["LMNA",]*mult.factor4["LMNA"], bw = dens.bw)
de_knnsmooth = density(knnsmooth.norm[,"LMNA" ]*mult.factor5["LMNA"], bw = dens.bw)
de_smurf = density(smurf.norm[,"LMNA" ]*mult.factor6["LMNA"], bw = dens.bw)

 
kld_base = function(x,y,...){
  integrand = function(x,y,t){
    de_fish = density(x, na.rm = TRUE)
    dens.bw = de_fish$bw
    de_impute = density(y, bw=dens.bw)
    f.x =  approx(de_fish$x, de_fish$y,t)$y
    f.y =  approx(de_impute$x,de_impute$y,t)$y
    tmpRatio = f.x *(log2(f.x) - log2(f.y))
    tmpRatio = ifelse(is.infinite(tmpRatio),0,ifelse(is.na(tmpRatio),0,tmpRatio))
    return(tmpRatio)
  }
  return(integrate(integrand,-Inf,Inf,x = x,y = y,stop.on.error=FALSE)$value)
}

set.seed(5)


ge = "BABAM1"
x = fish.norm[, ge]

y_dropseq = dropseq.norm[ge, ]*mult.factor[ge]
y_scimpute = scimpute.norm[ge,]*mult.factor4[ge]
y_magic =  magic.norm[ge, ]*mult.factor3[ge]
y_knnsmooth = knnsmooth.norm[,ge ]*mult.factor5[ge]
y_saver = Y.samp.norm[ge, ]
y_smurf = smurf.norm[,ge ]*mult.factor6[ge]


kld_base(x,y_dropseq)
kld_base(x,y_scimpute)
kld_base(x,y_magic)
kld_base(x,y_knnsmooth)
kld_base(x,y_saver)
kld_base(x,y_smurf)








plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.25), 
     xlim = c(-2, 60), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 60, by = 20))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("CCNA2", side = 3, line = -2, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "CCNA2"], na.rm = TRUE)$bw
lines(density(fish.norm[, "CCNA2"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["CCNA2", ]*mult.factor["CCNA2"], bw = dens.bw), lwd = 3)
lines(density(Y.samp.norm["CCNA2", ], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()

pdf("plots/suppfig3c_dens1.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.04), 
     xlim = c(-10, 500), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 500, by = 100))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("LMNA", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "LMNA"], na.rm = TRUE)$bw
lines(density(fish.norm[, "LMNA"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(magic.norm["LMNA", ]*mult.factor3["LMNA"], bw = dens.bw), lwd = 3)
lines(density(knnsmooth.norm[,"LMNA" ]*mult.factor5["LMNA"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "MAGIC", "KNN-smoothing"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.25), 
     xlim = c(-2, 100), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 100, by = 20))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("CCNA2", side = 3, line = 0, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "CCNA2"], na.rm = TRUE)$bw
lines(density(fish.norm[, "CCNA2"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(magic.norm["CCNA2", ]*mult.factor3["CCNA2"], bw = dens.bw), lwd = 3)
lines(density(knnsmooth.norm[,"CCNA2" ]*mult.factor5["CCNA2"], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()



pdf("plots/smurf1.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.04), 
     xlim = c(-10, 500), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 500, by = 100))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("LMNA", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "LMNA"], na.rm = TRUE)$bw
lines(density(fish.norm[, "LMNA"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["LMNA", ]*mult.factor["LMNA"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"LMNA" ]*mult.factor6["LMNA"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "Dropseq", "SCend"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.25), 
     xlim = c(-2, 100), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 100, by = 20))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("CCNA2", side = 3, line = 0, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "CCNA2"], na.rm = TRUE)$bw
lines(density(fish.norm[, "CCNA2"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["CCNA2", ]*mult.factor["CCNA2"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"CCNA2" ]*mult.factor6["CCNA2"], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()


pdf("plots/smurf2.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.7), 
     xlim = c(-2, 40), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 40, by = 8))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("BABAM1", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "BABAM1"], na.rm = TRUE)$bw
lines(density(fish.norm[, "BABAM1"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["BABAM1", ]*mult.factor["BABAM1"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"BABAM1" ]*mult.factor6["BABAM1"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "Dropseq", "SCend"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.06), 
     xlim = c(-10, 400), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 400, by = 80))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("SOX10", side = 3, line = 0, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "SOX10"], na.rm = TRUE)$bw
lines(density(fish.norm[, "SOX10"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["SOX10", ]*mult.factor["SOX10"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"SOX10" ]*mult.factor6["SOX10"], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()


pdf("plots/smurf3.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.5), 
     xlim = c(-2, 40), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 40, by = 8))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("KDM5A", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "KDM5A"], na.rm = TRUE)$bw
lines(density(fish.norm[, "KDM5A"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["KDM5A", ]*mult.factor["KDM5A"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"KDM5A" ]*mult.factor6["KDM5A"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "Dropseq", "SCend"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.7), 
     xlim = c(-2, 40), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 40, by = 8))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("KDM5B", side = 3, line = 0, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "KDM5B"], na.rm = TRUE)$bw
lines(density(fish.norm[, "KDM5B"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["KDM5B", ]*mult.factor["KDM5B"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"KDM5B" ]*mult.factor6["KDM5B"], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()

pdf("plots/smurf4.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.05), 
     xlim = c(-2, 400), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 400, by = 80))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("C1S", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "C1S"], na.rm = TRUE)$bw
lines(density(fish.norm[, "C1S"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["C1S", ]*mult.factor["C1S"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"C1S"]*mult.factor6["C1S"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "Dropseq", "SCend"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.5), 
     xlim = c(-2, 40), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 40, by = 8))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("FGFR1", side = 3, line = 0, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "FGFR1"], na.rm = TRUE)$bw
lines(density(fish.norm[, "FGFR1"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["FGFR1", ]*mult.factor["FGFR1"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"FGFR1" ]*mult.factor6["FGFR1"], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()


pdf("plots/smurf5.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 5), 
     xlim = c(-0.2, 3), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 3, by = 0.6))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("FOSL1", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "FOSL1"], na.rm = TRUE)$bw
lines(density(fish.norm[, "FOSL1"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["FOSL1", ]*mult.factor["FOSL1"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"FOSL1"]*mult.factor6["FOSL1"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "Dropseq", "SCend"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 1), 
     xlim = c(-2, 30), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 30, by = 6))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("JUN", side = 3, line = 0, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "JUN"], na.rm = TRUE)$bw
lines(density(fish.norm[, "JUN"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["JUN", ]*mult.factor["JUN"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"JUN" ]*mult.factor6["JUN"], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()

pdf("plots/smurf6.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.1), 
     xlim = c(-2, 200), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 200, by = 40))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("MITF", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "MITF"], na.rm = TRUE)$bw
lines(density(fish.norm[, "MITF"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["MITF", ]*mult.factor["MITF"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"MITF"]*mult.factor6["MITF"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "Dropseq", "SCend"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 3), 
     xlim = c(-0.6, 8), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 8, by = 1.6))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("RUNX2", side = 3, line = 0, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "RUNX2"], na.rm = TRUE)$bw
lines(density(fish.norm[, "RUNX2"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["RUNX2", ]*mult.factor["RUNX2"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"RUNX2" ]*mult.factor6["RUNX2"], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()


pdf("plots/smurf7.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.2), 
     xlim = c(-2, 200), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 200, by = 40))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("TXNRD1", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "TXNRD1"], na.rm = TRUE)$bw
lines(density(fish.norm[, "TXNRD1"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["TXNRD1", ]*mult.factor["TXNRD1"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"TXNRD1"]*mult.factor6["TXNRD1"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "Dropseq", "SCend"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.1), 
     xlim = c(-2, 300), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 300, by = 60))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("VCL", side = 3, line = 0, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "VCL"], na.rm = TRUE)$bw
lines(density(fish.norm[, "VCL"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["VCL", ]*mult.factor["VCL"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"VCL" ]*mult.factor6["VCL"], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()

pdf("plots/smurf8.pdf", 5.5, 4)
par(mfrow = c(1, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.3), 
     xlim = c(-2, 50), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 50, by = 10))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("VGF", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "VGF"], na.rm = TRUE)$bw
lines(density(fish.norm[, "VGF"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["VGF", ]*mult.factor["VGF"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"VGF"]*mult.factor6["VGF"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "Dropseq", "SCend"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)
dev.off()

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.1), 
     xlim = c(-2, 300), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 300, by = 60))
par(las = 1)
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("VCL", side = 3, line = 0, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, "VCL"], na.rm = TRUE)$bw
lines(density(fish.norm[, "VCL"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["VCL", ]*mult.factor["VCL"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"VCL" ]*mult.factor6["VCL"], bw = dens.bw), lwd = 3, lty = 2)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dev.off()


pdf("plot/fish_dropseq_smurf.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.05), 
     xlim = c(-2, 400), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 400, by = 80))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("C1S", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "C1S"], na.rm = TRUE)$bw
lines(density(fish.norm[, "C1S"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["C1S", ]*mult.factor["C1S"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"C1S"]*mult.factor6["C1S"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "Dropseq", "SCend"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.5), 
     xlim = c(-2, 40), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 40, by = 8))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("KDM5A", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "KDM5A"], na.rm = TRUE)$bw
lines(density(fish.norm[, "KDM5A"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(dropseq.norm["KDM5A", ]*mult.factor["KDM5A"], bw = dens.bw), lwd = 3)
lines(density(smurf.norm[,"KDM5A" ]*mult.factor6["KDM5A"], bw = dens.bw), lwd = 3, lty = 2)
dev.off()


pdf("plot/fish_magic_knnsmooth.pdf", 5.5, 8)
par(mfrow = c(2, 1), cex.main = 1.5, mar = c(3, 2, 4, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.05), 
     xlim = c(-2, 400), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 400, by = 80))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("C1S", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "C1S"], na.rm = TRUE)$bw
lines(density(fish.norm[, "C1S"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(magic.norm["C1S", ]*mult.factor3["C1S"], bw = dens.bw), lwd = 3)
lines(density(knnsmooth.norm[,"C1S"]*mult.factor5["C1S"], bw = dens.bw), lwd = 3, lty = 2)
legend("topright", c("RNA FISH", "MAGIC", "KNN-Smoothing"), lty = c(1, 1, 2),
       lwd = c(4, 3, 3), col = c("#bdbdbd", "black", "black"), box.lty = 0,
       cex = 1.1)

plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.5), 
     xlim = c(-2, 40), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 40, by = 8))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext("KDM5A", side = 3, line = 0.5, cex = 1.5, font = 3)
mtext("Expression Count", side = 1, line = 3, cex = 1.5)
dens.bw <- density(fish.norm[, "KDM5A"], na.rm = TRUE)$bw
lines(density(fish.norm[, "KDM5A"], na.rm = TRUE), lwd = 4, col = "#bdbdbd")
lines(density(magic.norm["KDM5A", ]*mult.factor3["KDM5A"], bw = dens.bw), lwd = 3)
lines(density(knnsmooth.norm[,"KDM5A" ]*mult.factor5["KDM5A"], bw = dens.bw), lwd = 3, lty = 2)
dev.off()
###############################################################################
## BABAM1 and LMNA correlation

library(reshape2)

adj.vec <- rep(0, 16)
for (i in 1:16) {
  adj.vec[i] <- 
    sqrt(var(saver$estimate[i, ]*sf, na.rm = TRUE)/
           (var(saver$estimate[i, ]*sf, na.rm = TRUE) + 
              mean(saver$se[i, ]^2*sf^2, na.rm = TRUE)))
}
scale.factor <- outer(adj.vec, adj.vec)


fish.cor <- cor(fish[, gene.names[genes]], use = "pairwise.complete.obs")
knnsmooth.cor <- cor(knnsmooth[, gene.names[genes]], use = "pairwise.complete.obs")
smurf.cor <- cor(smurf[, gene.names[genes]], use = "pairwise.complete.obs") 
dropseq.cor <- cor(t(dropseq[genes, ]))
saver.cor <- cor(t(sweep(saver$estimate, 2, sf, "*")))*scale.factor
magic.cor <- cor(t(sweep(magic, 2, sf, "*")))
scimpute.cor <- cor(t(scimpute))

fish.cor[lower.tri(fish.cor, diag = TRUE)] <- NA
fish.cor.melt <- melt(fish.cor, na.rm = TRUE, value.name = "FISH")
knnsmooth.cor[lower.tri(knnsmooth.cor, diag = TRUE)] <- NA
knnsmooth.cor.melt <- melt(knnsmooth.cor, na.rm = TRUE, value.name = "KNN-Smoothing")
smurf.cor[lower.tri(smurf.cor, diag = TRUE)] <- NA
smurf.cor.melt <- melt(smurf.cor, na.rm = TRUE, value.name = "SCend")
dropseq.cor[lower.tri(dropseq.cor, diag = TRUE)] <- NA
dropseq.cor.melt <- melt(dropseq.cor, na.rm = TRUE, value.name = "Dropseq")
saver.cor[lower.tri(saver.cor, diag = TRUE)] <- NA
saver.cor.melt <- melt(saver.cor, na.rm= TRUE, value.name = "SAVER")
magic.cor[lower.tri(magic.cor, diag = TRUE)] <- NA
magic.cor.melt <- melt(magic.cor, na.rm= TRUE, value.name = "MAGIC")
scimpute.cor[lower.tri(scimpute.cor, diag = TRUE)] <- NA
scimpute.cor.melt <- melt(scimpute.cor, na.rm= TRUE, value.name = "scImpute")

unnorm.cor <- Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2")), 
                     list(fish.cor.melt, dropseq.cor.melt, saver.cor.melt,
                          magic.cor.melt, scimpute.cor.melt))
unnorm.cor <- Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2")), 
                     list(fish.cor.melt, dropseq.cor.melt, saver.cor.melt,
                          magic.cor.melt, knnsmooth.cor.melt, smurf.cor.melt))

png("plot/fig1d_corplot.png", 8.5, 3.5, units = "in", res = 300)

par(mfrow = c(1, 3), cex.main = 1.5, mar = c(4, 4, 3, 0) + 0.1, oma = c(3, 4, 3, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
smoothScatter(fish[, "BABAM1"], fish[, "LMNA"], ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "", nrpoints = 0, bandwidth = c(1, 20))
axis(1)
axis(2)
par(las = 0)
mtext("LMNA", side = 2, line = 4.5, cex = 1.3, font = 3)
mtext("FISH", side = 3, cex = 1.3, line = 2)
text(52, 1150, paste("r = ", round(unnorm.cor[3, "FISH"], 2)),
     cex = 1.8)

smoothScatter(dropseq["BABAM1", ], bandwidth = c(0.1, 0.2),
              dropseq["LMNA", ], nrpoints = 0,
              ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "")
axis(1)
par(las = 1)
axis(2)
mtext("Drop-seq", side = 3, cex = 1.3, line = 2)
mtext("BABAM1", side = 1, cex = 1.3, font = 3, line = 4)
text(2.4, 15.8, paste("r = ", round(unnorm.cor[3, "Dropseq"], 2)),
     cex = 1.8)

library(scales)

babam1 <- matrix(0, ncol(saver$estimate), 5)
lmna <- matrix(0, ncol(saver$estimate), 5)
for (i in 1:ncol(saver$estimate)) {
  babam1[i, ] <- rgamma(5, saver$estimate["BABAM1", i]^2/saver$se["BABAM1", i]^2, 
                        saver$estimate["BABAM1", i]/saver$se["BABAM1", i]^2/sf[i])
  lmna[i, ] <- rgamma(5, saver$estimate["LMNA", i]^2/saver$se["LMNA", i]^2, 
                      saver$estimate["LMNA", i]/saver$se["LMNA", i]^2/sf[i])
}

smoothScatter(c(babam1), c(lmna), nrpoints = 0, bandwidth = c(0.08, 0.6), ylim = c(0, 15),
              xlim = c(0, 2), ylab = "", 
              xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "")
points(saver$estimate["BABAM1", ]*sf, saver$estimate["LMNA", ]*sf, pch = 20, cex = 0.8)

axis(1)
par(las = 1)
axis(2)
mtext("SAVER", side = 3, cex = 1.3, line = 2)
text(1.6, 14.8, paste("r = ", round(unnorm.cor[3, "SAVER"], 2)),
     cex = 1.8)

dev.off()

png("plots/suppfig3e_smooth.png", 8.5, 3.5, units = "in", res = 300)
par(mfrow = c(1, 3), cex.main = 1.5, mar = c(4, 4, 3, 0) + 0.1, oma = c(3, 4, 3, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
smoothScatter(fish[, "BABAM1"], fish[, "LMNA"], ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "", nrpoints = 0, bandwidth = c(1, 20))
axis(1)
axis(2)
par(las = 0)
mtext("LMNA", side = 2, line = 4.5, cex = 1.3, font = 3)
mtext("FISH", side = 3, cex = 1.3, line = 2)
text(52, 1150, paste("r = ", round(unnorm.cor[3, "FISH"], 2)),
     cex = 1.8)

smoothScatter(magic["BABAM1", ]*sf, magic["LMNA", ]*sf, ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "", nrpoints = 0, xlim = c(0, 20),
              ylim = c(0, 35))
axis(1)
axis(2)
par(las = 0)
mtext("MAGIC", side = 3, cex = 1.3, line = 2)
mtext("BABAM1", side = 1, cex = 1.3, font = 3, line = 4)
text(17, 34.3, paste("r = ", round(unnorm.cor[3, "MAGIC"], 2)),
     cex = 1.8)

smoothScatter(scimpute["BABAM1", ], scimpute["LMNA", ], ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "", nrpoints = 0, bandwidth = c(0.05, 0.2))
axis(1)
axis(2)
par(las = 0)
mtext("scImpute", side = 3, cex = 1.3, line = 2)
text(2.4, 15.6, paste("r = ", round(unnorm.cor[3, "scImpute"], 2)),
     cex = 1.8)
dev.off()



png("plot/fig1d_corplot.png", 8.5, 3.5, units = "in", res = 300)

par(mfrow = c(1, 3), cex.main = 1.5, mar = c(4, 4, 3, 0) + 0.1, oma = c(3, 4, 3, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
smoothScatter(fish[, "BABAM1"], fish[, "LMNA"], ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "", nrpoints = 0, bandwidth = c(1, 20))
axis(1)
axis(2)
par(las = 0)
mtext("LMNA", side = 2, line = 4.5, cex = 1.3, font = 3)
mtext("FISH", side = 3, cex = 1.3, line = 2)
text(52, 1150, paste("r = ", round(unnorm.cor[3, "FISH"], 2)),
     cex = 1.8)

smoothScatter(dropseq["BABAM1", ], bandwidth = c(0.1, 0.2),
              dropseq["LMNA", ], nrpoints = 0,
              ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "")
axis(1)
par(las = 1)
axis(2)
mtext("Drop-seq", side = 3, cex = 1.3, line = 2)
mtext("BABAM1", side = 1, cex = 1.3, font = 3, line = 4)
text(2.4, 15.8, paste("r = ", round(unnorm.cor[3, "Dropseq"], 2)),
     cex = 1.8)

library(scales)

babam1 <- matrix(0, ncol(saver$estimate), 5)
lmna <- matrix(0, ncol(saver$estimate), 5)
for (i in 1:ncol(saver$estimate)) {
  babam1[i, ] <- rgamma(5, saver$estimate["BABAM1", i]^2/saver$se["BABAM1", i]^2, 
                        saver$estimate["BABAM1", i]/saver$se["BABAM1", i]^2/sf[i])
  lmna[i, ] <- rgamma(5, saver$estimate["LMNA", i]^2/saver$se["LMNA", i]^2, 
                      saver$estimate["LMNA", i]/saver$se["LMNA", i]^2/sf[i])
}

smoothScatter(smurf[, "BABAM1"],smurf[, "LMNA"], ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "", nrpoints = 0, bandwidth = c(1, 20))

axis(1)
par(las = 1)
axis(2)
mtext("SAVER", side = 3, cex = 1.3, line = 2)
text(1.6, 14.8, paste("r = ", round(unnorm.cor[3, "SCend"], 2)),
     cex = 1.8)

dev.off()

png("plot/suppfig3e_smooth.png", 8.5, 3.5, units = "in", res = 300)
par(mfrow = c(1, 3), cex.main = 1.5, mar = c(4, 4, 3, 0) + 0.1, oma = c(3, 4, 3, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
smoothScatter(fish[, "BABAM1"], fish[, "LMNA"], ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "", nrpoints = 0, bandwidth = c(1, 20))
axis(1)
axis(2)
par(las = 0)
mtext("LMNA", side = 2, line = 4.5, cex = 1.3, font = 3)
mtext("FISH", side = 3, cex = 1.3, line = 2)
text(52, 1150, paste("r = ", round(unnorm.cor[3, "FISH"], 2)),
     cex = 1.8)

smoothScatter(magic["BABAM1", ]*sf, magic["LMNA", ]*sf, ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "", nrpoints = 0, xlim = c(0, 20),
              ylim = c(0, 35))
axis(1)
axis(2)
par(las = 0)
mtext("MAGIC", side = 3, cex = 1.3, line = 2)
mtext("BABAM1", side = 1, cex = 1.3, font = 3, line = 4)
text(17, 34.3, paste("r = ", round(unnorm.cor[3, "MAGIC"], 2)),
     cex = 1.8)

smoothScatter(scimpute["BABAM1", ], scimpute["LMNA", ], ylab = "", xlab = "", cex = 1.5,
              lwd = 2, axes = FALSE, main = "", nrpoints = 0, bandwidth = c(0.05, 0.2))
axis(1)
axis(2)
par(las = 0)
mtext("scImpute", side = 3, cex = 1.3, line = 2)
text(2.4, 15.6, paste("r = ", round(unnorm.cor[3, "scImpute"], 2)),
     cex = 1.8)
dev.off()

###############################################################################
## K-S

dropseq.ks <- sapply(1:16, function(i) 
  ks.test(dropseq.norm[i, ]*mult.factor[i], 
          fish.norm[, gene.names[genes[i]]])$statistic)

Y.ks <- sapply(1:16, function(i)
  ks.test(Y.samp.norm[i, ], fish.norm[, gene.names[genes[i]]])$statistic)


magic.ks <- sapply(1:16, function(i) ks.test(magic.norm[i, ]*mult.factor3[i],
                                             fish.norm[, gene.names[genes[i]]])$statistic)

scimpute.ks <- sapply(1:16, function(i) ks.test(scimpute.norm[i, ]*mult.factor4[i],
                                                fish.norm[, gene.names[genes[i]]])$statistic)

ks <- cbind(dropseq.ks, Y.ks)[-6, ]
rownames(ks) <- gene.names[genes[c(1:5, 7:16)]]


pdf("plots/suppfig3b_ks.pdf", 7, 4.5)
par(cex.main = 1.5, mar = c(5, 6, 4, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
x <- c(1.5, 2.5, 3.5, 4.5)
plot(x, c(-10, -10, -10, -10), type = "p", ylab = " ", xlab = " ", cex = 1.5, 
     ylim = c(0, 1), xlim = c(1, 5), lwd = 2, pch = 5, axes = FALSE, main = " ")
axis(1, at = c(1.5, 2.5, 3.5, 4.5), labels = c("Drop-seq", "MAGIC", "scImpute", 
                                               "SAVER"), cex.axis = 1.5)
par(las = 1)
axis(2, pos = 1.3)
par(las = 0)
mtext("Kolmogorov-Smirnov", side = 2, line = 0.5, cex = 1.5)

for (i in 1:length(dropseq.ks)) {
  if (i == 6) next
  lines(x, c(dropseq.ks[i], magic.ks[i], scimpute.ks[i], Y.ks[i]), 
        lwd = 1.5, type = "c")
}
points(rep(1.55, length(dropseq.ks)-1), dropseq.ks[-6], pch = 21, 
       bg = "#9ecae1", cex = 3)
points(rep(2.53, length(dropseq.ks)-1), magic.ks[-6], pch = 21, 
       bg = "#9ecae1", cex = 3)
points(rep(3.53, length(dropseq.ks)-1), scimpute.ks[-6], pch = 21, 
       bg = "#9ecae1", cex = 3)
points(rep(4.5, length(dropseq.ks)-1), Y.ks[-6], pch = 21, 
       bg = "#9ecae1", cex = 3)

dev.off()

###############################################################################
## Correlation

pdf("plots/suppfig3d_cor.pdf", 7.5, 4)
par(mfrow = c(2, 2), cex.main = 1.5, mar = c(0.5, 4, 0, 2) + 0.1, oma = c(4, 2, 1, 0), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(-0.2, 1), 
     xlim = c(-0.2, 1), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, labels = FALSE)
axis(2)
par(las = 0)
mtext("Drop-seq Cor", side = 2, line = 3.5, cex = 1.5)

abline(0, 1, lty = 2, col = "gray")
points(unnorm.cor[, "FISH"], unnorm.cor[, "Dropseq"], pch = 21, bg = "#9ecae1",
       cex = 3)
points(unnorm.cor[unnorm.cor$Var1 == "BABAM1" & unnorm.cor$Var2 == "LMNA", 
                  "FISH"],
       unnorm.cor[unnorm.cor$Var1 == "BABAM1" & unnorm.cor$Var2 == "LMNA", 
                  "Dropseq"], pch = 21, cex = 3,
       bg = "#e6550d")
text(0.93, 0.21, "BABAM1", xpd = TRUE, cex = 1.2)
text(0.93, 0.13, "LMNA", xpd = TRUE, cex = 1.2)


plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(-0.2, 1), 
     xlim = c(-0.2, 1), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, labels = FALSE)
par(las = 1)
axis(2)
par(las = 0)
mtext("SAVER Cor", side = 2, cex = 1.5, line = 3.5)
abline(0, 1, lty = 2, col = "gray")
points(unnorm.cor[, "FISH"], unnorm.cor[, "SAVER"], pch = 21, bg = "#9ecae1",
       cex = 3)
points(unnorm.cor[unnorm.cor$Var1 == "BABAM1" & unnorm.cor$Var2 == "LMNA", 
                  "FISH"],
       unnorm.cor[unnorm.cor$Var1 == "BABAM1" & unnorm.cor$Var2 == "LMNA", 
                  "SAVER"], pch = 21, cex = 3,
       bg = "#e6550d")
text(0.93, 0.74, "BABAM1", xpd = TRUE, cex = 1.2)
text(0.93, 0.66, "LMNA", xpd = TRUE, cex = 1.2)



plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(-0.2, 1), 
     xlim = c(-0.2, 1), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1)
par(las = 1)
axis(2)
par(las = 0)
mtext("MAGIC Cor", side = 2, line = 3.5, cex = 1.5)
mtext("FISH Cor", side = 1, line = 3, cex = 1.5)

abline(0, 1, lty = 2, col = "gray")
points(unnorm.cor[, "FISH"], unnorm.cor[, "MAGIC"], pch = 21, bg = "#9ecae1",
       cex = 3)
points(unnorm.cor[unnorm.cor$Var1 == "BABAM1" & unnorm.cor$Var2 == "LMNA", 
                  "FISH"],
       unnorm.cor[unnorm.cor$Var1 == "BABAM1" & unnorm.cor$Var2 == "LMNA", 
                  "MAGIC"], pch = 21, cex = 3,
       bg = "#e6550d")
text(0.93, 0.98, "BABAM1", xpd = TRUE, cex = 1.2)
text(0.93, 0.9, "LMNA", xpd = TRUE, cex = 1.2)


plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(-0.2, 1), 
     xlim = c(-0.2, 1), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1)
par(las = 1)
axis(2)
par(las = 0)
mtext("scImpute Cor", side = 2, cex = 1.5, line = 3.5)
mtext("FISH Cor", side = 1, line = 3, cex = 1.5)
abline(0, 1, lty = 2, col = "gray")
points(unnorm.cor[, "FISH"], unnorm.cor[, "scImpute"], pch = 21, bg = "#9ecae1",
       cex = 3)
points(unnorm.cor[unnorm.cor$Var1 == "BABAM1" & unnorm.cor$Var2 == "LMNA", 
                  "FISH"],
       unnorm.cor[unnorm.cor$Var1 == "BABAM1" & unnorm.cor$Var2 == "LMNA", 
                  "scImpute"], pch = 21, cex = 3,
       bg = "#e6550d")
text(0.90, 0.24, "BABAM1", xpd = TRUE, cex = 1.2)
text(0.90, 0.16, "LMNA", xpd = TRUE, cex = 1.2)
dev.off()









