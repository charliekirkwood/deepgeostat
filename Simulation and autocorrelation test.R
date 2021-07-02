library(gstat)
library(sp)
library(automap)
library(geoR)

dir.create(paste0(getwd(), "/plots/simulation"))

load_model_weights_hdf5(object = model, filepath = paste0(getwd(), "/models/modelweights"))

# National scale maps ####
rm(predimgs, predimgs_ann, predloc, predloc_ann)

gbelevcoarse <- aggregate(gbelev, fact=1)
predgrid <- as.data.table(raster::as.data.frame(gbelevcoarse, xy = TRUE, na.rm = TRUE))[x > 0 & y > 0]
setnames(predgrid, "alt", "elevation")
#predgrid <- predgrid[!elevation == 0]
nrow(predgrid)

ggplot(predgrid) + geom_raster(aes(x = x, y = y, fill = elevation)) + theme_bw() + scale_fill_viridis_c(option = "B") + coord_equal()

predimgs <- array(dim = c(nrow(predgrid), imagedim, imagedim))
str(predimgs)

time <- Sys.time()
for(xmeter in unique(cells$coordx)){
  for(ymeter in unique(cells$coordy)){
    predimgs[,xmeter,ymeter] <- extract(gbelev, method = "bilinear",
                                        cbind(predgrid$x + cells[coordx == xmeter & coordy == ymeter]$x, predgrid$y + cells[coordx == xmeter & coordy == ymeter]$y)) 
  }
}
Sys.time() - time

predimgs[is.na(predimgs)] = 0

predimgs_ann <- predimgs - predgrid$elevation
predimgs_ann <- (predimgs_ann)/sd
predimgs_ann[is.na(predimgs_ann)] = 0
#predimgs_ann <- predimgs_ann - rowMeans(predimgs_ann, na.rm = TRUE)

dim(predimgs_ann) <- c(nrow(predgrid), imagedim, imagedim, 1)

# Location information:
predloc <- predgrid[, c("x", "y", "elevation"), with = FALSE]
predlocmean <- apply(loc, 2, mean)
predlocsd <- apply(loc, 2, sd)
predloc_ann <- t(apply(predloc, 1, function(x) (x - predlocmean)/predlocsd))

nsim <- 9
dropratespat <- 0.5 #match what model was trained with
dropratedense <- 0.33 #match what model was trained with

i <- 1
for(i in 1:nsim){
  load_model_weights_hdf5(object = model, filepath = paste0(getwd(), "/models/modelweights"))
  weights <- model %>% get_weights()
  
  dropweights <- weights
  #str(dropweights)
  
  # spatial dropout masks
  maskl1 <- array(rep(rbinom(128*9, 1, 1-dropratespat), each = 9), dim = dim(dropweights[[1]])) # generate dropout mask for kernel and bias of layer 1
  maskl2 <- array(rep(rbinom(128*9, 1, 1-dropratespat), each = 9), dim = dim(dropweights[[3]])) # generate dropout mask for kernel and bias of layer 2
  maskl3 <- array(rep(rbinom(128*9, 1, 1-dropratespat), each = 9), dim = dim(dropweights[[5]])) # generate dropout mask for kernel and bias of layer 3
  maskl4 <- array(rep(rbinom(128*9, 1, 1-dropratespat), each = 9), dim = dim(dropweights[[7]])) # generate dropout mask for kernel and bias of layer 4
  
  # normal dropout masks
  maskl5 <- rbinom(dim(dropweights[[9]])[2], 1, 1-dropratedense) # generate dropout mask for kernel and bias of layer 5
  maskl6 <- rbinom(dim(dropweights[[11]])[2], 1, 1-dropratedense) # generate dropout mask for kernel and bias of layer 6
  maskl7 <- rbinom(dim(dropweights[[13]])[2], 1, 1-dropratedense) # generate dropout mask for kernel and bias of layer 7
  
  dropweights[[1]] <- (dropweights[[1]] * maskl1)*(1/(1-dropratespat))
  dropweights[[2]] <- (dropweights[[2]] * as.numeric(maskl1)[seq(1, length(as.numeric(maskl1)), 9)])*(1/(1-dropratespat))
  
  dropweights[[3]] <- (dropweights[[3]] * maskl2)*(1/(1-dropratespat))
  dropweights[[4]] <- (dropweights[[4]] * as.numeric(maskl1)[seq(1, length(as.numeric(maskl1)), 9)])*(1/(1-dropratespat))
  
  dropweights[[5]] <- (dropweights[[5]] * maskl3)*(1/(1-dropratespat))
  dropweights[[6]] <- (dropweights[[6]] * as.numeric(maskl1)[seq(1, length(as.numeric(maskl1)), 9)])*(1/(1-dropratespat))
  
  dropweights[[7]] <- (dropweights[[7]] * maskl4)*(1/(1-dropratespat))
  dropweights[[8]] <- (dropweights[[8]] * as.numeric(maskl1)[seq(1, length(as.numeric(maskl1)), 9)])*(1/(1-dropratespat))
  
  dropweights[[9]] <- t(t(dropweights[[9]]) * maskl5)*(1/(1-dropratedense))
  dropweights[[10]] <- (dropweights[[10]] * maskl5)*(1/(1-dropratedense))
  
  dropweights[[11]] <- t(t(dropweights[[11]]) * maskl6)*(1/(1-dropratedense))
  dropweights[[12]] <- (dropweights[[12]] * maskl6)*(1/(1-dropratedense))
  
  dropweights[[13]] <- t(t(dropweights[[13]]) * maskl7)*(1/(1-dropratedense))
  dropweights[[14]] <- (dropweights[[14]] * maskl7)*(1/(1-dropratedense))
  
  model %>% set_weights(dropweights)
  
  preds <- predgrid
  #preds[, simpred := as.numeric(predict(robust_model, list(predimgs_ann_tf, predloc_ann_tf), batch_size = 4096)), ] # predict with aleatoric uncertainty from robust_model
  preds[, simpred := as.numeric(predict(model, list(predimgs_ann, predloc_ann), batch_size = 4096)[,1]), ] # predict with aleatoric uncertainty from model
  preds[, run := i]
  preds
  
  
  if(i == 1){
    sim <- copy(preds)
  }else{
    sim <- rbindlist(list(sim, preds), use.names = FALSE, fill = FALSE, idcol = NULL)
  }
}

sim[, post_mean := mean(simpred), by = c("x", "y")]

postdat <- sim[run == sample(1:nsim, 1)]

ggplot(postdat) + geom_raster(aes(x = x, y = y, fill = simpred), interpolate = FALSE) + theme_bw(8) + coord_equal() +
  scale_fill_viridis_c(option = "B", name = paste("simulated\nlog(CaO)"), 
                       limits = c(-7.5, +7.5), oob = scales::squish) + 
  labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.background=element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), labels = paste(seq(0,1000000, length.out = 6)), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 670000))
ggsave(paste0("plots/simulation/", elem, "_national_simulation_qtr_", unique(postdat$run), ".png"), width = 45, height = 72.5, units = "mm", type = "cairo", dpi = 300, scale = 1.375)


# Analyse spatial-autocorrelation:
# Use held-out test data (compare autocorrelation of test set with autocorrelation of simulations at the same locations)

nsims <- 25
dropratespat <- 0.5 #match what model was trained with
dropratedense <- 0.33 #match what model was trained with
distbreaks <- 1500

#plot(x_test[[2]][,1:2])

spatdat <- as.data.table(cbind(loc[test], y_test))

variogramdat <- variog(coords = data.frame(spatdat[, c("Bng_easting", "Bng_northing"), with = FALSE]), data = spatdat$y_test, 
                       max.dist = 102000, breaks = seq(0,102000, distbreaks))
variogramdat

datvar <- as.data.table(data.frame(dist = variogramdat$u, gamma = variogramdat$v))
datvar[, dist := dist/1000, ]
datvar[, type := "observations", ]
datvar[, rep := 1, ]
datvar

ggplot(datvar) + 
  geom_point(aes(x = dist, y = gamma, group = interaction(type, rep), col = type), alpha = 0.5, shape = 16) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  theme_bw() + scale_y_continuous(limits=c(0, max(datvar$gamma)*1.2), expand=c(0,0)) +
  scale_x_continuous(limits=c(0, 302), expand=c(0,0)) +
  labs(x = "distance (km)", y = "semivariance") +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background=element_blank(), legend.title=element_blank(), 
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5))

rm(variogramdata)
time <- Sys.time()
for(i in 1:nsims){
  load_model_weights_hdf5(object = model, filepath = paste0(getwd(), "/models/modelweights"))
  weights <- model %>% get_weights()
  
  dropweights <- weights
  #str(dropweights)
  
  # spatial dropout masks
  maskl1 <- array(rep(rbinom(128*9, 1, 1-dropratespat), each = 9), dim = dim(dropweights[[1]])) # generate dropout mask for kernel and bias of layer 1
  maskl2 <- array(rep(rbinom(128*9, 1, 1-dropratespat), each = 9), dim = dim(dropweights[[3]])) # generate dropout mask for kernel and bias of layer 2
  maskl3 <- array(rep(rbinom(128*9, 1, 1-dropratespat), each = 9), dim = dim(dropweights[[5]])) # generate dropout mask for kernel and bias of layer 3
  maskl4 <- array(rep(rbinom(128*9, 1, 1-dropratespat), each = 9), dim = dim(dropweights[[7]])) # generate dropout mask for kernel and bias of layer 4
  
  # normal dropout masks
  maskl5 <- rbinom(dim(dropweights[[9]])[2], 1, 1-dropratedense) # generate dropout mask for kernel and bias of layer 5
  maskl6 <- rbinom(dim(dropweights[[11]])[2], 1, 1-dropratedense) # generate dropout mask for kernel and bias of layer 6
  maskl7 <- rbinom(dim(dropweights[[13]])[2], 1, 1-dropratedense) # generate dropout mask for kernel and bias of layer 7
  
  dropweights[[1]] <- (dropweights[[1]] * maskl1)*(1/(1-dropratespat))
  dropweights[[2]] <- (dropweights[[2]] * as.numeric(maskl1)[seq(1, length(as.numeric(maskl1)), 9)])*(1/(1-dropratespat))
  
  dropweights[[3]] <- (dropweights[[3]] * maskl2)*(1/(1-dropratespat))
  dropweights[[4]] <- (dropweights[[4]] * as.numeric(maskl1)[seq(1, length(as.numeric(maskl1)), 9)])*(1/(1-dropratespat))
  
  dropweights[[5]] <- (dropweights[[5]] * maskl3)*(1/(1-dropratespat))
  dropweights[[6]] <- (dropweights[[6]] * as.numeric(maskl1)[seq(1, length(as.numeric(maskl1)), 9)])*(1/(1-dropratespat))
  
  dropweights[[7]] <- (dropweights[[7]] * maskl4)*(1/(1-dropratespat))
  dropweights[[8]] <- (dropweights[[8]] * as.numeric(maskl1)[seq(1, length(as.numeric(maskl1)), 9)])*(1/(1-dropratespat))
  
  dropweights[[9]] <- t(t(dropweights[[9]]) * maskl5)*(1/(1-dropratedense))
  dropweights[[10]] <- (dropweights[[10]] * maskl5)*(1/(1-dropratedense))
  
  dropweights[[11]] <- t(t(dropweights[[11]]) * maskl6)*(1/(1-dropratedense))
  dropweights[[12]] <- (dropweights[[12]] * maskl6)*(1/(1-dropratedense))
  
  dropweights[[13]] <- t(t(dropweights[[13]]) * maskl7)*(1/(1-dropratedense))
  dropweights[[14]] <- (dropweights[[14]] * maskl7)*(1/(1-dropratedense))
  
  model %>% set_weights(dropweights)
  
  simdat <- as.data.table(cbind(loc[test], y_test))
  #preds[, simpred := as.numeric(predict(robust_model, list(predimgs_ann_tf, predloc_ann_tf), batch_size = 4096)), ] # predict with aleatoric uncertainty from robust_model
  simdat[, simpred := as.numeric(predict(model, x_test, batch_size = 4096)[,1]), ] # predict with aleatoric uncertainty from model
  
  variogramsim <- variog(coords = data.frame(spatdat[, c("Bng_easting", "Bng_northing"), with = FALSE]), data = simdat$simpred, 
                         max.dist = 102000, breaks = seq(0,102000, distbreaks))
  
  simvar <- as.data.table(data.frame(dist = variogramsim$u, gamma = variogramsim$v))
  simvar[, dist := dist/1000, ]
  simvar[, type := "simulations", ]
  simvar[, rep := i, ]
  simvar
  
  if(exists("variogramdata")){
    variogramdata <- rbind(variogramdata, simvar)
  }else{
    variogramdata <- rbind(datvar, simvar)
  }
}
Sys.time() - time

ggplot() + 
  #geom_point(aes(x = dist, y = gamma, group = interaction(type, run), col = type), alpha = 0.5, shape = 16) +
  geom_line(data = variogramdata[type == "simulations"], aes(x = dist, y = gamma, group = interaction(type, rep), col = rep), alpha = 0.33, size = 0.5) +
  geom_line(data = variogramdata[type == "observations"], aes(x = dist, y = gamma, group = interaction(type, rep)), col = "#e41a1c", size = 1) +
  #scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  theme_bw() + scale_y_continuous(limits=c(0, max(variogramdata$gamma)*1.2), expand=c(0,0)) +
  scale_x_continuous(limits=c(0, 102), expand=c(0,0)) +
  labs(x = "distance (km)", y = "semivariance") +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background=element_blank(), 
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(paste0("plots/simulation/autocorrelation_test_6.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

variogramdata[dist == min(dist)]

ggplot() + stat_density(data = variogramdata[type == "simulations"][dist == min(dist)], aes(x = gamma), col = "darkblue", geom = "line", size = 1.5) +
  geom_vline(data = variogramdata[type == "observations"][dist == min(dist)], aes(xintercept = gamma), col = "red") + theme_bw()

plot(ecdf(variogramdata[type == "simulations"][dist == min(dist)]$gamma))
ecdf(variogramdata[type == "simulations"][dist == min(dist)]$gamma)(variogramdata[type == "observations"][dist == min(dist)]$gamma)

ecdf(variogramdata[type == "simulations"][dist == min(dist)]$gamma)(0.27)
quantile(variogramdata[type == "simulations"][dist == min(dist)]$gamma, 0.05)

ggplot() + 
  #geom_point(aes(x = dist, y = gamma, group = interaction(type, run), col = type), alpha = 0.5, shape = 16) +
  geom_line(data = variogramdata[type == "simulations"], aes(x = dist, y = gamma, group = interaction(type, run), col = run), alpha = 0.33, size = 1) +
  geom_ribbon(data = variogramdata[type == "simulations"][, list(upper = quantile(gamma, 0.975), lower = quantile(gamma, 0.025)), by = "dist"],
              aes(x = dist, ymin = lower, ymax = upper), fill = "black", alpha = 0.33) +
  geom_line(data = variogramdata[type == "observations"], aes(x = dist, y = gamma, group = interaction(type, run)), col = "#e41a1c", size = 1.5) +
  #scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  theme_bw() + scale_y_continuous(limits=c(0, max(variogramdata$gamma)*1.2), expand=c(0,0)) +
  scale_x_continuous(limits=c(0, 302), expand=c(0,0)) +
  labs(x = "distance (km)", y = "semivariance") +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background=element_blank(), 
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(paste0("plots/simulation/autocorrelation_test_with_95.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)


