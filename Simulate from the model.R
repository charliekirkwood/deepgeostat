library(gstat)
library(sp)
library(automap)
library(geoR)

dir.create(paste0(getwd(), "/plots/simulation"))

load_model_weights_hdf5(object = model, filepath = paste0(getwd(), "/models/modelweights"))

# National scale maps ####
rm(predimgs, predimgs_ann, predloc, predloc_ann)

gbelevcoarse <- aggregate(gbelev, fact=2)
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

nsim <- 100
dropratespat <- 0.5 #match what model was trained with
dropratedense <- 0.125 #match what model was trained with

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


spatdat <- SpatialPointsDataFrame(gchem[, c("Bng_easting", "Bng_northing"), with = FALSE], data.frame(CaO = log(gchem[,get(elem)])))
spatdat$CaO <- log(spatdat$CaO)

variogramdat <- autofitVariogram(as.formula(paste(elem, "~", 1)), spatdat[sample(1:nrow(spatdat), 10000), ], model="Mat", verbose = TRUE)

ggplot(variogramdat$exp_var) + geom_point(aes(x = dist, y = gamma, shape='binned data', linetype='binned data')) + theme_bw(7) +
  geom_line(data = variogramLine(variogramdat$var_model, maxdist = variogramdat$exp_var[12,"dist"]), 
            aes(x = dist, y = gamma, shape='fitted model', linetype='fitted model')) +
  scale_y_continuous(limits=c(0, max(variogramdat$exp_var$gamma)*1.1), expand=c(0,0)) +
  scale_x_continuous(limits=c(0, max(variogramdat$exp_var$dist)), expand=c(0,0)) +
  xlab(paste("Distance (metres)")) + ylab(paste("Semivariance")) +
  scale_shape_manual(paste(i, "variogramdat"), values=c(19, NA)) + 
  scale_linetype_manual(paste(i, "variogramdat"), values=c(0, 1)) +
  theme(legend.position = "none", legend.justification=c(1,1)) +
  theme(legend.background = element_blank()) +
  annotate("text", size=as.numeric(grid::convertX(grid::unit(7, "points"), "mm")), 
           fontface="plain", y= -Inf, x =max(variogramdat$exp_var$dist)*0.83,hjust=0, vjust=-0.1,
           label=paste("Model parameters:", "\nType: Matern","\nKappa:", as.numeric(format(variogramdat$var_model$kappa[2], digits=3)),
                       "\nNugget:", as.numeric(format(variogramdat$var_model$psill[1], digits=3)), 
                       "\nSill:", as.numeric(format((variogramdat$var_model$psill[1]+variogramdat$var_model$psill[2]), digits=3)), 
                       "\nRange:", as.numeric(format(variogramdat$var_model$range[2], digits=3))))

simdat <- sim[run == sample(1:max(run), 1)][sample(1:.N, 10000)]
simdat
simdatsp <- SpatialPointsDataFrame(simdat[, c("x", "y"), with = FALSE], data.frame(CaO = simdat$simpred))

variogramsim <- autofitVariogram(as.formula(paste(elem, "~", 1)), simdatsp, model="Mat", verbose = TRUE)

ggplot(variogramsim$exp_var) + geom_point(aes(x = dist, y = gamma, shape='binned data', linetype='binned data')) + theme_bw(7) +
  geom_line(data = variogramLine(variogramsim$var_model, maxdist = variogramsim$exp_var[12,"dist"]), 
            aes(x = dist, y = gamma, shape='fitted model', linetype='fitted model')) +
  scale_y_continuous(limits=c(0, max(variogramsim$exp_var$gamma)*1.1), expand=c(0,0)) +
  scale_x_continuous(limits=c(0, max(variogramsim$exp_var$dist)), expand=c(0,0)) +
  xlab(paste("Distance (metres)")) + ylab(paste("Semivariance")) +
  scale_shape_manual(paste(i, "variogramsim"), values=c(19, NA)) + 
  scale_linetype_manual(paste(i, "variogramsim"), values=c(0, 1)) +
  theme(legend.position = "none", legend.justification=c(1,1)) +
  theme(legend.background = element_blank()) +
  annotate("text", size=as.numeric(grid::convertX(grid::unit(7, "points"), "mm")), 
           fontface="plain", y= -Inf, x =max(variogramsim$exp_var$dist)*0.83,hjust=0, vjust=-0.1,
           label=paste("Model parameters:", "\nType: Matern","\nKappa:", as.numeric(format(variogramsim$var_model$kappa[2], digits=3)),
                       "\nNugget:", as.numeric(format(variogramsim$var_model$psill[1], digits=3)), 
                       "\nSill:", as.numeric(format((variogramsim$var_model$psill[1]+variogramsim$var_model$psill[2]), digits=3)), 
                       "\nRange:", as.numeric(format(variogramsim$var_model$range[2], digits=3))))

# Both on one plot

spatdat <- SpatialPointsDataFrame(gchem[, c("Bng_easting", "Bng_northing"), with = FALSE], data.frame(CaO = log(gchem[,get(elem)])))
spatdat$CaO <- log(spatdat$CaO)

variogramdat <- autofitVariogram(as.formula(paste(elem, "~", 1)), spatdat[sample(1:nrow(spatdat), 40000), ], model="Mat", verbose = TRUE)


simdat <- sim[run == sample(1:max(run), 1)][sample(1:.N, 40000)]
simdat
simdatsp <- SpatialPointsDataFrame(simdat[, c("x", "y"), with = FALSE], data.frame(CaO = simdat$simpred))

variogramsim <- autofitVariogram(as.formula(paste(elem, "~", 1)), simdatsp, model="Mat", verbose = TRUE)

ggplot(variogramsim$exp_var) + geom_point(aes(x = dist, y = gamma, shape='binned data', linetype='binned data'), col = "orange") + theme_bw(7) +
  geom_line(data = variogramLine(variogramsim$var_model, maxdist = variogramsim$exp_var[12,"dist"]), 
            aes(x = dist, y = gamma, shape='fitted model', linetype='fitted model'), col = "orange") +
  geom_point(data = variogramdat$exp_var, aes(x = dist, y = gamma, shape='binned data', linetype='binned data'), col = "blue") +
  geom_line(data = variogramLine(variogramdat$var_model, maxdist = variogramsim$exp_var[12,"dist"]), 
            aes(x = dist, y = gamma, shape='fitted model', linetype='fitted model'), col = "blue") +
  scale_y_continuous(limits=c(0, max(variogramsim$exp_var$gamma)*1.1), expand=c(0,0)) +
  scale_x_continuous(limits=c(0, max(variogramsim$exp_var$dist)), expand=c(0,0)) +
  xlab(paste("Distance (metres)")) + ylab(paste("Semivariance")) +
  scale_shape_manual(paste(i, "variogramsim"), values=c(19, NA)) + 
  scale_linetype_manual(paste(i, "variogramsim"), values=c(0, 1)) +
  theme(legend.position = "none", legend.justification=c(1,1)) +
  theme(legend.background = element_blank()) +
  annotate("text", size=as.numeric(grid::convertX(grid::unit(7, "points"), "mm")), 
           fontface="plain", y= -Inf, x =max(variogramsim$exp_var$dist)*0.83,hjust=0, vjust=-0.1,
           label=paste("Model parameters:", "\nType: Matern","\nKappa:", as.numeric(format(variogramsim$var_model$kappa[2], digits=3)),
                       "\nNugget:", as.numeric(format(variogramsim$var_model$psill[1], digits=3)), 
                       "\nSill:", as.numeric(format((variogramsim$var_model$psill[1]+variogramsim$var_model$psill[2]), digits=3)), 
                       "\nRange:", as.numeric(format(variogramsim$var_model$range[2], digits=3))))

# Multiple runs on one plot

nvars <- 100
nsamples <- 1000

rm(variogramdata)

for(i in 1:nvars){
  spatdat <- SpatialPointsDataFrame(gchem[, c("Bng_easting", "Bng_northing"), with = FALSE], data.frame(CaO = log(gchem[,get(elem)])))
  variogramdat <- autofitVariogram(as.formula(paste(elem, "~", 1)), na.omit(spatdat)[sample(1:nrow(na.omit(spatdat)), nsamples), ], model="Exp", verbose = FALSE)
  
  datvar <- as.data.table(variogramLine(variogramdat$var_model, maxdist = variogramsim$exp_var[12,"dist"]))
  datvar[, type := "observations", ]
  datvar[, run := i, ]
  datvar
  
  simdat <- sim[run == sample(1:max(run), 1)][sample(1:.N, nsamples)]
  simdat
  simdatsp <- SpatialPointsDataFrame(simdat[, c("x", "y"), with = FALSE], data.frame(CaO = simdat$simpred))
  
  variogramsim <- autofitVariogram(as.formula(paste(elem, "~", 1)), simdatsp, model="Exp", verbose = FALSE)
  
  simvar <- as.data.table(variogramLine(variogramsim$var_model, maxdist = variogramsim$exp_var[12,"dist"]))
  simvar[, type := "simulations", ]
  simvar[, run := i, ]
  simvar
  
  if(exists("variogramdata")){
    variogramdata <- rbind(variogramdata, rbind(datvar, simvar))
  }else{
      variogramdata <- rbind(datvar, simvar)
      }
}

variogramdata

# ggplot(variogramdata) + geom_line(aes(x = dist, y = gamma, group = interaction(type, run), col = type, alpha = 1)) +
#   theme_bw() + scale_y_continuous(limits=c(0, max(variogramdata$gamma)), expand=c(0,0)) +
#   scale_x_continuous(limits=c(0, max(variogramdata$dist)), expand=c(0,0))

variogramsummary <- copy(variogramdata)
variogramsummary[, dist := round(dist/4400, 0)*4.4, ]
variogramsummary
variogramsummary <- unique(variogramsummary[, list(upper = quantile(gamma, 0.975), lower = quantile(gamma, 0.025), mean = mean(gamma)), by = c("dist", "type")])

ggplot() + 
  geom_ribbon(data = variogramsummary[type == "simulations"], aes(x = dist, ymin = lower, ymax = upper, group = type, col = type, fill = type), alpha = 0.33) +
  geom_line(data = variogramsummary[type == "simulations"], aes(x = dist, y = mean, group = type, col = type), size = 1.25, se = FALSE) +
  geom_ribbon(data = variogramsummary[type == "observations"], aes(x = dist, ymin = lower, ymax = upper, group = type, col = type, fill = type), alpha = 0.33) +
  geom_line(data = variogramsummary[type == "observations"], aes(x = dist, y = mean, group = type, col = type, style = type), size = 1.25, se = FALSE) +
  theme_bw() + scale_y_continuous(limits=c(0, max(variogramdata$gamma)), expand=c(0,0)) +
  scale_x_continuous(limits=c(0, 420), expand=c(0,0)) +
  labs(x = "distance (km)", y = "semivariance") +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background=element_blank(), legend.title=element_blank(), 
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(paste0("plots/simulation/autocorrelation_comparison.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)


# Multiple runs on one - data without model fit

nvars <- 10
nsamples <- 5000
distbreaks <- 1000

rm(variogramdata)

for(i in 1:nvars){
  #spatdat <- SpatialPointsDataFrame(gchem[, c("Bng_easting", "Bng_northing"), with = FALSE], data.frame(CaO = log(gchem[,get(elem)])))
  #variogramdat <- autofitVariogram(as.formula(paste(elem, "~", 1)), na.omit(spatdat)[sample(1:nrow(na.omit(spatdat)), nsamples), ], model="Exp", verbose = FALSE)
  
  spatdat <- copy(gchem[, c("Bng_easting", "Bng_northing", "CaO"), with = FALSE])
  spatdat[, CaO := log(CaO), ]
  spatdat <- spatdat[sample(1:nrow(spatdat), nsamples)]
  
  variogramdat <- variog(coords = spatdat[, c("Bng_easting", "Bng_northing"), with = FALSE], data = as.numeric(spatdat$CaO), max.dist = 350000, breaks = seq(0,350000, distbreaks))
  variogramdat
  
  #datvar <- as.data.table(variogramdat$exp_var)
  datvar <- as.data.table(data.frame(dist = variogramdat$u, gamma = variogramdat$v))
  datvar[, type := "observations", ]
  datvar[, run := i, ]
  datvar
  
  simdat <- sim[run == sample(1:max(run), 1)][sample(1:.N, nsamples)]
  simdat
  #simdatsp <- SpatialPointsDataFrame(simdat[, c("x", "y"), with = FALSE], data.frame(CaO = simdat$simpred))
  
  #variogramsim <- autofitVariogram(as.formula(paste(elem, "~", 1)), simdatsp, model="Exp", verbose = FALSE)
  #simvar <- as.data.table(variogramsim$exp_var)
  
  variogramsim <- variog(coords = simdat[, c("x", "y"), with = FALSE], data = simdat$simpred, max.dist = 350000, breaks = seq(0,350000, distbreaks))
  variogramsim
  
  simvar <- as.data.table(data.frame(dist = variogramsim$u, gamma = variogramsim$v))
  simvar[, type := "simulations", ]
  simvar[, run := i, ]
  simvar
  
  if(exists("variogramdata")){
    variogramdata <- rbind(variogramdata, rbind(datvar, simvar))
  }else{
    variogramdata <- rbind(datvar, simvar)
  }
}

variogramdata[, dist := dist/1000, ]
variogramdata

ggplot(variogramdata) + 
  geom_smooth(aes(x = dist, y = gamma, group = interaction(type, run), col = type), alpha = 0.5, method = "gam", se = FALSE) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  theme_bw() + scale_y_continuous(limits=c(0, max(variogramdata$gamma)*1.2), expand=c(0,0)) +
  scale_x_continuous(limits=c(0, 350), expand=c(0,0)) +
  labs(x = "distance (km)", y = "semivariance") +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background=element_blank(), legend.title=element_blank(), 
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5))

variogramsummary <- copy(variogramdata)
variogramsummary[, dist := round(dist/4400, 0)*4.4, ]
variogramsummary
variogramsummary <- unique(variogramsummary[, list(upper = quantile(gamma, 0.975), lower = quantile(gamma, 0.025), mean = mean(gamma)), by = c("dist", "type")])

ggplot() + 
  geom_ribbon(data = variogramsummary[type == "simulations"], aes(x = dist, ymin = lower, ymax = upper, group = type, col = type, fill = type), alpha = 0.33) +
  geom_line(data = variogramsummary[type == "simulations"], aes(x = dist, y = mean, group = type, col = type), size = 1.25, se = FALSE) +
  geom_ribbon(data = variogramsummary[type == "observations"], aes(x = dist, ymin = lower, ymax = upper, group = type, col = type, fill = type), alpha = 0.33) +
  geom_line(data = variogramsummary[type == "observations"], aes(x = dist, y = mean, group = type, col = type, style = type), size = 1.25, se = FALSE) +
  theme_bw() + scale_y_continuous(limits=c(0, max(variogramdata$gamma)), expand=c(0,0)) +
  scale_x_continuous(limits=c(0, 420), expand=c(0,0)) +
  labs(x = "distance (km)", y = "semivariance") +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background=element_blank(), legend.title=element_blank(), 
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5))
ggsave(paste0("plots/simulation/autocorrelation_comparison.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

