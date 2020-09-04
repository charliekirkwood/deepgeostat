library(data.table)
library(ggplot2)
library(ggnewscale)
library(plotly)
library(tidyverse)
library(rgdal)
library(raster)
library(keras)
library(tfprobability)
library(scoringRules)
library(scoringutils)

options(scipen = 6000)

dir.create(paste0(getwd(), "/plots_loc"))
dir.create(paste0(getwd(), "/models_loc"))

# Download UK elevation data ####
gbelev <- projectRaster(getData("alt", country = "GBR", mask = FALSE), crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")
names(gbelev) <- "alt"
elevdat <- as.data.table(raster::as.data.frame(gbelev, xy = TRUE, na.rm = TRUE))
elevdat

# Plot of UK elevation data
# ggplot(elevdat) + geom_raster(aes(x = x, y = y, fill = alt)) + theme_bw() +
#   scale_fill_viridis_c(option = "cividis", begin = 0.2, end = 1) + coord_equal() +
#   theme(legend.position = "none") + labs(x = "Easting (BNG)", y = "Northing (BNG)") +
#   scale_x_continuous(limits = c(-10000,660000)) + scale_y_continuous(limits = c(0,1200000))
# ggsave(paste0("plots_loc/UKelevation.png"), width = 68, height = 100, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Download UK outline data ####
gbout <- spTransform(getData("GADM", country = "GBR", level = 0), crs("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")) 
gbout <- fortify(gbout, region = "id")

# Read in stream sediment geochemistry data ####
gchem <- fread("data/OpenGeoScience2016 national geochemical map release stream sediment data points.csv")

# Plot of geochem sampe sites within UK outline
# ggplot(gbout) + geom_path(aes(x = long, y = lat, group = group)) + theme_bw() +
#   geom_point(data = gchem, aes(x = Bng_easting, y = Bng_northing), 
#              col = "black", shape = 16, size = 0.01, alpha = 0.2) +
#   coord_equal() + 
#   scale_x_continuous(limits = c(-10000,660000)) + scale_y_continuous(limits = c(0,1200000)) +
#   theme(legend.position = "none") + labs(x = "Easting (BNG)", y = "Northing (BNG)")
# ggsave(paste0("plots_loc/UKgeochem.png"), width = 68, height = 100, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Choose an element name to continue with, out of available names:
names(gchem)
elem <- "CaO"

# Histogram of the element
hist(gchem[,get(elem)], breaks = 20)
hist(gchem[,log(get(elem))], breaks = 20)
# Do we wish to model the element on the log scale?
logtrans <- TRUE

gchem[, elevation := extract(gbelev, methods = "bilinear", gchem[, c("Bng_easting", "Bng_northing"), with = FALSE]), ]
gchem <- gchem[!is.na(gchem[,get(elem)]) & !is.na(gchem$elevation)]
#gchem <- gchem[sample(1:nrow(gchem), 50000)][order(-As )]

# Plot point data for chosen element
ggplot(na.omit(gchem[, c("Bng_easting", "Bng_northing", paste(elem)), with = FALSE])) + 
  geom_point(aes(x = Bng_easting, y = Bng_northing, col = if(logtrans == TRUE){log(get(elem))}else{get(elem)})) +
  scale_colour_viridis_c(option = "B", name = if(logtrans == TRUE){paste0("log(",elem,")")}else{paste(elem)}) + labs(x = "Easting (BNG)", y = "Northing (BNG)") +
  theme(legend.justification = c(1, 0.5), legend.position = c(0.98, 0.7), legend.background=element_blank())

# Extract sample-centred terrain images, of specified dimension and resolution ####
imagedim <- 32
imageres <- 250

imgs <- array(dim = c(nrow(gchem), imagedim, imagedim))

seq.int(from = imageres/2, by = imageres, length.out = imagedim)-(imageres*imagedim)/2

cells <- as.data.table(expand.grid(x = seq.int(from = imageres/2, by = imageres, length.out = imagedim)-(imageres*imagedim)/2, 
                                   y = seq.int(from = imageres/2, by = imageres, length.out = imagedim)-(imageres*imagedim)/2))
cells[ ,coordx := rep(1:imagedim, imagedim),]
cells[ ,coordy := rep(1:imagedim, each = imagedim),]

time <- Sys.time()
for(xmeter in unique(cells$coordx)){
  for(ymeter in unique(cells$coordy)){
    imgs[,xmeter,ymeter] <- extract(gbelev, method = "bilinear",
                                    cbind(gchem$Bng_easting + cells[coordx == xmeter & coordy == ymeter]$x, gchem$Bng_northing + cells[coordx == xmeter & coordy == ymeter]$y)) 
  }
}
Sys.time() - time

imgs[is.na(imgs)] = 0

# Plot an image at random:
ggplot(reshape2::melt(imgs[sample(1:nrow(gchem), 1),,]), aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() + scale_fill_viridis_c(option = "B", name = "elevation") +
  coord_fixed()

# Scale sample-centered terrain images for use in neural network ####

mean = mean(as.data.frame(gbelev, na.rm = TRUE)[,1])
sd =  sd(as.data.frame(gbelev, na.rm = TRUE)[,1])

imgs_ann <- imgs - gchem$elevation
imgs_ann <- (imgs_ann)/sd
imgs_ann[is.na(imgs_ann)] = 0
#imgs_ann <- imgs_ann - rowMeans(imgs_ann, na.rm = TRUE)

# Plot a normalised image at random
ggplot(reshape2::melt(imgs_ann[sample(1:nrow(gchem), 1),,]), aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() + scale_fill_viridis_c(option = "B", name = "normalised \n elevation") +
  coord_fixed()

# Extract easting, northing and elevation as location variables

loc <- gchem[, c("Bng_easting", "Bng_northing", "elevation"), with = FALSE]

locmean <- apply(loc, 2, mean)
locsd <- apply(loc, 2, sd)

loc_ann <- t(apply(loc, 1, function(x) (x - locmean)/locsd))

# Plot and save three random terrain images to use as figure in paper
set.seed(789)
img1 <- reshape2::melt(imgs_ann[sample(1:nrow(gchem), 1),,])
img2 <- reshape2::melt(imgs_ann[sample(1:nrow(gchem), 1),,])
img3 <- reshape2::melt(imgs_ann[sample(1:nrow(gchem), 1),,])

ggplot(img1, aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() + 
  scale_fill_viridis_c(option = "cividis", limits = c(min(rbind(img1,img2,img3)$value),max(rbind(img1,img2,img3)$value))) +
  coord_fixed() + labs(x = "", y = "") + theme(legend.position = "none")
ggsave(paste0("plots_loc/img1.png"), width = 45, height = 45, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

ggplot(img2, aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() +
  scale_fill_viridis_c(option = "cividis", limits = c(min(rbind(img1,img2,img3)$value),max(rbind(img1,img2,img3)$value))) +
  coord_fixed() + labs(x = "", y = "") + theme(legend.position = "none")
ggsave(paste0("plots_loc/img2.png"), width = 45, height = 45, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

ggplot(img3, aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() + 
  scale_fill_viridis_c(option = "cividis", limits = c(min(rbind(img1,img2,img3)$value),max(rbind(img1,img2,img3)$value))) +
  coord_fixed() + labs(x = "", y = "") + theme(legend.position = "none")
ggsave(paste0("plots_loc/img3.png"), width = 45, height = 45, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Use neural network to learn relationship between terrain features and geochemistry ####

set.seed(321)
fold_size = nrow(gchem)/10
test <- sample(1:nrow(gchem), fold_size)
val <- sample(which(!1:nrow(gchem) %in% test), fold_size)
#train <- sample(which(!1:nrow(gchem) %in% c(test, val)), length(which(!1:nrow(gchem) %in% c(test, val)))*0.1)
train <- which(!1:nrow(gchem) %in% c(test, val))

# Parameters --------------------------------------------------------------

batch_size <- 256
epochs <- 500

# Data Preparation --------------------------------------------------------

x_train <- list(imgs_ann[train,,], loc_ann[train, ])
x_val <- list(imgs_ann[val,,], loc_ann[val, ])
x_test <- list(imgs_ann[test,,], loc_ann[test, ])

dim(x_train[[1]]) <- c(nrow(x_train[[1]]), imagedim, imagedim, 1)
dim(x_val[[1]]) <- c(nrow(x_val[[1]]), imagedim, imagedim, 1)
dim(x_test[[1]]) <- c(nrow(x_test[[1]]), imagedim, imagedim, 1)

if(logtrans == TRUE){
  y_train <- as.numeric(log(gchem[,get(elem)]))[train]
  y_val <- as.numeric(log(gchem[,get(elem)]))[val]
  y_test <- as.numeric(log(gchem[,get(elem)]))[test]
} else {
  y_train <- as.numeric(gchem[,get(elem)])[train]
  y_val <- as.numeric(gchem[,get(elem)])[val]
  y_test <- as.numeric(gchem[,get(elem)])[test]
}

# Defining Model ----------------------------------------------------------
time <- Sys.time()

# Convolutional stack:
conv_input <- layer_input(shape = c(imagedim, imagedim, 1), name = 'conv_input')

conv_output <- conv_input %>%
  #layer_average_pooling_2d(pool_size = c(3,3), strides = c(1,1)) %>%
  
  layer_conv_2d(filter = 128, kernel_size = c(3,3), dilation_rate = 1, strides = 3) %>%
  layer_activation("relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_spatial_dropout_2d(rate = 0.25) %>%
  
  layer_conv_2d(filter = 128, kernel_size = c(3,3), dilation_rate = 1, strides = 1) %>%
  layer_activation("relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_spatial_dropout_2d(rate = 0.25) %>%
  
  layer_conv_2d(filter = 128, kernel_size = c(3,3), dilation_rate = 1, strides = 1) %>%
  layer_activation("relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_spatial_dropout_2d(rate = 0.25) %>%
  
  layer_conv_2d(filter = 128, kernel_size = c(3,3), dilation_rate = 1, strides = 1) %>%
  layer_activation("relu") %>%
  layer_dropout(rate = 0.2) %>% 
  layer_spatial_dropout_2d(rate = 0.25) %>%
  
  layer_average_pooling_2d(pool_size = c(3,3), strides = c(1,1)) %>%
  layer_flatten()

# Auxiliary input:
auxiliary_input <- layer_input(shape = c(3), name = 'aux_input')

auxiliary_output <- auxiliary_input %>%
  layer_dense(512) %>%
  layer_activation("relu") %>%
  layer_dropout(rate = 0.1) %>%
  layer_flatten()

# Main output:
main_output <- layer_concatenate(c(conv_output, auxiliary_output)) %>%  
  layer_dense(256) %>%
  layer_activation("relu") %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(128) %>%
  layer_activation("relu") %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 2, activation = "linear", name = "dist_param") %>%
  layer_distribution_lambda(function(x)
    tfd_normal(loc = x[, 1, drop = FALSE],
               scale = 1e-3 + tf$math$softplus(0.1 * x[, 2, drop = FALSE])
    )
  )

rm(model)
model <- keras_model(
  inputs = c(conv_input, auxiliary_input), 
  outputs = main_output
)

summary(model)

negloglik <- function(y, model) - (model %>% tfd_log_prob(y))

opt <- optimizer_adam(lr = 0.001, decay = 1e-6)

model %>% compile(
  loss = negloglik,
  optimizer = opt
)

# Training ----------------------------------------------------------------

history <- model %>% fit(
  x_train, y_train,
  batch_size = batch_size,
  epochs = epochs,
  validation_data = list(x_val, y_val),
  shuffle = TRUE,
  callbacks = callback_early_stopping(monitor = "val_loss", patience = 50, restore_best_weights = TRUE)
)
min(history$metrics$val_loss)
Sys.time() - time


holdout <- data.frame(obs = y_test, preds = as.numeric(tfd_mean(model(x_test))))
#holdout <- readRDS(paste0(elem, "_holdout_data_",63, ".rds"))
ggplot(holdout) + geom_point(aes(x = obs, y = preds), shape = 16, size = 1, alpha = 0.1) + 
  coord_equal() + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) +
  labs(x = if(logtrans == TRUE){paste0("observed log(", elem, ")")}else{paste("observed", elem)},
       y = if(logtrans == TRUE){paste0("predicted log(", elem, ")")}else{paste("predicted", elem)}, 
       subtitle = paste0("R\U00B2 = ", round(cor(holdout$preds, holdout$obs)^2, 2), "     RMSE = ", round(sqrt(mean((holdout$preds-holdout$obs)^2)), 2)),
       tag = expression(bold("a")))
ggsave(paste0("plots_loc/", elem, "_mean_holdout_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# holdoutsim <- data.frame(obs = y_test, preds = replicate(10, predict(model, x_test, batch_size = 4096)))
# holdout <- melt(holdoutsim, id.vars = "obs", value.name = "preds")
# ggplot(holdout) + geom_point(aes(x = obs, y = preds), shape = 16, size = 0.4, alpha = 0.2) + 
#   coord_equal() + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) +
#   labs(x = if(logtrans == TRUE){paste0("observed log(", elem, ")")}else{paste("observed", elem)},
#        y = if(logtrans == TRUE){paste0("predicted log(", elem, ")")}else{paste("predicted", elem)}, 
#        subtitle = paste0("R\U00B2 = ", round(cor(holdout$preds, holdout$obs)^2, 2), "     RMSE = ", round(sqrt(mean((holdout$preds-holdout$obs)^2)), 2)))
# cor(holdout$preds, holdout$obs)^2
# ggsave(paste0("plots_loc/", elem, "_stochastic_holdout_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 45, height = 45, units = "mm", type = "cairo", dpi = 300, scale = 1.375)


#saveRDS(holdout, paste0(elem, "_holdout_data_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".rds"))
model %>% save_model_weights_hdf5(paste0(getwd(), "/models_loc/", elem, "_dist_terrain_cnn_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".hdf5"))
# model %>% load_model_weights_hdf5(paste0(getwd(), "/models_loc/", elem, "_dist_terrain_cnn_72.hdf5"))
#k_set_learning_phase(1)

trainhist <- as.data.table(history$metrics)[,c("loss", "val_loss"), with = FALSE]
names(trainhist) <- c("training", "testing")
trainhist[, epoch := 1:.N, ]

ggplot(melt(trainhist, id.vars = "epoch", value.name = "MSE", variable.name = "dataset")) + 
  geom_line(aes(x = epoch, y = MSE, col = dataset)) + theme_bw() + scale_y_continuous(limits = c(0,quantile(c(trainhist$testing, trainhist$testing), 0.999))) +
  geom_vline(xintercept = trainhist[testing == min(testing)]$epoch, linetype = "dashed")
ggsave(paste0("plots_loc/", elem, "_training_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 136, height = 60, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Evaluate calibration and sharpness

# calib <- verifRankHist(replicate(10, predict(model, x_test, batch_size = 4096)), y_test)
# 
# calib <- apply(cbind(y_test, as.matrix(replicate(20, predict(model, x_test, batch_size = 4096)))),1,function(x) rank(x,ties="random")[1])
# str(calib)
# hist(calib, breaks = 20, prob = TRUE, xaxt = "n",
#      xlab = "", ylab = "", main = "Verification Rank Histogram")

k_set_learning_phase(1)

holdoutcdf <- data.frame(draw = seq(1,10,1), obs = rep(y_test, 10), preds = as.numeric(replicate(10, predict(model, x_test, batch_size = 4096))))
str(melt(holdoutcdf, id.vars = "draw"))

# decdf <- function(x, pred, obs)  ecdf(pred)(x) - ecdf(obs)(x)
# curve(decdf(x,holdoutcdf$preds,holdoutcdf$obs), from=min(holdoutcdf$preds,holdoutcdf$obs), to=max(holdoutcdf$preds,holdoutcdf$obs))
# curve(decdf(x,holdoutcdf$preds,holdoutcdf$obs), from=-6, to=6)

# ggplot() +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = as.numeric(replicate(1, predict(model, x_test, batch_size = 4096)))), aes(x = x), geom = "line", col = "red", alpha = 0.5) +
#   stat_density(data = data.frame(x = y_test), aes(x = x), geom = "line", size = 1) +
#   scale_x_continuous(limits = c(-5,5)) + theme_bw()

testcrps <- crps_sample(y_test, drop(replicate(50, predict(model, x_test, batch_size = 4096))))
mean(testcrps)

testlogs <- logs_sample(y_test, drop(replicate(50, predict(model, x_test, batch_size = 4096))))
mean(testlogs)

ggplot() +
  stat_density(data = data.frame(x = y_test, data = "observed"), aes(x = x, col = data, linetype = data), geom = "line") +
  stat_density(data = data.frame(x = as.numeric(replicate(50, predict(model, x_test, batch_size = 4096))), data = "predicted"), 
               aes(x = x, col = data, linetype = data), geom = "line") +
  scale_x_continuous(limits = c(-6,6)) + theme_bw() + theme(legend.justification = c(1,1), legend.position = c(0.95, 0.95)) +
  labs(x = if(logtrans == TRUE){paste0("log(", elem, ")")}else{paste(elem)},
       subtitle = paste0("CRPS = ", round(mean(testcrps), 2), "     Log Score = ", round(mean(testlogs), 2)),
       tag = expression(bold("b")))
ggsave(paste0("plots_loc/distribution_comparison.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Coverage:

covtest <- melt(as.data.table(data.frame(y = y_test, x = drop(replicate(100, predict(model, x_test, batch_size = 4096))))), id.vars = "y")
covtest
paste("95% =", 1-round(nrow(covtest[y > quantile(value, 0.95)])/nrow(covtest), 3))
paste("70% =", 1-round(nrow(covtest[y > quantile(value, 0.70)])/nrow(covtest), 3))
paste("50% =", 1-round(nrow(covtest[y > quantile(value, 0.50)])/nrow(covtest), 3))

# Create covariate maps using the trained neural network

meanmodel <- keras_model(
  inputs = model$input, 
  outputs = get_layer(model, "dist_param")$output[,1],
  trainable = TRUE
)

save_model_tf(meanmodel, "meanmodelCaO.tf")
meanmodel <- load_model_tf("meanmodelCaO.tf")

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

predgrid[, paste(elem) := predict(meanmodel, list(predimgs_ann, predloc_ann), batch_size = 4096), ]
#predgrid[, paste(elem) := predict(model, list(predimgs_ann, predloc_ann), batch_size = 4096), ]

ggplot(predgrid) + geom_raster(aes(x = x, y = y, fill = get(elem)), interpolate = FALSE) + theme_bw() + coord_equal() +
  scale_fill_viridis_c(option = "B", name = if(logtrans == TRUE){paste0("predicted\nlog(",elem,")")}else{paste("predicted\n", elem)}) + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.90), legend.background=element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), labels = paste(seq(0,1000000, length.out = 6)), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 670000))
ggsave(paste0("plots_loc/", elem, "_national_map_full_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 183, height = 225, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
ggsave(paste0("plots_loc/", elem, "_national_map_half_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 89, height = 145, units = "mm", type = "cairo", dpi = 300, scale = 1.375)


# Local scale maps ####
rm(predgrid, predimgs, predimgs_ann, predloc, predloc_ann)

gbelevcoarse <- aggregate(gbelev, fact=1)
#predgrid <- as.data.table(raster::as.data.frame(gbelevcoarse, xy = TRUE, na.rm = TRUE))[x > 350000 & x < 700000 & y > 75000 & y < 225000]
#predgrid <- as.data.table(raster::as.data.frame(gbelevcoarse, xy = TRUE, na.rm = TRUE))[x > 200000 & x < 510000 & y > 75000 & y < 225000]
predgrid <- as.data.table(raster::as.data.frame(gbelevcoarse, xy = TRUE, na.rm = TRUE))[x > 350000 & x < 500000 & y > 300000 & y < 500000]
gchemsub <- gchem[Bng_easting >  350000 & Bng_easting < 500000 & Bng_northing > 300000 & Bng_northing < 500000]
#gchemsub <- gchemsub[sample(1:nrow(gchemsub), 500)]
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

predgrid[, paste(elem) := predict(meanmodel, list(predimgs_ann, predloc_ann), batch_size = 4096), ]
#predgrid[, paste(elem) := as.numeric(tfd_mean(model(list(predimgs_ann, predloc_ann)))), ]
#

ggplot(predgrid) + geom_raster(aes(x = x, y = y, fill = get(elem)), interpolate = FALSE) + theme_bw() + coord_equal() +
  scale_fill_viridis_c(option = "B", name = if(logtrans == TRUE){paste0("log(",elem,")")}else{paste(elem)}) + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)", tag = expression(bold("b"))) +
  theme(legend.justification = c(0, 0.5), legend.position = "none", legend.background=element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(350000,450000, length.out = 3), labels = paste(seq(350000,450000, length.out = 3)), expand = c(0,0))
ggsave(paste0("plots_loc/", elem, "_mean_local_map_third_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 61, height = 72, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
ggsave(paste0("plots_loc/", elem, "_mean_local_map_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 136, height = 136, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
ggsave(paste0("plots_loc/", elem, "_mean_local_map_double_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 272, height = 272, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

k_set_learning_phase(1)
predgrid[, paste(elem) := predict(model, list(predimgs_ann, predloc_ann), batch_size = 4096), ]

ggplot(predgrid) + geom_raster(aes(x = x, y = y, fill = get(elem)), interpolate = FALSE) + theme_bw() + coord_equal() +
  scale_fill_viridis_c(option = "B", name = if(logtrans == TRUE){paste0("log(",elem,")")}else{paste(elem)}) + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)", tag = expression(bold("c"))) +
  theme(legend.justification = c(0, 0.5), legend.position = "none", legend.background=element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(350000,450000, length.out = 3), labels = paste(seq(350000,450000, length.out = 3)), expand = c(0,0))
ggsave(paste0("plots_loc/", elem, "_stoch_local_map_third_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 61, height = 72, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
ggsave(paste0("plots_loc/", elem, "_stoch_local_map_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 136, height = 136, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
ggsave(paste0("plots_loc/", elem, "_stoch_local_map_double_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 272, height = 272, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
#

# uncertainty map:

k_set_learning_phase(1)
predgrid[, paste(elem) := apply(drop(replicate(50, predict(model, list(predimgs_ann, predloc_ann), batch_size = 4096))), 1, function(x) sd(x)), ]

ggplot(predgrid) + geom_raster(aes(x = x, y = y, fill = get(elem)), interpolate = FALSE) + theme_bw() + coord_equal() +
  scale_fill_viridis_c(option = "B", name = if(logtrans == TRUE){paste0("log(",elem,")")}else{paste(elem)}) + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)", tag = expression(bold("c"))) +
  theme(legend.justification = c(0, 0.5), legend.position = "none", legend.background=element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(350000,450000, length.out = 3), labels = paste(seq(350000,450000, length.out = 3)), expand = c(0,0))
ggsave(paste0("plots_loc/", elem, "_unc_local_map_third_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 61, height = 72, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
ggsave(paste0("plots_loc/", elem, "_unc_local_map_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 136, height = 136, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
ggsave(paste0("plots_loc/", elem, "_unc_local_map_double_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 272, height = 272, units = "mm", type = "cairo", dpi = 300, scale = 1.375)


ggplot(predgrid) + geom_raster(aes(x = x, y = y, fill = get(elem)), interpolate = FALSE) +
  labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") + theme_bw() + coord_equal() +
  theme(legend.justification = c(0, 0.5), legend.position = "none", legend.background=element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(350000,450000, length.out = 3), labels = paste(seq(350000,450000, length.out = 3)), expand = c(0,0)) +
  geom_point(data = gchemsub, aes(x = Bng_easting, y = Bng_northing, fill = if(logtrans == TRUE){log(get(elem))}else{get(elem)}), shape = 21, col = "grey", size = 0.8) +
  scale_fill_viridis_c(option = "B", name = if(logtrans == TRUE){paste0("log(",elem,")")}else{paste(elem)},
                       limits = c(min(c(predgrid[,get(elem)], gchemsub[, if(logtrans == TRUE){log(get(elem))}else{get(elem)}])), 
                                  max(c(predgrid[,get(elem)], gchemsub[, if(logtrans == TRUE){log(get(elem))}else{get(elem)}]))))
ggsave(paste0("plots_loc/", elem, "_local_map_with_points_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 136, height = 60, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

#

ggplot(predgrid) + geom_raster(aes(x = x, y = y, fill = elevation), interpolate = FALSE) + theme_bw() + coord_equal() +
  scale_fill_viridis_c(option = "cividis") + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)", tag = expression(bold("a"))) +
  theme(legend.justification = c(0, 0.5), legend.position = "none", legend.background=element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(350000,450000, length.out = 3), labels = paste(seq(350000,450000, length.out = 3)), expand = c(0,0))
ggsave(paste0("plots_loc/", elem, "_local_elevation_civid_third_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 61, height = 72, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
ggsave(paste0("plots_loc/", elem, "_local_elevation_civid_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 136, height = 136, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Cross section probabilities ####
rm(predgrid, predimgs, predimgs_ann, predloc, predloc_ann)
gc()

predgrid <- as.data.table(data.frame(x = 400000, y = seq(000000, 700000, by = 2000)))
predgrid <- na.omit(predgrid[, elevation := extract(gbelev, methods = "bilinear", predgrid[, c("x", "y"), with = FALSE]), ])
predgrid

gchembuffer <- 500
gchemsub <- gchem[Bng_easting >  min(predgrid$x) - gchembuffer & 
                    Bng_easting < max(predgrid$x) + gchembuffer & 
                    Bng_northing > min(predgrid$y) - gchembuffer & 
                    Bng_northing < max(predgrid$y) + gchembuffer]
nrow(gchemsub)

ggplot(predgrid) + geom_line(aes(x = y, y = elevation*50, col = elevation)) + theme_bw() + scale_colour_viridis_c(option = "cividis") + coord_equal()

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

# predict(meanmodel, list(predimgs_ann, predloc_ann), batch_size = 4096)
# 
# plot(as.numeric(tfd_mean(model(list(predimgs_ann, predloc_ann)))))
# plot(as.numeric(tfd_variance(model(list(predimgs_ann, predloc_ann)))))
# plot(as.numeric(tfd_sample(model(list(predimgs_ann, predloc_ann)))))
# points(as.numeric(predict(model, (list(predimgs_ann, predloc_ann)), batch_size = 4096)), col = "red")
# plot(as.matrix(tfd_quantile(model(list(predimgs_ann, predloc_ann)), c(0.05, 0.5, 0.95))))

k_set_learning_phase(1)

reps <- 250

quantconf <- replicate(reps, as.matrix(tfd_quantile(model(list(predimgs_ann, predloc_ann)), c(0.05, 0.5, 0.95))))
quantconf <- lapply(seq_len(dim(quantconf)[3]), function(i) as.data.frame(quantconf[,,i]))
quantconf <- rbindlist(quantconf, idcol = TRUE)
quantconf
names(quantconf) <- c("draw", "lower", "mean", "upper")
quantconf <- cbind(do.call("rbind", replicate(reps, predgrid, simplify = FALSE)), quantconf)
quantconf_melt <- melt(quantconf, id.vars = c("x", "y", "elevation", "draw"), value.name = paste(elem))
#quantconf_melt[, group := paste(draw, variable), ]
quantconf_melt

meanconf <- cbind(predgrid, as.data.frame(replicate(reps, as.numeric(tfd_mean(model(list(predimgs_ann, predloc_ann)))))))
meanconf_melt <- melt(meanconf, id.vars = c("x", "y", "elevation"), value.name = paste(elem))
meanconf_melt[, confup := quantile(get(elem), 0.95), by = y]
meanconf_melt[, conflow := quantile(get(elem), 0.05), by = y]
meanconf_melt[, mean := mean(get(elem)), by = y]

postsim <- cbind(predgrid, as.data.frame(replicate(reps, as.numeric(tfd_sample(model(list(predimgs_ann, predloc_ann)))))))
postsim_melt <- melt(postsim, id.vars = c("x", "y", "elevation"), value.name = paste(elem))
postsim_melt[, confup := quantile(get(elem), 0.95), by = y]
postsim_melt[, conflow := quantile(get(elem), 0.05), by = y]
postsim_melt[, mean := mean(get(elem)), by = y]

ggplot() +
  geom_jitter(data = postsim_melt, aes (x = y, y = get(elem), col = get(elem)), alpha = 0.25, shape = 16, size = 1, width = 1000) +
  scale_colour_viridis_c(option = "B", name = "samples from\npredictive\ndistribution") + 
  new_scale_colour() + 
  geom_ribbon(data = meanconf_melt, aes(x = y, ymin = conflow, ymax = confup, group = variable, fill = "")) +
  geom_point(data = gchemsub, aes(x = Bng_northing, y = log(get(elem)), col = ""), shape = 16, size = 1.5) +
  scale_fill_manual(name = "predicted mean\n(95% credible interval)", values = "lightgrey") + scale_color_manual(name = "observations", values = "black") + theme_bw() +
  theme(legend.position = c(0.05, 0.05), legend.justification = c(0,0), legend.direction = "horizontal", legend.box = "horizontal") +
  labs(x = "Northing (metres BNG)", y = if(logtrans == TRUE){paste0("log(", elem, ")")}else{paste(elem)}) +
  guides(colour = guide_legend(title.vjust = 0.5), fill = guide_legend(title.vjust = 0.5))
ggsave(paste0("plots_loc/", elem, "_cross_section_2_", round(cor(holdout$preds, holdout$obs)^2, 2)*100, ".png"), width = 183, height = 72, units = "mm", dpi = 300, scale = 1.375)