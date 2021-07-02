library(data.table)
library(ggplot2)
library(stringr)
library(compositions)
library(rgl)

options(scipen = 6)

# Load and wrangle data ####
# data from https://mrdata.usgs.gov/geochem/
alldata <- fread("data/geochem.csv")
str(alldata)
sort(alldata[, lapply(.SD, function(x) sum(is.na(x))), ], decreasing = TRUE)

sort(table(alldata$DATASET), decreasing = TRUE)
ggplot(alldata[DATASET == "NURE2002"]) + geom_point(aes(x = LONGITUDE, y = LATITUDE, col = DATASET))

# Crop to mainland USA and select only the analysis method with the least NAs (ICP40)
data <- alldata[LATITUDE < 50 & LONGITUDE > -150 & LONGITUDE < 0, c("LONGITUDE", "LATITUDE", names(alldata)[grep("_ICP40", names(alldata))]), with = FALSE]
names(data) <- c("Longitude", "Latitude", str_to_title(str_replace(names(data[,-(1:2)]), "\\_.*", "")))
data

data <- data[,lapply(.SD, as.numeric), ]
str(data)

sort(data[, lapply(.SD, function(x) sum(is.na(x))), ], decreasing = FALSE)

#elements <- names(data)[-(1:2)]
elements <- names(sort(data[, lapply(.SD, function(x) sum(is.na(x))), ], decreasing = FALSE))[3:25]
elements

# Plot to get our eye in
ggplot(data) + geom_point(aes(x = Longitude, y = Latitude, col = log(Ti)))

# Convert wt% elements into ppm like all the rest, to allow compositional analysis
data[, c("Al", "Ca", "Fe", "K", "Mg", "Na", "P", "Ti") := lapply(.SD, function(x) x*10000), .SDcols = c("Al", "Ca", "Fe", "K", "Mg", "Na", "P", "Ti")]
data

data <- na.omit(data[, c("Longitude", "Latitude", elements), with = FALSE])

# Set negative values (measurement errors) to zero
for(j in elements){set(data, i= which(data[[j]]<0), j= j, value=0)}
for(j in elements){set(data, i= which(data[[j]]==0), j= j, value=0.001)}
#for(j in elements){set(data, i= which(is.na(data[[j]])), j= j, value=0)}

# Add an unmeasured column
data[, R := 1000000 - rowSums(.SD, na.rm = TRUE), .SDcols = names(data)%in%elements]
#data[, Total := rowSums(.SD, na.rm = TRUE), .SDcols = names(data)%in%c(elements, "R")]
data

hist(data$R, breaks = 200)

# Perform log-ratio transformations ####

datacomp <- acomp(data, parts = c(elements, "R"))
datailr <- ilr(datacomp)

datailrpca <- prcomp(datailr, center = TRUE, scale = FALSE)
str(datailrpca)
ggplot(data.frame(Component = 1:length(datailrpca$sdev), SD = as.numeric(datailrpca$sdev))) + geom_point(aes(x = Component, y = SD)) + theme_bw() +
  geom_hline(yintercept = 1) + scale_x_continuous(limits = c(0, length(datailrpca$sdev)+1), expand = c(0,0))

clrloadings <- as.data.frame(ilr2clr(datailrpca$rotation))
names(clrloadings) <- c(elements, "R")

scores <- as.data.frame(predict(datailrpca, datailr))
names(scores) <- paste0("Comp.", 1:ncol(scores))

clrscores <- as.data.frame(ilr2clr(scores))
names(clrscores) <- c(elements, "R")

# datailrpca <- princomp(datailr)
# str(datailrpca)
# ggplot(data.frame(Component = 1:length(datailrpca$sdev), SD = as.numeric(datailrpca$sdev))) + geom_point(aes(x = Component, y = SD)) + theme_bw() +
#   geom_hline(yintercept = 1) + scale_x_continuous(limits = c(0, length(datailrpca$sdev)+1), expand = c(0,0))
# 
# clrloadings <- as.data.frame(ilr2clr(datailrpca$loadings))
# names(clrloadings) <- c(elements, "R")
# 
# clrscores <- as.data.frame(ilr2clr(datailrpca$scores))
# names(clrscores) <- c(elements, "R")
# 
# scores <- as.data.frame(datailrpca$scores)
# names(scores) <- paste0("Comp.", 1:ncol(scores))

colcomp <- c(1,2,3)
PCAcol <- rgb((scores[,colcomp[1]]-min(scores[,colcomp[1]]))/(max(scores[,colcomp[1]])-min(scores[,colcomp[1]])), 
              (scores[,colcomp[2]]-min(scores[,colcomp[2]]))/(max(scores[,colcomp[2]])-min(scores[,colcomp[2]])), 
              (scores[,colcomp[3]]-min(scores[,colcomp[3]]))/(max(scores[,colcomp[3]])-min(scores[,colcomp[3]])))
ggplot(data) + geom_point(aes(x = Longitude, y = Latitude, col = scores$Comp.1), alpha = 0.5)

ss <- sample(1:nrow(data), 10000)

ggplot(scores) + geom_point(aes(x = Comp.4, y = Comp.5), col = PCAcol)

pairs(scores[ss, 1:8], col = PCAcol[ss])

plot3d(scores[ss, 1:3], col = PCAcol[ss])

ss <- sample(1:nrow(data), 1000)

ggplot(data) + geom_point(aes(x = Longitude, y = Latitude), col = PCAcol)

head(clrscores)
head(datailrpca$clrscores)

data.frame(datailrpca$scores, names =)
