# This file contains the function for freezing the dropout mask of the deep neural network #
# Calling this function allows us to evaluate single posterior samples at multiple locations #
# Useful for generating maps, cross sections, etc that are coherent simulations #
# Contact c.kirkwood@exeter.ac.uk with any queries
# Code also available at https://github.com/charliekirkwood/deepgeostat

dropout_sample <- function(...){
  drop_model <- keras_model(
    inputs = c(conv_input, auxiliary_input), 
    outputs = main_output
  )
  
  load_model_weights_hdf5(object = drop_model, filepath = paste0(getwd(), "/models/modelweights.hdf5"))
  weights <- drop_model %>% get_weights()
  
  dropweights <- weights
  #str(dropweights)
  
  # spatial dropout masks
  sampl1 <- rep(rbinom(dim(dropweights[[1]])[4], 1, 1-dropratespat))
  sampl2 <- rep(rbinom(dim(dropweights[[3]])[4], 1, 1-dropratespat))
  sampl3 <- rep(rbinom(dim(dropweights[[5]])[4], 1, 1-dropratespat))
  sampl4 <- rep(rbinom(dim(dropweights[[7]])[4], 1, 1-dropratespat))
  
  maskl1 <- array(rep(sampl1, each = 9), dim = dim(dropweights[[1]])) # generate dropout mask for kernel and bias of layer 1
  maskl2 <- array(rep(sampl1, each = 9), dim = dim(dropweights[[3]])) # generate dropout mask for kernel and bias of layer 2
  maskl3 <- array(rep(sampl1, each = 9), dim = dim(dropweights[[5]])) # generate dropout mask for kernel and bias of layer 3
  maskl4 <- array(rep(sampl1, each = 9), dim = dim(dropweights[[7]])) # generate dropout mask for kernel and bias of layer 4

  # normal dropout masks
  maskl5 <- rbinom(dim(dropweights[[9]])[2], 1, 1-dropratedense) # generate dropout mask for kernel and bias of layer 5
  maskl6 <- rbinom(dim(dropweights[[11]])[2], 1, 1-dropratedense) # generate dropout mask for kernel and bias of layer 6
  maskl7 <- rbinom(dim(dropweights[[13]])[2], 1, 1-dropratedense) # generate dropout mask for kernel and bias of layer 7
  
  dropweights[[1]] <- (dropweights[[1]] * maskl1)*(1/(1-dropratespat))
  dropweights[[2]] <- (dropweights[[2]] * sampl1)*(1/(1-dropratespat))
  
  dropweights[[3]] <- (dropweights[[3]] * maskl2)*(1/(1-dropratespat))
  dropweights[[4]] <- (dropweights[[4]] * sampl2)*(1/(1-dropratespat))
  
  dropweights[[5]] <- (dropweights[[5]] * maskl3)*(1/(1-dropratespat))
  dropweights[[6]] <- (dropweights[[6]] * sampl3)*(1/(1-dropratespat))
  
  dropweights[[7]] <- (dropweights[[7]] * maskl4)*(1/(1-dropratespat))
  dropweights[[8]] <- (dropweights[[8]] * sampl4)*(1/(1-dropratespat))
  
  dropweights[[9]] <- t(t(dropweights[[9]]) * maskl5)*(1/(1-dropratedense))
  dropweights[[10]] <- (dropweights[[10]] * maskl5)*(1/(1-dropratedense))
  
  dropweights[[11]] <- t(t(dropweights[[11]]) * maskl6)*(1/(1-dropratedense))
  dropweights[[12]] <- (dropweights[[12]] * maskl6)*(1/(1-dropratedense))
  
  dropweights[[13]] <- t(t(dropweights[[13]]) * maskl7)*(1/(1-dropratedense))
  dropweights[[14]] <- (dropweights[[14]] * maskl7)*(1/(1-dropratedense))
  
  drop_model %>% set_weights(dropweights)
  
  return(predict(drop_model, ...))
}