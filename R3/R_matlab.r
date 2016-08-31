rm(list=ls())
require(R.matlab)
data <- readMat('~/git/abbg/R3/StephaneNew/data_hermite_cons2.mat')
save(data, file='~/git/abbg/R3/mat_new.dat')