#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(ChIPQC)

configFile=args[1]
annotationName=args[2]

cat("configFile=", configFile, "\n")
cat("annotationName=", annotationName, "\n")

samples <- read.table(configFile, sep="\t", header=T)

experiment = ChIPQC(samples, annotaiton = annotationName)

ChIPQCreport(experiment)
