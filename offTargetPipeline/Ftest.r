#!/usr/bin/env Rscript

library(stats)
args <- commandArgs(trailingOnly=TRUE)

df <- read.table(paste(args[1], "/txt/indel_info.txt", sep=""), header = TRUE, sep = " ")

for (i in 1:length(df$site)) {
    indel <- df$indels[[i]]
    c_indel <- df$c_indels[[i]]
    regular <- df$regular[[i]]
    c_regular <- df$c_regular[[i]]
    df$freq[[i]] <- (indel / (indel + regular))*100
    df$c_freq[[i]] <- (c_indel / (c_indel + c_regular))*100
    tdf <- data.frame("indels" = c(indel, c_indel), "regular" = c(regular, c_regular), row.names = c("cas", "control"))
    df$Pvalue[[i]] <- fisher.test(tdf)[1][[1]]
}

df = as.matrix(df)
write.table(df, paste(args[1], "/txt/cor_indels_info_stat.txt", sep=""), dec=",")
