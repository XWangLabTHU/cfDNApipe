
data <- readRDS("CNV_0.5.RDS")

v.concordance <- c()

for (sample in names(data)) {
    sample.data <- data[[sample]]
    sample.data <- sample.data[which(sample.data$group %in% c("1000", "3000", "7000", "auto")), ]
    for (bin in unique(sample.data$group)) {
        this.data <- sample.data[which(sample.data$group == bin), ]
        
        # prediction for precision
        loess.precision <- loess(m.concordance ~ cov, this.data)
        v.concordance <- c(v.concordance, predict(loess.precision, 5))
    }
}

mean(v.concordance)



