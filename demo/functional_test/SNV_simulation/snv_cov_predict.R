
data <- readRDS("SNV_4.RDS")

v.concordance <- c()

for (sample in names(data)) {
    sample.data <- data[[sample]]

    # prediction for precision
    loess.precision <- loess(m.concordance ~ cov, sample.data)
    v.concordance <- c(v.concordance, predict(loess.precision, 15))
}


mean(v.concordance)


