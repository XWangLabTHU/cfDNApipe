
data <- readRDS("SNV.RDS")

v.precision <- c()
v.recall <- c()

for (sample in names(data)) {
    sample.data <- data[[sample]]

    # prediction for precision
    loess.precision <- loess(m.precision ~ cov, sample.data)
    v.precision <- c(v.precision, predict(loess.precision, 15))
    
    # prediction for recall
    loess.recall <- loess(m.recall ~ cov, sample.data)
    v.recall <- c(v.recall, predict(loess.recall, 15))
}





