
data <- readRDS("CNV.RDS")

v.precision <- c()
v.recall <- c()

for (sample in names(data)) {
    sample.data <- data[[sample]]
    for (bin in unique(sample.data$group)) {
        this.data <- sample.data[which(sample.data$group == bin), ]
        
        # prediction for precision
        loess.precision <- loess(m.precision ~ cov, this.data)
        v.precision <- c(v.precision, predict(loess.precision, 5))
        
        # prediction for recall
        loess.recall <- loess(m.recall ~ cov, this.data)
        v.recall <- c(v.recall, predict(loess.recall, 5))
    }
}





