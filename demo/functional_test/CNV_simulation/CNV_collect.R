library(ggplot2)

chr_filter <- function (file) {
    data <- read.table(file = file, header = TRUE)
    data <- data[which(data$chromosome %in% c("chr20")), ]
    
    return(data$gene)
}

samples <- c("IC15", "IC17", "IC20", "IC35", "IC37")
covs <- c(29.77, 42.08, 23.38, 18.22, 38.22)
bin_type <- c("1000", "3000", "5000", "7000", "10000", "auto")



df.summary <- list()

for (ll in seq(5)) {
    sample <- samples[ll]
    cov <- covs[ll]
    
    df.sample <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(df.sample) <- c("group", "cov", 
                             "m.concordance", "s.concordance")
    
    for (bin_length in bin_type) {
        path <- paste0("./", sample, "/CNV/Bin_", bin_length, "/")
        chr20_positive <- chr_filter(file = file.path(path, paste0("chr20/intermediate_result/step_CNV03_cnvTable/", sample, "_chr20_genemetrics_cnr.txt")))
        
        count.concordance <- list()
        
        for (i in c("01", "05", seq(9))) {
            this.concordance <- c()
            for (j in seq(10)) {
                print(paste0("Now, processing ", sample, "_", i, "_", j))
                this.file <- file.path(path, paste0("ds_", i, "_", j, "/intermediate_result/step_CNV03_cnvTable/", sample, "_", i, "_", j, "_genemetrics_cnr.txt"))
                if (file.exists(this.file)) {
                    this.data <- chr_filter(file = this.file)
                    # TP + FP
                    this.allcount <- length(this.data)
                    # TP
                    this.TP <- length(intersect(chr20_positive, this.data))
                    # FN
                    this.FN <- length(chr20_positive) - this.TP
                    # FP
                    this.FP <- this.allcount - this.TP
                    # concordance
                    tmp_concordance <- this.TP / (this.TP + this.FP + this.FN)
                    this.concordance <- c(this.concordance, tmp_concordance)
                } else {
                    print(paste0(sample, "_", i, "_", j, " do not exist!"))
                    # concordance
                    this.concordance <- c(this.concordance, 0)
                }
            }
            name <- paste0("ds_", i)
            count.concordance[[name]] <- this.concordance
        }
        
        m.concordance <- unlist(lapply(count.concordance, mean))
        s.concordance <- unlist(lapply(count.concordance, sd))
        
        # concordance rate
        df <- data.frame(m.concordance = m.concordance, s.concordance = s.concordance,
                         cov = c(0.01, 0.05, seq(0.1, 0.9, 0.1)) * cov, group = rep(bin_length, 11))
        
        df.sample <- rbind(df.sample, df)
    }
    
    df.summary[[sample]] <- df.sample

}

saveRDS(object = df.summary, file = "CNV_0.5.RDS")














