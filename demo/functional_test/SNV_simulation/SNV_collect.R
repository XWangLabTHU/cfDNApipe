library(vcfR)
library(tidyr)
library(GenomicRanges)
library(gridExtra)
library(ggplot2)


TLOD_filter <- function (vcf.file, threshold, makeGR = FALSE) {
    vcf.data <- read.vcfR(vcf.file)
    vcf.data <- as.data.frame(vcf.data@fix)
    vcf.data <- separate(data = vcf.data, col = INFO, into = c("INFO", "TLOD"), sep = ";TLOD=")
    vcf.data <- vcf.data[which(as.numeric(vcf.data$TLOD) > threshold), ]
    
    if (makeGR == TRUE) {
        vcf.data <- vcf.data[c("CHROM", "POS", "POS")]
        colnames(vcf.data) <- c("chr", "start", "end")
        vcf.data <- makeGRangesFromDataFrame(vcf.data)
    }
    
    return(vcf.data)
}


samples <- c("IC15", "IC17", "IC20", "IC35", "IC37")
covs <- c(29.77, 42.08, 23.38, 18.22, 38.22)

threshold <- 4

df.summary <- list()

for (ll in seq(5)) {
    sample <- samples[ll]
    cov <- covs[ll]
    
    path <- paste0("./", sample, "/SNV/")
    
    # create granges
    chr20 <- TLOD_filter(vcf.file = paste0(path, "chr20/intermediate_result/step_somatic_bcftoolsVCF/", sample, "_chr20.somatic.vcf.gz"), 
                         threshold = threshold, makeGR = TRUE)
    
    count.FN <- list()
    count.TP <- list()
    count.FP <- list()
    
    for (i in c("01", "05", seq(9))) {
        this.FP <- c()
        this.TP <- c()
        this.FN <- c()
        for (j in seq(10)) {
            vcf.file <- paste0(path, "ds_", i, "_", j, "/intermediate_result/step_somatic_bcftoolsVCF/", sample, "_", i, "_", j, ".somatic.vcf.gz")
            print(vcf.file)
            this.vcf <- TLOD_filter(vcf.file = vcf.file, threshold = threshold, makeGR = TRUE)
            this.allcount <- length(this.vcf)
            this.validcount <- sum(countOverlaps(query = chr20, subject = this.vcf))
            this.misscount <- length(chr20) - this.validcount
            this.fakecount <- this.allcount - this.validcount
            this.FP <- c(this.FP, this.fakecount)
            this.TP <- c(this.TP, this.validcount)
            this.FN <- c(this.FN, this.misscount)
        }
        name <- paste0("ds_", i)
        count.FN[[name]] <- this.FP
        count.TP[[name]] <- this.TP
        count.FP[[name]] <- this.FN
    }
    
    m.FN <- unlist(lapply(count.FN, mean))
    s.FN <- unlist(lapply(count.FN, sd))
    
    m.TP <- unlist(lapply(count.TP, mean))
    s.TP <- unlist(lapply(count.TP, sd))
    
    m.FP <- unlist(lapply(count.FP, mean))
    s.FP <- unlist(lapply(count.FP, sd))
    
    # precision
    m.precision <- m.TP / (m.TP + m.FP)
    s.precision <- s.TP / (m.TP + m.FP)
    
    # recall
    m.recall <- m.TP / (m.TP + m.FN)
    s.recall <- s.TP / (m.TP + m.FN)
    
    df.sample <- data.frame(m.precision = m.precision, s.precision = s.precision, 
                     m.recall = m.recall, s.recall = s.recall,
                     cov = c(0.01, 0.05, seq(0.1, 0.9, 0.1)) * cov)
    
    df.summary[[sample]] <- df.sample
    
}


saveRDS(object = df.summary, file = "SNV.RDS")







