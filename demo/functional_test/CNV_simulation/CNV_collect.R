library(ggplot2)

chr_filter <- function (file) {
    data <- read.table(file = file, header = TRUE)
    data <- data[which(data$chromosome %in% c("chr20")), ]
    
    return(data$gene)
}

samples <- c("IC15", "IC17", "IC20", "IC35", "IC37")
covs <- c(29.77, 42.08, 23.38, 18.22, 38.22)
bin_type <- c("auto", "1000", "2000", "5000", "7000", "10000")

df.summary <- list()

for (ll in seq(5)) {
    sample <- samples[ll]
    cov <- covs[ll]
    
    df.sample <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(df.sample) <- c("group", "cov", "m.precision", "s.precision", "m.recall", "s.recall")
    
    for (bin_length in bin_type) {
        path <- paste0("./", sample, "/CNV/Bin_", bin_length, "/")
        
        # get chr20 info
        data_chr20 <- chr_filter(file = file.path(path, paste0("chr20/intermediate_result/step_CNV02_cnvTable/", sample, "_chr20_genemetrics_cnr.txt")))
        
        count.FN <- list()
        count.TP <- list()
        count.FP <- list()
        
        for (i in c("01", "05", seq(9))) {
            this.FP <- c()
            this.TP <- c()
            this.FN <- c()
            for (j in seq(10)) {
                print(paste0("Now, processing ", sample, "_", i, "_", j))
                this.file <- file.path(path, paste0("ds_", i, "_", j, "/intermediate_result/step_CNV02_cnvTable/", sample, "_", i, "_", j, "_genemetrics_cnr.txt"))
                if (file.exists(this.file)) {
                    this.data <- chr_filter(file = this.file)
                    this.allcount <- length(this.data)
                    this.validcount <- length(intersect(data_chr20, this.data))
                    this.misscount <- length(data_chr20) - this.validcount
                    this.fakecount <- this.allcount - this.validcount
                    this.FP <- c(this.FP, this.fakecount)
                    this.TP <- c(this.TP, this.validcount)
                    this.FN <- c(this.FN, this.misscount)
                } else {
                    print(paste0(sample, "_", i, "_", j, " do not exist!"))
                }
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
        
        
        df <- data.frame(m.precision = m.precision, s.precision = s.precision, 
                         m.recall = m.recall, s.recall = s.recall,
                         cov = c(0.01, 0.05, seq(0.1, 0.9, 0.1)) * cov, group = rep(bin_length, 11))
        
        df.sample <- rbind(df.sample, df)
    }
    
    df.summary[[sample]] <- df.sample
    
    ## plot precision
    # p1 <- ggplot(df.sample, aes(x=cov, y=m.precision, group=group, color=group)) + 
    #     geom_errorbar(aes(ymin=m.precision-s.precision, ymax=m.precision+s.precision), width=0.5) +
    #     geom_line(lwd=1) +
    #     geom_hline(yintercept=0.95, linetype="dashed", color = "red") +
    #     geom_point() +
    #     scale_color_brewer(palette="Paired") + 
    #     theme_minimal() +
    #     xlim(-0.5, max(cov) + 1) +                                                              
    #     ylim(0, 1) +
    #     xlab("Sequence Coverage") +
    #     ylab("Precision") +
    #     theme(axis.text.x = element_text(size = 13),
    #           axis.text.y = element_text(size = 13),  
    #           axis.title.x = element_text(size = 18),
    #           axis.title.y = element_text(size = 18),
    #           legend.text=element_text(size=12))
    
    # ggsave(filename = fig1, plot = p1, width = 7, height = 4.5)
    
    ## plot recall
    # p2 <- ggplot(df.sample, aes(x=cov, y=m.recall, group=group, color=group)) + 
    #     geom_errorbar(aes(ymin=m.recall-s.recall, ymax=m.recall+s.recall), width=0.5) +
    #     geom_line(lwd=1) +
    #     geom_hline(yintercept=0.95, linetype="dashed", color = "red") +
    #     geom_point() +
    #     scale_color_brewer(palette="Paired") + 
    #     theme_minimal() +
    #     xlim(-0.5, max(cov) + 1) +
    #     ylim(0, 1) +
    #     xlab("Sequence Coverage") +
    #     ylab("Recall") +
    #     theme(axis.text.x = element_text(size = 13),
    #           axis.text.y = element_text(size = 13),  
    #           axis.title.x = element_text(size = 18),
    #           axis.title.y = element_text(size = 18),
    #           legend.text=element_text(size=12))
    
    # ggsave(filename = fig2, plot = p2, width = 7, height = 4.5)
}

saveRDS(object = df.summary, file = "CNV.RDS")

