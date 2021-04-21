library(ggplot2)
library(gridExtra)
library(tidyverse)
library(wesanderson)
library(dplyr)

data <- readRDS("SNV_4.RDS")

plot_lst <- vector("list")

for (i in seq(length(data))) {
    name <- names(data)[i]
    df.sample <- data[[name]]
    
    p1 <- ggplot(df.sample, aes(x=cov, y=m.concordance)) +
        geom_line(size=1) +
        geom_point(size=1)+
        geom_errorbar(aes(ymin=m.concordance-s.concordance, ymax=m.concordance+s.concordance), 
                      width=0.6, size=1,
                      position=position_dodge(0)) +
        theme_minimal() +
        xlim(-0.5, max(df.sample$cov) + 1) +
        ylim(0, 1) +
        xlab("Sequence Coverage") +
        ylab("Concordance") +
        theme(axis.text.x = element_text(size = 13),
              axis.text.y = element_text(size = 13),
              axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text = element_text(size=12))
    
    plot_lst[[i]] <- p1
}

fi.fig <- marrangeGrob(plot_lst, nrow = 2, ncol = 3)

# fi.fig

ggsave(filename = "SNV_plot.pdf", plot = fi.fig, width = 750, height = 350, units = "mm")

