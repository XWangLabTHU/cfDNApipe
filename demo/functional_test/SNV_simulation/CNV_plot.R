library(ggplot2)
library(gridExtra)
library(tidyverse)
library(wesanderson)
library(dplyr)
library(jcolors)

data <- readRDS("CNV_0.5.RDS")

plot_lst <- vector("list")

for (i in seq(length(data))) {
    name <- names(data)[i]
    df.sample <- data[[name]]
    df.sample <- df.sample[which(df.sample$group %in% c("1000", "3000", "7000", "auto")), ]
    
    p1 <- ggplot(df.sample, aes(x=cov, y=m.concordance, group=group, color=group)) +
        geom_line(size=1) +
        geom_point(size=1)+
        geom_errorbar(aes(ymin=m.concordance-s.concordance, ymax=m.concordance+s.concordance), 
                      width=2.8, size=1,
                      position=position_dodge(0)) +
        scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
        theme_minimal() +
        xlim(-1, max(df.sample$cov) + 1) +
        ylim(0, 1.1) +
        xlab("Sequence Coverage") +
        ylab("Concordance") +
        theme(axis.text.x = element_text(size = 13),
              axis.text.y = element_text(size = 13),
              axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text=element_text(size=12))
    
    plot_lst[[i]] <- p1
}

fi.fig <- marrangeGrob(plot_lst, layout_matrix = matrix(seq(6), nrow = 2, byrow = TRUE),nrow = 2, ncol = 3)

fi.fig

ggsave(filename = "CNV_plot.pdf", plot = fi.fig, width = 750, height = 350, units = "mm")

