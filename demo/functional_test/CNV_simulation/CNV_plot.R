library(ggplot2)
library(gridExtra)
library(tidyverse)
library(wesanderson)
library(dplyr)

data <- readRDS("CNV.RDS")

plot_lst <- vector("list", length = 10)

for (i in seq(length(data))) {
    name <- names(data)[i]
    df.sample <- data[[name]]
    
    p1 <- ggplot(df.sample, aes(x=cov, y=m.precision, group=desc(group), color=group)) +
        geom_errorbar(aes(ymin=m.precision-s.precision, ymax=m.precision+s.precision)) +
        geom_point() +
        geom_line() +
        theme_minimal() +
        xlim(-0.5, max(df.sample$cov) + 1) +
        ylim(0, 1) +
        xlab("Sequence Coverage") +
        ylab("Precision") +
        scale_color_manual(values=wes_palette(n=4, name="Darjeeling1")) + 
        theme(axis.text.x = element_text(size = 13),
              axis.text.y = element_text(size = 13),
              axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text=element_text(size=12))
    
    p2 <- ggplot(df.sample, aes(x=cov, y=m.recall, group=desc(group), color=group)) +
        geom_errorbar(aes(ymin=m.recall-s.recall, ymax=m.recall+s.recall)) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        xlim(-0.5, max(df.sample$cov) + 1) +
        ylim(0, 1) +
        xlab("Sequence Coverage") +
        ylab("Recall") +
        scale_color_manual(values=wes_palette(n=4, name="Darjeeling1")) + 
        theme(axis.text.x = element_text(size = 13),
              axis.text.y = element_text(size = 13),
              axis.title.x = element_text(size = 18),
              axis.title.y = element_text(size = 18),
              legend.text=element_text(size=12))
    
    plot_lst[[i]] <- p1
    plot_lst[[i + 5]] <- p2
}

fi.fig <- marrangeGrob(plot_lst, nrow = 5, ncol = 2)

fi.fig

ggsave(filename = "CNV_plot.pdf", plot = fi.fig, width = 210, height = 297, units = "mm")

