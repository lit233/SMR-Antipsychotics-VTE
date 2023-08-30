library(forestplot)
library(forestploter)
library(readxl)
library(xlsx)
dt <- read.csv("result.csv",header = TRUE) 
dt$` ` <- paste(rep(" ", 20), collapse = " ")
dt$mean <- signif(dt$mean,digits = 4)
dt$lower <- signif(dt$lower,digits = 4)
dt$upper <- signif(dt$upper,digits = 4)
colnames(dt)[4]="mean"
colnames(dt)[7]=""
colnames(dt)[1]="Drug target gene"
colnames(dt)[2]="P value"
colnames(dt)[3]="OR(95% CI)"
p <- forest(dt[,c(1, 2, 7)],
            est = dt$mean,     
            lower = dt$lower,   
            upper = dt$upper,   
            ci_column = 3,   
            ref_line = 1,
            xlim = c(0.75, 1.25))
p
g <- insert_text(p,
                 text = "plot_name",
                 col = 2:3,
                 part = "header",
                 gp = gpar(fontface = "bold"))
g
g <- edit_plot(g,
               row = c(12),
               col = 3,
               which = "ci",
               gp = gpar(col = "orange"))
g
pdf(file = "forest.pdf",width = 10,height = 20)
plot(g)
dev.off()