

load("/Users/gthb3/Dropbox/PhD students/Shahriari/My code/fig-fact-cost.RData")
#----------------------heatmap code--------------------------------------------------
Function.PATH = "/Users/gthb3/Dropbox/PhD students/Shahriari/My code"
library(ggplot2)
library(viridis)
library(reshape2)

library (patchwork)

cormat =  Aa[, , 1]
rownames(cormat) = paste("X", as.character(1:55),sep = " ")
colnames(cormat) = paste(expression(lambda), as.character(1:11),sep = "")

melted_cormat = melt(cormat)
F1 = ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_viridis(alpha = 1, begin = 0.1, end = 0.97) +
  theme(legend.title = element_blank()) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 0, size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title="First Group Factor loadings") +
  scale_x_discrete(labels=c(expression(lambda[1]), expression(lambda[2]),
                            expression(lambda[3]), expression(lambda[4]),
                            expression(lambda[5]), expression(lambda[6]),
                            expression(lambda[7]), expression(lambda[8]),
                            expression(lambda[9]), expression(lambda[10]),
                            expression(lambda[11])))

#library(RColorBrewer)
#coul <- colorRampPalette(brewer.pal(8, "Blues"))(50)

#heatmap(cormat, Colv = NA, Rowv = NA,
#        scale="column", col = coul, xlab="Factors",
#        ylab="Variables", main="First Group Factor loadings")
#


#----------------------heatmap code--------------------------------------------------


cormat <-  Aa[,, 2]
rownames(cormat) = paste("X", as.character(1:55),sep = " ")
colnames(cormat) = paste(expression(lambda), as.character(1:11),sep = "")

melted_cormat <- melt(cormat)
F2 = ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_viridis(alpha = 1, begin = 0.1, end = 0.97) +
  theme(legend.title = element_blank()) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 0, size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title="Second Group Factor loadings") +
  scale_x_discrete(labels=c(expression(lambda[1]), expression(lambda[2]),
                            expression(lambda[3]), expression(lambda[4]),
                            expression(lambda[5]), expression(lambda[6]),
                            expression(lambda[7]), expression(lambda[8]),
                            expression(lambda[9]), expression(lambda[10]),
                            expression(lambda[11])))

#----------------------heatmap code--------------------------------------------------

cormat <-  Aa[,,3]
rownames(cormat) = paste("X", as.character(1:55),sep = " ")
colnames(cormat) = paste(expression(lambda), as.character(1:11),sep = "")

melted_cormat <- melt(cormat)
F3 = ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_viridis(alpha = 1, begin = 0.1, end = 0.97) +
  theme(legend.title = element_blank()) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 0, size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title="Third Group Factor loadings") +
  scale_x_discrete(labels=c(expression(lambda[1]), expression(lambda[2]),
                            expression(lambda[3]), expression(lambda[4]),
                            expression(lambda[5]), expression(lambda[6]),
                            expression(lambda[7]), expression(lambda[8]),
                            expression(lambda[9]), expression(lambda[10]),
                            expression(lambda[11])))

#----------------------heatmap code--------------------------------------------------



cormat <- Aa[,,4]
rownames(cormat) = paste("X", as.character(1:55),sep = " ")
colnames(cormat) = paste(expression(lambda), as.character(1:11),sep = "")

melted_cormat <- melt(cormat)
F4 = ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_viridis(alpha = 1, begin = 0.1, end = 0.97) +
  theme(legend.title = element_blank()) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 0, size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title="Fourth Group Factor loadings") +
  scale_x_discrete(labels=c(expression(lambda[1]), expression(lambda[2]),
                            expression(lambda[3]), expression(lambda[4]),
                            expression(lambda[5]), expression(lambda[6]),
                            expression(lambda[7]), expression(lambda[8]),
                            expression(lambda[9]), expression(lambda[10]),
                            expression(lambda[11])))


postscript(paste(Function.PATH, '/F1.eps', sep=''), width=27, height=27)
F1 | F2
dev.off()

postscript(paste(Function.PATH, '/F2.eps', sep=''), width=27, height=27)
F3 | F4
dev.off()



ggsave(filename = "F1.eps",
       plot = P1,
       device = "eps",
       dpi = 1200,
       width = 15,
       height = 10,
       units = "cm")
