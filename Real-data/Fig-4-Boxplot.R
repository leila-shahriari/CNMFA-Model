library(readr)
library(tidyverse)

cost_of_living = read_csv("Dropbox/PhD students/Shahriari/My code/cost-of-living.csv")

cost_of_living = cost_of_living %>% 
  select(4:58) %>%
  mutate_at(1:55, ~(scale(.) %>% as.vector)) %>% 
  pivot_longer(cols = 1:55, names_to = "Variable", values_to = "Amount")

cost_of_living$Variable <- factor(cost_of_living$Variable, levels = unique(cost_of_living$Variable))


Function.PATH = "/Users/gthb3/Dropbox/PhD students/Shahriari/My code"

postscript(paste(Function.PATH, '/Box.eps', sep=''), width=27, height=27)
ggplot(cost_of_living, aes(x = Variable, y = Amount, fill = Variable)) + 
  geom_boxplot() +
  coord_cartesian(ylim=c(-3, 20))+
  theme(legend.position = "")
dev.off()
