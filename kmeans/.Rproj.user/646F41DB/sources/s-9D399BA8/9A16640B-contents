library(ggplot2)
library(reshape2)
library(tidyr)
ab <- read.table("Bacteria+Archaea.txt_summary.txt", header = T, stringsAsFactors = F, sep = '\t')
ab <- melt(ab, id.vars = c("k", "MeanDistance"))
ggplot(ab,aes(x = factor(k), y = value, color = variable, group = variable)) +
  geom_line() +
  scale_x_discrete() +
  geom_point() +
  ggtitle("Bacteria+Archaea")

a <- read.table("Archaea.txt_summary.txt", header = T, stringsAsFactors = F, sep = '\t')
a <- melt(a, id.vars = c("k", "MeanDistance")) 
ggplot(a,aes(x = factor(k), y = value, color = variable, group = variable)) +
  geom_line() +
  scale_x_discrete() +
  geom_point() +
  ggtitle("Archaea")

b <- read.table("Bacteria.txt_summary.txt", header = T, stringsAsFactors = F, sep = '\t')
b <- melt(b, id.vars = c("k", "MeanDistance")) 
ggplot(b,aes(x = factor(k), y = value, color = variable, group = variable)) +
  geom_line() +
  scale_x_discrete() +
  geom_point() +
  ggtitle("Bacteria")

k5table <- read.csv("k5table.csv", stringsAsFactors = F, header = T) %>% tidyr::separate(col = info, into = c("species","optimal_growth_temperature","GCcontent"), sep = '\\|')
ggplot(k5table, aes(x = as.numeric(GCcontent), y = as.numeric(optimal_growth_temperature), color = factor(cluster), shape = factor(cluster))) +
  geom_point() +
  xlab("GC content") + ylab("Temperature") +
  xlim(0,100) + ylim(-5,100)
