library(ggplot2)
library(tidyverse)

bench <- read.csv("data_output/benchmark_run_3.csv")
bench <- read.csv("data_output/benchmark_run_pc2.csv")

bench_summary <- bench %>% 
  mutate("dimension" = cols*rows) %>% 
  filter(expr != 'pcaMethods_wraper(matrix, "bpca")') %>% 
  group_by(dimension, expr) %>% 
  summarise("mean_time" = mean(time/(1e9)))

bench_summary
class(bench)
View(bench_summary)
bench$time
sum(bench$time)/(1e9*60)
ggplot(bench_summary, aes(y = mean_time, x = dimension, color = expr))+geom_line(size = 2)+
  xlab("Number of values treated (dimension)")+
  ylab("Time in seconds")+
  labs(color = "Algorithm")+
  theme(legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1,"cm"),
    legend.text = element_text(size=15),
    legend.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 17))

ggsave("figures/benchmark.png", width = 18, height = 11)
