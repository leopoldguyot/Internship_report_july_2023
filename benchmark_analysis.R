library(ggplot2)
library(tidyverse)

bench2 <- read.csv("data_output/benchmark_run_pc2.csv")
bench2$nPCs <- 2
bench10 <- read.csv("data_output/benchmark_run_pc10.csv")
bench10$nPCs <- 10
bench25 <- read.csv("data_output/benchmark_run_pc25.csv")
bench25$nPCs <- 25

bench <- rbind(bench2, bench10, bench25)

bench_summary <- bench %>% 
  mutate("dimension" = cols*rows) %>% 
  group_by(dimension, expr, nPCs) %>% 
  summarise("mean_time" = mean(time/(1e9)))
sum(bench$time)/(1e9*60)
bench_summary$nPCs <- as.factor(bench_summary$nPCs)
ggplot(bench_summary, aes(y = mean_time, x = dimension, color = expr))+geom_line(size = 2)+
  facet_wrap(vars(nPCs))+
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
