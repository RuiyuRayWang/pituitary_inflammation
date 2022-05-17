library(ggplot2)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

qpcr <- read.csv("../data/qpcr/qpcr_before_ratio_raw.csv", row.names = 1)

cols_two = brewer.pal(8, "Set1")

df <- tibble(
  gene = character(),
  tissue = character(),
  treat = character(),
  rep = integer(),
  value = double()  ## Relative mRNA multipled by 1000
)

for (tissue in rownames(qpcr)){
  for (gene_treat_rep in colnames(qpcr)){
    gene = unlist(strsplit(gene_treat_rep, split = "_"))[1]
    treat = unlist(strsplit(gene_treat_rep, split = "_"))[2]
    rep = as.integer(unlist(strsplit(gene_treat_rep, split = "_"))[3])
    df <- df %>% tibble::add_row(gene = gene, tissue = tissue, treat = treat, rep = rep, value = qpcr[tissue,gene_treat_rep])
  }
}

## Plot order
df <- df %>% dplyr::mutate(tissue = forcats::fct_relevel(tissue, "gut", "muscle", "fat", "adrenal gland", "kidney", "spleen", "lung",
                                                         "cerebellum", "midbrain", "cortex", "liver", "heart", "pituitary")) %>%
  dplyr::mutate(treat = forcats::fct_relevel(treat, "LPS", "Saline"))

## Boxplots of relative mRNA level w.r.t. Actin, split by treatment
genes_of_interest <- unique(df$gene)
cols = c(brewer.pal(9, name = "Pastel1"), brewer.pal(4, name = "Pastel2"))

### Cartpt
g = "Cartpt"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(2,5)],
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(angle = 90, size = 12),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black")
    )
ggsave(filename = paste0(g,"_rel_exp.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 3, height = 10, dpi = 300, family = "Arial")
  
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    # axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Fgg
g = "Fgg"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value/1000) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Fndc4
g = "Fndc4"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Il1r2
g = "Il1r2"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value/10) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Insl6
g = "Insl6"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value/10) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Isg15
g = "Isg15"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value/10) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Lbp
g = "Lbp"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value/10) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0), breaks = c(0,10,20,30)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Lgals9
g = "Lgals9"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value/10) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Lox
g = "Lox"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Nptx2
g = "Nptx2"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Nrtn
g = "Nrtn"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.04,0), breaks = c(0,0.8)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Serpina3n
g = "Serpina3n"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value/1000) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")

### Vgf
g = "Vgf"
df %>% dplyr::filter(gene == g) %>% dplyr::mutate(value = value) %>%
  ggboxplot(
    x = "tissue", y = "value", fill = "treat",
    ylab = g,
    width = 0.7,
    orientation = "horizontal",
    palette = cols_two[c(5,2)],
    
  ) +
  scale_y_continuous(expand = c(0.03,0)) +
  theme_replace() +
  theme(
    axis.line = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    # panel.grid.major.x = element_line(color = "#CCCCCC"),
    panel.grid.major.y = element_line(color = "#CCCCCC"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none"
  )
ggsave(filename = paste0(g,"_rel_exp.clean.eps"), device = "eps", path = '../figures/Fig3/qpcr/',
       width = 1.3, height = 8, dpi = 300, family = "Arial")