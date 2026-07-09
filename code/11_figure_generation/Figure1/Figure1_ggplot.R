library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(TDbook)
library(castor)
library(ape)

tree0 <- read.tree("mammal302.anno.BL_support.nwk")
node_support <- read.delim("mammal.branch_support.txt")
output_file <- "Figure1A.main_framework.three_category.pdf"

info <- read.delim("mammal_clades_info.txt")
rownames(info) <- info$Species
length(unique(info$Hightlighten))
info0 <- info[,c(1,2)];head(info0)
info0 <- na.omit(info0)

hightlight_nodes <- unique(info0$Hightlighten)
hightlight_nodes
tips <- data.frame(tree0$tip.label)

ancestor_nodes <- data.frame()
for (i in hightlight_nodes){
  info_species <- as.character(info0[info0$Hightlighten == i,1])
  anc <- get_mrca_of_set(tree0, info_species)
  anc_f <- data.frame(node = anc, label = i)
  ancestor_nodes <- rbind(ancestor_nodes, anc_f)
}
ancestor_nodes

ancestor_nodes$image <- paste(ancestor_nodes$label, ".png", sep = "")

p1 <- ggtree(tree0, layout="fan", size=0.15, open.angle=5) +
  geom_hilight(data=ancestor_nodes, mapping=aes(node=node),
               extendto=0.05, alpha=0.3, fill="grey", color="grey50",
               size=0.05) +
  geom_cladelab(data=ancestor_nodes, 
                mapping=aes(node=node, 
                            label=label),
                hjust=0.5,
                angle="auto",
                barsize=NA,
                horizontal=FALSE, 
                fontsize=2,
                fontface="italic")
traits <- read.delim("Mammal_ecotype.txt");head(traits)
traits$Trait <- factor(
  traits$Trait,
  levels = c("Marine", "Non-marine aquatic", "Semi-aquatic")
)

p2 <- p1 %<+% traits +
      geom_tippoint(mapping=aes(color=Trait, shape=Trait), 
            size=2,
            show.legend=TRUE) +
      scale_color_manual(
        name = "Trait",
        values = c(
          "Marine" = "#2C7FB8",
          "Non-marine aquatic" = "#D95F0E",
          "Semi-aquatic" = "#F1A340"
        ),
        drop = FALSE,
        na.translate = FALSE
      ) +
      scale_shape_manual(
        name = "Trait",
        values = c(
          "Marine" = 16,
          "Non-marine aquatic" = 15,
          "Semi-aquatic" = 17
        ),
        drop = FALSE,
        na.translate = FALSE
      ) +
      guides(
        color = guide_legend(override.aes = list(size = 3)),
        shape = guide_legend(override.aes = list(size = 3))
      ) +
      theme(
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
      #+
      #scale_color_manual(values=c("#1f78b4", "#e31a1c"))+
      #scale_size_manual(values = c(0.5))

ggsave(
  filename = output_file,
  plot = p2,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

message("Figure saved to: ", normalizePath(output_file, mustWork = FALSE))
