# generic colors to use for divergent and sequential continuous color scales
library(viridis)
library(colorjam)

# generate ggplot colors
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


# colors for heatmaps, as well as the FDR on the original knee plot for libraries
heatmap_colors <- viridis(1000)

# colors for scatter intensity plots (count depth, marker gene expression, etc.)
scatter_intensity_colors <- viridis(1000)

# colors for dot plots
dotplot_colors <- viridis(1000)

# colors for scatter plots of feature intensity
feature_intensity_colors <- c("grey90", "blue")

# need discrete colors that are divergent but distinct from viridis
# colors for cell and gene filters
filter_colors <- c("#F8766D", # fail (red)
                   "dodgerblue") # pass (blue)

# metadata colors
# colors for plCoA regions
region_colors <- c("#FF0080", # ASt (magenta)
                   "#FFCE00", # CeA (yellow)
                   "#48EDED", # DS (cyan)
                   "#00FF00") # TS (green)

# colors for individual library batches
ASt_colors <- c("#FF4DA6", "#B3005A") # hot magenta/pink
CeA_colors <- c("#FF8E00", "#FFB900", "#FFDD4D") # warm yellow/orange
DS_colors <- c("#48EDC4", "#15D4D4", "#8EF4F4", "#48C4ED") # cool cyan/blue
TS_colors <- c("#00B300","#4DFF4D", "#00E600") # neon/green
batch_colors <- c(ASt_colors, CeA_colors, DS_colors, TS_colors)

# cell type colors
# neurons vs glia
neuron_color <- "dodgerblue"
glia_color <- "forestgreen"

tissue_colors <- c(neuron_color, # neurons
                   glia_color) # glia

# first-level cell classes (glut, GABA, astrocytes, microglia, etc)
# for run 1
class_colors <- c("#FF6666", # GABA-D1
                   "#93F998", # GABA-D2
                   "#053061", # non-D1/2 neuron
                   "forestgreen", # astrocytes
                   "wheat3", # microglia
                   "orange", # OPC
                   "orangered", # NFOL
                   "darkorange2", # OLG
                   "gold3", # endothelial
                   "gold1", # mural
                   "gold4") # ependymal

class_colors2 <- c("#FF6666", # GABA-D1
                  "#93F998", # GABA-D2
                  "#053061", # non-D1/2 neuron
                  "forestgreen", # astrocytes
                  "wheat3", # microglia
                  "orange", # OPC
                  "darkorange2", # OLG
                  "gold3", # endothelial
                  "gold1", # mural
                  "gold4", # ependymal
                  "lightgrey") # low depth neuron

class_colors3 <- c("forestgreen", # astrocytes
                   "gold3", # endothelial
                   "gold4", # ependymal
                   "#FF6666", # GABA-D1
                   "#93F998", # GABA-D2
                   "#053061", # non-D1/2 neuron
                   "wheat3", # microglia
                   "gold1", # mural
                   "orangered", # NFOL
                   "darkorange2", # OLG
                   "orange") # OPC

nond12_color <- "#053061"
nond12_colors <- c("#3944BC", # Blue
                   "#59788E", # Stone
                   "#489090", # Sky
                   "#0492C2", # Arctic
                   "#016064", # Ocean
                   "#1F456E", # Aegean
                   "#663046", # Grape
                   "#2C041C") # Wine

# for subtypes of neurons
d1_color <- "#FF6666" # D1, shades of light green
d1_colors <- c("#9DC183", # Sage Green
               "#D0F0C0", # Tea Green
               "#98FB98") # Mint Green

d2_color <- "#93F998" # D2, shades of red
d2_colors <- c("#FA8072", # Salmon
               "#ED2939", # Imperial red
               "#E0115F") # Ruby

# for subtypes of glia
astro_color <- "forestgreen"
astro_colors <- c("#50C878", # Emerald Green
                  "#3F704D", # Hunter Green
                  "#043927") # Sacramento Green

micro_color <- "wheat3"
micro_colors <- c("#CDBA96", # Wheat3
                  "#9A7B4F", # Tortilla
                  "#4A3728") # Cedar

type_colors <- c(astro_colors, # astrocyte colors
                 "gold3", # endothelial
                 "gold4", # ependymal
                 d1_colors, # D1
                 d2_colors, # D2
                 micro_colors, # Micro
                 "gold1", # mural
                 nond12_colors, # nonD12 neurons
                 "orangered", # NFOL
                 "darkorange2", # OLG
                 "orange") # OPC

type_colors_ast <- c(astro_colors, # astrocyte colors
                     "gold3", # endothelial
                     "gold4", # ependymal
                     d1_colors, # D1
                     d2_colors[2:3], # D2
                     micro_colors, # Micro
                     "gold1", # mural
                     nond12_colors[1:7], # nonD12 neurons
                     "orangered", # NFOL
                     "darkorange2", # OLG
                     "orange") # OPC

type_colors_cea <- c(astro_colors, # astrocyte colors
                     "gold3", # endothelial
                     "gold4", # ependymal
                     d1_colors, # D1
                     d2_colors, # D2
                     micro_colors, # Micro
                     "gold1", # mural
                     nond12_colors, # nonD12 neurons
                     "orangered", # NFOL
                     "darkorange2", # OLG
                     "orange") # OPC

region_class_colors <- NULL
for (i in 1:length(class_colors)){
  for (j in 1:length(region_colors)){
    region_class_colors <- c(region_class_colors, substr(blend_colors(c(class_colors[i], region_colors[j])), 1, 7))
  }
}

region_type_colors <- NULL
for (i in 1:length(type_colors)){
  for (j in 1:length(region_colors)){
    region_type_colors <- c(region_type_colors, substr(blend_colors(c(type_colors[i], region_colors[j])), 1, 7))
  }
}

region_d1_colors <- NULL
for (i in 1:length(d1_colors)){
  for (j in 1:length(region_colors)){
    region_d1_colors <- c(region_d1_colors, substr(blend_colors(c(d1_colors[i], region_colors[j])), 1, 7))
  }
}

region_d2_colors <- NULL
for (i in 1:length(d2_colors)){
  for (j in 1:length(region_colors)){
    region_d2_colors <- c(region_d2_colors, substr(blend_colors(c(d2_colors[i], region_colors[j])), 1, 7))
  }
}

region_nond12_colors <- NULL
for (i in 1:length(nond12_colors)){
  for (j in 1:length(region_colors)){
    region_nond12_colors <- c(region_nond12_colors, substr(blend_colors(c(nond12_colors[i], region_colors[j])), 1, 7))
  }
}

region_astro_colors <- NULL
for (i in 1:length(astro_colors)){
  for (j in 1:length(region_colors)){
    region_astro_colors <- c(region_astro_colors, substr(blend_colors(c(astro_colors[i], region_colors[j])), 1, 7))
  }
}

region_micro_colors <- NULL
for (i in 1:length(micro_colors)){
  for (j in 1:length(region_colors)){
    region_micro_colors <- c(region_micro_colors, substr(blend_colors(c(micro_colors[i], region_colors[j])), 1, 7))
  }
}

d1_d2_colors_level1 <- c(d1_color, d2_color, nond12_color)
#d1_d2_colors_complete <- c(d1_colors, d2_colors, other_colors)

# glutamatergic, GABAergic, and non-neuronal
#tissue_neuron_colors <- c(d1_d2_colors,
#                          glia_color) # non-neuronal

batch_color_list <- c("batch_colors", "region_colors")
tissue_color_list <- c("tissue_colors", "class_colors")

ast_color_list <- c("tissue_colors", "class_colors3", "type_colors_ast")
cea_color_list <- c("tissue_colors", "class_colors3", "type_colors_cea")
