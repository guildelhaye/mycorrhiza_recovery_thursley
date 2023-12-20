#### Script for analyses of ErM communities between N and Control treatment
#### 3 blocs and 4 plots in each block (2N and 2C)
#### Code written by Guillaume Delhaye - Last edited on 29-09-2023
library(tidyverse)
library(FactoMineR)
library(lme4)
library(vegan)
library(ggpubr)
library(indicspecies)
library(plotly)
library(ggfortify)

# ### Data preparation ----
load("data_thursley.RData")

# soil chemistry
soil_chem <- soil_chem |>
  mutate(block = as.factor(block),
         plot = as.factor(plot),
         treatment = as.factor(treatment),
         rep = as.factor(rep), 
         Soil_density = as.numeric(Soil_density), 
         Conduct = as.numeric(Conduct))

# leaf chemistry
leaf_chem <- leaf_chem |> 
  mutate(block = as.factor(block),
         plot = as.factor(plot), 
         treatment = as.factor(treatment),
         Tot_lf_NP = Tot_lf_N/T_lf_P)

# removed units in stem_ring file. ring_max, pith, all_rings_fb 
# and max_first_5yr_rings in µm
stem_ring <- stem_ring |>
  mutate(block = as.factor(block), 
         plot = as.factor(plot), 
         treatment = as.factor(treatment))

# change in stem ring size through time
stem_ring_change <- stem_ring %>%
  pivot_longer(c('2007':'2021'), 
               names_to = "year", values_to = "ring_size") |>
  mutate(block = as.factor(block),
         plot = as.factor(plot),
         stem = as.factor(stem),
         treatment = as.factor(treatment),
         year = as.numeric(year))

# Survey of the vegetation cover
veg_survey <- veg_survey |>
  mutate(block = as.factor(block),
         plot = as.factor(plot),
         treatment = as.factor(treatment),
         subplot = as.factor(subplot)) |>
  mutate_if(is.integer, as.numeric) |>
  rename(Ul_pt_20cm = Ul_.pt_20cm)

# Root fungal communities
root_com_raw <- root_com_raw %>%
  mutate(putative_function = as.factor(putative_function))

# Root colonisation 
root_colonisation

# Soil fungal communities
soil_com_raw

# Lichens
lichen_cover
lichen_abundance

#### PCA exploring soil chemistry ------------
pca_res <- prcomp(soil_chem |> select(-c(block, plot, rep, treatment)), 
                  scale. = TRUE)

p <- autoplot(pca_res, data = soil_chem, colour = 'treatment',
              loadings = TRUE, loadings.colour = 'black',
              loadings.label = TRUE, loadings.label.size = 4, 
              loadings.label.colour = "black") +
  stat_ellipse(aes(colour = treatment),  position = "identity") +
  theme_minimal() 
ggplotly(p)

setwd("../figures-tables")
ggsave("PCA_soil.pdf")

#### PCA exploring leaf chemistry ----
leaf_ch <- leaf_chem %>% select(-c(block, plot, col_prop))
leaf_pca <- PCA(leaf_ch, quali.sup = 1)
coord_leaf <- data.frame(leaf_chem[,1:3], leaf_pca$ind$coord)

(pca_leaf_var <- plot(leaf_pca, choix = "var", graph.type = "ggplot", title = ""))

(pca_leaf_ind <- ggplot(coord_leaf, aes(y = Dim.2, x = Dim.1, colour = treatment)) +
  geom_point() +
  theme_bw() + 
  stat_ellipse() +
  scale_color_manual(values=c("blue", "red"), name = "Treatment", labels = c('Control', 'N-addition')) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
)

#### Comparison soil chemistry between treatments ----
# between C and N for soil variables with a model nesting samples in plots in blocs
mod.soilchem <- function(y) {
  summary(lmerTest::lmer(y ~ treatment + (1 | block/plot), 
                         data = soil_chem))}

mod.soilchem(soil_chem$pH) 
mod.soilchem(soil_chem$H20)
mod.soilchem(soil_chem$Tot_soil_N)
mod.soilchem(soil_chem$Tot_soil_C)
mod.soilchem(soil_chem$CN_ratio)
mod.soilchem(soil_chem$P)
mod.soilchem(soil_chem$K)
mod.soilchem(soil_chem$Mg) 
mod.soilchem(soil_chem$Nitrate)
mod.soilchem(soil_chem$Ammonium)
mod.soilchem(soil_chem$Dry_matter)
mod.soilchem(soil_chem$Soil_density)
mod.soilchem(soil_chem$Na)
mod.soilchem(soil_chem$Ca)
mod.soilchem(soil_chem$Conduct)
mod.soilchem(soil_chem$CEC) 


#### Comparison leaf Chemistry between treatments ----
mod.leafchem <- function(y) {summary(lmerTest::lmer(y ~ treatment + (1 | block), data = leaf_chem))}

mod.leafchem(leaf_chem$col_prop)
mod.leafchem(leaf_chem$Tot_lf_N)
mod.leafchem(leaf_chem$T_lf_P) 
mod.leafchem(leaf_chem$Tot_lf_C)
mod.leafchem(leaf_chem$Tot_lf_CN)
mod.leafchem(leaf_chem$Tot_lf_NP)

#### Comparison stem amatomy between treatments ----
mod.anat <- function(y){
  summary(lmerTest::lmer(y ~ treatment + (1 | block/plot), data = stem_ring))}

mod.anat(stem_ring$n_rings_max)
mod.anat(stem_ring$rings_max)
mod.anat(stem_ring$diam_mm)
mod.anat(stem_ring$circ_mm)
mod.anat(stem_ring$pith) 
mod.anat(stem_ring$all_rings_fb)
mod.anat(stem_ring$max_first_5yr_rings)

#### Vegetation physiognomy analyses #----
# Total density of Calluna
veg_survey <- veg_survey %>% 
  mutate(Cal_pt_tot = Cal_pt_1cm + Cal_pt_20cm) 

# Helper function to model the response of vegetation physiognomy
veg_mod <- function(y) {
  summary(lmerTest::lmer(y ~ treatment + (1 | block), data = veg_survey))}

veg_mod(veg_survey$Cal_pt_1cm)
veg_mod(veg_survey$Cal_pt_20cm)
veg_mod(veg_survey$Cal_hi_cm) 
veg_mod(veg_survey$Cal_rand_ht) 
veg_mod(veg_survey$Camp_int_1cm) 
veg_mod(veg_survey$Hyp_jut_1cm) 
veg_mod(veg_survey$Cal_pt_tot)

#### Analyse of root colonisation between treatments ----
## This is a bit more complicated because the response variable is a 
## proportion of cells colonised by ErM out of a total number of roots.
## The appartenance of each root to a specific plant was not noted and therefore
## we have roots per plot per bloc crossed with the two treatments.
## Use a GLMM with a binomial response variable and root nested in plot nested in bloc
root_colonisation <- root_colonisation |>
  select(-plant) |># not use the plant level because not complete
  mutate(bloc = as.factor(bloc), 
         plot = as.factor(plot), 
         root = as.factor(root), 
         treatment = as.factor(treatment))
str(root_colonisation)

# significant difference, but wrong because considers all points as independant. 
summary(glm(cbind(colonised, not_colonised) ~ treatment,
                 data = root_colonisation,
                 family = binomial, 
                 weights = colonised + not_colonised))

# right organisation of the data with a random level for plotand bloc
summary(glmer(cbind(colonised, not_colonised) ~ treatment + (1 | bloc/plot),
                    data = root_colonisation,
                    family = binomial, 
                    weights = colonised + not_colonised))

## Figure
root_colonisation |> 
  rename(Block = bloc) |>
  mutate(colonisation_proportion = colonised/(colonised + not_colonised)) |>
  ggplot(aes(y = colonisation_proportion, x = treatment)) +
  #geom_jitter() + 
  geom_boxplot(aes(colour = Block)) +
  ylab("Proporion of colonisation") +
  theme_minimal()
ggsave("colonisation_proportion.jpeg")

## correlation between soil chem and colonisation
tab_soil_col <- soil_chem |> 
  select(-c(rep, treatment)) |>
  unite(plot, block, plot,  sep = "")|>
  group_by(plot) |>
  summarise_all(mean) |>
  ungroup() |>
  mutate(plot = as.factor(plot))

prop_col <- root_colonisation |> 
  mutate(prop = colonised / (colonised + not_colonised)) |>
  select(plot, prop) |>
  group_by(plot) |>
  summarise_all(mean) |> 
  ungroup()

tab_soil_col <- prop_col|>
  left_join(tab_soil_col)
  

cor <- round(cor(tab_soil_col[,-c(1)]),2)
write.csv(cor, "../figures-tables/correl_test_colonisation.csv")

####  Fungal Communities Analyses ------------------------
## Root fungal communities ####
design <- data.frame(sample = c("TH5A", "TH5B", "TH5C" , "TH5D", "TH6A", "TH6B", "TH6C", "TH6D", "TH7A", "TH7B", "TH7C", "TH7D"),
                     block = c(5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7),
                     plot = c("A", "B", "C", "D", "A", "B", "C", "D","A", "B", "C", "D"),
                     treatment = c("C", "C", "N", "N", "N", "N", "C", "C", "N", "C", "N", "C"))

root_com <- root_com_raw %>%
  select(uniqueSHs, tax_ID, putative_function, 
         TH5A, TH5B, TH5C, TH5D, 
         TH6A, TH6B, TH6C, TH6D,
         TH7A, TH7B, TH7C, TH7D) %>%
  group_by(uniqueSHs, tax_ID, putative_function) %>%
  summarize_all(sum) %>%
  ungroup()

## Root ErM communities
root_erm <- root_com %>% 
  filter(putative_function == "ErM") 

# Root potential ErM communities
root_pot_erm <- root_com %>% 
  filter(putative_function %in% c("ErM","potential ErM"))

## Compare proportions of ErM relative to others in both treatments 
erm_sum_root <- root_com %>% 
  select(-c(uniqueSHs, tax_ID)) %>%
  group_by(putative_function) %>%
  summarize_all(sum) %>%
  ungroup() 

prop_ErM_roots = cbind(design, 
                       prop_erm =  t(erm_sum_root[1,-1] / colSums(erm_sum_root[,-1])),
                       prop_pot_erm = colSums(erm_sum_root[c(1, 3),-1]) / colSums(erm_sum_root[,-1]))

## test for differences between treatments
(summary(lmerTest::lmer(prop_erm ~ treatment + (1 | block), data = prop_ErM_roots))) #ns
# No difference in the proportion of ErM in C and N
(summary(lmerTest::lmer(prop_pot_erm ~ treatment + (1 | block), data = prop_ErM_roots))) #ns
# No difference in the proportion of potential ErM in C and N

# Multivariate comparison of community composition between treatment
# Test permanova and draw graph
y = root_com# change root_com for root_erm or root_pot_erm depending on which community we want to use

# prepare data
x =  y |>
  select(-c(uniqueSHs, tax_ID, putative_function)) |> 
  t() |> 
  tibble() |> 
  mutate(sample = design$sample)|> 
  filter(rowSums(across(where(is.numeric)))!=0)

# keep the plots with non zero communities
design_filt = design[design$sample %in% x$sample,]

# Hellinger transform
x_hel <- decostand(x |> select(-sample), method = "hellinger")

# permanova 
ad <- adonis2(x_hel ~ treatment, 
              data = design_filt, 
              permutation = 999)
R2 = round(ad$R2[1], 2); f = round(ad$F[1], 3); p = round(ad$`Pr(>F)`[1], 3)

# NMDS
NMDS <- metaMDS(x_hel , distance = "bray", k = 2)
d <- data.frame(NMDS$points[,]) %>% cbind(design_filt)

# Plot
ggplot(data = d, aes(y = MDS2, x = MDS1, colour = treatment)) +
  geom_point() +
  stat_ellipse() +
  theme_classic() +
  labs(subtitle = paste0("PERMANOVA: F = ", f, ", R2 = ", R2, ", p = ", p)) +
  scale_color_manual(values=c("blue", "red"),
                     name = "Treatment",
                     labels = c('Control', 'N-addition'))


## soil fungal communities ####
#  filtered at 100 reads beforehand
# design (two samples per plot instead of one for roots)
env <- data.frame(
  block = rep(c(5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7), each = 2),
  plot = rep(c("A", "B", "C", "D", "A", "B", "C", "D","A", "B", "C", "D"), each = 2),
  treatment = rep(c("C", "C", "N", "N", "N", "N", "C", "C", "N", "C", "N", "C"), each = 2),
  rep = rep(c(1, 2)))

soil_com <- soil_com_raw %>%
  select(uniqueSHs, tax_ID, putative_function, 
         TH5A1, TH5A2, TH5B1, TH5B2, TH5C1, TH5C2, TH5D1, TH5D2,
         TH6A1, TH6A2, TH6B1, TH6B2, TH6C1, TH6C2, TH6D1, TH6D2,
         TH7A1, TH7A2, TH7B1, TH7B2, TH7C1, TH7C2, TH7D1, TH7D2) %>%
  group_by(uniqueSHs, tax_ID, putative_function) %>%
  summarize_all(sum) %>%
   ungroup() 

## soil erm communities
soil_erm <- soil_com %>% 
  filter(putative_function == "ErM")

# soil potential ErM communities
soil_pot_erm <- soil_com %>% 
  filter(putative_function %in% c("ErM","potential ErM"))

## Compare proportions of ErM relative to others in both treatments
erm_sum_soil <- soil_com %>% 
  select(-c(uniqueSHs, tax_ID)) %>%
  group_by(putative_function) %>%
  summarize_all(sum) %>%
  ungroup() 

prop_ErM_soil = data.frame(cbind(treatment = as.factor(rep(design$treatment, each = 2)), 
                                  block = rep(design$block, each = 2),
                                  prop_erm =  t(erm_sum_soil[1,-1] / colSums(erm_sum_soil[,-1])),
                                  prop_pot_erm = colSums(erm_sum_soil[c(1, 3),-1]) / colSums(erm_sum_soil[,-1]))) |>
  rename(prop_erm = V3)

(summary(lmerTest::lmer(prop_erm ~ treatment + (1 | block), data = prop_ErM_soil))) # NS
(summary(lmerTest::lmer(prop_pot_erm ~ treatment + (1 | block), data = prop_ErM_soil))) # NS
# No differences in proportions of ERM or potential ERM between both treatments 

# Multivariate comparison of community composition between treatment
# Test permanova and draw graph
x = soil_pot_erm |> 
  select(-c(uniqueSHs, tax_ID, putative_function)) |> 
  t() 
#colnames(x) <- soil_pot_erm$uniqueSHs

# Hellinger transform
x_hel <- decostand(x, method = "hellinger")

# permanova 
ad <- adonis2(x_hel ~ treatment, data = env)
R2 = round(ad$R2[1], 2); f = round(ad$F[1], 3); p = round(ad$`Pr(>F)`[1], 3)

# NMDS
NMDS <- metaMDS(x_hel , distance = "bray", k = 6)
d <- data.frame(NMDS$points[,]) %>% cbind(env) |>
  # add line to identify both replicates of one site
  unite(replic, c(block, plot), sep = "_")

# Plot
soil_pot_erm_g <- ggplot(data = d, aes(y = MDS2, x = MDS1, colour = treatment)) +
  geom_point(size = 4) +
  stat_ellipse() +
  theme_classic() +
  labs(subtitle = paste0("PERMANOVA: F = ", f, ", R2 = ", R2, ", p = ", p)) 

################################################################################
#### Analyses of lichen communities ----
lichen_diversity <- data.frame(
  lichen_abundance[1:4], 
  shannon = diversity(lichen_abundance[,c(5:length(lichen_abundance))], 
                      index = "shannon"),
  inv_simpson = diversity(lichen_abundance[,c(5:length(lichen_abundance))], 
                          index = "invsimpson")) |>
  left_join(lichen_cover)

## Test for differences between treaments
summary(lmerTest::lmer(Abundance ~ treatment + (1 | Block), 
                       data = lichen_diversity))
summary(lmerTest::lmer(Richness ~ treatment + (1 | Block), 
                       data = lichen_diversity))

summary(lmerTest::lmer(shannon ~ treatment + (1 | Block), 
                       data = lichen_diversity))
summary(lmerTest::lmer(inv_simpson ~ treatment + (1 | Block), 
                       data = lichen_diversity))

## Permanova testing the effect of treatment on communities composition
x = lichen_abundance[,5:length(lichen_abundance)] 

# Hellinger transform
x_hel <- x# decostand(x, method = "hellinger")

# permanova 
ad <- adonis2(x_hel ~ treatment, data = lichen_abundance)
R2 = round(ad$R2[1], 2); f = round(ad$F[1], 3); p = round(ad$`Pr(>F)`[1], 3)

# NMDS
NMDS <- metaMDS(x_hel , distance = "bray", k = 6)
d <- data.frame(NMDS$points[,]) %>% 
  cbind(lichen_abundance |>  select(Block, Subplot, treatment)) #|>
  # add line to identify both replicates of one site
 # unite(replic, c(block, plot), sep = "_")

# Plot
lichen_com <- ggplot(data = d, 
                         aes(y = MDS2, x = MDS1, colour = treatment)) +
  geom_point(size = 4) +
  stat_ellipse() +
  theme_classic() +
  labs(subtitle = paste0("PERMANOVA: F = ", f, ", R2 = ", R2, ", p = ", p)) + 
  theme(legend.position = "none")
ggsave("nmds_lichens.pdf")

## Indic species anaysis
is_lichen <- multipatt(x = lichen_abundance[,5:length(lichen_abundance)], 
                       cluster = lichen_abundance$treatment, 
                       func = "IndVal.g",
                       control = how(nperm=999))

summary(is_lichen)

lichen_abundance |> 
  select(treatment, fimbriata:ciliata) |>
  pivot_longer(!treatment, values_to = "abundance") |>
  ggplot(aes(y = log(abundance + 1), x = treatment, colour = treatment)) +
  geom_violin() +
  geom_jitter(height = 0.1, width = 0.1) +
  facet_wrap(~name) +
  ylab("abundance (log)") + 
  xlab("") +
  theme_minimal() + 
  theme(legend.position = "none")
ggsave("boxplot_lichen_specificity.jpeg")

### Test of specificity using the Clamtest of Chezon et al 2011.
library(vegan)
CT <- clamtest(comm = lichen_abundance[,5:length(lichen_abundance)],
               groups = lichen_abundance[,4],
               specialization = 2/3, 
               alpha = 0.05/21) # divide by the number of tests. Can be changed
CT
### END ####
