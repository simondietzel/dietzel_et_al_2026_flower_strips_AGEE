# This R-script includes the R code to download and analyze the data of the manuscript

# "Synergistic effects of flower strips and landscape complexity buffer arthropods against warming and drought - A meta-analysis". https://doi.org/xxx.yyy

# Dietzel et al. 2026, Agriculture, Ecosystems & Environment
# Version 1.0


# Before you start, please have a look at the README file.
# Before you start, please have a look at the metadata file.

# The orchaRd package can be downloaded from the github repository of Daniel Noble: https://daniel1noble.github.io/orchaRd/
install.packages("pacman")
rm(list = ls())
devtools::install_github("daniel1noble/orchaRd", ref = "main", force = TRUE)
pacman::p_load(
  devtools, tidyverse, metafor, patchwork, R.rsp, orchaRd, emmeans, ape, phytools, flextable
  )

# Restore package versions: uncomment if needed
# library(renv)
# renv::restore()

# libraries----
library(zen4R)
library(tidyverse)
library(scales)
library(conflicted)
library(janitor)
library(metafor)
library(orchaRd)
library(patchwork)
library(ggimage)
library(png)
library(performance)
library(MuMIn)
library(grid)


# solve conflicting functions----
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("update", "stats")

# data----
## download from Zenodo----
zenodo <- ZenodoManager$new()
download_zenodo(path = "data/",
                "10.5281/zenodo.18016073",
                timeout = 3600)

## read data----
beetlespider <- read.csv2("data/beetlespider_shared_data.csv", dec = ".") %>%
  mutate_if(is.character, as.factor)
str(beetlespider)
beetlespider

## project colors (color-blind friendly)----
colors_beetlespider <- c("#26496B",
                         "#605E78",
                         "#8A616D",
                         "#B16364",
                         "#EB8867",
                         "#E9CBA3")
ramp_beetlespider <- colorRampPalette(colors_beetlespider)
ramp_beetlespider(100)
barplot(1:100, col = ramp_beetlespider(100))


#________________________________----
# >> BEETLE ABUNDANCE----
# create subset for beetle abundance
beetle_abund <-
  beetlespider %>%
  filter(tax_group == "beetles") %>% # please make sure to use the function dplyr::filter throughout the code
  filter(!is.na(t_abund)) %>%
  select(-t_rich, -t_rich_sd, -c_rich, -c_rich_sd) %>% # remove richness data
  droplevels() %>% # drop spider taxa levels
  as.data.frame()

## 0 SMD for abundance-----
beetle_abund <-
  escalc(
    measure = "SMD",# SMD: the positive bias in the standardized mean difference automatically corrected within the function, yielding Hedges' g (Hedges, 1981)
    n1i = t_sites,# sample size treatment
    n2i = c_sites,# sample size control
    m1i = t_abund,# mean treatment
    m2i = c_abund,# mean control
    sd1i = t_abund_sd,# SD treatment
    sd2i = c_abund_sd, # SD control
    data = beetle_abund,
    slab = lit_id,
    vtype = "LS",
    correct = TRUE, # bias correction for small sample sizes
    var.names = c("yi_abund_beetles", "vi_abund_beetles"),
    add.measure = FALSE
  )

beetle_abund <-
  beetle_abund %>%
  filter(!is.na(yi_abund_beetles)) # exclude NAs in the data
beetle_abund

# 1 Nullmodel rma----
m_beetle_abund_null_0 <-
  rma(
    yi_abund_beetles,
    vi_abund_beetles,
    method = "REML",
    data = beetle_abund,
    slab = lit_id
  )
summary(m_beetle_abund_null_0) # AIC 292.9347

# Forest plot
forest(
  m_beetle_abund_null_0,
  cex = .7,
  main = "Beetle abundance",
  showweights = TRUE,
  slab = paste(
    beetle_abund$lit_id,
    beetle_abund$authors,
    beetle_abund$taxon,
    sep = " - "
  ),
  order = "obs",
  pch = 10
)
# Study 301.3 (Boetzl et al., carabidae) extreme far on the negative side, otherwise fine

# QQ-Plot
res <- rstandard(m_beetle_abund_null_0)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)
# looks a little out of line in the lower left corner

# Baujat plot to identify extreme outliers
baujat(m_beetle_abund_null_0, symbol = "slab") # Extreme outlier: effect 301.3 = Boetzl. et al. (Carabidae)

## 1.1 exclude outlier ID 301.3----
beetle_abund_outl <-
  beetle_abund %>%
  filter(yi_abund_beetles > -5, yi_abund_beetles < 5)

## 1.2 compare models----
# Model without outlier
m_beetle_abund_null_00 <-
  rma(
    yi_abund_beetles,
    vi_abund_beetles,
    method = "REML",
    data = beetle_abund_outl,
    slab = lit_id
  )
summary(m_beetle_abund_null_00)

# Forest plot
forest(
  m_beetle_abund_null_00,
  cex = .7,
  main = "Beetle abundance",
  showweights = TRUE,
  order = "obs",
  slab = paste(
    beetle_abund_outl$lit_id,
    beetle_abund_outl$authors,
    beetle_abund_outl$taxon,
    sep = " "
  ),
  pch = 10
)
# Cleaner distribution of effects on negative and positive side

# QQ-Plot
res <- rstandard(m_beetle_abund_null_00)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)
# Looks fine now

# Baujat plot
baujat(m_beetle_abund_null_00, symbol = "slab")
# Still extremes to be addressed in sensitivity analysis

# Compare models
compare_performance(m_beetle_abund_null_0, m_beetle_abund_null_00)
summary(m_beetle_abund_null_0)
summary(m_beetle_abund_null_00)
# Minimal change in estimate, variances (I², tau²), and p-value after exclusion of outlier
# >> Decision: Exclude effect 301.3. Reasons:
## 1. Extreme in Baujat plot on x axis compared to next values (110, 66.6, 66.8)
## 2. Extreme low SMD in forest plot < -5
## 3. Contortion in QQ plot
## 4. Exclusion does not affect overall model fit, variance and outcome


# 2 Sensitivity analysis----
## Aggregate data for sensitivity tests (can only be done with rma-objects) by considering study-level variance
beetle_abund_agg <-
  aggregate(beetle_abund_outl, cluster = lit_id, struct = "ID")

# rma Model with aggregated data
m_ba_agg <-
  rma(yi_abund_beetles,
      vi_abund_beetles,
      data = beetle_abund_agg,
      slab = lit_id)
summary(m_ba_agg)


# ## 2.1 Failsafe number----
# Rosenthal's FS N
fsn(yi_abund_beetles,
    vi_abund_beetles,
    method = "Rosenthal",
    data = beetle_abund_agg)
# 375
# > minimum = n x 5 + 10 (here = 160)


## 2.2 Baujat plot----
baujat(m_ba_agg, symbol = "slab")
# IDs 66 and 110 need to be checked in leave 1 out

## 2.3 Leave 1 Out----
ba_l1o <- leave1out(m_ba_agg)

ba_l1o %>%
  as.data.frame() %>%
  rownames_to_column(var = "study") %>%
  ggplot(aes(x = reorder(study, -estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.2) +
  geom_hline(
    yintercept = coef(m_ba_agg),
    linetype = "dashed",
    color = "red"
  ) +
  coord_flip() +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = .05)) +
  labs(title = "Leave-one-out Analysis", x = "Removed study", y = "Estimate without study") +
  theme_minimal() +
  theme(
    title = element_text(size = 18),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15)
  )

## 2.4 Cook´s D, Hat, rstudent----
inf_ba <- influence(m_ba_agg)

inf_ba$inf %>%
  as.data.frame() %>%
  rownames_to_column(var = "lit_id") %>%
  select(lit_id, rstudent, cook.d, hat) %>%
  pivot_longer(cols = 2:4,
               names_to = "variable",
               values_to = "value") %>%
  ggplot(aes(x = lit_id, y = value, group = variable)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ variable, scales = "free", ncol = 1) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

inf_ba$inf %>%
  as.data.frame() %>%
  rownames_to_column("lit_id") %>%
  filter(abs(rstudent) > 2 | cook.d > 0.2 | hat > 2 / nrow(.))
# Extremes: Studies 110 and 66

## 2.5 Model without outliers----
m_ba_agg_update <-
  update(m_ba_agg, subset = !(lit_id %in% c("66", "110")))
summary(m_ba_agg_update)

## 2.6 Compare models----
compare_performance(m_ba_agg, m_ba_agg_update)
summary(m_ba_agg_update)
summary(m_ba_agg)
# IDs 66 and 110 have notable influence on the model estimate, indicated by the sensitivity checks.
# Decision: Only small change in estimates (0.01) indicates that studies 66 and 110 compensate each other due to their opposite effect directions, which speaks for the robustness of the model itself. Therefore, they are not excluded.


# 3 Publication bias----
## 3.1 Egger´s test----
funnel(m_ba_agg)
regtest(m_ba_agg)
# Test for Funnel Plot Asymmetry: z = 0.3865, p = 0.6991
## no significant assymetry

## 3.2 Trimfill----
trim_ba <-
  trimfill(m_ba_agg)
summary(trim_ba)
# Estimated number of missing studies on the left side: 1 (SE = 3.4475)
funnel(trim_ba,
       bg = c("black", "red"),
       pch = 21,
       main = "Beetle abundance\nTrim-and-Fill")
legend(
  "topright",
  legend = c("Original", "Trim-and-Fill"),
  pt.bg = c("black", "red"),
  pch = 21,
  col = "black",
  cex = 1.2,
  # inset = c(-0, -0.005),
  x.intersp = .4,
  y.intersp = .8
)
# One missing study in the lower left part of the funnel means one study with a negative effect size of around -2 and low replication (high variance, SE = 1.3) is missing -> No indication for strong publication bias.


# 4 Multivariate models----
### 4.0 Variance structure----
ba_variance <-
  rma.mv(
    yi_abund_beetles,
    vi_abund_beetles,
    random = ~ 1 | lit_id / year / taxon / effect_id,
    data = beetle_abund_outl,
    method = "REML"
  )
summary(ba_variance)
# Q-test significant: Q(df = 94) = 270.7218, p-val < .0001. Based on the total weighted deviation of effect sizes from the model mean. Tests whether observed effect sizes vary more than expected by sampling error alone, ignoring the random structure.
# Intercept (0.3654) still significant (p = 0.0203) -> Total variance is not fully explained by random structure.

# Random level τ² (estim)	Interpretation:
# lit_id	0.4621: high variance among studies
# lit_id/year	0.000: no variance among years in studies
# lit_id/year/taxon	0.2149: moderate variance among taxa within studies
# lit_id/taxon/effect_id 	0.0082: minimal residual heterogeneity (within taxon–study combinations)
# “Q” < 0.0001 = based on the total weighted deviation of effect sizes from the model mean. tests whether observed effect sizes vary more than expected by sampling error alone, ignoring the random structure

# create and export table for appendix (A.Table 1)
ba_var_table <-
  data.frame(
    Model = "Beetle abundance",
    Factor_level = unname(ba_variance$s.names),
    Variance = "σ ^2",
    Estimate = round(ba_variance$sigma2, 3),
    SD = round(sqrt(ba_variance$sigma2), 3),
    Levels = ba_variance$s.nlevels,
    k = ba_variance$k,
    Model_inclusion = c("mandatory", "excluded", "in fixed model term", "excluded")
  )
ba_var_table

write.table(
  ba_var_table,
  "output/tables/ba_var_table.csv",
  dec = ".",
  sep = ";",
  row.names = FALSE
)




## 4.1 Nullmodel rma.mv----
m_beetle_abund_null <-
  rma.mv(
    yi_abund_beetles,
    vi_abund_beetles,
    random = ~ 1 | lit_id,
    method = "REML",
    data = beetle_abund_outl,
    slab = lit_id
  )
summary(m_beetle_abund_null)

### 4.1.1 Robust variance estimation (Cluster-robust test)----
robust(m_beetle_abund_null,
       cluster = lit_id,
       clubSandwich = TRUE)
# Still gives comparable results. Robust model.

###  4.1.2 Forest plot----
forest(
  m_beetle_abund_null,
  cex = .7,
  main = "Beetle abundance",
  showweights = TRUE,
  order = "obs",
  slab = paste(
    beetle_abund_outl$lit_id,
    beetle_abund_outl$authors,
    beetle_abund_outl$taxon,
    sep = " "
  ),
  pch = 10
)
# Looks fine

### 4.2.2 Model performance----
res_ba <-
  residuals(m_beetle_abund_null, type = "response", cluster = lit_id)
qqnorm(res_ba, cex = 2)
qqline(res_ba)
hist(res_ba, breaks = 20)
# looks fine

## 4.2 Moderator model----
# Full model
m_beetle_abund_full <-
  rma.mv(
    yi_abund_beetles,
    vi_abund_beetles,
    mods = ~
      plant_rich + # 2 NAs in plant_richness
      aes_age + # 21 NAs in age of treatments.
      pr_sum +
      pr_winter +
      tas_mean +
      tas_winter +
      land_div +
      land_div * tas_mean +
      land_div * tas_winter +
      land_div * pr_winter +
      land_div * pr_sum +
      taxon,
    random = ~ 1 | lit_id,
    method = "REML",
    data = beetle_abund_outl,
    slab = lit_id
  )
summary(m_beetle_abund_full)

# final model
m_beetle_abund <-
  rma.mv(
    yi_abund_beetles,
    vi_abund_beetles,
    mods = ~
      # plant_rich + # 2 NAs in plant_richness
      aes_age + # 21 NAs in age of treatments
      pr_sum +
      pr_winter +
      tas_mean +
      tas_winter +
      land_div +
      land_div * tas_mean +
      land_div * tas_winter +
      # land_div * pr_winter +
      # land_div * pr_sum +
      taxon,
    random = ~ 1 | lit_id,
    method = "REML",
    data = beetle_abund_outl,
    slab = lit_id
  )
summary(m_beetle_abund)

# Compare nullmodel with final model
compare_performance(m_beetle_abund, m_beetle_abund_null)
AICc(m_beetle_abund, m_beetle_abund_null)


### 4.2.1 Forest plot----
forest(
  m_beetle_abund,
  cex = .7,
  main = "Beetle abundance",
  showweights = TRUE,
  order = "obs",
  slab = paste(
    beetle_abund_outl$lit_id,
    beetle_abund_outl$authors,
    beetle_abund_outl$taxon,
    sep = " - "
  ),
  pch = 10
)


### 4.2.2 VIF----
# Tested full model without interaction terms
vif(
  update(
    m_beetle_abund_full,
    mods = ~
      plant_rich + # 2 NAs in plant_richness
      aes_age + # 21 NAs in age of treatments.
      pr_sum +
      pr_winter +
      tas_mean +
      tas_winter +
      land_div +
      taxon
  )
)
# Fine. All < 5

# 5 Plot results----
## 5.1 Overall effect----
png_beetle <- readPNG("img/beetle_img.png")
g <- rasterGrob(
  png_beetle,
  interpolate = TRUE,
  x = unit(0.93, "npc"),
  # % from left
  y = unit(0.90, "npc"),
  # % from bottom
  width = unit(0.10, "npc")
)

ba <-
  orchard_plot(
    m_beetle_abund_null,
    mod = "1",
    group = "lit_id",
    xlab = "Hedges' g",
    transfm = "none",
    trunk.size = 2.5,
    # central point
    twig.size = "none",
    # prediction interval
    branch.size = 2,
    # confidence intervals
    alpha = .6,
    flip = TRUE,
    fill = TRUE,
    k = FALSE,
    # number of effect sizes: 95
    g = FALSE
  ) + # number of studies: 30
  labs(title = "") +
  scale_fill_manual(values =  rev(ramp_beetlespider(1))) +
  scale_color_manual(values =  "white") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 20),
    plot.subtitle = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.ticks.y = element_blank()
  ) +
  annotation_custom(g) +
  annotate(
    "text",
    x = 0.6,
    y = -2,
    label = "k = 95 (30)",
    size = 7,
    hjust = 0,
    vjust = 1
  )
ba


## 5.2 grouping = taxon----
# remove taxa with insufficient data
ba_tax_reduced <-
  mod_results(
    m_beetle_abund,
    mod = "taxon",
    group = "lit_id",
    at = list(taxon = c(
      "coccinellidae", "staphylinidae", "carabidae"
    )),
    subset = TRUE
  )

ba_tax <-
  orchard_plot(
    ba_tax_reduced,
    mod = "taxon",
    group = "lit_id",
    xlab = "Hedges' g",
    transfm = "none",
    trunk.size = 2.5,
    # central point
    twig.size = "none",
    # prediction interval
    branch.size = 2,
    # confidence intervals
    alpha = .5,
    flip = TRUE,
    fill = TRUE,
    k = FALSE,
    # number of effect sizes
    g = FALSE
  ) + # number of studies
  labs(title = "") +
  scale_fill_manual(values =  c("#B16364", "#605E78", "#26496B")) +
  scale_color_manual(values =  c("white", "white", "white", "white")) +
  # ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 20),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  ) +
  annotation_custom(g) +
  annotate(
    "text",
    x = c(2.8, 1.8, 0.8),
    y = c(-2.7, -2.7, -2.7),
    label = c("k = 41 (20)", "k = 20 (8)", "k = 8 (5)"),
    size = 7,
    hjust = 0,
    vjust = 0
  )
ba_tax


## 5.3 Interaction: land_div * tas_mean----
# Prepare data for plotting
# separate land_div in quartiles
quant_vals <- quantile(beetle_abund_outl$land_div,
                       probs = c(0.25, 0.50, 0.75),
                       na.rm = TRUE)
names(quant_vals) <- c("low", "medium", "high")
quant_vals

# pred_grid for all combinations of tas_mean and the three fixed quant_vals, rest is fixed to mean or 0
pred_grid_ba <-
  expand.grid(
    tas_mean   = seq(
      min(beetle_abund_outl$tas_mean, na.rm = TRUE),
      max(beetle_abund_outl$tas_mean, na.rm = TRUE),
      length.out = 50
    ),
    land_div   = quant_vals,
    plant_rich = mean(beetle_abund_outl$plant_rich, na.rm = TRUE),
    aes_age    = mean(beetle_abund_outl$aes_age, na.rm = TRUE),
    pr_sum = mean(beetle_abund_outl$pr_sum, na.rm = TRUE),
    pr_winter = mean(beetle_abund_outl$pr_winter, na.rm = TRUE),
    tas_winter = mean(beetle_abund_outl$tas_winter, na.rm = TRUE),
    taxoncoccinellidae = 0,
    taxonothers        = 0,
    taxonstaphylinidae = 0
  )
pred_grid_ba # 150 rows (50 × 3 land_div groups)

# add interaction terms
pred_grid_ba <-
  pred_grid_ba %>%
  mutate(
    `tas_mean:land_div`   = tas_mean * land_div,
    `tas_winter:land_div` = tas_winter * land_div
  )
pred_grid_ba

# order variables the same as in model object
coef_names <- names(coef(m_beetle_abund))
newmods_ba <- as.matrix(pred_grid_ba[, coef_names[coef_names != "intrcpt"]])
newmods_ba


# Now the prediction
# addx = TRUE will produce result table with all variables
pred_ba <- predict(m_beetle_abund, newmods = newmods_ba, addx = TRUE)
pred_ba

# Then bind results to grid
pred_grid_ba$pred <- pred_ba$pred
pred_grid_ba$ci_low  <- pred_ba$ci.lb
pred_grid_ba$ci_up  <- pred_ba$ci.ub

# take over grouping names
pred_grid_ba$group <- rep(names(quant_vals), each = 50)

### 5.3.1 Plot interaction: land_div * tas_mean (finally!)----
ba_landtas_1plot <-
  pred_grid_ba %>%
  mutate(group = factor(group, levels = c("high", "medium", "low"))) %>%
  ggplot(aes(
    x = tas_mean,
    y = pred,
    color = group,
    fill = group
  )) +
  geom_point(
    data = beetle_abund_outl %>%
      mutate(
        group = cut(
          land_div,
          breaks = quantile(land_div, probs = 0:3 /
                              3),
          include.lowest = TRUE,
          labels = c("high", "medium", "low")
        )
      ),
    aes(
      x = tas_mean,
      y = yi_abund_beetles,
      size = 1 / sqrt(vi_abund_beetles),
      color = group
    ),
    position = position_jitter(width = .2, height = 0),
    alpha = .5,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_up),
              alpha = 0.3,
              color = NA) +
  geom_line(linewidth = 1.5) +
  scale_y_continuous(limits = c(-3, 4.2)) +
  scale_x_continuous(limits = c(6, 13)) +
  labs(title = "", x = "Annual temperature [°C]", y = "Hedges´g") +
  scale_fill_manual(
    name = "Landscape complexity [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  scale_color_manual(
    name = "Landscape complexity [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 16),
    title = element_text(size = 20),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = c(6.05, 6.05),
    y = c(-2.5, -3),
    vjust = 0,
    hjust = 0,
    label = c("k = 74 (25)", "LC = [0.74] [0.52] [0.37]"),
    size = 7
  )
ba_landtas_1plot

# > quant_vals
# low    medium      high
# 0.3714592 0.5160853 0.7325880

# Find out k and g from model
bubble_plot(
  m_beetle_abund,
  mod = "land_div",
  group = "lit_id",
  k = TRUE,
  g = TRUE
)
# k = 74 (25)

## 5.4 land_div * tas_winter----
# Three fixed values for land_div in quartiles
quant_vals <- quantile(beetle_abund_outl$land_div,
                       probs = c(0.25, 0.50, 0.75),
                       na.rm = TRUE)
names(quant_vals) <- c("low", "medium", "high")
quant_vals

# pred_grid for all combinations of tas_winter and land_div categories
pred_grid_ba2 <-
  expand.grid(
    tas_winter   = seq(
      min(beetle_abund_outl$tas_winter, na.rm = TRUE),
      max(beetle_abund_outl$tas_winter, na.rm = TRUE),
      length.out = 50
    ),
    land_div   = quant_vals,
    plant_rich = mean(beetle_abund_outl$plant_rich, na.rm = TRUE),
    aes_age    = mean(beetle_abund_outl$aes_age, na.rm = TRUE),
    tas_mean = mean(beetle_abund_outl$tas_mean, na.rm = TRUE),
    pr_sum = mean(beetle_abund_outl$pr_sum, na.rm = TRUE),
    pr_winter = mean(beetle_abund_outl$pr_winter, na.rm = TRUE),
    taxoncoccinellidae = 0,
    taxonothers        = 0,
    taxonstaphylinidae = 0
  )
pred_grid_ba2 # 150 rows (50 × 3 groups).

# add interaction terms
pred_grid_ba2 <-
  pred_grid_ba2 %>%
  mutate(
    `tas_mean:land_div`   = tas_mean * land_div,
    `tas_winter:land_div` = tas_winter * land_div
  )
pred_grid_ba2

# order of columns like in original model
coef_names <- names(coef(m_beetle_abund))
newmods_ba2 <- as.matrix(pred_grid_ba2[, coef_names[coef_names != "intrcpt"]])
newmods_ba2


# predict regressions with beta of original model
pred_ba2 <- predict(m_beetle_abund, newmods = newmods_ba2, addx = TRUE)
pred_ba2

# add results to pred_grid
pred_grid_ba2$pred <- pred_ba2$pred
pred_grid_ba2$ci_low  <- pred_ba2$ci.lb
pred_grid_ba2$ci_up  <- pred_ba2$ci.ub

# add group names for landscape diversity
pred_grid_ba2$group <- rep(names(quant_vals), each = 50)

### 5.4.1 Plot interaction: land_div * tas_winter----
ba_landwint_1plot <-
  pred_grid_ba2 %>%
  mutate(group = factor(group, levels = c("high", "medium", "low"))) %>%
  ggplot(aes(
    x = tas_winter,
    y = pred,
    color = group,
    fill = group
  )) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_up),
              alpha = 0.3,
              color = NA) +
  geom_point(
    data = beetle_abund_outl %>%
      mutate(
        group = cut(
          land_div,
          breaks = quantile(land_div, probs = 0:3 /
                              3),
          include.lowest = TRUE,
          labels = c("high", "medium", "low")
        )
      ),
    aes(
      x = tas_winter,
      y = yi_abund_beetles,
      size = 1 / sqrt(vi_abund_beetles),
      color = group
    ),
    position = position_jitter(width = .2, height = 0),
    alpha = .5,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  geom_line(linewidth = 1.5) +
  labs(title = "", x = "Mean winter temperature [°C]", y = "") +
  scale_y_continuous(limits = c(-3, 4.2)) +
  scale_fill_manual(
    name = "Landscape complexity  [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  scale_color_manual(
    name = "Landscape complexity  [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 16),
    title = element_text(size = 20),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = c(-3.2, -3.2),
    y = c(-2.5, -3),
    hjust = 0,
    vjust = 0,
    label = c("k = 74 (25)", "LC = [0.74] [0.52] [0.37]"),
    size = 7
  )

ba_landwint_1plot

# k and g
bubble_plot(
  m_beetle_abund,
  mod = "tas_winter",
  group = "lit_id",
  g = TRUE,
  k = TRUE
)
# k = 74 (25)

# 5.5 site age----
ba_age <-
  bubble_plot(
    m_beetle_abund,
    group = "lit_id",
    mod = "aes_age",
    xlab = "Site age [a]",
    ci.lwd = .7,
    ci.col = "black",
    pi.lwd = NA,
    alpha = .7,
    k = FALSE,
    g = FALSE
  ) +
  labs(title = "", y = "Hedges´ g") +
  scale_y_continuous(limits = c(-2.5, 4)) +
  scale_fill_manual(values = "#26496B") +
  scale_color_manual(values = "white") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = 1.7,
    y = -2,
    label = "k = 74 (25)",
    size = 6,
    hjust = 0,
    vjust = 1
  )

ba_age$layers[[1]]$aes_params$colour <- "white"

ba_age




#________________________________----
# >> BEETLE RICHNESS----
beetle_rich <-
  beetlespider %>%
  filter(tax_group == "beetles") %>%
  filter(!is.na(t_rich)) %>%
  select(-t_abund, -t_abund_sd, -c_abund, -c_abund_sd) %>%
  droplevels() %>%
  as.data.frame()
beetle_rich
levels(beetle_rich$taxon)

# 0 SMD for richness----
beetle_rich <-
  escalc(
    measure = "SMD",
    n1i = t_sites, # sample size treatment
    n2i = c_sites, # sample size control
    m1i = t_rich, # mean treatment
    m2i = c_rich, # mean control
    sd1i = t_rich_sd,# SD treatment
    sd2i = c_rich_sd, # SD control
    data = beetle_rich,
    slab = lit_id,
    vtype = "LS",
    correct = TRUE,
    var.names = c("yi_rich_beetles", "vi_rich_beetles"),
    add.measure = FALSE
  )
beetle_rich

# remove NAs
beetle_rich <-
  beetle_rich %>%
  filter(!is.na(yi_rich_beetles))


# 1 Nullmodel rma----
m_beetle_rich_null_0 <-
  rma(
    yi_rich_beetles,
    vi_rich_beetles,
    method = "REML",
    data = beetle_rich,
    slab = lit_id
  )
summary(m_beetle_rich_null_0)

# Forest plot
forest(
  m_beetle_rich_null_0,
  cex = .7,
  main = "Beetle richness",
  showweights = TRUE,
  order = "obs",
  pch = 10
)
# No extreme outlier visible

# QQ-Plot
res <- rstandard(m_beetle_rich_null_0)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)
# looks a little out of line in the lower and upper corners

# Baujat plot to identify extreme outliers
baujat(m_beetle_rich_null_0, symbol = "slab")
# Decision: No extreme outlier to exclude, but need to address effects 223.2,  223.1 and others in sensitivity analysis


# 2 Sensitivity analysis----
## Aggregate data for sensitivity tests
beetle_rich_agg <-
  aggregate(beetle_rich, cluster = lit_id, struct = "ID")

# rma Model with aggregated data
m_br_agg <-
  rma(
    yi_rich_beetles,
    vi_rich_beetles,
    # mods = ~ vi_abund_beetles,
    data = beetle_rich_agg,
    slab = lit_id
  )
summary(m_br_agg)

# Forest plot
forest(
  m_br_agg,
  cex = .7,
  main = "Beetle richness",
  showweights = TRUE,
  order = "obs",
  slab = paste(beetle_rich_agg$lit_id, beetle_rich_agg$authors, sep = " - "),
  pch = 10
)
# No extreme outlier visible

# QQ-Plot
res <- rstandard(m_br_agg)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)


## 2.1 Failsafe number----
# Rosenthal's FSN
fsn(yi_rich_beetles,
    vi_rich_beetles,
    method = "Rosenthal",
    data = beetle_rich_agg)
# 652
# > minimum of n x 5 + 10 (here = 120)


## 2.2 Baujat plot----
baujat(m_br_agg, symbol = "slab")
# IDs 108, 115 and 49 need to be checked in leave 1 out

## 2.3 Leave 1 Out----
ba_l1o <- leave1out(m_br_agg)

ba_l1o %>%
  as.data.frame() %>%
  rownames_to_column(var = "study") %>%
  ggplot(aes(x = reorder(study, -estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.2) +
  geom_hline(
    yintercept = coef(m_br_agg),
    linetype = "dashed",
    color = "red"
  ) +
  coord_flip() +
  scale_y_continuous(breaks = seq(from = 0.4, to = 1.4, by = .05)) +
  labs(title = "Leave-one-out Analysis", x = "Removed study", y = "Estimate without study") +
  theme_minimal() +
  theme(
    title = element_text(size = 18),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15)
  )
# Extremes on negative (ID 108, 49) and positive (115) side. Leaving out does not seem to overly change the model estimate (~ 0.1 each)

## 2.4 Cook´s D, Hat, rstudent----
inf_br <- influence(m_br_agg)

inf_br$inf %>%
  as.data.frame() %>%
  rownames_to_column(var = "lit_id") %>%
  select(lit_id, rstudent, cook.d, hat) %>%
  pivot_longer(cols = 2:4,
               names_to = "variable",
               values_to = "value") %>%
  ggplot(aes(x = lit_id, y = value, group = variable)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ variable, scales = "free", ncol = 1) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
# Extremes: Studies 108, 115, 49

## 2.5 Model without extremes----
m_br_agg_update <-
  update(m_br_agg, subset = !(lit_id %in% c("108", "115", "49")))
summary(m_br_agg_update)
summary(m_br_agg)

## 2.6 Compare models----
compare_performance(m_br_agg, m_br_agg_update)

# IDs 49, 108 and 115 have notable influence on the model estimate. However, they compensate each other partly due to their opposite effect directions, which speaks for the robustness of the model itself. Excluding the effects strongly reduces model variance (I², tau²) and brings the estimate down by 0.13 from 0.81 to 0.68. p-value stays very low.
# Decision: Leave the effects in the data, which represent the overall heterogeneity between the studies. This can later help in the rma.mv models to explain variation by moderators. Excluding the effects cannot be justified in these cases.

# 3 Publication bias----
## 3.1 Egger´s test----
funnel(m_br_agg)
regtest(m_br_agg)
# Test for Funnel Plot Asymmetry: z = 1.5225, p = 0.1279
## No significant asymmetry

## 3.2 Trimfill----
trim_br <-
  trimfill(m_br_agg)
summary(trim_br)
# Estimated number of missing studies on the left side: 0 (SE = 2.8268)
# Looks asymmetric, but not detected



# 4 Multivariate models----
### 4.0 Variance structure----
br_variance <-
  rma.mv(
    yi_rich_beetles,
    vi_rich_beetles,
    random = ~ 1 | lit_id / year / taxon / effect_id,
    data = beetle_rich,
    method = "REML"
  )
summary(br_variance)
# Q-test significant
# Intercept still positive (p = 0.0003). Total variance is not fully explained by random structure

# Random level τ² (estim)	Interpretation
# lit_id	0.815: strong variance among studies
# lit_id/year 0.000: no variance between years in studies
# lit_id/year/taxon 0.000: no	variance among taxa within studies and years
# lit_id/year/taxon/effect_id	0.4773:	moderate residual heterogeneity (within taxon–study combinations). non-independence of multiple effects from the same study/taxon/year. keep in the moderator model

# Create and export table for appendix (A.Table 1)
### Beetle richness
br_var_table <-
  data.frame(
    Model = "Beetle species richness",
    Factor_level = unname(br_variance$s.names),
    Variance = "σ ^2",
    Estimate = round(br_variance$sigma2, 3),
    SD = round(sqrt(br_variance$sigma2), 3),
    Levels = br_variance$s.nlevels,
    k = br_variance$k,
    Model_inclusion = c("mandatory", "excluded", "in fixed model term", "included")
  )

br_var_table

write.table(
  br_var_table,
  "output/tables/br_var_table.csv",
  dec = ".",
  sep = ";",
  row.names = FALSE
)


## 4.1 Nullmodel----
m_beetle_rich_null <-
  rma.mv(
    yi_rich_beetles,
    vi_rich_beetles,
    random = ~ 1 | lit_id / effect_id,
    method = "REML",
    data = beetle_rich,
    slab = lit_id
  )
summary(m_beetle_rich_null)

### 4.1.1 Robust variance estimation (Cluster-robust test)----
robust(m_beetle_rich_null,
       cluster = lit_id,
       clubSandwich = TRUE)
# Estimate stays, p-value increases, but still very positive result -> robust model

### 4.1.2 Forest plot----
forest(
  m_beetle_rich_null,
  cex = 1,
  main = "Beetle species richness",
  showweights = TRUE,
  order = "obs",
  # slab = paste(beetle_rich$lit_id, beetle_rich$authors, sep = " "),
  pch = 12
)
# IDs 223.1 and 223.2 very strong (Pfiffner et al.). These are two effects with simulated SD values. Check model outcome without IDs

# Model without Pfiffner et al.
m_beetle_rich_null_update <-
  update(m_beetle_rich_null, subset = !(lit_id %in% c(223)))
summary(m_beetle_rich_null_update)
summary(m_beetle_rich_null)

# Compare with original
compare_performance(m_beetle_rich_null_update, m_beetle_rich_null)
# Excluding does not change model estimate strongly (estimate reduced by 0.08). Decision: Leave 223 in the data set.


### 4.1.3 Model performance----
res_br <-
  residuals(m_beetle_rich_null, type = "response", cluster = lit_id)
qqnorm(res_br, cex = 2)
qqline(res_br)
hist(res_br, breaks = 20)
# looks okay


## 4.2 Moderator model----
# Full model
m_beetle_rich_full <-
  rma.mv(
    yi_rich_beetles,
    vi_rich_beetles,
    mods = ~
      plant_rich +
      aes_age +
      pr_sum +
      pr_winter +
      tas_mean +
      tas_winter +
      land_div +
      land_div * tas_mean +
      land_div * tas_winter +
      land_div * pr_winter +
      land_div * pr_sum +
      taxon,
    random = ~ 1 | lit_id / effect_id,
    #
    method = "REML",
    data = beetle_rich
  )
summary(m_beetle_rich_full)

# Final model
m_beetle_rich <-
  rma.mv(
    yi_rich_beetles,
    vi_rich_beetles,
    mods = ~
      # plant_rich +
      # aes_age +
      pr_sum +
      pr_winter +
      # tas_mean +
      tas_winter +
      land_div +
      # land_div * tas_mean +
      land_div * tas_winter +
      # land_div * pr_winter +
      land_div * pr_sum +
      taxon,
    random = ~ 1 | lit_id / effect_id,
    #
    method = "REML",
    data = beetle_rich
  )
summary(m_beetle_rich)


### 4.2.1 Forest plot----
forest(
  m_beetle_rich,
  cex = 1,
  main = "Beetle species richness",
  showweights = TRUE,
  order = "obs",
  slab = paste(
    beetle_rich$lit_id,
    beetle_rich$authors,
    beetle_rich$taxon,
    sep = " "
  ),
  pch = 12
)

### 4.2.2 VIF----
# With all moderators "on", without interactions
vif(
  update(
    m_beetle_rich_full,
    mods = ~ plant_rich +
      aes_age +
      pr_sum +
      pr_winter +
      tas_mean +
      tas_winter +
      land_div +
      taxon,
  )
)
# All < 5

# 5 Plot results----
## 5.1 Overall effect----
br <-
  orchard_plot(
    m_beetle_rich_null,
    mod = "1",
    group = "lit_id",
    xlab = "Hedges' g",
    transfm = "none",
    trunk.size = 2.5,
    # central point
    twig.size = "none",
    # prediction interval
    branch.size = 2,
    # confidence intervals
    alpha = .6,
    flip = TRUE,
    fill = TRUE,
    k = FALSE,
    # number of effect sizes
    g = FALSE
  ) + # number of studies
  labs(title = "") +
  scale_fill_manual(values =  rev(ramp_beetlespider(1))) +
  scale_color_manual(values =  "white") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 20),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  annotation_custom(g) +
  annotate(
    "text",
    x = 0.55,
    y = -2,
    label = "k = 57 (22)",
    size = 7,
    hjust = 0,
    vjust = 0
  )
br

## 5.2 grouping = taxon----
# remove taxa with insufficient replicates
br_tax_reduced <-
  mod_results(
    m_beetle_rich,
    mod = "taxon",
    group = "lit_id",
    at = list(taxon = c("staphylinidae", "carabidae")),
    subset = TRUE
  )

br_tax <-
  orchard_plot(
    br_tax_reduced,
    mod = "taxon",
    group = "lit_id",
    xlab = "Hedges' g",
    transfm = "none",
    trunk.size = 2.5,# central point
    twig.size = "none",# prediction interval
    branch.size = 2,# confidence intervals
    alpha = .5,
    flip = TRUE,
    fill = TRUE,
    k = FALSE,# number of effect sizes
    g = FALSE # number of studies
  ) + 
  labs(title = "") +
  scale_fill_manual(values =  c("#605E78", "#26496B")) +
  scale_color_manual(values =  c("white", "white", "white", "white")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  annotation_custom(g) +
  annotate(
    "text",
    x = c(1.8, 0.8),
    y = c(-5, -5),
    label = c("k = 39 (20)", "k = 11 (7)"),
    size = 7,
    hjust = 0,
    vjust = 0
  )
br_tax


## 5.3 pr_winter----
br_prw <-
  bubble_plot(
    m_beetle_rich,
    group = "lit_id",
    mod = "pr_winter",
    xlab = "Winter precipitation [mm]",
    ci.lwd = .7,
    ci.col = "black",
    pi.lwd = NA,
    alpha = .7,
    k = FALSE,
    g = FALSE
  ) +
  labs(title = "", y = "Hedges´ g") + # scale_y_continuous(limits = c(-2.5, 5)) +
  scale_fill_manual(values = "#26496B") +
  scale_color_manual(values = "white") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = c(170),
    y = c(-1.3),
    label = "k = 57 (22)",
    size = 6,
    hjust = 0,
    vjust = 1
  )

br_prw$layers[[1]]$aes_params$colour <- "white"
br_prw


## 5.4 land_div * tas_winter----
# land_div quartiles
quant_vals_r <- quantile(beetle_rich$land_div,
                         probs = c(0.25, 0.50, 0.75),
                         na.rm = TRUE)
names(quant_vals_r) <- c("low", "medium", "high")
quant_vals_r

# pred_grid for all combinations of tas_winter and landscape categories
pred_grid_br <-
  expand.grid(
    tas_winter = seq(
      min(beetle_rich$tas_winter, na.rm = TRUE),
      max(beetle_rich$tas_winter, na.rm = TRUE),
      length.out = 50
    ),
    land_div  = quant_vals_r,
    pr_sum = mean(beetle_rich$pr_sum, na.rm = TRUE),
    pr_winter = mean(beetle_rich$pr_winter, na.rm = TRUE),
    taxoncoccinellidae = 0,
    taxonothers        = 0,
    taxonstaphylinidae = 0
  )
pred_grid_br # 150 rows (50 × 3 groups).

# add interaction terms
pred_grid_br <-
  pred_grid_br %>%
  mutate(
    `pr_sum:land_div`   = pr_sum * land_div,
    `tas_winter:land_div` = tas_winter * land_div
  )
pred_grid_br

# bring column order as in nmodel
coef_names <- names(coef(m_beetle_rich))
newmods_br <- as.matrix(pred_grid_br[, coef_names[coef_names != "intrcpt"]])
newmods_br


# Predict
# addx = TRUE brings all variables into the result table = data frame with variables and calculated predictions
pred_br <- predict(m_beetle_rich, newmods = newmods_br, addx = TRUE)
pred_br

# include results in prediction grid
pred_grid_br$pred <- pred_br$pred
pred_grid_br$ci_low  <- pred_br$ci.lb
pred_grid_br$ci_up  <- pred_br$ci.ub

# add landscape groups
pred_grid_br$group <- rep(names(quant_vals_r), each = 50)

### 5.3.1 Plot interaction: land_div * tas_winter----
br_landwint_1plot <-
  pred_grid_br %>%
  mutate(group = factor(group, levels = c("high", "medium", "low"))) %>%
  ggplot(aes(
    x = tas_winter,
    y = pred,
    color = group,
    fill = group
  )) +
  geom_point(
    data = beetle_rich %>%
      mutate(
        group = cut(
          land_div,
          breaks = quantile(land_div, probs = 0:3 /
                              3),
          include.lowest = TRUE,
          labels = c("high", "medium", "low")
        )
      ),
    aes(
      x = tas_winter,
      y = yi_rich_beetles,
      size = 1 / sqrt(vi_rich_beetles),
      color = group
    ),
    position = position_jitter(width = .2, height = 0),
    alpha = .5,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  scale_y_continuous(limits = c(-2.6, 5)) +
  scale_x_continuous(limits = c(1, 6)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_up),
              alpha = 0.3,
              color = NA) +
  geom_line(linewidth = 1.5) +
  labs(title = "", x = "Mean winter temperature [°C]", y = "Hedges´g") +
  scale_fill_manual(
    name = "Landscape complexity [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  scale_color_manual(
    name = "Landscape complexity [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 16),
    title = element_text(size = 20),
    axis.text = element_text(size = 16),
    # axis.text.y = element_blank(),
    # axis.title.y = element_blank(),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = c(1.1, 1.1),
    y = c(-2.1, -2.6),
    label = c("k = 57 (22)", "LC = [0.57] [0.46] [0.38]"),
    hjust = 0,
    vjust = 0,
    size = 7
  )
br_landwint_1plot

# low       medium    high
# 0.3807558 0.4579463 0.5697144

# k and g for interactions
bubble_plot(
  m_beetle_rich,
  mod = "tas_winter",
  group = "lit_id",
  k = TRUE,
  g = TRUE
)
# k = 57 (22)


## 5.5 land_div * pr_sum----
# land_div quartiles
quant_vals_r <- quantile(beetle_rich$land_div,
                         probs = c(0.25, 0.50, 0.75),
                         na.rm = TRUE)
names(quant_vals_r) <- c("low", "medium", "high")
quant_vals_r

# pred_grid for all combinations of tas_winter and landscape categories
pred_grid_br2 <-
  expand.grid(
    pr_sum = seq(
      min(beetle_rich$pr_sum, na.rm = TRUE),
      max(beetle_rich$pr_sum, na.rm = TRUE),
      length.out = 50
    ),
    land_div  = quant_vals_r,
    tas_winter = mean(beetle_rich$tas_winter, na.rm = TRUE),
    pr_winter = mean(beetle_rich$pr_winter, na.rm = TRUE),
    taxoncoccinellidae = 0,
    taxonothers        = 0,
    taxonstaphylinidae = 0
  )
pred_grid_br2 # 150 rows (50 × 3 groups).

# add interaction terms
pred_grid_br2 <-
  pred_grid_br2 %>%
  mutate(
    `pr_sum:land_div`   = pr_sum * land_div,
    `tas_winter:land_div` = tas_winter * land_div
  )
pred_grid_br2

# bring into column order as in nmodel
coef_names <- names(coef(m_beetle_rich))
newmods_br2 <- as.matrix(pred_grid_br2[, coef_names[coef_names != "intrcpt"]])
newmods_br2


# Predict
# addx = TRUE brings all variables into the result = data frame with variables and calculated predictions
pred_br2 <- predict(m_beetle_rich, newmods = newmods_br2, addx = TRUE)
pred_br2

# include in prediction grid
pred_grid_br2$pred <- pred_br2$pred
pred_grid_br2$ci_low  <- pred_br2$ci.lb
pred_grid_br2$ci_up  <- pred_br2$ci.ub

# add landscape groups
pred_grid_br2$group <- rep(names(quant_vals_r), each = 50)

### 5.3.1 Plot interaction: land_div * pr_sum----
br_landpr_1plot <-
  pred_grid_br2 %>%
  mutate(group = factor(group, levels = c("high", "medium", "low"))) %>%
  ggplot(aes(
    x = pr_sum,
    y = pred,
    color = group,
    fill = group
  )) +
  geom_point(
    data = beetle_rich %>%
      mutate(
        group = cut(
          land_div,
          breaks = quantile(land_div, probs = 0:3 /
                              3),
          include.lowest = TRUE,
          labels = c("high", "medium", "low")
        )
      ),
    aes(
      x = pr_sum,
      y = yi_rich_beetles,
      size = 1 / sqrt(vi_rich_beetles),
      color = group
    ),
    position = position_jitter(width = .2, height = 0),
    alpha = .5,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  scale_y_continuous(limits = c(-2.6, 5)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_up),
              alpha = 0.3,
              color = NA) +
  geom_line(linewidth = 1.5) +
  labs(title = "", x = "Annual precipitation [mm]", y = "Hedges´g") +
  scale_fill_manual(
    name = "Landscape diversity",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  scale_color_manual(
    name = "Landscape diversity",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    # c(.85, .85),
    legend.text = element_text(size = 16),
    title = element_text(size = 20),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = c(415, 415),
    y = c(-2.1, -2.6),
    label = c("k = 57 (22)", "LC = [0.57] [0.46] [0.38]"),
    hjust = 0,
    vjust = 0,
    size = 7
  )

br_landpr_1plot

# > quant_vals_r
# low    medium      high
# 0.3807558 0.4579463 0.5697144



#________________________________----
# >> SPIDER ABUNDANCE----
spider_abund <-
  beetlespider %>%
  filter(tax_group == "spiders") %>%
  select(-t_rich, -t_rich_sd, -c_rich, -c_rich_sd) %>%
  droplevels() %>%
  as.data.frame()
spider_abund
levels(spider_abund$taxon)

# 0 SMD for spider abundance-----
spider_abund <-
  escalc(
    measure = "SMD",
    n1i = t_sites, # sample size treatment
    n2i = c_sites, # sample size control
    m1i = t_abund, # mean treatment
    m2i = c_abund, # mean control
    sd1i = t_abund_sd, # SD treatment
    sd2i = c_abund_sd, # SD control
    data = spider_abund,
    slab = lit_id,
    vtype = "LS",
    correct = TRUE,
    var.names = c("yi_abund_spiders", "vi_abund_spiders"),
    add.measure = FALSE
  )
spider_abund

spider_abund <-
  spider_abund %>%
  filter(!is.na(yi_abund_spiders))
# remove NAs

# 1 Nullmodel rma----
m_spider_abund_null_0 <-
  rma(
    yi_abund_spiders,
    vi_abund_spiders,
    method = "REML",
    data = spider_abund,
    slab = lit_id
  )
summary(m_spider_abund_null_0)

# Forest plot
forest(
  m_spider_abund_null_0,
  cex = .7,
  main = "Spider abundance",
  showweights = TRUE,
  order = "obs",
  pch = 10,
  slab = paste(
    spider_abund$lit_id,
    spider_abund$authors,
    spider_abund$taxon
  )
)
# Extreme outlier: 253.1 (Marshall et al. "linyphiidae"), and 139 (Mei et al. "others")

# QQ-Plot
res <- rstandard(m_spider_abund_null_0)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)
# Extreme outliers at both ends of the axes, otherwise fine

# Baujat plot to identify extreme outliers
baujat(m_spider_abund_null_0, symbol = "slab")
#  Effects 253.1 and 139 are extreme here, too.

## 1.1 Exclude outliers----
m_spider_abund_null_0_update <-
  update(m_spider_abund_null_0, subset = !(
    lit_id == "139" | # Mei et al.
      lit_id == "253"
  )) # Marshall et al.
m_spider_abund_null_0_update

# Forest plot
forest(
  m_spider_abund_null_0_update,
  cex = .7,
  main = "Spider abundance",
  showweights = TRUE,
  order = "obs",
  pch = 10,
  slab = paste(
    spider_abund$lit_id,
    spider_abund$authors,
    spider_abund$taxon
  )
)
# Range of yi now less extreme. Looks fine.

# QQ-Plot
res <- rstandard(m_spider_abund_null_0_update)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)
# Improved a lot

# Baujat plot to identify extreme outliers
baujat(m_spider_abund_null_0_update, symbol = "slab")
#  Several extreme effects to be addressed in the Sensitivity analysis


## 1.2 Compare models----
compare_performance(m_spider_abund_null_0_update, m_spider_abund_null_0)
# Variance I² and tau² reduced.

summary(m_spider_abund_null_0_update)
# estimate      se    zval    pval   ci.lb   ci.ub
# 0.5278  0.1007  5.2437  <.0001  0.3305  0.7251  ***

summary(m_spider_abund_null_0)
# estimate      se    zval    pval   ci.lb   ci.ub
# 0.5395  0.1546  3.4902  0.0005  0.2365  0.8425  ***

# Variance measures reduced in updated model, estimates, SE and p stay comparable.
## >> Decision: Exclude outliers to improve model stability and fit. Reasons:
## 1. Extreme low SMD in forest plot < +/-10
## 2. Contortion in QQ plot
## 3. Exclusion improves model roobustness

spider_abund_outl <-
  spider_abund %>%
  filter(!lit_id %in% c("139", "253")) # Mei et al. + Marshal et al.

# 2 Sensitivity analysis----
## Aggregate data for sensitivity tests
spider_abund_agg <-
  aggregate(spider_abund_outl, cluster = lit_id, struct = "ID")
spider_abund_agg

# rma Model with aggregated data
m_sa_agg <-
  rma(yi_abund_spiders,
      vi_abund_spiders,
      data = spider_abund_agg,
      slab = lit_id)
summary(m_sa_agg)

# Forest plot
forest(
  m_sa_agg,
  cex = .7,
  main = "Spider abundance",
  showweights = TRUE,
  order = "obs",
  pch = 10
)

# QQ-Plot
res <- rstandard(m_sa_agg)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)
# looks okay. Somehow extreme in the upper right

## 2.1 Failsafe number----
# #  Lower boundary of the model CI
# fsn(yi_abund_spiders,
#     vi_abund_spiders,
#     alpha = 0.05,
#     target = m_sa_agg$ci.lb,
#     type = "General",
#     method = "REML",
#     data = spider_abund_agg)

# Rosenthal's FSN
fsn(yi_abund_spiders,
    vi_abund_spiders,
    method = "Rosenthal",
    data = spider_abund_agg)
# 339 > minimum of n x 5 + 10 (here = 120)


## 2.2 Baujat plot----
baujat(m_sa_agg, symbol = "slab")
# IDs 66, 169, 62, 301 and 223 need to be checked in leave 1 out

## 2.3 Leave 1 Out----
sa_l1o <- leave1out(m_sa_agg)

sa_l1o %>%
  as.data.frame() %>%
  rownames_to_column(var = "study") %>%
  ggplot(aes(x = reorder(study, -estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.2) +
  geom_hline(
    yintercept = coef(m_sa_agg),
    linetype = "dashed",
    color = "red"
  ) +
  coord_flip() +
  scale_y_continuous(breaks = seq(from = 0.2, to = 1.8, by = .05)) +
  labs(title = "Leave-one-out Analysis", x = "Removed study", y = "Estimate without study") +
  theme_minimal() +
  theme(
    title = element_text(size = 18),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15)
  )
# Extreme on negative (ID 66) side. Leaving out reduces the model estimate  by 0.1

## 2.4 Cook´s D, Hat, rstudent----
inf_sa <- influence(m_sa_agg)

inf_sa$inf %>%
  as.data.frame() %>%
  rownames_to_column(var = "lit_id") %>%
  select(lit_id, rstudent, cook.d, hat) %>%
  pivot_longer(cols = 2:4,
               names_to = "variable",
               values_to = "value") %>%
  ggplot(aes(x = lit_id, y = value, group = variable)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ variable, scales = "free", ncol = 1) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
# Extreme only ID 66 in Cook's D and rstudent

## 2.5 Model without outliers----
m_sa_agg_update <-
  update(m_sa_agg, subset = !(lit_id %in% c("66")))
summary(m_sa_agg_update)
# estimate      se    zval    pval   ci.lb   ci.ub
# 0.3817  0.1315  2.9036  0.0037  0.1241  0.6394  **

summary(m_sa_agg)
# estimate      se    zval    pval   ci.lb   ci.ub
# 0.4909  0.1782  2.7548  0.0059  0.1416  0.8402  **

## 2.6 Compare models----
compare_performance(m_sa_agg, m_sa_agg_update)

forest(
  m_sa_agg,
  cex = .7,
  main = "Spider abundance",
  showweights = TRUE,
  order = "obs",
  pch = 10,
  slab = paste(spider_abund_agg$authors)
)

# ID 66 has notable influence on the model estimate and p-value. Excluding the effect reduces model variance (I², tau²) and brings the estimate down by 0.1 from 0.49 to 0.38, however, removing also reduces the p-value to the next significance level.
## Decision: Leave the effects in the data, which represents the overall heterogeneity between the studies. This can later help in the rma.mv models to explain variation by moderators. ID 66 (Ganser et al) was a study without flaws and high replication, also indicated by the small CI band in the forest plot. There, it is also not extreme in yi. Excluding the effects cannot be justified.

# 3 Publication bias----
## 3.1 Egger´s test----
funnel(m_sa_agg)
regtest(m_sa_agg)
# Test for Funnel Plot Asymmetry: z = -0.0290, p = 0.9769
## No significant asymmetry, although it looks asymmetric in the lower part


## 3.2 Trimfill----
trim_sa <-
  trimfill(m_sa_agg)
summary(trim_sa)
# Estimated number of missing studies on the right side: 0 (SE = 2.9495)

funnel(
  trim_sa,
  bg = c("black", "green"),
  pch = 21,
  main = "Spider abundance\nFunnelplot with Trim-and-Fill"
)
legend(
  "topright",
  legend = c("Original", "Trim-and-Fill"),
  pt.bg = c("black", "green"),
  pch = 21,
  col = "black"
)
# Trim-and-Fill did not detect potentially missing studies, and Egger's test was n.s., too. However, visual inspection of the funnel suggests a possible under-representation of positive effects in the lower right corner. If there was a publication bias, then it woudl be a underestimation of the overall positive effect on spider abundance.



# 4 Multivariate models----
## 4.0 Variance structure----
sa_variance <-
  rma.mv(
    yi_abund_spiders,
    vi_abund_spiders,
    random = ~ 1 | lit_id / year / taxon / effect_id,
    data = spider_abund_outl,
    method = "REML"
  )
summary(sa_variance)
# Q-test significant
# Intercept still positive (p = 0.004). Total variance is not fully explained by random structure

# Random level τ² (estim)	Interpretation
# lit_id	0.3795: moderate variance among studies
# lit_id/year: can be excluded
# lit_id/year/taxon 0.0992  :	minimal	variance among taxa: will be used as fixed factor
# lit_id/year/taxon/effect_id	0.0857: minimal residual heterogeneity (within taxon–study combinations). Non-independence of multiple effects from the same study/taxon/year. Leave out to retain variance in the fixed model part

# create and export summary table
sa_var_table <-
  data.frame(
    Model = "Spider abundance",
    Factor_level = unname(sa_variance$s.names),
    Variance = "σ ^2",
    Estimate = round(sa_variance$sigma2, 3),
    SD = round(sqrt(sa_variance$sigma2), 3),
    Levels = sa_variance$s.nlevels,
    k = sa_variance$k,
    Model_inclusion = c("mandatory", "excluded", "in fixed model term", "excluded")
  )

sa_var_table

write.table(
  sa_var_table,
  "output/tables/sa_var_table.csv",
  dec = ".",
  sep = ";",
  row.names = FALSE
)


## 4.1 Nullmodel----
m_spider_abund_null <-
  rma.mv(
    yi_abund_spiders,
    vi_abund_spiders,
    random = ~ 1 | lit_id,
    method = "REML",
    data = spider_abund_outl,
    slab = lit_id
  )
summary(m_spider_abund_null)

### 4.1.1 Robust variance estimation (Cluster-robust test)----
robust(m_spider_abund_null,
       cluster = lit_id,
       clubSandwich = TRUE)
# Estimate stays, p-value increases, but still positive result

### 4.1.2 Forest plot----
forest(
  m_spider_abund_null,
  cex = 1,
  main = "Spider abundance",
  showweights = TRUE,
  order = "obs",
  slab = paste(
    spider_abund_outl$lit_id,
    spider_abund_outl$authors,
    spider_abund_outl$taxon,
    sep = " - "
  ),
  pch = 12
)
# ID 66 (Ganser et al.) was influential in sensitivity analysis. Check model outcome without these effects

# Model without Ganser et al.
m_spider_abund_null_update <-
  update(m_spider_abund_null, subset = !(lit_id %in% c(66)))
summary(m_spider_abund_null_update)
summary(m_spider_abund_null)

# Compare with original
compare_performance(m_spider_abund_null_update, m_spider_abund_null)
# Excluding changes model estimate by ~ 0.1, like in the leave-one-out analysis (which speaks for its robustness). Let's see if moderators can explain the very positive effect.

# forest plot without Ganser et al.
forest(
  m_spider_abund_null_update,
  cex = 1,
  main = "Spider abundance",
  showweights = TRUE,
  order = "obs",
  slab = paste(spider_abund_outl$lit_id, spider_abund_outl$authors, sep = " "),
  pch = 12
)


### 4.1.3 Model performance----
res_sa <-
  residuals(m_spider_abund_null, type = "response", cluster = lit_id)
qqnorm(res_sa, cex = 2)
qqline(res_sa)
hist(res_sa, breaks = 20)
# looks fine


# 4.2 Moderator model----
# Final model
m_spider_abund_full <-
  rma.mv(
    yi_abund_spiders,
    vi_abund_spiders,
    mods = ~
      plant_rich + # 2 NA
      aes_age + # 28 NA
      pr_sum +
      pr_winter +
      tas_mean +
      tas_winter +
      land_div +
      land_div * tas_mean +
      land_div * tas_winter +
      land_div * pr_winter +
      land_div * pr_sum +
      taxon,
    random = ~ 1 | lit_id,
    method = "REML",
    data = spider_abund_outl
  )
summary(m_spider_abund_full)


# Final model
m_spider_abund <-
  rma.mv(
    yi_abund_spiders,
    vi_abund_spiders,
    mods = ~
      # plant_rich + # 2 NA
      # aes_age + # 28 NA
      # pr_sum +
      pr_winter +
      # tas_mean +
      # tas_winter +
      land_div +
      # land_div * tas_mean +
      # land_div * tas_winter +
      land_div * pr_winter +
      # land_div * pr_sum +
      taxon,
    random = ~ 1 | lit_id,
    method = "REML",
    data = spider_abund_outl
  )
summary(m_spider_abund)
# Test of Moderators (coefficients 2:7):
# QM(df = 6) = 10.4746, p-val = 0.1060
# moderators together are not significant (p = 0.1060). This means adding them doesn’t strongly explain heterogeneity, but the interaction is still meaningful. Report QM in MS!

### 4.2.1 Variance Inflation----
vif(
  update(
    m_spider_abund,
    mods = ~ # with all moderators except interactions
      plant_rich + # 2 NA
      aes_age + # 28 NA
      pr_sum +
      pr_winter +
      tas_mean +
      tas_winter +
      land_div +
      taxon
  )
)
# Factors all < 5

## 4.2.2 Model performance----
res_sa <- residuals(m_spider_abund, type = "response", parallel = "multicore")
qqnorm(res_sa, cex = 2)
qqline(res_sa)
hist(res_sa, breaks = 20)
# fine



# 5 Plot results----
# Load spider image
png_spider <- readPNG("img/spider_img.png")
h <- rasterGrob(
  png_spider,
  interpolate = TRUE,
  x = unit(0.93, "npc"),
  y = unit(0.90, "npc"),
  width = unit(0.10, "npc")
)

## 5.1 Overall effect----
sa <-
  orchard_plot(
    m_spider_abund_null,
    mod = "1",
    group = "lit_id",
    xlab = "Hedges' g",
    transfm = "none",
    trunk.size = 2.5, # central point
    twig.size = "none", # prediction interval
    branch.size = 2, # confidence intervals
    alpha = .6,
    flip = TRUE,
    fill = TRUE,
    k = FALSE, # number of effect sizes
    g = FALSE  # number of studies
  ) + 
  labs(title = "") +
  scale_fill_manual(values =  "#26496B") +
  scale_color_manual(values =  "white") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 20),
    axis.text = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 18),
    axis.ticks.y = element_blank()
  ) +
  annotation_custom(h) +
  annotate(
    "text",
    x = .55,
    y = -2,
    label = "k = 71 (23)",
    size = 7,
    hjust = 0,
    vjust = 0
  )
sa

## 5.2 grouping = taxon----
# remove insufficient data
sa_tax_reduced <-
  mod_results(
    m_spider_abund,
    mod = "taxon",
    group = "lit_id",
    at = list(taxon = c("others", "araneae")),
    subset = TRUE
  )
sa_tax <-
  orchard_plot(
    sa_tax_reduced,
    mod = "taxon",
    group = "lit_id",
    xlab = "Hedges' g",
    transfm = "none",
    trunk.size = 2.5, # central point
    twig.size = "none", # prediction interval
    branch.size = 2, # confidence intervals
    alpha = .5,
    flip = TRUE,
    fill = TRUE,
    k = FALSE, # number of effect sizes
    g = FALSE  # number of studies
  ) + 
  labs(title = "") +
  scale_fill_manual(values =  c("#E9B38D", "#26496B")) +
  scale_color_manual(values =  c("white", "white", "white", "white")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  annotation_custom(h) +
  annotate(
    "text",
    x = c(1.8, 0.8),
    y = c(-2, -2),
    label = c("k = 23 (10)", "k = 34 (15)"),
    size = 7,
    hjust = 0,
    vjust = 0
  )
sa_tax


## 5.3 Interaction land_div * pr_winter----
# land_div quartiles
quant_vals_sa <- quantile(spider_abund_outl$land_div,
                          probs = c(0.25, 0.50, 0.75),
                          na.rm = TRUE)
names(quant_vals_sa) <- c("low", "medium", "high")
quant_vals_sa

# pred_grid for all combinations of tas_winter and landscape categories
pred_grid_sa <-
  expand.grid(
    pr_winter = seq(
      min(spider_abund_outl$pr_winter, na.rm = TRUE),
      max(spider_abund_outl$pr_winter, na.rm = TRUE),
      length.out = 50
    ),
    land_div  = quant_vals_sa,
    taxonlinyphiidae   = 0,
    taxonlycosidae     = 0,
    taxonothers        = 0
  )
pred_grid_sa # 150 rows (50 × 3 groups).

# add interaction terms
pred_grid_sa <-
  pred_grid_sa %>%
  mutate(`pr_winter:land_div` = pr_winter * land_div)
pred_grid_sa

# bring column order as in nmodel
coef_names_sa <- names(coef(m_spider_abund))
newmods_sa <- as.matrix(pred_grid_sa[, coef_names_sa[coef_names_sa != "intrcpt"]])
newmods_sa

# Predict
# addx = TRUE brings all variables into the result = data frame with variables and calculated predictions
pred_sa <- predict(m_spider_abund, newmods = newmods_sa, addx = TRUE)
pred_sa

# include in prediction grid
pred_grid_sa$pred <- pred_sa$pred
pred_grid_sa$ci_low  <- pred_sa$ci.lb
pred_grid_sa$ci_up  <- pred_sa$ci.ub

# add landscape groups
pred_grid_sa$group <- rep(names(quant_vals_sa), each = 50)

### 5.3.1 Plot interaction: land_div * tas_winter----
sa_landwintpr_1plot <-
  pred_grid_sa %>%
  mutate(group = factor(group, levels = c("high", "medium", "low"))) %>%
  ggplot(aes(
    x = pr_winter,
    y = pred,
    color = group,
    fill = group
  )) +
  geom_point(
    data = spider_abund_outl %>%
      mutate(
        group = ntile(land_div, 3),
        group = factor(
          group,
          levels = 1:3,
          labels = c("low", "medium", "high")
        )
      ),
    aes(
      x = pr_winter,
      y = yi_abund_spiders,
      size = 1 / sqrt(vi_abund_spiders),
      color = group
    ),
    position = position_jitter(width = .2, height = 0),
    alpha = .5,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  # scale_y_continuous(limits = c(-2.5, 5)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_up),
              alpha = 0.3,
              color = NA) +
  geom_line(linewidth = 1.5) +
  labs(title = "", x = "Winter precipitation [mm]", y = "Hedges´g") +
  scale_fill_manual(
    name = "Landscape complexity [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  scale_color_manual(
    name = "Landscape complexity [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 16),
    title = element_text(size = 20),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = c(155, 155),
    y = c(-2.05, -2.4),
    label = c("k = 71 (23)", "LC = [0.68] [0.41] [0.00]"),
    hjust = 0,
    vjust = 0,
    size = 7
  )
sa_landwintpr_1plot

# > quant_vals_sa
# low    medium      high
# 0.0000000 0.4130294 0.6799131

# check k and g
bubble_plot(
  m_spider_abund,
  mod = "land_div",
  group = "lit_id",
  k = TRUE,
  g = TRUE
)


#________________________________----
## >> SPIDER RICHNESS----
spider_rich <-
  beetlespider %>%
  filter(tax_group == "spiders") %>%
  filter(!is.na(t_rich)) %>%
  select(-t_abund, -t_abund_sd, -c_abund, -c_abund_sd) %>%
  mutate(c_rich_sd = if_else(c_rich_sd < 0, 0, c_rich_sd)) %>%
  droplevels() %>%
  as.data.frame()
spider_rich
levels(spider_rich$taxon)

### 0 SMD for spider richness-----
spider_rich <-
  escalc(
    measure = "SMD",
    n1i = t_sites, # sample size treatment
    n2i = c_sites, # sample size control
    m1i = t_rich,  # mean treatment
    m2i = c_rich,  # mean control
    sd1i = t_rich_sd, # SD treatment
    sd2i = c_rich_sd, # SD control
    data = spider_rich,
    slab = lit_id,
    vtype = "LS",
    correct = TRUE,
    var.names = c("yi_rich_spiders", "vi_rich_spiders"),
    add.measure = FALSE
  )

spider_rich <-
  spider_rich %>%
  filter(!is.na(yi_rich_spiders))


# 1 Nullmodel rma----
m_spider_rich_null_0 <-
  rma(
    yi_rich_spiders,
    vi_rich_spiders,
    method = "REML",
    data = spider_rich,
    slab = lit_id
  )
summary(m_spider_rich_null_0)

# Forest plot
forest(
  m_spider_rich_null_0,
  cex = .7,
  main = "Spider species richness",
  showweights = TRUE,
  order = "obs",
  pch = 10,
  slab = paste(spider_rich$lit_id, spider_rich$authors, spider_rich$taxon)
)
# Extremes: 169 (Rosas-Ramos et al. "others"), and 233 (Pfiffner et al. "araneae")

# QQ-Plot
res <- rstandard(m_spider_rich_null_0)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)
# Extreme outliers at both ends of the axes, otherwise fine

# Baujat plot to identify extreme outliers
baujat(m_spider_rich_null_0, symbol = "slab")
#  Effects 223, 108.3, 49, 169 are extreme here, too.

## 1.1 Exclude outliers----
### Identify extremes from forest plot (Pfiffner, Rosas-Ramos)
spider_rich %>%
  filter(lit_id %in% c("223", "169"))
# filter(yi_rich_spiders < -3 |
#        yi_rich_spiders > 5)

## 1.2 Compare models----
# Update model
m_spider_rich_null_0_update <-
  update(m_spider_rich_null_0, subset = !(lit_id %in% c("223", "169")))
summary(m_spider_rich_null_0_update)

# Forest plot
forest(
  m_spider_rich_null_0_update,
  cex = .7,
  main = "Spider species richness",
  showweights = TRUE,
  order = "obs",
  pch = 10,
  slab = paste(spider_rich$lit_id, spider_rich$authors, spider_rich$taxon)
)


# QQ-Plot
res <- rstandard(m_spider_rich_null_0_update)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)
# Extreme outliers top right corner, otherwise fine

# Baujat plot to identify extreme outliers
baujat(m_spider_rich_null_0_update, symbol = "slab")
#  Effects 108.3, 108.4, 277.7, 49, are extreme here. Need to be addressed in sensitivity analysis.

compare_performance(m_spider_rich_null_0_update, m_spider_rich_null_0)
summary(m_spider_rich_null_0_update)
# estimate      se    zval    pval    ci.lb   ci.ub
# 0.3535  0.1875  1.8852  0.0594  -0.0140  0.7209  .

summary(m_spider_rich_null_0)
# estimate      se    zval    pval    ci.lb   ci.ub
# 0.3859  0.2094  1.8434  0.0653  -0.0244  0.7963  .

# Extreme outliers detected
# >> Decision: Exclude extremes Rosas-Ramos and Pfiffner et al. Reasons:
# 1 smooths out the forest and QQ plots
# 2 significantly reduces heterogeneity I² and tau² to moderate levels.

spider_rich_outl <-
  spider_rich %>%
  filter(!lit_id %in% c("223", "169"))

# 2 Sensitivity analysis----
## Aggregate data for sensitivity tests
spider_rich_agg <-
  aggregate(spider_rich_outl, cluster = lit_id, struct = "ID")
spider_rich_agg

# rma Model with aggregated data
m_sr_agg <-
  rma(yi_rich_spiders,
      vi_rich_spiders,
      data = spider_rich_agg,
      slab = lit_id)
summary(m_sr_agg)

# Forest plot
forest(
  m_sr_agg,
  cex = .7,
  main = "Spider abundance",
  showweights = TRUE,
  order = "obs",
  pch = 10
)
# ID 49 on the far right

# QQ-Plot
res <- rstandard(m_sr_agg)
qqnorm(res$z, main = "Normal Q-Q Plot of standardized residuals")
qqline(res$z, col = "red", lwd = 2)
# looks okay. Somehow extreme in the upper right

## 2.1 Failsafe number----
# Rosenthal's FSN
fsn(yi_rich_spiders,
    vi_rich_spiders,
    method = "Rosenthal",
    data = spider_rich_agg)
# 58 < minimum of n x 5 + 10 (here = 85)
## Unfortunately low FS N. Needs to be reported in MS!

## 2.2 Baujat plot----
baujat(m_sr_agg, symbol = "slab")
# IDs 108 and 49 need to be checked in leave 1 out

## 2.3 Leave 1 Out----
sr_l1o <- leave1out(m_sr_agg)

sr_l1o %>%
  as.data.frame() %>%
  rownames_to_column(var = "study") %>%
  ggplot(aes(x = reorder(study, -estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.2) +
  geom_hline(
    yintercept = coef(m_sr_agg),
    linetype = "dashed",
    color = "red"
  ) +
  coord_flip() +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = .05)) +
  labs(title = "Leave-one-out Analysis", x = "Removed study", y = "Estimate without study") +
  theme_minimal() +
  theme(
    title = element_text(size = 18),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15)
  )
# Leaving out Id 49 or 108 reduces the model estimate by 0.25 or 0.2 respectively. Both are positive effects.

## 2.4 Cook´s D, Hat, rstudent----
inf_sr <- influence(m_sr_agg)

inf_sr$inf %>%
  as.data.frame() %>%
  rownames_to_column(var = "lit_id") %>%
  select(lit_id, rstudent, cook.d, hat) %>%
  pivot_longer(cols = 2:4,
               names_to = "variable",
               values_to = "value") %>%
  ggplot(aes(x = lit_id, y = value, group = variable)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ variable, scales = "free", ncol = 1) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
# Extreme IDs 49 and 108 in in all indicator variables

## 2.5 Model without outliers----
m_sr_agg_update <-
  update(m_sr_agg, subset = !(lit_id %in% c("49", "108")))
summary(m_sr_agg_update)
# estimate      se    zval    pval    ci.lb   ci.ub
# 0.0760  0.1531  0.4961  0.6198  -0.2241  0.3760

summary(m_sr_agg)
# estimate      se    zval    pval    ci.lb   ci.ub
# 0.4058  0.2323  1.7472  0.0806  -0.0494  0.8611  .

## 2.6 Compare models----
compare_performance(m_sr_agg, m_sr_agg_update)

# IDs 49 and 108 have notable influence on the model estimate and p-value. Excluding the effects reduces model variance (I², tau²) and brings the estimate down by 0.35 from 0.41 to 0.06. The p-value increases strongly. ID 108 (Jeanneret et al.) has extensive sampling over several years (1997-2003) and therefore can be considered a very reliable result. ID 49 (Dittner et al.) had high landscape diversity (0.8) with a mean cover of SMH of ~22 %. This could be the reason for the strong positive effect of flowering strips, considering the intermediate level of sampling effort n_sites = 7), represented by the longer CI band.
## Decision: leave them in the data, as they are reliable sources. Check possible publication bias, and the influence of the effects in the multivariate models.


# 3 Publication bias----
## 3.1 Egger´s test----
funnel(m_sr_agg)
regtest(m_sr_agg)
# Barely not significant test for Funnel Plot Asymmetry: z =  1.9550, p = 0.0506


## 3.2 Trimfill----
trim_sr <-
  trimfill(m_sr_agg)
summary(trim_sr)
# Estimated number of missing studies on the left side: 0 (SE = 2.5174)
# Egger's test indicates marginally significant funnel asymmetry (p = 0.051). Visual inspection of the funnel plot suggests that min. 2–3 effects might be missing in the lower left, indicating a possible under-representation of negative and small effects. However, the trim-and-fill method did not identify any missing studies, which is not uncommon when heterogeneity is high (I² > 80%) or when outliers are present. Therefore, indications of publication bias should be interpreted with caution. Outlier IDs need to be checked for their influence in the rma.mv models. Overall a minor publication bias is indicated and needs to be reported in the MS



# 4 Multivariate models----
## 4.0 Variance structure----
sr_variance <-
  rma.mv(
    yi_rich_spiders,
    vi_rich_spiders,
    random = ~ 1 | lit_id / year / taxon / effect_id,
    data = spider_rich_outl,
    method = "REML"
  )
summary(sr_variance)
# Q-test significant
# Intercept still positive (p = 0.0931). Total variance is not fully explained by random structure

# Random level τ² (estim)	Interpretation
# lit_id	0.4533: moderate variance among studies
# lit_id/year: can be excluded
# lit_id/taxon:	can be excluded
# lit_id/year/taxon/effect_id	0.4925: high residual heterogeneity (within taxon–study combinations). Non-independence of multiple effects from the same study/taxon/year.

# “Q” < 0.0001 = observed effect sizes vary more than expected by sampling error alone (by ignoring the random structure)



## 4.1 Nullmodel----
m_spider_rich_null <-
  rma.mv(
    yi_rich_spiders,
    vi_rich_spiders,
    random = ~ 1 |
      lit_id / effect_id,
    method = "REML",
    data = spider_rich_outl,
    slab = lit_id
  )
summary(m_spider_rich_null)

### 4.1.1 Robust variance estimation (Cluster-robust test)----
robust(m_spider_rich_null,
       cluster = lit_id,
       clubSandwich = TRUE)
# No dramatic changes

### 4.1.2 Forest plot----
forest(
  m_spider_rich_null,
  cex = 1,
  main = "Spider abundance",
  showweights = TRUE,
  order = "obs",
  slab = paste(spider_rich_outl$lit_id, spider_rich_outl$authors, sep = " "),
  pch = 12
)
# ID 49 and 108 were considered very influential. Check model outcome without effects

# Model without influential studies
m_spider_rich_null_update <-
  update(m_spider_rich_null, subset = !(yi_rich_spiders > 3))
summary(m_spider_rich_null_update)
summary(m_spider_rich_null)
# Estimate drops from 0.41 to 0.21, and p increases from 0.0989 to 0.2556
# Let's see the residuals and if moderators can explain the effects.


### 4.1.3 Model performance----
res_sr <-
  residuals(m_spider_rich_null, type = "response", cluster = lit_id)
qqnorm(res_sr, cex = 2)
qqline(res_sr)
hist(res_sr, breaks = 20)
# Extremes in the upper right corner


## 4.2 Moderator model----
# Full model
m_spider_rich_full <-
  rma.mv(
    yi_rich_spiders,
    vi_rich_spiders,
    mods = ~
      plant_rich + # 1 NA
      aes_age + # 11 NA
      pr_sum +
      pr_winter +
      tas_mean +
      tas_winter +
      land_div +
      land_div * tas_mean +
      land_div * tas_winter +
      land_div * pr_winter +
      land_div * pr_sum  +
      taxon,
    random = ~ 1 | lit_id / effect_id,
    method = "REML",
    data = spider_rich_outl
  )
summary(m_spider_rich_full)

# Final model
m_spider_rich <-
  rma.mv(
    yi_rich_spiders,
    vi_rich_spiders,
    mods = ~
      plant_rich + # 1 NA
      aes_age + # 11 NA
      # pr_sum +
      # pr_winter +
      # tas_mean +
      # tas_winter +
      # land_div +
      # land_div * tas_mean +
      land_div * tas_winter +
      # land_div * pr_winter +
      # land_div * pr_sum  +
      taxon,
    random = ~ 1 | lit_id / effect_id,
    method = "REML",
    data = spider_rich_outl
  )
summary(m_spider_rich)

AICc(m_spider_rich, m_spider_rich_null)

#### 4.2.1.1 Variance Inflation----
vif(
  update(
    m_spider_rich_full,
    # without interactions
    mods = ~
      plant_rich + # 1 NA
      aes_age + # 11 NA
      pr_sum +
      pr_winter +
      tas_mean +
      tas_winter +
      land_div +
      taxon
  )
)
# all < 5

# Forest plot
forest(
  m_spider_rich,
  cex = .7,
  main = "Spider species richness",
  showweights = TRUE,
  order = "obs",
  pch = 10,
  slab = paste(
    spider_rich_outl$lit_id,
    spider_rich_outl$authors,
    spider_rich_outl$taxon,
    sep = " - "
  )
)


## 4.3 Check influential studies from sensitivity analysis----
# ID 49 Ditner et al., ID 108 Jeanneret et al.
### 4.3.1 Winter temperature----
bubble_plot(
  m_spider_rich,
  group = "lit_id",
  mod = "tas_winter",
  xlab = "Winter temperature [°C]",
  alpha = .8
) +
  geom_text(
    data = spider_rich_outl,
    aes(x = tas_winter, y = yi_rich_spiders, label = authors),
    size = 3,
    vjust = -0.7,
    check_overlap = TRUE
  ) +
  theme(legend.position = "none")
# Jeanneret et al. with high winter temperature (> 4.5 °C)

### 4.3.2 Age----
bubble_plot(
  m_spider_rich,
  group = "lit_id",
  mod = "aes_age",
  xlab = "Years after sowing",
  alpha = .8
) +
  geom_text(
    data = spider_rich_outl,
    aes(x = aes_age, y = yi_rich_spiders, label = authors),
    size = 3,
    vjust = -0.7,
    check_overlap = TRUE
  ) +
  theme(legend.position = "none")

### 4.3.2 Landscape----
bubble_plot(
  m_spider_rich,
  group = "lit_id",
  mod = "land_div",
  xlab = "Landscape diversity",
  alpha = .8
) +
  geom_text(
    data = spider_rich_outl,
    aes(x = land_div, y = yi_rich_spiders, label = authors),
    size = 3,
    vjust = -0.7,
    check_overlap = TRUE
  ) +
  theme(legend.position = "none")

### 4.3.2 Plant rich----
bubble_plot(
  m_spider_rich,
  group = "lit_id",
  mod = "plant_rich",
  xlab = "Seed mix richnes",
  alpha = .8
) +
  geom_text(
    data = spider_rich_outl,
    aes(x = plant_rich, y = yi_rich_spiders, label = authors),
    size = 3,
    vjust = -0.7,
    check_overlap = TRUE
  ) +
  theme(legend.position = "none")

# Check influential IDs from the plots
# Frank et al. and Jeanneret et al. (272, 108) were long-time studies over several years around the year 2000 with site-age ranges 1-7 years. The influential effect of ID 108 in the sensitivity analysis can be explained by combination of high age of the sites and winter temperatures.
# Ditner et al. (49) has low site age and moderate winter temperature, but very high landscape diversity and low seed mix richness. This combination explains the high effect size
# The studies can be included in the analysis




# 5 Plot results----
## 5.1 Overall effect----
sr <-
  orchard_plot(
    m_spider_rich_null,
    mod = "1",
    group = "lit_id",
    xlab = "Hedges' g",
    transfm = "none",
    trunk.size = 2.5, # central point
    twig.size = "none", # prediction interval
    branch.size = 2, # confidence intervals
    alpha = .6,
    flip = TRUE,
    fill = TRUE,
    k = FALSE, # number of effect sizes
    g = FALSE  # number of studies
  ) + 
  labs(title = "") +
  scale_fill_manual(values =  ramp_beetlespider(1)) +
  scale_color_manual(values = "white") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    axis.ticks.y = element_blank()
  ) +
  annotation_custom(h) +
  annotate(
    "text",
    x = .55,
    y = -2,
    label = "k = 36 (16)",
    size = 7,
    hjust = 0,
    vjust = 0
  )
sr

## 5.2 grouping = taxon----
sr_tax <-
  orchard_plot(
    m_spider_rich,
    mod = "taxon",
    group = "lit_id",
    xlab = "Hedges' g",
    transfm = "none",
    trunk.size = 2.5, # central point
    twig.size = "none", # prediction interval
    branch.size = 2, # confidence intervals
    tree.order = c("Others", "Araneae"),
    alpha = .5,
    flip = TRUE,
    fill = TRUE,
    k = FALSE, # number of effect sizes
    g = FALSE  # number of studies
  ) + 
  labs(title = "") +
  scale_fill_manual(values =  rev(ramp_beetlespider(2))) +
  scale_color_manual(values =  c("white", "white", "white", "white")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  annotation_custom(h) +
  annotate(
    "text",
    x = c(1.9, 0.9),
    y = c(-2, -2),
    label = c("k = 11 (4)", "k = 13 (8)"),
    size = 7,
    hjust = 0,
    vjust = 0
  )
sr_tax

## 5.3 Continuous = plant richness----
sr_plr <-
  bubble_plot(
    m_spider_rich,
    group = "lit_id",
    mod = "plant_rich",
    xlab = "Seed mixture species richness",
    ci.lwd = .7,
    ci.col = "black",
    pi.lwd = NA,
    alpha = 1,
    k = FALSE,
    g = FALSE
  ) +
  labs(title = "", y = "Hedges' g") +
  scale_y_continuous(limits = c(-2.5, 4.3), breaks = c(-2, 0, 2, 4)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = 0.8,
    y = -2,
    label = "k = 24 (12)",
    size = 7,
    hjust = 0,
    vjust = 1
  )

sr_plr$layers[[1]]$aes_params$colour <- "white"
sr_plr$layers[[1]]$aes_params$fill <- "#EB8867"
sr_plr


## 5.4 Continuous = Site age----
sr_age <-
  bubble_plot(
    m_spider_rich,
    group = "lit_id",
    mod = "aes_age",
    xlab = "Site age [a]",
    ci.lwd = .7,
    ci.col = "black",
    pi.lwd = NA,
    alpha = 1,
    k = FALSE,
    g = FALSE
  ) +
  labs(title = "", y = "Hedges´g") +
  scale_y_continuous(limits = c(-2.5, 4.3), breaks = c(-2, 0, 2, 4)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = .8,
    y = -2,
    label = "k = 24 (12)",
    size = 7,
    hjust = 0,
    vjust = 1
  )

sr_age$layers[[1]]$aes_params$colour <- "white"
sr_age$layers[[1]]$aes_params$fill <- "#EB8867"
sr_age


## 5.5 Interaction land_div * tas_winter----
# land_div quartiles
quant_vals_sr <- quantile(spider_rich_outl$land_div,
                          probs = c(0.25, 0.50, 0.75),
                          na.rm = TRUE)
names(quant_vals_sr) <- c("low", "medium", "high")
quant_vals_sr

# pred_grid for all combinations of tas_winter and landscape categories
pred_grid_sr <-
  expand.grid(
    tas_winter = seq(
      min(spider_rich_outl$tas_winter, na.rm = TRUE),
      max(spider_rich_outl$tas_winter, na.rm = TRUE),
      length.out = 50
    ),
    land_div  = quant_vals_sr,
    plant_rich = mean(spider_rich_outl$plant_rich, na.rm = TRUE),
    aes_age = mean(spider_rich_outl$aes_age, na.rm = TRUE),
    taxonothers = 0
  )
pred_grid_sr # 150 rows (50 × 3 groups).

# add interaction terms
pred_grid_sr <-
  pred_grid_sr %>%
  mutate(`land_div:tas_winter` = tas_winter * land_div)
pred_grid_sr

# bring column order as in model
coef_names_sr <- names(coef(m_spider_rich))
newmods_sr <- as.matrix(pred_grid_sr[, coef_names_sr[coef_names_sr != "intrcpt"]])
newmods_sr


# Predict
# addx = TRUE brings all variables into the result = data frame with variables and calculated predictions
pred_sr <- predict(m_spider_rich, newmods = newmods_sr, addx = TRUE)
pred_sr

# include in prediction grid
pred_grid_sr$pred <- pred_sr$pred
pred_grid_sr$ci_low  <- pred_sr$ci.lb
pred_grid_sr$ci_up  <- pred_sr$ci.ub

# add landscape groups
pred_grid_sr$group <- rep(names(quant_vals_sr), each = 50)

### 5.5.1 Plot interaction: land_div * tas_winter----
sr_landwinttas_1plot <-
  pred_grid_sr %>%
  mutate(group = factor(group, levels = c("high", "medium", "low"))) %>%
  ggplot(aes(
    x = tas_winter,
    y = pred,
    color = group,
    fill = group
  )) +
  geom_point(
    data = spider_rich_outl %>%
      mutate(
        group = cut(
          land_div,
          breaks = quantile(land_div, probs = 0:3 /
                              3),
          include.lowest = TRUE,
          labels = c("high", "medium", "low")
        )
      ),
    aes(
      x = tas_winter,
      y = yi_rich_spiders,
      size = 1 / sqrt(vi_rich_spiders),
      color = group
    ),
    alpha = .5,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2, 8)) +
  scale_x_continuous(limits = c(1, 6)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_up),
              alpha = 0.3,
              color = NA) +
  geom_line(linewidth = 1.5) +
  labs(title = "", x = "Mean winter temperature [°C]", y = "") +
  scale_fill_manual(
    name = "Landscape complexity [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  scale_color_manual(
    name = "Landscape complexity [LC]",
    values = c("#26496B", "#8A616D", "#EAA581"),
    labels = c("High", "Medium", "Low")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 16),
    title = element_text(size = 20),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = c(1.05, 1.05),
    y = c(-3.9, -4.4),
    label = c("k = 24 (12)", "LC = [0.52] [0.42] [0.23]"),
    hjust = 0,
    vjust = 0,
    size = 7
  )
sr_landwinttas_1plot

# > quant_vals_sr
# low    medium      high
# 0.2290441 0.4219131 0.5163950

# k and g for interactions
bubble_plot(
  m_spider_rich,
  mod = "tas_winter",
  group = "lit_id",
  k = TRUE,
  g = TRUE
)


#________________________________----
# Final plots----

## 01 Overall effects----
final_plot_overall <-
  (((ba + labs(title = "Abundance")) + (br + labs(title = "Species richness")))) /
  
  (((sa + labs(title = "Abundance")) + (sr + labs(title = "Species richness"))) &
     scale_fill_manual(values = "#EB8867")) &
  
  theme(
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 22),
    axis.text  = element_text(size = 20),
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 22),
    plot.tag = element_text(size = 26, face = "bold")
  ) &
  plot_annotation(title = "Overall effects",
                  tag_levels = "a",
                  tag_suffix = ")")
final_plot_overall

# export
png(
  "output/figures/fig2_overall_effects.png",
  height = 55 * 2 / 4,
  width = 45 * 3 / 4,
  units = "cm",
  res = 320
)
final_plot_overall
dev.off()



## 02 Taxon effects----
final_plot_taxa <-
  (((ba_tax + labs(title = "Abundance")) + (br_tax + labs(title = "Species richness"))) /
     
     ((sa_tax + labs(title = "Abundance")) + (sr_tax + labs(title = "Species richness")))) &
  
  theme(
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 22),
    axis.text  = element_text(size = 20),
    legend.title = element_text(size = 24),
    legend.text  = element_text(size = 22),
    plot.tag = element_text(size = 26, face = "bold")
  ) &
  plot_annotation(title = "Taxon-specific effects",
                  tag_levels = "a",
                  tag_suffix = ")")
final_plot_taxa

# export
png(
  "output/figures/fig3_taxa.png",
  height = 55 * 2 / 4 + 10,
  width = 45,
  units = "cm",
  res = 320
)
final_plot_taxa
dev.off()



### 3.2 Main effects----
main_effects <-
  ((ba_age & labs(y = "Beetle abundance\nHedges' g", x = "Site age [a]")) +
     
     (br_prw & labs(y = "Beetle species richness\nHedges' g")) &
     annotation_custom(g)) /
  
  
  ((sr_age & labs(y = "Spider species richness\nHedges'g")) +
     (sr_plr & labs(y = "\n")) &
     annotation_custom(h)) &
  
  theme(
    plot.title = element_text(size = 34),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 28),
    legend.title = element_text(size = 32),
    legend.text  = element_text(size = 30),
    plot.tag = element_text(size = 34, face = "bold")
  ) &
  plot_annotation(title = "Main effects",
                  tag_levels = "a",
                  tag_suffix = ")")
main_effects

png(
  "output/figures/fig4_main_effects.png",
  height = 55 * 2 / 3,
  width = 45,
  units = "cm",
  res = 320
)
main_effects
dev.off()



## 3.3 Interactions----
interactions <-
  
  # Beetle abundance
  (((
    ba_landtas_1plot & labs(y = "Beetle abundance\nHedges'g")
  ) +
    
    (
      ba_landwint_1plot & labs(x = "Mean winter temperature [°C]") &
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    )) /
    
    
    # Beetle richness
    ((
      br_landpr_1plot & labs(y = "Beetle species richness\nHedges'g")
    ) +
      
      (
        br_landwint_1plot & labs(y = "", x = "Mean winter temperature [°C]") &
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      )) &
    annotation_custom(g)
  ) /
  
  
  # Spider abundance & richness
  (((
    sa_landwintpr_1plot & labs(y = "Spider abundance\nHedges'g") &
      theme(legend.position = "none")
  ) +
    
    (
      sr_landwinttas_1plot & labs(y = "Spider species richness\nHedges'g", x = "Mean winter temperature [°C]") &
        theme(legend.position = "none")
    ) &
    annotation_custom(h)
  )) &
  
  theme(
    plot.title = element_text(size = 34),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 28),
    legend.title = element_text(size = 32),
    legend.text  = element_text(size = 30),
    plot.tag = element_text(size = 34, face = "bold")
  ) &
  plot_annotation(title = "Interactions",
                  tag_levels = "a",
                  tag_suffix = ")")
interactions


png(
  "output/figures/fig5_interactions.png",
  height = 55,
  width = 45,
  units = "cm",
  res = 320
)
interactions
dev.off()

# Note: In Fig. 6, the legend showing the landscape complexity levels is missing. There was no way to plot it into the original png file (or better: I did not find one). Therefore, the legend was manually integrated afterwards with other software.


## A. Table 2----
# Final Model results

### Nullmodel beetle abundance----
## A. Table 2----
# Final Model results

### Nullmodel beetle abundance----
summary(m_beetle_abund_null)
ba_summary_null <-
  data.frame(
    Model = "Beetle abundance - Overall effect",
    Factor_level = unname(m_beetle_abund_null$s.names),
    σ_2 = m_beetle_abund_null$sigma2,
    Q = round(m_beetle_abund_null$QE, 3),# Test for heterogeneity
    Q_df = m_beetle_abund_null$QEdf,
    Q_pval = "<0.0001",
    k = m_beetle_abund_null$k,# Number of studies
    AICc = round(AICc(m_beetle_abund_null), 1), # AICc value, if metafor::AICc available
    Factor = colnames(m_beetle_abund_null$X),
    Estimate = round(m_beetle_abund_null$b, 3),# Coefficients
    SE = m_beetle_abund_null$se,
    z = round(m_beetle_abund_null$zval, 3),
    ci.lb = round(m_beetle_abund_null$ci.lb, 3),
    ci.ub = round(m_beetle_abund_null$ci.ub, 3),
    pval = round(m_beetle_abund_null$pval, 3),
    row.names = NULL
  )
ba_summary_null

### Minimal Model beetle abundance----
summary(m_beetle_abund)
ba_summary <-
  data.frame(
    Model = "Beetle abundance - Minimal model",
    Factor_level = unname(m_beetle_abund$s.names),
    σ_2 = m_beetle_abund$sigma2,
    Q = round(m_beetle_abund$QE, 3),# Test for heterogeneity
    Q_df = m_beetle_abund$QEdf,
    Q_pval = "<0.0001",
    k = m_beetle_abund$k,# Number of studies
    AICc = round(AICc(m_beetle_abund), 1),# AICc value, if metafor::AICc available
    Factor = colnames(m_beetle_abund$X),
    Estimate = round(m_beetle_abund$b, 3),# Coefficients
    SE = m_beetle_abund$se,
    z = round(m_beetle_abund$zval, 3),
    ci.lb = round(m_beetle_abund$ci.lb, 3),
    ci.ub = round(m_beetle_abund$ci.ub, 3),
    pval = round(m_beetle_abund$pval, 3),
    row.names = NULL
  )
ba_summary


### Nullmodel beetle species richness----
summary(m_beetle_rich_null)
br_summary_null <-
  data.frame(
    Model = "Beetle species richness - Overall effect",
    Factor_level = unname(m_beetle_rich_null$s.names),
    σ_2 = m_beetle_rich_null$sigma2,
    Q = round(m_beetle_rich_null$QE, 3), # Test for heterogeneity
    Q_df = m_beetle_rich_null$QEdf,
    Q_pval = "<0.0001",
    k = m_beetle_rich_null$k,# Number of studies
    AICc = round(AICc(m_beetle_rich_null), 1),# AICc value, if metafor::AICc available
    Factor = colnames(m_beetle_rich_null$X),
    Estimate = round(m_beetle_rich_null$b, 3),# Coefficients
    SE = m_beetle_rich_null$se,
    z = round(m_beetle_rich_null$zval, 3),
    ci.lb = round(m_beetle_rich_null$ci.lb, 3),
    ci.ub = round(m_beetle_rich_null$ci.ub, 3),
    pval = m_beetle_rich_null$pval,
    row.names = NULL
  )
br_summary_null

### Minimal Model beetle species richness----
summary(m_beetle_rich)
br_summary <-
  data.frame(
    Model = "Beetle species richness - Minimal model",
    Factor_level = unname(m_beetle_rich$s.names),
    σ_2 = m_beetle_rich$sigma2,
    Q = round(m_beetle_rich$QE, 3),# Test for heterogeneity
    Q_df = m_beetle_rich$QEdf,
    Q_pval = "<0.0001",
    k = m_beetle_rich$k,# Number of studies
    AICc = round(AICc(m_beetle_rich), 1),# AICc value, if metafor::AICc available
    Factor = colnames(m_beetle_rich$X),
    Estimate = round(m_beetle_rich$b, 3), # Coefficients
    SE = m_beetle_rich$se,
    z = round(m_beetle_rich$zval, 3),
    ci.lb = round(m_beetle_rich$ci.lb, 3),
    ci.ub = round(m_beetle_rich$ci.ub, 3),
    pval = round(m_beetle_rich$pval, 3),
    row.names = NULL
  )
br_summary


### Nullmodel spider abundance----
summary(m_spider_abund_null)
sa_summary_null <-
  data.frame(
    Model = "Spider abundance - Overall effect",
    Factor_level = unname(m_spider_abund_null$s.names),
    σ_2 = m_spider_abund_null$sigma2,
    Q = round(m_spider_abund_null$QE, 3),# Test for heterogeneity
    Q_df = m_spider_abund_null$QEdf,
    Q_pval = "<0.0001",
    k = m_spider_abund_null$k,# Number of studies
    AICc = round(AICc(m_spider_abund_null), 1),# AICc value, if metafor::AICc available
    Factor = colnames(m_spider_abund_null$X),
    Estimate = round(m_spider_abund_null$b, 3),# Coefficients
    SE = m_spider_abund_null$se,
    z = round(m_spider_abund_null$zval, 3),
    ci.lb = round(m_spider_abund_null$ci.lb, 3),
    ci.ub = round(m_spider_abund_null$ci.ub, 3),
    pval = round(m_spider_abund_null$pval, 3),
    row.names = NULL
  )
sa_summary_null

### Minimal Model spider abundance----
summary(m_spider_abund)
sa_summary <-
  data.frame(
    Model = "spider abundance - Minimal model",
    Factor_level = unname(m_spider_abund$s.names),
    σ_2 = m_spider_abund$sigma2,
    Q = round(m_spider_abund$QE, 3), # Test for heterogeneity
    Q_df = m_spider_abund$QEdf,
    Q_pval = "<0.0001",
    k = m_spider_abund$k, # Number of studies
    AICc = round(AICc(m_spider_abund), 1),
    Factor = colnames(m_spider_abund$X),
    Estimate = round(m_spider_abund$b, 3),# Coefficients
    SE = m_spider_abund$se,
    z = round(m_spider_abund$zval, 3),
    ci.lb = round(m_spider_abund$ci.lb, 3),
    ci.ub = round(m_spider_abund$ci.ub, 3),
    pval = round(m_spider_abund$pval, 3),
    row.names = NULL
  )
sa_summary



### Nullmodel spider species richness----
sr_summary_null <-
  data.frame(
    Model = "Spider species richness - Overall effect",
    Factor_level = unname(m_spider_rich_null$s.names),
    σ_2 = m_spider_rich_null$sigma2,
    Q = round(m_spider_rich_null$QE, 3), # Test for heterogeneity
    Q_df = m_spider_rich_null$QEdf,
    Q_pval = "<0.0001",
    k = m_spider_rich_null$k,# Number of studies
    AICc = round(AICc(m_spider_rich_null), 1),# AICc Wert, falls metafor::AICc verfügbar
    Factor = colnames(m_spider_rich_null$X),
    Estimate = round(m_spider_rich_null$b, 3),# Koeffizienten
    SE = m_spider_rich_null$se,
    z = round(m_spider_rich_null$zval, 3),
    ci.lb = round(m_spider_rich_null$ci.lb, 3),
    ci.ub = round(m_spider_rich_null$ci.ub, 3),
    pval = round(m_spider_rich_null$pval, 3),
    row.names = NULL
  )
sr_summary_null

summary(m_spider_rich_null)


### Minimal Model spider species richness----
summary(m_spider_rich)
sr_summary <-
  data.frame(
    Model = "Spider species richness - Minimal model",
    Factor_level = unname(m_spider_rich$s.names),
    σ_2 = m_spider_rich$sigma2,
    Q = round(m_spider_rich$QE, 3),# Test for heterogeneity
    Q_df = m_spider_rich$QEdf,
    Q_pval = "0.0005",
    k = m_spider_rich$k,# Number of studies
    AICc = round(AICc(m_spider_rich), 0),
    Factor = c(colnames(m_spider_rich$X), "space_to_fill"),# Spacer to be deleted afterwards
    Estimate = c(round(m_spider_rich$b, 3), NA),# Coefficient
    SE = c(round(m_spider_rich$se, 3), NA),
    z = c(round(m_spider_rich$zval, 3), NA),
    ci.lb = c(round(m_spider_rich$ci.lb, 3), NA),
    ci.ub = c(round(m_spider_rich$ci.ub, 3), NA),
    pval = c(round(m_spider_rich$pval, 3), NA),
    row.names = NULL
  )
sr_summary

### Export model summary table
ba_summary_null %>%
  bind_rows(
    ba_summary,
    br_summary_null,
    br_summary,
    sa_summary_null,
    sa_summary,
    sr_summary_null,
    sr_summary
  ) %>%
  write.table(
    "output/tables/model_results_A_Table2.csv",
    dec = ".",
    sep = ";",
    row.names = FALSE
  )