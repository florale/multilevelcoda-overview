# ## Notes for running on Ubuntu = Within Ubuntu 24.04 terminal:
# 
# ## get some dependencies for R and related packages
# sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev cmake
# sudo apt install libfontconfig1-dev libcairo2-dev libharfbuzz-dev libfribidi-dev
# 
# ## get some fonts for graphs & make sure to accept the EULA or you must reinstall it
# sudo apt install ttf-mscorefonts-installer
# ## update the font cache
# sudo fc-cache -fv
# ## make sure that Times New Roman font has been installed successfully
# fc-match TimesNewRoman
# 
# ## get the latest version of R for Noble Ubuntu need to setup repo access
# wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/"
# 
# ## update packages then install R base and dev
# sudo apt update && sudo apt upgrade -y
# sudo apt-get install r-base r-base-dev -y
# 
# ## now open R and within R, get relevant packages
# 
# install.packages(c("data.table", "future", "foreach", "ggplot2", "scales", "ggpubr", "extrafont"))
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# 
# ## installing cmdstanr is not enough as we need cmdstan installed on the OS as well
# library(cmdstanr)
# install_cmdstan(cores=8)
# 
# ## now we can get multilevelcoda and its dependencies installed
# install.packages("multilevelcoda", dependencies = TRUE)
# 
# ## fonts are installed on the OS but we need R to be able to access and know about them
# library(extrafont)
# font_import()
# 
# ## now close and restart R

library(data.table)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(scales)
library(extrafont)
library(brms)

library(multilevelcoda)

# a working install of cmdstanr is required for the replication script to work
# you can install by running the following code
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# next the actual cmdstan CLI must be installed, which can be done by:
# library(cmdstanr)
# install_cmdstan()

# helper functions -----------------------------
theme_multilevelcoda <- function() {
  hrbrthemes::theme_ipsum() +
    theme(
      axis.ticks         = element_blank(),
      panel.background   = element_blank(),
      panel.border       = element_blank(),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.background    = element_rect(fill = "transparent", colour = NA),
      strip.background   = element_rect(fill = "transparent", colour = NA),
      strip.text         = element_text(size = 12, hjust = .5),
      strip.placement    = "outside",
      axis.title.x       = element_blank(),
      axis.title.y       = element_text(size = 13, hjust = .5),
      axis.title.y.right = element_text(size = 13, hjust = .5, angle = 270),
      # plot.title         = element_text(size = 12, face = "bold"),
      # plot.title.position= "plot",
      plot.margin        = margin(.5, .5, .5, .5, "cm"),
      legend.position    = "none"
    )
}
fill_multilevelcoda   <- c(`LPA` = "#BBA9A7", `MVPA` = "#FAF7F3", `SB` = "#CFDAE2", `TST` = "#647F9A", `WAKE` = "#8399AE")
colour_multilevelcoda <- c(`LPA` = "#978787", `MVPA` = "#BBA9A7", `SB` = "#A1B2C2", `TST` = "#647F9A", `WAKE` = "#8399AE")
fill2_multilevelcoda  <- c(`LPA` = "#978787", `MVPA` = "#BBA9A7", `SB` = "#A1B2C2", `TST` = "#647F9A", `WAKE` = "#8399AE")

# Illustrations ----------------
data(mcompd)
data(psub)

head(mcompd)

# Transforming Multilevel Compositional Data ----------------------------------------------------
cilr <- complr(
  data = mcompd,
  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
  idvar = "ID",
  total = 1440
)

# a summary of complr object
summary(cilr)

# Fitting Bayesian Multilevel Models with Compositional Predictors ------------------------------------
m <- brmcoda(
  complr = cilr,
  formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
    wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
  warmup = 1000, iter = 2000, seed = 123,
  chains = 4, cores = 4, backend = "cmdstanr"
)

## Model Summary
summary(m)

## NOTE: results produced might not be identical as
# the seed and chain identifier determine the behavior of the underlying random number generator.
# For complete reproducibility, every aspect of the environment needs to be locked down from 
# the OS and version to the C++ compiler and version to the version of Stan and all dependent libraries.
# For further information, please see: https://mc-stan.org/docs/2_28/reference-manual/reproducibility.html

## Estimating and Interpreting Pivot Coordinate Coefficients
m_coordinates <- pivot_coord(m, method = "rotate")
summary(m_coordinates)

# Running Multilevel Compositional Substitution Analysis ---------------------------------------
sub_simple <- substitution(
  object = m,
  delta = 1:10,
  ref = "grandmean",
  level = c("between", "within")
)

# check base pairwise substitutionå
psub

# obtain substitution model summary for 10-min reallocation
print(summary(sub_simple, delta = 10, digits = 2), row.names = FALSE, class = FALSE)

# 5.4. Presenting substitution model results ------------------------------------------------------
plot_between <- plot(sub_simple, to = "WAKE", ref = "grandmean", level = "between") +
  scale_fill_manual(values = fill2_multilevelcoda) +
  scale_colour_manual(values = colour_multilevelcoda) +
  scale_x_continuous(limits = c(-10, 10),
                     breaks = c(-10, 0, 10),
                     labels = c(10, 0, 10)
  ) +
  scale_y_continuous(limits = c(-0.08, 0.10),
                     breaks = c(-0.05, 0, 0.05, 0.10),
  ) +
  labs(x = "Difference in time in bed at between level",
       y = "Difference in stress") +
  facet_wrap(ggplot2::vars(From, To),
             labeller = label_bquote(cols = .(as.character(From)) %<-% minutes %->% .(as.character(To))),
             strip.position = "bottom", ncol = 4) +
  theme_multilevelcoda()

plot_within <- plot(sub_simple, to = "WAKE", ref = "grandmean", level = "within") +
  scale_fill_manual(values = fill2_multilevelcoda) +
  scale_colour_manual(values = colour_multilevelcoda) +
  scale_x_continuous(limits = c(-10, 10),
                     breaks = c(-10, 0, 10),
                     labels = c(10, 0, 10)
  ) +
  scale_y_continuous(limits = c(-0.075, 0.10),
                     breaks = c(-0.05, 0, 0.05, 0.10),
  ) +
  labs(x = "Difference in time in bed at within level",
       y = "Difference in stress") +
  facet_wrap(ggplot2::vars(From, To),
             labeller = label_bquote(cols = .(as.character(From)) %<-% minutes %->% .(as.character(To))),
             strip.position = "bottom", ncol = 4) +
  theme_multilevelcoda()

grDevices::cairo_pdf(
  file = paste0("plot_sub", ".pdf"),
  width = 9,
  height = 7,
)
ggarrange(plot_between, plot_within, 
          nrow = 2, legend = "none",
          labels = c("A. Awake in Bed Reallocations at Between-person level", 
                     "B. Awake in Bed Reallocations at Within-person level"),
          hjust = 0,
          font.label = list(size = 13, color = "black", family = "Arial Narrow")
) #10x9
dev.off()
