
#########################################################
# Script 4 - Atnarko Camera Traps
# Multinomial regression 
#########################################################

# Load Packages
list.of.packages <- c("mclogit", "lme4", "lmerTest", "nlme", "MuMIn", "performance", "formattable", "sjPlot", "dplyr", "ggiraphExtra", "DescTools", "arm")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

data <- cam.multi[which(cam.multi$AgeSex != "Unknown adult"),]
data$Deployment.Location.ID <- as.factor(data$Deployment.Location.ID)
data$AgeSex <- as.factor(data$AgeSex)
data$AgeSex <- factor(data$AgeSex, levels = c("Adult male", "Adult female", 
                                              "Female with young", "Sub-adult"))
# Join co-variates
data <- data %>%
  dplyr::select(year, week, site_name, Deployment.Location.ID, Date_Time.Captured, DarkPoint, Treatment, riv.section.weekid, siteweek_id,TimeDiff,AgeSex) # need to add a treatweek_id

data <- data %>%
  group_by(year, week) %>%
  mutate(yearweek_id = paste(year, week, sep = "")) %>%
  mutate(treatweek_id = paste(Treatment,yearweek_id, sep = ""))

# join_by's per dataframe:
# Salmon... dataframe: salmon; column: riv.section.weekid
# Berries... dataframe: berries; column: yearweek_id
# Humans... dataframe: humans; column: treatweek_id
# Water level... dataframe: hydro; column: yearweek_id

master_multi <- left_join(data, salmon, by = "riv.section.weekid") %>%
  dplyr::select(-year.x, -yearweek_id.x)
master_multi$yearweek_id <- as.character(master_multi$yearweek_id)
master_multi <- left_join(master_multi, berries, by = "yearweek_id")
master_multi$year <- substr(master_multi$yearweek_id,1,4)
master_multi <- left_join(master_multi, humans, by = "treatweek_id") 
master_multi$yearweek_id <- master_multi$yearweek_id.x
master_multi <- left_join(master_multi, hydro, by = "yearweek_id") 

master_multi <- master_multi %>%
  dplyr::select(week.x, site_name.x, Deployment.Location.ID, Date_Time.Captured, DarkPoint, Treatment, riv.section.weekid, siteweek_id, TimeDiff, AgeSex, treatweek_id, total.biomass, number_spp_ripe, total_estimated_humans, year, weekly_mean_waterlevel.m, Treatment)
master_multi <- master_multi %>%
  rename(week= week.x,
         site_name = site_name.x)

# Scale by 2sd as recommended by Gelman et al
master_multi$total.biomass.rescale <- rescale(master_multi$total.biomass)
master_multi$estimated.humans.rescale <- rescale(master_multi$total_estimated_humans)
master_multi$number_spp_ripe <- as.numeric(master_multi$number_spp_ripe)
master_multi$num.ripe.rescale <- rescale(master_multi$number_spp_ripe)
master_multi$riv_lev.rescale <- rescale(master_multi$weekly_mean_waterlevel.m)
master_multi$TimeDiff.rescale <- rescale(master_multi$TimeDiff)

# Build multinomial models
master_multi$AgeSex <- as.factor(master_multi$AgeSex)

# Subsample so that each deployment location ID has at least one of each AgeSex
master_subsample <- master_multi %>%
  group_by(Deployment.Location.ID) %>%
  filter(all(c("Female with young", "Adult female", "Adult male", "Sub-adult") %in% AgeSex)) %>%
  ungroup
# Complete cases only
master_subsample <- na.omit(master_subsample)

# Null
m1 <- mblogit(formula = AgeSex ~ 1,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

# Environmental
m2 <- mblogit(AgeSex ~ total.biomass.rescale, 
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

m3 <- mblogit(AgeSex ~ total.biomass.rescale +
                num.ripe.rescale,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

m4 <- mblogit(AgeSex ~ total.biomass.rescale +
              TimeDiff.rescale,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

m5 <- mblogit(AgeSex ~ total.biomass.rescale +
                TimeDiff.rescale +
                num.ripe.rescale,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

# Environmental and humans
m6 <- mblogit(AgeSex ~ total.biomass.rescale +
              estimated.humans.rescale,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

m7 <- mblogit(AgeSex ~ total.biomass.rescale +
                Treatment,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

m8 <- mblogit(AgeSex ~ total.biomass.rescale +
                estimated.humans.rescale +
                Treatment,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

m9 <- mblogit(AgeSex ~ total.biomass.rescale +
                estimated.humans.rescale +
                num.ripe.rescale,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

m10 <- mblogit(AgeSex ~ total.biomass.rescale +
                Treatment +
                num.ripe.rescale,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

m11 <- mblogit(AgeSex ~ total.biomass.rescale +
                estimated.humans.rescale +
                Treatment +
                num.ripe.rescale,
              random = ~1 | Deployment.Location.ID,
              data = master_subsample)

m12 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 estimated.humans.rescale:total.biomass.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m13 <- mblogit(AgeSex ~ total.biomass.rescale +
                 Treatment +
                 Treatment:total.biomass.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m14 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 estimated.humans.rescale:total.biomass.rescale +
                 Treatment +
                 Treatment:total.biomass.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m15 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 estimated.humans.rescale:total.biomass.rescale +
                 num.ripe.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m16 <- mblogit(AgeSex ~ total.biomass.rescale +
                 Treatment +
                 Treatment:total.biomass.rescale +
                 num.ripe.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m17 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 estimated.humans.rescale:total.biomass.rescale +
                 Treatment +
                 Treatment:total.biomass.rescale +
                 num.ripe.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m18 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m19 <- mblogit(AgeSex ~ total.biomass.rescale +
                 Treatment +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m20 <- mblogit(AgeSex ~ total.biomass.rescale +
                 Treatment +
                 estimated.humans.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m21 <- mblogit(AgeSex ~ total.biomass.rescale +
                 Treatment +
                 num.ripe.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m22 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 num.ripe.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m23 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 Treatment +
                 num.ripe.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m24 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 estimated.humans.rescale:total.biomass.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m25 <- mblogit(AgeSex ~ total.biomass.rescale +
                 Treatment +
                 Treatment:total.biomass.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m26 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 estimated.humans.rescale:total.biomass.rescale +
                 Treatment +
                 Treatment:total.biomass.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m27 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 estimated.humans.rescale:total.biomass.rescale +
                 num.ripe.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m28 <- mblogit(AgeSex ~ total.biomass.rescale +
                 Treatment +
                 Treatment:total.biomass.rescale +
                 num.ripe.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

m29 <- mblogit(AgeSex ~ total.biomass.rescale +
                 estimated.humans.rescale +
                 estimated.humans.rescale:total.biomass.rescale +
                 Treatment +
                 Treatment:total.biomass.rescale +
                 num.ripe.rescale +
                 TimeDiff.rescale,
               random = ~1 | Deployment.Location.ID,
               data = master_subsample)

# Define a function to calculate AICc 
calculate_AICc <- function(model) {
  n <- nobs(model)
  k <- df.residual(model) + 1
  AICc <- AIC(model) + 2 * k * (k + 1) / (n - k - 1)
  return(AICc)
}

# Calculate AICc for each model
AICc_values <- sapply(list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23, m24, m25, m26, m27, m28, m29), calculate_AICc)
# Find the model with the minimum AICc
min_AICc <- min(AICc_values)
# Calculate Delta AIC
delta <- AICc_values - min_AICc
# Extract log-likelihood from model summary
log_likelihood <- -0.5 * AICc_values

# Calculate Akaike Weights
weights <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))
results <- data.frame(Model = c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12"," m13", "m14", "m15", "m16", "m17", "m18", "m19", "m20", "m21", "m22", "m23", "m24", "m25", "m26", "m27", "m28", "m29"),
                      AICc = AICc_values,
                      Delta = delta,
                      Log_Likelihood = log_likelihood,
                      Weight = weights)
results <- results[order(results$AICc),]
results$cumulative.weight = cumsum(results$Weight)
results$Weight <- as.numeric(results$Weight)
print(results)

# Calculate RVI of salmon, dirun, humans, and berries

# Salmon
rviSalmon <- sum(results$Weight) # All models had salmon
Salmon <- rviSalmon/nrow(results)
# Diurn
diurn_models <- c("m4","m5","m18","m19","m20","m21","m22","m23","m24","m25","m26","m27","m28","m29")
rviDirun <- sum(sapply(diurn_models, function(model) {
  subset(results, model == Model, select = Weight)$Weight
}))
Diurn <- rviDirun/14
# Humans
human_models <- c("m6","m8","m9","m11","m12","m14","m15","m17","m18","m20","m22","m23","m24","m26", "m27","m29")
rvihuman <- sum(sapply(human_models, function(model) {
  subset(results, model == Model, select = Weight)$Weight
}))
LandHum <- rvihuman/16
# Berries
berry_models <- c("m3","m5","m9","m10","m11","m15","m16","m17","m21","m22","m23","m27","m28","m29")
rviberry <- sum(sapply(berry_models, function(model) {
  subset(results, model == Model, select = Weight)$Weight
}))
Berry <- rviberry/14

# Coef plot
cp3 <- broom.mixed::tidy(m22, conf.int = TRUE)

cp3$model <- "m22"

cp3 <- filter (cp3, term != "Adult female~(Intercept)")
cp3 <- filter (cp3, term != "Female with young~(Intercept)")
cp3 <- filter (cp3, term != "Sub-adult~(Intercept)")

cp3 <- cp3 %>% arrange(factor(term, levels = c(
  "Female with young~total.biomass.rescale" ,
  "Female with young~TimeDiff.rescale" ,
  "Female with young~estimated.humans.rescale",
  "Female with young~num.ripe.rescale" ,
  "Female with young~total.biomass.rescale:estimated.humans.rescale" ,
  "Adult female~total.biomass.rescale",
  "Adult female~TimeDiff.rescale" ,
  "Adult female~estimated.humans.rescale",
  "Adult female~num.ripe.rescale",
  "Adult female~total.biomass.rescale:estimated.humans.rescale",
  "Sub-adult~total.biomass.rescale" ,
  "Sub-adult~TimeDiff.rescale" ,
  "Sub-adult~estimated.humans.rescale" ,
  "Sub-adult~num.ripe.rescale",
  "Sub-adult~total.biomass.rescale:estimated.humans.rescale")))

#Relative variable importance
cp3$RVIscale <- c(Salmon, Diurn, LandHum, Berry, Salmon, Diurn, LandHum, Berry, Salmon, Diurn, LandHum, Berry)
cp3$RVIscale <- round(cp3$RVIscale, digits = 3)
dotwhisker::dwplot(cp3)

# MS Figure
cp3 <-dotwhisker::dwplot(cp3,
                         dot_args = list(color = "red"), # color for the dot
                         whisker_args = list(color = "black")) + 
  scale_y_discrete(breaks = c( "Female with young~total.biomass.rescale" ,
                               "Female with young~TimeDiff.rescale" ,
                               "Female with young~estimated.humans.rescale",
                               "Female with young~num.ripe.rescale" ,
                               "Adult female~total.biomass.rescale",
                               "Adult female~TimeDiff.rescale" ,
                               "Adult female~estimated.humans.rescale",
                               "Adult female~num.ripe.rescale",
                               "Sub-adult~total.biomass.rescale" ,
                               "Sub-adult~TimeDiff.rescale" ,
                               "Sub-adult~estimated.humans.rescale" ,
                               "Sub-adult~num.ripe.rescale"), 
                   labels=c( "Female with young ~ Salmon biomass" ,
                             "Female with young ~ Diurnality" ,
                             "Female with young ~ Visitors",
                             "Female with young ~ Berries" ,
                             "Adult female ~ Salmon biomass",
                             "Adult female ~ Diurnality" ,
                             "Adult female ~ Visitors",
                             "Adult female ~ Berries",
                             "Sub-adult ~ Salmon biomass" ,
                             "Sub-adult ~ Diurnality" ,
                             "Sub-adult ~ Visitors" ,
                             "Sub-adult ~ Berries")) +
  geom_text(label = cp3$RVIscale, x = 2.6, y = 12:1, size = 3) +
  geom_text(label = "RVI", x = 2.6, y = 12.6, size = 2.5) +
  theme_bw() + xlab("Coefficient") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  coord_cartesian(clip = "off", xlim = c(-1.5, 2.6)) +
  geom_segment(aes(x = -3.6, xend = 2.8,
                   y = 4.4, yend = 4.4),
               color = "black", size = 0.4) +
  geom_segment(aes(x = -3.6, xend = 2.8,
                   y = 8.4, yend = 8.4),
               color = "black", size = 0.4) +
  geom_segment(aes(x = -3.6, xend = 2.8,
                   y = 0, yend = 0),
               color = "black", size = 0.4) +
  geom_segment(aes(x = -3.6, xend = -3.6,
                   y = 12.85, yend = 0),
  color = "black", size = 0.4) +
  geom_segment(aes(x = -3.6, xend = -1.7,
                   y = 12.85, yend = 12.85),
               color = "black", size = 0.4) +
  theme(plot.title = element_text(face="bold", hjust = -0.7), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cp3
ggsave(file="multino-coef.png", width=8, height=3.6, dpi=300)

library(ggeffects)

# Nocturnality plot
noc.plotm22 <- ggpredict(m22, terms = c("TimeDiff.rescale"))

# Un-center and scale the days since data by multiplying by sd and adding mean
# Mean and SD of TimeDiff

mean(master_subsample$TimeDiff)
sd(master_subsample$TimeDiff)

noc.plotm22$x <- (noc.plotm22$x*2.908812) + 6.349798 

plotm22.1 <- ggplot(noc.plotm22, aes(x = x, y = predicted[,"Female with young"])) + 
  geom_line(aes(color = "Female with young"), show.legend = FALSE) +
  geom_line(aes(y = predicted[,"Adult male"], color = "Adult male"), show.legend = FALSE) +
  geom_line(aes(y = predicted[,"Adult female"], color = "Adult female"), show.legend = FALSE) +
  geom_line(aes(y = predicted[,"Sub-adult"], color = "Sub-adult"), show.legend = FALSE) +
  xlab("Hours since darkest point") + ylab("Age sex probability") +
  geom_ribbon(aes(ymin = conf.low[,"Adult male"], ymax = conf.high[,"Adult male"], 
                  fill = "Adult male", color = NULL), alpha = .15) +
  geom_ribbon(aes(ymin = conf.low[,"Adult female"], ymax = conf.high[,"Adult female"], 
                  fill = "Adult female", color = NULL), alpha = .15) +
  geom_ribbon(aes(ymin = conf.low[,"Female with young"], ymax = conf.high[,"Female with young"], 
                  fill = "Female with young", color = NULL), alpha = .15) +
  geom_ribbon(aes(ymin = conf.low[,"Sub-adult"], ymax = conf.high[,"Sub-adult"], 
                  fill = "Sub-adult", color = NULL), alpha = .15) +
  theme_classic() + 
  guides(fill = guide_legend(override.aes = list(color = NA)), 
                           color = FALSE)   +
  theme(legend.title=element_blank(),
        text = element_text(size = 24))

plotm22.1

# Land human index
humansm22 <- ggpredict(m22, terms = c("estimated.humans.rescale"))
humansm22$x <- (humansm22$x*705.2613) + 336.6244 

plotm22.3 <- ggplot(humansm22, aes(x = x, y = predicted[,"Sub-adult"])) + 
  geom_line(aes(color = "Sub-adult"), show.legend = FALSE) +
    geom_line(aes(y = predicted[,"Adult male"], color = "Adult male"), show.legend = FALSE) +
    geom_line(aes(y = predicted[,"Adult female"], color = "Adult female"), show.legend = FALSE) +
    geom_line(aes(y = predicted[,"Female with young"], color = "Female with young"), show.legend = FALSE) +
  xlab("Visitors") + ylab("Age sex probability") +
    geom_ribbon(aes(ymin = conf.low[,"Adult male"], ymax = conf.high[,"Adult male"], 
                    fill = "Adult male", color = NULL), alpha = .15) +
    geom_ribbon(aes(ymin = conf.low[,"Adult female"], ymax = conf.high[,"Adult female"], 
                    fill = "Adult female", color = NULL), alpha = .15) +
  geom_ribbon(aes(ymin = conf.low[,"Sub-adult"], ymax = conf.high[,"Sub-adult"], 
                  fill = "Sub-adult", color = NULL), alpha = .15) +
    geom_ribbon(aes(ymin = conf.low[,"Female with young"], ymax = conf.high[,"Female with young"], 
                    fill = "Female with young", color = NULL), alpha = .15) +
  theme_classic() + 
  guides(fill = guide_legend(override.aes = list(color = NA)), 
         color = FALSE)  +
  theme(legend.title=element_blank(),
        text = element_text(size = 24))
plotm22.3

# MS Figure
ggarrange(plotm22.3,plotm22.1, 
          common.legend = TRUE, 
          legend = "left", 
          labels = c("A", "B"))
ggsave(file="multinom.analy.plot.png", width=19, height=9, dpi=300)