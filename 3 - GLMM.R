
#############################################
# Script 3 - Atnarko Camera Traps
# Weekly Detection Analysis 
#############################################

library(optimx)
library(lme4)
library(AICcmodavg)
library(arm)
library(MuMIn) 
library(coefplot)
library(car)
library(ggpubr)
library(corrplot)
library(spdep)
library(sp)
library(broom.mixed)
library(dplyr)

# setwd() 

grizz <- cam.sum %>%
  dplyr::select(year, week, site_name, Treatment, riv.section.weekid, siteweek_id, detections.per.siteweek) 

grizz <- grizz %>%
  group_by(year, week) %>%
  mutate(yearweek_id = paste(year, week, sep = "")) %>%
  mutate(treatweek_id = paste(Treatment,yearweek_id, sep = ""))

# join_by's per dataframe:
# Salmon... dataframe: salmon; column: riv.section.weekid
# Berries... dataframe: berries; column: yearweek_id
# Humans... dataframe: humans; column: treatweek_id
# Water level... dataframe: hydro; column: yearweek_id

salmon <- salmon %>%
  rename(riv.section.weekid = riversection_weekid)
master <- left_join(grizz, salmon, by = "riv.section.weekid") %>%
  dplyr::select(-year.x, -yearweek_id.x)
master$yearweek_id <- as.character(master$yearweek_id)

berries$yearweek_id <- as.character(berries$yearweek_id)
master <- left_join(master, berries, by = "yearweek_id") %>%
  dplyr::select(-year.y)

master <- left_join(master, humans, by = "treatweek_id") %>%
  dplyr::rename(year = year.x,
         site_name = site_name.x,
         yearweek_id = yearweek_id.x) %>%
  dplyr::select(-c(year.y, week, site_name.y, yearweek_id.y))

master <- left_join(master, hydro, by = "yearweek_id") %>%
  dplyr::rename(year = year.x) %>%
  dplyr::select(-c(year.y, week.y))

master <- master %>%
  dplyr::select(-c(week.x, year))

# Correlation between variables 
# salmon, humans, berries, and water level
corr.test <- master %>%
  dplyr::select(total.biomass, total_estimated_humans, weekly_mean_waterlevel.m, number_spp_ripe)
c <- cor(data.matrix(corr.test), use = "complete.obs")
c
corrplot(c, method = "color")

## GLMM ##

master[substr(master$yearweek_id,1,4) =="2019", "year"] <- "2019"
master[substr(master$yearweek_id,1,4) =="2020", "year"] <- "2020"
master[substr(master$yearweek_id,1,4) =="2021", "year"] <- "2021"

ggplot(master) +
  geom_point(aes(x = week, y = detections.per.siteweek)) +
  facet_grid(row = vars(year), scales = "free_x") +
  labs(y = "Detections per site-week", x = "Week") +
  theme_classic() # Each data point is a camera. 

ggplot(master, aes(x = detections.per.siteweek)) +
  geom_histogram(position = "identity",binwidth=1, fill="grey", color="black") +
  labs(x = "Detections per site-week", y = "Frequency") +
  theme_classic()
ggsave(file="frequency of number of detections.png", width=6, height=4, dpi=300)

# Scale by 2sd as recommended by Gelman et al
master$total.biomass.rescale <- rescale(master$total.biomass)
master$estimated.humans.rescale <- rescale(master$total_estimated_humans)
master$number_spp_ripe <- as.numeric(master$number_spp_ripe)
master$num.ripe.rescale <- rescale(master$number_spp_ripe)
master$riv_lev.rescale <- rescale(master$weekly_mean_waterlevel.m)

# Null
dr.m1 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name),
  data = master,
  family = poisson) 

# Food only models 
dr.m2 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale,
  data = master,
  family = poisson)

dr.m3 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale,
  data = master,
  family = poisson) 

dr.m4 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale,
  data = master,
  family = poisson) 

dr.m5 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    num.ripe.rescale,
  data = master,
  family = poisson) 

# Food and humans models

dr.m6 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale,
  data = master,
  family = poisson)

dr.m7 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    Treatment,
  data = master,
  family = poisson) 

dr.m8 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale +
    Treatment,
  data = master,
  family = poisson)

dr.m9 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale +
    num.ripe.rescale,
  data = master,
  family = poisson) 

dr.m10 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    Treatment,
  data = master,
  family = poisson) 

dr.m11 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    estimated.humans.rescale +
    Treatment,
  data = master,
  family = poisson) 

dr.m12 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master,
  family = poisson) 

dr.m13 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    Treatment +
    total.biomass.rescale:Treatment,
  data = master,
  family = poisson) 

dr.m14 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale +
    Treatment +
    total.biomass.rescale:Treatment,
  data = master,
  family = poisson,
  control = glmerControl(
    optimizer ='optimx', optCtrl=list(method='nlminb')))

dr.m15 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master,
  family = poisson)

dr.m16 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    Treatment +
    total.biomass.rescale:Treatment,
  data = master,
  family = poisson)

dr.m17 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master,
  family = poisson)

dr.m18 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    estimated.humans.rescale,
  data = master,
  family = poisson)

dr.m19 <- glmer(
  detections.per.siteweek ~
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment,
  data = master,
  family = poisson)

dr.m20 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    estimated.humans.rescale,
  data = master,
  family = poisson)

dr.m21 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    num.ripe.rescale +
    estimated.humans.rescale,
  data = master,
  family = poisson)

dr.m22 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    num.ripe.rescale +
    Treatment,
  data = master,
  family = poisson)

dr.m23 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    num.ripe.rescale +
    Treatment +
    estimated.humans.rescale,
  data = master,
  family = poisson)

dr.m24 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master,
  family = poisson)

dr.m25 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    total.biomass.rescale:Treatment,
  data = master,
  family = poisson)

dr.m26 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    total.biomass.rescale:Treatment +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master,
  family = poisson,
  control = glmerControl(
    optimizer ='optimx', optCtrl=list(method='nlminb')))

dr.m27 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale +
    num.ripe.rescale,
  data = master,
  family = poisson)

dr.m28 <- glmer(
  detections.per.siteweek ~
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    total.biomass.rescale:Treatment +
    num.ripe.rescale,
  data = master,
  family = poisson)

dr.m29 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    total.biomass.rescale:Treatment +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale +
    num.ripe.rescale,
  data = master,
  family = poisson)

AICc(dr.m1, dr.m2, dr.m3, dr.m4, dr.m5, dr.m6, dr.m7, dr.m8, dr.m9, dr.m10, dr.m11, dr.m12, dr.m13, dr.m14, dr.m15, dr.m16, dr.m17, dr.m18, dr.m19, dr.m20, dr.m21, dr.m22, dr.m23, dr.m24, dr.m25, dr.m26, dr.m27, dr.m28, dr.m29)

# Model selection table with AICc 
model_sel.dr.AICc <- model.sel(dr.m1, dr.m2, dr.m3, dr.m4, dr.m5, dr.m6, dr.m7, dr.m8, dr.m9, dr.m10, dr.m11, dr.m12, dr.m13, dr.m14, dr.m15, dr.m16, dr.m17, dr.m18, dr.m19, dr.m20, dr.m21, dr.m22, dr.m23, dr.m24, dr.m25, dr.m26, dr.m27, dr.m28, dr.m29)
model_sel.dr.AICc$cumulative.weight = cumsum(model_sel.dr.AICc$weight)

# spatial autocorrelation in the residuals
res <- resid(dr.m29)
log_transformed_residuals <- log(abs(res) + 1)
plot(log_transformed_residuals, ylab = "Log-Transformed Residuals")
hist(log_transformed_residuals)
qqnorm(log_transformed_residuals)
qqline(log_transformed_residuals)

#Moran's I

#Prepare coordinate data
df <- as.data.frame(dr.m29@frame)
sa <- df[, c("detections.per.siteweek", "year", "site_name")]
sta2 <- sta[, c("Year", "site_name", "Deployment.Location.ID", "Latitude", "Longitude")]
sta2$year <- sta2$Year

coordinates_agg <- sta2 %>%
  group_by(Deployment.Location.ID, year, site_name) %>%
  summarize(Latitude = first(Latitude),
            Longitude = first(Longitude))
sa <- merge(sa, coordinates_agg, by = c("year", "site_name"))
coordinates(sa) <- c("Longitude", "Latitude")
w <- knn2nb(knearneigh(coordinates(sa), k = 5))
lw <- nb2listw(w)
moran.test(log_transformed_residuals, lw)

# Moran I test under randomisation
# 
# data:  log_transformed_residuals  
# weights: lw    
# 
# Moran I statistic standard deviate = 1.2323, p-value = 0.1089
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.0275811201     -0.0021321962      0.0005813663 

# R2
r2 <- r.squaredGLMM(dr.m29)

#Relative variable importance
model_sel.dr.rvi <- as.data.frame(model_sel.dr.AICc)
utils::View(model_sel.dr.rvi)

SalmonTable.dr <- subset(model_sel.dr.rvi, (!is.na(model_sel.dr.rvi[,2])))
rviSalmon.dr <- sum(SalmonTable.dr$weight)
Salmon.dr <- rviSalmon.dr/nrow(SalmonTable.dr)

BerryTable.dr <- subset(model_sel.dr.rvi, (!is.na(model_sel.dr.rvi[,3])))
rviBerry.dr <- sum(BerryTable.dr$weight)
Berry.dr <- rviBerry.dr/nrow(BerryTable.dr)

WaterTable.dr <- subset(model_sel.dr.rvi, (!is.na(model_sel.dr.rvi[,4])))
rviWater.dr <- sum(WaterTable.dr$weight)
Water.dr <- rviWater.dr/nrow(WaterTable.dr)

WaterSalmonTable.dr <- subset(model_sel.dr.rvi, (!is.na(model_sel.dr.rvi[,5])))
rviWaterSalmon.dr <- sum(WaterSalmonTable.dr$weight)
WaterSalmon.dr <- rviWaterSalmon.dr/nrow(WaterSalmonTable.dr)

LandHumTable.dr <- subset(model_sel.dr.rvi, (!is.na(model_sel.dr.rvi[,6])))
rviLandHum.dr <- sum(LandHumTable.dr$weight)
LandHum.dr <- rviLandHum.dr/nrow(LandHumTable.dr)

Treat.dr <- subset(model_sel.dr.rvi, (!is.na(model_sel.dr.rvi[,7])))
rviTreat.dr <- sum(Treat.dr$weight)
Treat.dr <- rviTreat.dr/nrow(Treat.dr)

LandHumSalmonTable.dr <- subset(model_sel.dr.rvi, (!is.na(model_sel.dr.rvi[,8])))
rviLandHumSalmon.dr <- sum(LandHumSalmonTable.dr$weight)
LandHumSalmon.dr <- rviLandHumSalmon.dr/nrow(LandHumSalmonTable.dr)

TreatSalmon.dr <- subset(model_sel.dr.rvi, (!is.na(model_sel.dr.rvi[,9])))
rviTreatSalmon.dr <- sum(TreatSalmon.dr$weight)
TreatSalmon.dr <- rviTreatSalmon.dr/nrow(TreatSalmon.dr)

# Coefficient plots
coef.p <- broom.mixed::tidy(dr.m29, conf.int = TRUE)
coef.p$model <- "dr.m29"
coef.p <- filter (coef.p, term != "sd__(Intercept)")
coef.p <- filter (coef.p, term != "(Intercept)")

coef.p <- coef.p %>% arrange(factor(term, levels = c(
  "total.biomass.rescale:estimated.humans.rescale",
  "total.biomass.rescale:TreatmentLandbased",
  "total.biomass.rescale:TreatmentReference",
  "TreatmentReference",
  "TreatmentLandbased",
  "riv_lev.rescale",
  "num.ripe.rescale",
  "estimated.humans.rescale",
  "total.biomass.rescale:riv_lev.rescale",
  "total.biomass.rescale")))

#Relative variable importance
coef.p$RVIscale <- c(LandHumSalmon.dr,
                    TreatSalmon.dr,
                    TreatSalmon.dr,
                    Treat.dr,
                    Treat.dr,
                    Water.dr,
                    Berry.dr,
                    LandHum.dr,
                    WaterSalmon.dr,
                    Salmon.dr)
coef.p$RVIscale <- round(coef.p$RVIscale, digits = 3)
dotwhisker::dwplot(coef.p)
coef.p <-dotwhisker::dwplot(coef.p,
                         dot_args = list(color = "red"), # color for the dot
                         whisker_args = list(color = "black")) + 
  scale_y_discrete(labels=c( "Salmon biomass" ,
                             "Salmon biomass : Water level" ,
                             "Visitors",
                             "Number species with ripe berries" ,
                             "Water level",
                             "Land-based tour area" ,
                             "No tour area",
                             "Salmon biomass : No tour area",
                             "Salmon biomass : Land-based tour area" ,
                             "Salmon biomass : Visitors")) +
  geom_text(label = coef.p$RVIscale, x = 2.3, y = 10:1, size = 3.5) +
  geom_text(label = "RVI", x = 2.4, y = 10.5, size = 3) +
  theme_bw() + xlab("Coefficient") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  coord_cartesian(clip = "off", ylim = c(1, 10.2), xlim = c(-2.5, 2.7)) +
  theme(plot.title = element_text(face="bold", hjust = -0.7), 
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

coef.p

require(ggiraph)
require(ggiraphExtra)
library(ggeffects)

# plot predicted values
# Land Human Index
hum.p <- ggpredict(dr.m29, "estimated.humans.rescale") 

# Un-center and scale human data by multiplying by sd and adding mean
# Mean and SD of human land index
getsd.and.mean <- master[!is.na(master$total_estimated_humans),]
mean(getsd.and.mean$total_estimated_humans)
sd(getsd.and.mean$total_estimated_humans)
hum.p$x <- (hum.p$x*705.2613) + 336.6244 

# What is the average land human index in each treatment? 
human.stats <- 
  getsd.and.mean %>%
  group_by(Treatment, year, yearweek_id) %>%
  summarise("total.hum.per.treat.per.week" = total_estimated_humans)

human.stats <- human.stats %>%
  distinct(yearweek_id, .keep_all = TRUE)

# Below gets the mean weekly human numbers per treatment across years
human.stats.by.treatment.across.years <- human.stats %>%
  group_by(Treatment) %>%
  mutate("weekly.mean.hum.per.treatment" = mean(total.hum.per.treat.per.week)) %>%
  mutate("weekly.max.hum.per.treatment" = max(total.hum.per.treat.per.week)) %>%
  mutate("weekly.min.hum.per.treatment" = min(total.hum.per.treat.per.week))

# 1363 is the mean number of humans in the platform drift area across years
# 11 is the mean in the reference

hum.p <- hum.p %>% plot()+
  labs(y = "Grizzly bear detections per site-week", x = "Visitors", title = "") +
  theme_classic() +
  theme(plot.title = element_text(face="bold", hjust = -0.7), 
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) 
hum.p

# Land hum and index and salmon biomass interaction 
humsal.p <- ggpredict(dr.m29, c("estimated.humans.rescale", "total.biomass.rescale [meansd]")) 
# [meansd] returns predictions for the values one SD below and above mean, as well as mean.

humsal.p$x <- (humsal.p$x*705.2613) + 336.6244 

roc_humsal <- humsal.p

humsal.p <- humsal.p %>%
  plot() + 
  labs(y = "Grizzly bear detections per site-week", x = "Visitors", colour  = "Salmon biomass (kg)", title = "") +
  scale_color_manual(values=c('#ff7f00','#4daf4a', '#984ea3'))+
  scale_color_manual(values=c('#ff7f00','#4daf4a', '#984ea3'),labels = c("-1 SD", "Mean", "+1 SD"))+
  scale_fill_manual(values=c('#ff7f00','#4daf4a', '#984ea3'), name="fill") +
#  geom_vline(xintercept=1363, linetype='dotted', col = 'red', size = 1)+
  theme_classic()+
  theme(legend.position = c(0.3, .85),
        legend.background = element_rect(fill = "transparent"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

humsal.p

# Rate of change per x humans under the different salmon conditions
roc_humsal <- roc_humsal %>%
  group_by(group) %>%
  arrange(group, x) %>%
  mutate(percent.change.in.detect.per.141.humans = 100 * (predicted - lag(predicted))/lag(predicted)) %>% # This is the rate of change per 141.0523 humans
  ungroup()

# In high salmon scenarios, the rate of change is -18.67946 per 141.0523 increase in humans. 
# (-18.67946/141.0523) = (x/100) to get the rate of change per 100 people

roc_humsal <- roc_humsal %>%
  mutate(per.500.people = 500*(roc_humsal$percent.change.in.detect.per.141.humans/141.0523)) %>%
  mutate(per.250.people = 250*(roc_humsal$percent.change.in.detect.per.141.humans/141.0523)) %>%
  mutate(per.100.people = 100*(roc_humsal$percent.change.in.detect.per.141.humans/141.0523)) %>%
  mutate(per.50.people = 50*(roc_humsal$percent.change.in.detect.per.141.humans/141.0523)) %>%
  mutate(per.1.person = 1*(roc_humsal$percent.change.in.detect.per.141.humans/141.0523)) %>%
  ungroup()

# Treatment and salmon biomass interaction
treatsal.p <- ggpredict(dr.m29, c("Treatment", "total.biomass.rescale")) 
treatsal.p <- treatsal.p %>% 
  plot(dodge = 0.4) + 
  scale_x_discrete(limits = c("Platform/Drift", "Landbased", "Reference"), labels = c("Land- & boat-based", "Land-based", "No tour")) +
  scale_color_manual(values=c('#ff7f00','#4daf4a', '#984ea3'),labels = c("-1 SD", "Mean", "+1 SD"))+
  labs(y = "Grizzly bear detections per site-week", x = "Spatial treatment", colour= "Salmon biomass (kg)", title = "") +
  theme_classic() +
  theme(legend.position = c(0.80, .85),
        legend.background = element_rect(fill = "transparent"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) 
treatsal.p

# Use ggpredict for predicted values
treatsal.p2 <- ggpredict(dr.m29, c("Treatment", "total.biomass.rescale"))

# Plot the interaction effect
treatsal.plot <- plot(treatsal.p)

# Customize plot
treatsal.plot <- treatsal.plot +
  scale_x_discrete(limits = c("Platform/Drift", "Landbased", "Reference"),
                   labels = c("Land- & boat-based", "Land-based", "No tour")) +
  scale_color_manual(values = c('#ff7f00', '#4daf4a', '#984ea3'),
                     labels = c("-1 SD", "Mean", "+1 SD")) +
  labs(y = "Grizzly bear detections per site-week",
       x = "Spatial treatment",
       colour = "Salmon biomass (kg)",
       title = "") +
  theme_classic() +
  theme(legend.position = c(0.80, .85),
        legend.background = element_rect(fill = "transparent"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# MS Figure
ggarrange(coef.p, hum.p, treatsal.p,humsal.p, labels = c("A", "B", "C", "D"))
ggsave(file="Figure2.png", width=12, height=8, dpi=300)

# Assess how our methods affected the general patterns in our bear activity data by running models with and without imputed 2019 values. 
master_omit2019 <- subset(master, master$year != "2019") 

# Null
dr.omit2019.m1 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name),
  data = master_omit2019 ,
  family = poisson) 

# Food only models 
dr.omit2019.m2 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m3 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale,
  data = master_omit2019 ,
  family = poisson) 

dr.omit2019.m4 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale,
  data = master_omit2019 ,
  family = poisson) 

dr.omit2019.m5 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    num.ripe.rescale,
  data = master_omit2019 ,
  family = poisson) 

# Food and humans models

dr.omit2019.m6 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m7 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    Treatment,
  data = master_omit2019 ,
  family = poisson) 

dr.omit2019.m8 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale +
    Treatment,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m9 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale +
    num.ripe.rescale,
  data = master_omit2019 ,
  family = poisson) 

dr.omit2019.m10 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    Treatment,
  data = master_omit2019 ,
  family = poisson) 

dr.omit2019.m11 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    estimated.humans.rescale +
    Treatment,
  data = master_omit2019 ,
  family = poisson) 

dr.omit2019.m12 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson) 

dr.omit2019.m13 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    Treatment +
    total.biomass.rescale:Treatment,
  data = master_omit2019 ,
  family = poisson) 

dr.omit2019.m14 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale +
    Treatment +
    total.biomass.rescale:Treatment,
  data = master_omit2019 ,
  family = poisson,
  control = glmerControl(
    optimizer ='optimx', optCtrl=list(method='nlminb')))

dr.omit2019.m15 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m16 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    Treatment +
    total.biomass.rescale:Treatment,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m17 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    num.ripe.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m18 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m19 <- glmer(
  detections.per.siteweek ~
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m20 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m21 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    num.ripe.rescale +
    estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m22 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    num.ripe.rescale +
    Treatment,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m23 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    num.ripe.rescale +
    Treatment +
    estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m24 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m25 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    total.biomass.rescale:Treatment,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m26 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    total.biomass.rescale:Treatment +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale,
  data = master_omit2019 ,
  family = poisson,
  control = glmerControl(
    optimizer ='optimx', optCtrl=list(method='nlminb')))

dr.omit2019.m27 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale +
    num.ripe.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m28 <- glmer(
  detections.per.siteweek ~
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    total.biomass.rescale:Treatment +
    num.ripe.rescale,
  data = master_omit2019 ,
  family = poisson)

dr.omit2019.m29 <- glmer(
  detections.per.siteweek ~ 
    (1|year) +
    (1|site_name) +
    total.biomass.rescale +
    riv_lev.rescale +
    total.biomass.rescale:riv_lev.rescale +
    Treatment +
    total.biomass.rescale:Treatment +
    estimated.humans.rescale +
    total.biomass.rescale:estimated.humans.rescale +
    num.ripe.rescale,
  data = master_omit2019 ,
  family = poisson,
  control = glmerControl(
    optimizer ='optimx', optCtrl=list(method='nlminb')))

# AICc
AICc(dr.omit2019.m1, dr.omit2019.m2, dr.omit2019.m3, dr.omit2019.m4, dr.omit2019.m5, dr.omit2019.m6, dr.omit2019.m7, dr.omit2019.m8, dr.omit2019.m9, dr.omit2019.m10, dr.omit2019.m11, dr.omit2019.m12, dr.omit2019.m13, dr.omit2019.m14, dr.omit2019.m15, dr.omit2019.m16, dr.omit2019.m17, dr.omit2019.m18, dr.omit2019.m19, dr.omit2019.m20, dr.omit2019.m21, dr.omit2019.m22, dr.omit2019.m23, dr.omit2019.m24, dr.omit2019.m25, dr.omit2019.m26, dr.omit2019.m27, dr.omit2019.m28, dr.omit2019.m29)
# Model selection table with AICc 
model_sel.dr_omit2019_AICc <- model.sel(dr.omit2019.m1, dr.omit2019.m2, dr.omit2019.m3, dr.omit2019.m4, dr.omit2019.m5, dr.omit2019.m6, dr.omit2019.m7, dr.omit2019.m8, dr.omit2019.m9, dr.omit2019.m10, dr.omit2019.m11, dr.omit2019.m12, dr.omit2019.m13, dr.omit2019.m14, dr.omit2019.m15, dr.omit2019.m16, dr.omit2019.m17, dr.omit2019.m18, dr.omit2019.m19, dr.omit2019.m20, dr.omit2019.m21, dr.omit2019.m22, dr.omit2019.m23, dr.omit2019.m24, dr.omit2019.m25, dr.omit2019.m26, dr.omit2019.m27, dr.omit2019.m28, dr.omit2019.m29)
model_sel.dr_omit2019_AICc$cumulative.weight = cumsum(model_sel.dr_omit2019_AICc$weight)

# Coefficient plots
# Model dr.omit2019.m29

coef.p_omit2019.m29 <- broom.mixed::tidy(dr.omit2019.m29, conf.int = TRUE)

coef.p_omit2019.m29$model <- "dr.omit2019.m29"

coef.p_omit2019.m29 <- filter (coef.p_omit2019.m29, term != "sd__(Intercept)")
coef.p_omit2019.m29 <- filter (coef.p_omit2019.m29, term != "(Intercept)")

coef.p_omit2019.m29 <- coef.p_omit2019.m29 %>% arrange(factor(term, levels = c(
  "total.biomass.rescale:estimated.humans.rescale",
  "total.biomass.rescale:TreatmentLandbased",
  "total.biomass.rescale:TreatmentReference",
  "TreatmentReference",
  "TreatmentLandbased",
  "riv_lev.rescale",
  "num.ripe.rescale",
  "estimated.humans.rescale",
  "total.biomass.rescale:riv_lev.rescale",
  "total.biomass.rescale")))

dotwhisker::dwplot(coef.p_omit2019.m29)

coef.p_omit2019.m29 <-dotwhisker::dwplot(coef.p_omit2019.m29,
                            dot_args = list(color = "red"), # color for the dot
                            whisker_args = list(color = "black")) + 
  scale_y_discrete(labels=c( "Salmon biomass" ,
                             "Salmon biomass : Water level" ,
                             "Visitors",
                             "Number species with ripe berries" ,
                             "Water level",
                             "Land-based tour area" ,
                             "No tour area",
                             "Salmon biomass : No tour area",
                             "Salmon biomass : Land-based tour area" ,
                             "Salmon biomass : Visitors")) +
  theme_bw() + xlab("Coefficient") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  coord_cartesian(clip = "off", ylim = c(1, 10.2), xlim = c(-3, 2.5)) +
  theme(plot.title = element_text(face="bold", hjust = -0.7), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Build coefficient plot that compares the two models with and without 2019 data. 

df_dr.omit2019.m29 <- data.frame(Variable = rownames(summary(dr.omit2019.m29)$coef),
                                 Coefficient = summary(dr.omit2019.m29)$coef[, 1],
                                 SE = summary(dr.omit2019.m29)$coef[, 2],
                                 modelName = "Model 29: 2019 data omitted")

df_dr.m29 <- data.frame(Variable = rownames(summary(dr.m29)$coef),
                        Coefficient = summary(dr.m29)$coef[, 1],
                        SE = summary(dr.m29)$coef[, 2],
                        modelName = "Model 29: All years")

allModelFrame <- data.frame(rbind(df_dr.omit2019.m29, df_dr.m29))
interval <- -qnorm((1-0.95)/2)  # 95% multiplier

allModelFrame <- filter (allModelFrame, Variable != "(Intercept)")
allModelFrame$term <- as.factor(allModelFrame$Variable)

allModelFrame$term <- factor(x = allModelFrame$term, 
                             levels = c(
                               "total.biomass.rescale",
                               "total.biomass.rescale:riv_lev.rescale",
                               "estimated.humans.rescale",
                               "num.ripe.rescale",
                               "riv_lev.rescale",
                               "TreatmentLandbased",
                               "TreatmentReference",
                               "total.biomass.rescale:TreatmentReference",
                               "total.biomass.rescale:TreatmentLandbased",
                               "total.biomass.rescale:estimated.humans.rescale"))
# SI Figure
compare <- ggplot(allModelFrame, aes(colour = modelName))
compare <- compare + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
compare <- compare + geom_linerange(aes(x = term, ymin = Coefficient - SE*interval,
                                        ymax = Coefficient + SE*interval),
                                    lwd = 1, position = position_dodge(width = 1/2))
compare <- compare + geom_pointrange(aes(x = term, y = Coefficient, ymin = Coefficient - SE*interval,
                                         ymax = Coefficient + SE*interval),
                                     lwd = 1/2, position = position_dodge(width = 1/2),
                                     shape = 21, fill = "WHITE") +
  xlab("") +
  scale_x_discrete(labels=c("total.biomass.rescale:estimated.humans.rescale" = "Salmon biomass : Visitors",
                            "total.biomass.rescale:TreatmentLandbased" = "Salmon biomass : Land-based tour area" ,
                            "total.biomass.rescale:TreatmentReference" = "Salmon biomass : No tour area",
                            "TreatmentReference" = "No tour area",
                            "TreatmentLandbased" = "Land-based tour area" ,
                            "riv_lev.rescale" = "Water level",
                            "num.ripe.rescale" = "Number species with ripe berries" ,
                            "estimated.humans.rescale" = "Visitors",
                            "total.biomass.rescale:riv_lev.rescale" = "Salmon biomass : Water level" ,
                            "total.biomass.rescale" = "Salmon biomass")) +
  scale_colour_brewer(palette = "Dark2", name ="") 
compare <- compare + coord_flip() + theme_classic() +
  theme(legend.position = c(0.85, .70)) 
print(compare) 
ggsave(file="model-compare.png", width=9, height=3.5, dpi=300)

# Management table
# Create a management table. Explore the predicted bear detections from top model output under each salmon scenario. 

# What are the following bear predictions under each salmon scenario? 
# 50% of mean predicted bears
# mean predicted bears
# 200% of the mean predicted bears

roc_humsal <- roc_humsal %>%
  group_by(group) %>%
  arrange(group, x) %>%
  group_by(group) %>%
  mutate("mean.predicted.bears" = mean(predicted)) %>%
  mutate("50%.mean.predicted.bears" = (mean(predicted)*.5)) %>%
  mutate("200%.mean.predicted.bears" = (mean(predicted)*2)) 

# What is the human output from the model for each of the following bear scenarios, and under each salmon scenario (high,mean,low)?
# 50% of mean predicted bears
# mean predicted bears
# 200% of the mean predicted bears

# Use the following table to identify how many humans are associated with predicted grizzly bears under each of the scenarios above. 
management.table <- ggpredict(dr.m29, terms = c("estimated.humans.rescale [-0.5:5 by =.001]", "total.biomass.rescale [meansd]"))
management.table$x <- (management.table$x*705.2613) + 336.6244 

management.table <- management.table %>%
  group_by(group) %>%
  arrange(group, x)

write.csv(management.table, "management.table.csv")
