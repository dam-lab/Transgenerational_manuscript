
surv.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Survival/"
fit.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Fitness/"
epr.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/EPR/"
dev.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Development_time/"
sex.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Sex_ratio/"

library(popbio)
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(broom)
library(sjPlot)
library(car)
library(emmeans)
library(itsadug)
library(mgcv)
library(lme4)
library(itsadug)
library(lm.beta)
library(mnormt)
library(lavaan)
library(semPlot)
library(glmmTMB)
library(ggpubr)


#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################

######################### SURVIVAL DATA ###############################
SurvTot <- fread(paste(surv.directory, "SurvDataFiles/Survival_data_total.txt", sep = ""))



##### Creating survival curves for each generation and treatment #####

### Original: Pop_growth_epr_and_survival.R
Surv.files <- filter(SurvTot, nx > 0)

Surv.files$nx <- as.numeric(Surv.files$nx)

Surv.files <- Surv.files[rep(seq(nrow(Surv.files)), 
                             Surv.files$nx),#the column which holds the value you wish to have everything repeated by
]


Surv.files <- Surv.files %>%
  group_by(Generation, Treatment, Rep, Beak) %>%
  mutate(
    survival = case_when(
      # Food == 0 ~ 1,# all starved treatments did not survive
      nx == min(nx) ~ 0, #need to use the "==" and VALUE OF 1 EQUALS DEATH HAPPENING!!!!!!
      nx > min(nx) ~ 1 )) %>% 
  group_by(Generation, #have to add Generation when you are separating by generation!!!!
           Treatment, 
           Rep, 
           Beak) %>%
  mutate(keep = max(nx) - nx) %>% 
  group_by(Generation, 
           Treatment, 
           Rep, 
           Beak, 
           keep) %>% #need to group by unique and keep to get all the animals that died in between the start and the lowest amount of dead animals
  filter(survival == 1 & nx == min(nx) |
           survival == 0) %>% 
  filter(row_number() %in% 1:unique(keep) |
           survival == 0) %>% 
  select(-keep) %>% 
  ungroup()

#get rid of starting point since they all started alive

Surv.files <- Surv.files %>% 
  group_by(Generation, Treatment, Rep, Beak) %>% 
  filter(nx != max(nx))


#Surv.files <- Surv.files[!Surv.files$nx == 25,]
Surv.files <- as.data.frame(Surv.files)

# put the table in order
Surv.files$Generation <- as.numeric(Surv.files$Generation)
Surv.files <- Surv.files[order(Surv.files$Generation),]



Surv.files <- unite(Surv.files,#the data frame
                    unique, #the name of the new column
                    c(Generation, Treatment, Rep, Beak), #the existing columns you want to combine
                    remove = FALSE)



# must be a data frame for the grouping variable to be sorted appropriately
Surv.files <- as.data.frame(Surv.files)


#Surv.files.split <- split(Surv.files, f = Surv.files$Gen.Treat)




# create survival lists for all the data grouped by Gen, Treat, Rep, and Beak
Surv.object <- Surv.files %>%
  surv_group_by(grouping.vars = "unique") %>%
  separate(data = ., unique, into = c("Generation", "Treatment", "Rep", "Beak"))

Surv.object


# create the survival fit curves for statistical analysis 

Surv.fits <- surv_fit(Surv(time, #the time indicator
                           survival # the event indicator
)~Temp+pH, # the variables to test survival by
data = Surv.files, # the DATA FRAME with data you want to fit
group.by = "Generation") 

## survMisc package



# stats for each generation by temp and pH
Surv.stats <- surv_pvalue(Surv.fits, method = "log", test.for.trend = TRUE)
Surv.stats


Surv.plots <-ggsurvplot(Surv.fits, 
                        color = "strata", 
                        conf.int = TRUE, 
                        palette = c("Green", "Blue", "Red", "Orange"),
                        test.for.trend = TRUE,
                        pval = TRUE)

Surv.plots


## create an accurate plot for Generation 12 which is missing AH data

Surv.files.Gen12 <- Surv.files[ which(Surv.files$Generation==12),]

Surv.fits.Gen12 <- surv_fit(Surv(time, #the time indicator
                                 
                                 survival) # the event indicator
                            
                            ~Temp+pH, # the variables to test survival by
                            
                            data = Surv.files.Gen12) # the DATA FRAME with data you want to fit



Surv.stats.Gen12 <- surv_pvalue(Surv.fits.Gen12)

Surv.plots.Gen12 <- ggsurvplot(Surv.fits.Gen12, color = "strata", conf.int = TRUE, palette = c("Blue", "Red", "Orange"))

Surv.plots.Gen12



##### Plots for mean survival across generations #####
SurvTot1 <- filter(SurvTot, SurvTot$lx > 0)



SurvTot.agg <- aggregate(SurvTot1$lx, 
                         list(SurvTot1$Treatment, SurvTot1$Rep, SurvTot1$Beak, SurvTot1$Generation), #use a list so you don't have to use "interaction"
                         min)#use the minimum value for each of the individual beakers which should be the final survival

colnames(SurvTot.agg) <- c("Treatment", "Rep", "Beak", "Generation", "lx")

is.even <- function(x) x %% 2 == 0



if (wd == "C:/Users/james/Documents/Grad_school/OA_Project/Survival/SurvDataFiles/") {
  
  SurvTot.agg$Temp <- if_else(SurvTot.agg$Treatment < 3, 18, 22)
  
  SurvTot.agg$pH <- if_else(is.even(SurvTot.agg$Treatment), 7.5, 8.2)
  
} else if (wd == "C:/Users/james/Documents/Grad_school/OA_hudsonica/Survival/SurvDataFiles/") {
  
  SurvTot.agg$Temp <- if_else(SurvTot.agg$Treatment < 3, 13, 15)
  
  SurvTot.agg$pH <- if_else(is.even(SurvTot.agg$Treatment), 7.8, 8.2)
  
  
}



SurvTot.agg <- unite(SurvTot.agg,
                     Treat.Rep,
                     c(Treatment, Rep),
                     remove = FALSE)

SurvTot.mean <- SurvTot.agg %>%
  group_by(Treatment, Generation) %>%
  summarise(mean = mean(lx, na.rm = TRUE),
            sd = sd(lx, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

SurvTot.agg$Generation.c <- as.numeric(SurvTot.agg$Generation)
SurvTot.agg$Generation <- as.factor(as.numeric(SurvTot.agg$Generation))
SurvTot.agg$Treatment <- as.factor(as.numeric(SurvTot.agg$Treatment))

SurvTot.mean$Generation.c <- as.numeric(SurvTot.mean$Generation)
SurvTot.mean$Generation <- as.factor(as.numeric(SurvTot.mean$Generation))

SurvTot.mean[is.na(SurvTot.mean)] <- 0

survPlotTotal <- ggplot(data = SurvTot.mean, aes(Generation.c, mean, color=factor(Treatment)))+
  geom_line(size=1,
            position = position_dodge(width = 2))+
  geom_point(size = 1, position = position_dodge(width = 2))+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width = 0,
                position = position_dodge(width = 2))+
  theme(legend.title = element_text(colour = "black", size=12))+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
                                )
                     )+ 
  theme_classic()+
  labs(y="Survival", x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00))+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  geom_text(x=15.5, y=0.51, size = 5, label = "*", color = "black")+
  geom_text(x=16, y=0.6, size = 5, label = "*", color = "black")+
  geom_text(x=26, y=0.3, size = 5, label = "*", color = "black")

survPlotTotal

ggsave(filename = paste0(surv.directory,"Survival_plot.pdf"), plot = survPlotTotal, height = 101, width = 180, units = "mm")


### Boxplot for supplementary information
SurvTot.agg$Generation <- factor(SurvTot.agg$Generation, levels = c("0","3","6", "9","12","15","25"), ordered = T, 
                                 labels = c("0","3","6", "9","12","15","25"))
SurvTot.mean$Generation <- factor(SurvTot.mean$Generation, levels = c("0","3","6", "9","12","15","25"), ordered = T, 
                                  labels = c("0","3","6", "9","12","15","25"))

survBoxplot <- ggplot()+
  geom_boxplot(data = SurvTot.agg, aes(x = Generation, y = lx, fill = factor(Treatment)), lwd = 1.1,
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = SurvTot.mean, aes(Generation, mean, color=factor(Treatment)),
             size = 3, 
             shape = 18,
             position = position_dodge(width = 0.75)
  ) +
  
  
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  ))+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  )
  )+
  theme_classic()+
  labs(y="Survival", x="Generation")+
  scale_x_discrete(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


survBoxplot

ggsave(filename = paste0(surv.directory,"Surv_boxplot.pdf"), plot = survBoxplot, height = 101, width = 180, units = "mm")

## Create dataframes for fitness landscapes

SurvTot.agg.0 <- filter(SurvTot.agg, Generation.c == 0)

last.gen <- max(SurvTot.agg$Generation.c)

SurvTot.agg.last <- filter(SurvTot.agg, Generation == last.gen)
SurvTot.agg.0.last <- rbind(SurvTot.agg.0, SurvTot.agg.last)

SurvTot.agg.0.last$Treatment2 <- case_when(SurvTot.agg.0.last$Treatment == 1 ~ "AM",
                                           SurvTot.agg.0.last$Treatment == 2 ~ "OA",
                                           SurvTot.agg.0.last$Treatment == 3 ~ "OW",
                                           SurvTot.agg.0.last$Treatment == 4 ~ "OWA")

# change the order of the factor to be in the order you want. Help from: https://stackoverflow.com/questions/5490638/how-to-change-the-order-of-facet-labels-in-ggplot-custom-facet-wrap-labels
SurvTot.agg.0.last <- within(SurvTot.agg.0.last, 
                             Treatment2 <- factor(Treatment2, 
                                                  
                                                  # put the specific levels in the order you want
                                                  
                                                  levels = c("AM", # 1
                                                             "OA", # 2
                                                             "OW", # 3
                                                             "OWA") # 4 
                                                  
                             ))

##### Statistics for Survival data #####


## model for pairwise tukey comparisons

SurvTot.lm <- lmer(lx~Generation*Treatment+(1|Treat.Rep), 
                   REML = FALSE,
                   #family = gaussian,
                   data = SurvTot.agg)
summary(SurvTot.lm)

plot_model(SurvTot.lm, type = 'int')

tab_model(SurvTot.lm)

surv.emm <- emmeans(SurvTot.lm, pairwise ~ Generation | Treatment)

pairs(surv.emm)

Surv.pairwise <- tidy(pairwise.t.test(SurvTot.agg$lx, SurvTot.agg$Generation:SurvTot.agg$Treatment, p.adjust.method = "none"))


## continuous model for a linear mixed model anova
SurvTot.lm.c <- lmer(lx~Generation.c*Treatment+(1|Treat.Rep), 
                     REML = FALSE,
                     #family = gaussian,
                     data = SurvTot.agg)

tab_model(SurvTot.lm.c)


## gam model

# s() indicates a smooth function
gam1 <- gam(lx ~ s(Generation.c, by = Treatment, k = 3), data = SurvTot.agg)

library(itsadug)
plot(gam1, pages = 1)
plot_smooth(gam1, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
gam.stats <- summary(gam1)



## three way anova
Surv.lm2 <- aov(lx ~ Generation.c*Temp*pH, data = SurvTot.agg)

Surv.3way.anova <- Anova(Surv.lm2)

Surv.3way.anova$factors <- rownames(Surv.3way.anova)

fwrite(Surv.3way.anova, file = paste(surv.directory, "Statistics/Survival_3way_anova.txt", sep = ""), sep = "\t")






##### Development Time #####


SurvTot.Cdev <- filter(SurvTot, SurvTot$Cdev > 0)


SurvTot.Cdev <- SurvTot.Cdev[rep(seq(nrow(SurvTot.Cdev)), SurvTot.Cdev$Cdev),]

SurvTot.Cdev <- unite(SurvTot.Cdev,
                      Treat.Rep,
                      c(Treatment, Rep),
                      remove = FALSE)

Cdev.sum <- SurvTot.Cdev %>%
  group_by(Treatment, Generation) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            sd = sd(time, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

fwrite(Cdev.sum, file = paste(dev.directory, "cdev_sum.txt", sep = ""), sep = "\t")

Cdev.sum$Generation.c <- as.numeric(Cdev.sum$Generation)


CdevPlotTotal <- ggplot(data = Cdev.sum, aes(Generation.c, mean, color=factor(Treatment)))+
  geom_line(size=1,
            position = position_dodge(width = 2))+
  geom_point(size = 1, position = position_dodge(width = 2))+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width=0, position = position_dodge(width = 2))+
  theme(legend.title = element_text(colour = "black", size=12))+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  )
  )+ 
  theme_classic()+
  labs(y="Total Development Time (days)", x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(breaks = c(10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15))+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))

CdevPlotTotal


ggsave(filename = paste0(dev.directory,"Dev_plot.pdf"), plot = CdevPlotTotal, height = 101, width = 180, units = "mm")

Cdev.sum$Generation <- factor(Cdev.sum$Generation, levels = c("0","3","6", "9","12","15","25"), ordered = T)
SurvTot.Cdev$Generation <- factor(SurvTot.Cdev$Generation, levels = c("0","3","6", "9","12","15","25"), ordered = T)

CdevBoxplot <- ggplot()+
  geom_boxplot(data = SurvTot.Cdev, aes(x = Generation, y = time, fill = factor(Treatment)), lwd = 1.1,
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = Cdev.sum, aes(Generation, mean, color=factor(Treatment)),
             size = 3, 
             shape = 18,
             position = position_dodge(width = 0.75)
  ) +
  
  #geom_errorbar(data = SurvTot.mean, aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width = 0,
  #             position = position_dodge(width = 2))+
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#0c7bdc", #AM
                               "#009e73", #OA
                               "#ffa500", #OW
                               "#d55e00" #OWA
  ))+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  )
  )+
  theme_classic()+
  labs(y="Total Development time (days)", x="Generation")+
  scale_x_discrete(breaks = c(0,3,6,9,12,15,25))+
  #scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


CdevBoxplot

ggsave(filename = paste0(dev.directory,"Cdev_boxplot.pdf"), plot = survBoxplot, height = 101, width = 180, units = "mm")

Cdev.plot.complete <- ggarrange(CdevPlotTotal, CdevBoxplot, ncol = 1, nrow = 2)

ggsave(filename = paste0(dev.directory, "Cdev_plots_complete.pdf"), plot = Cdev.plot.complete, height = 202, width = 180, units = "mm")


SurvTot.Cdev$Generation.c <- as.numeric(SurvTot.Cdev$Generation)

SurvTot.Cdev$Generation <- as.factor(as.numeric(SurvTot.Cdev$Generation))

SurvTot.Cdev$time <- as.numeric(SurvTot.Cdev$time)

SurvTot.Cdev$Treatment <- as.factor(SurvTot.Cdev$Treatment)



## plots to see how reps are different from each other

ggplot(data = SurvTot.Cdev, aes(Generation.c, time, color = factor(Rep))) +
  geom_point()+
  geom_smooth(method = "lm", formula = y~x)+
  facet_wrap(~as.factor(Treatment))

ggplot(data = SurvTot.Cdev, aes(Generation, time, color = factor(Rep))) +
  geom_boxplot()+
  facet_wrap(~as.factor(Treatment))



#SurvTot.Cdev$Generation[SurvTot.Cdev$Generation == 10] <- 9

SurvTot.Cdev$pH <- as.factor(as.numeric(SurvTot.Cdev$pH))
SurvTot.Cdev$Temp <- as.factor(as.numeric(SurvTot.Cdev$Temp))

##### Statistics for Development Time #####


## linear mixed models for pairwise comparisons

Cdev.lmm <- lmer(time ~ Generation * Treatment + (1|Treat.Rep),
                 REML = FALSE,
                 #family = gaussian,
                 data = SurvTot.Cdev)

tab_model(Cdev.lmm)


Cdev.emm <- emmeans(Cdev.lmm, pairwise ~ Generation | Treatment)

pairs(Cdev.emm)



Cdev.lmm2 <- lmer(time ~ Generation.c * Treatment + (1|Treat.Rep),
                  REML = FALSE,
                  #family = gaussian,
                  data = SurvTot.Cdev)

plot_model(Cdev.lmm2, 'int')



Cdev.pairwise <- tidy(pairwise.t.test(SurvTot.Cdev$time, SurvTot.Cdev$Generation:SurvTot.Cdev$Treatment, p.adjust.method = "none"))

Cdev.pairwise

#Cdev.pairwise <- filter(Cdev.pairwise, Cdev.pairwise$p.value < 0.05)
fwrite(Cdev.pairwise, file = paste(dev.directory,"Statistics/Dev_time_pairwise.txt", sep = ""), sep = "\t")

## model for three-way anova

Cdev.lm.2 <-  lm(time~Generation.c*Temp*pH, data = SurvTot.Cdev)

Cdev.anova.2 <- Anova(Cdev.lm.2)

Cdev.anova.2$factors <- rownames(Cdev.anova.2)

fwrite(Cdev.anova.2, file = paste(dev.directory,"Statistics/Dev_time_anova_ph_temp.txt", sep = ""), sep = "\t")

plot_model(Cdev.lm.2, 'int')

summary(Cdev.lm.2)


## gam model for analyzing development time across generations

Cdev.gam <- gam(time ~ s(Generation.c, by = Treatment, k = 3), data = SurvTot.Cdev)

plot(Cdev.gam, pages = 1)
plot_smooth(Cdev.gam, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
summary(Cdev.gam)


##### Sex ratio #####


## Extract the sex ratio data for scaling fecundity
SurvTot$M.Ratio <- as.numeric(SurvTot$M.Ratio)
SurvTot$F.Ratio <- as.numeric(SurvTot$F.Ratio)

SurvTotSex <- SurvTot %>%
  filter(F.Ratio >= 0 | M.Ratio >= 0) %>%
  group_by(Generation, Treatment, Rep, Beak) %>%
  summarise(F.Ratio = last(F.Ratio))

SurvTotSex$Generation.c <- as.numeric(SurvTotSex$Generation)
SurvTotSex$Generation <- as.factor(as.numeric(SurvTotSex$Generation))

gam.sexratio <- gam(F.Ratio ~ s(Generation.c, by = Treatment, k = 3), data = SurvTotSex)
summary(gam.sexratio)


plot(gam.sexratio, pages = 1)
plot_smooth(gam.sexratio, view="Generation.c", plot_all="Treatment", rug=FALSE) ## SAVE THIS FIGURE AS 6 X 4.5 INCHES IN PORTRAIT MODE
gam.stats.sexratio <- summary(gam.sexratio)

gam.stats.sexratio

# an HTML table of the results
tab_model(gam.sexratio)


SurvTotSex <- unite(SurvTotSex,#the data frame
                    Treat.Rep, #the name of the new column
                    c(Treatment, Rep), #the existing columns you want to combine
                    remove = FALSE)

SurvTotSex$Treatment <- as.factor(SurvTotSex$Treatment)


sexratio.lmm <- lmer(F.Ratio ~ Generation * Treatment + (1|Treat.Rep),
                     REML = FALSE,
                     #family = gaussian,
                     data = SurvTotSex)

sexratio.emm <- emmeans(sexratio.lmm, pairwise ~ Generation | Treatment)



pairs(sexratio.emm)




survPlotSex <- ggplot(data = SurvTotSex, aes(x = Generation.c,
                              y = F.Ratio/0.5,
                              color = factor(Treatment)))+
  
  geom_point() + 
  
  geom_smooth(method = "lm", se = F) +
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  )
  )+
  theme_classic()+
  labs(y="Female:Male", x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))

survPlotSex

ggsave(filename = paste0(sex.directory,"Sex_ratio_plot.pdf"), plot = survPlotSex, height = 112, width = 180, units = "mm")



######################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################

######################### EPR DATA ###############################


EPRtot <- fread(paste(epr.directory,"Data_frames/EPR_HF_data_total.txt", sep = ""))


#create linear models that are tested against each other grouped by generation and analyzes EPR/HF by temp, pH, and temp*pH
EPRtot <- EPRtot[!is.na(EPRtot),]
EPRtot <- EPRtot[!is.na(HFtot),]


#EPRtot <- na.omit(EPRtot)
EPRtot$Generation.c <- as.numeric(EPRtot$Generation)
EPRtot$Generation <- as.factor(as.numeric(EPRtot$Generation))
EPRtot$Treatment <- as.factor(EPRtot$Treatment)


EPRtot <- unite(EPRtot,
                Gen.treat,
                c(Generation, Treatment),
                remove = FALSE)


EPRtot <- subset(EPRtot, Gen.treat != "10_HA")


if (EPRtot$pH == 8.13) {
  
  EPRtot$pH[EPRtot$pH == '8.13'] <- 8.20
  
}


eprStatsAll <- EPRtot %>% 
  group_by(Treatment, Generation.c) %>% 
  summarise(Rate = mean(EPRtot, na.rm = TRUE),
            sd = sd(EPRtot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = Rate - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = Rate + qt(1 - (0.05 / 2), n.count-1)*se)


hfStatsAll <- EPRtot %>% 
  group_by(Treatment, Generation.c) %>% 
  summarise(HF = mean(HFtot, na.rm = TRUE),
            sd = sd(HFtot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = HF - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = HF + qt(1 - (0.05 / 2), n.count-1)*se)




###combine all plots with HF on the same plot

colnames(hfStatsAll) <- c("Treatment2", "Generation2", "HF", "sd2", "n.count2", "se2", "lower.ci2", "upper.ci2")

eprSumComplete <- cbind(eprStatsAll, hfStatsAll)
eprSumComplete$Treatment <- as.factor(eprSumComplete$Treatment)


## have to make the HF values on the same scale as EPR values
eprSumComplete$HF <- eprSumComplete$HF*max(eprSumComplete$Rate)
eprSumComplete$lower.ci2 <- eprSumComplete$lower.ci2*max(eprSumComplete$Rate)
eprSumComplete$upper.ci2 <- eprSumComplete$upper.ci2*max(eprSumComplete$Rate)

#eprSumComplete$Treatment[eprSumComplete$Treatment=="AA"] <- expression("Control")

eprSumComplete$Rate[eprSumComplete$Rate==0] <- NA
eprSumComplete$HF[eprSumComplete$HF==0] <- NA

eprSumComplete <- eprSumComplete[complete.cases(eprSumComplete),]

eprSumAA <- subset(eprSumComplete, Treatment == "AA")
eprSumAH <- subset(eprSumComplete, Treatment == "AH")
eprSumHA <- subset(eprSumComplete, Treatment == "HA")
eprSumHH <- subset(eprSumComplete, Treatment == "HH")


#### separate bar graphs for each treatment

## AA
epr.bar.graphs.AA <- ggplot(data = eprSumAA, aes(Generation.c, Rate))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8,
           fill = "#0c7bdc")+
  geom_errorbar(aes(ymin=lower.ci, 
                    ymax=upper.ci), 
                size=1.3, 
                width=1.8,
                position = position_dodge())+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), 
       x="Generation")+
  geom_line(aes(y=HF), 
            size=1.3,
            colour="purple")+
  geom_errorbar(aes(ymin=lower.ci2, ymax=upper.ci2), 
                size=1.3,
                width = 1.8,
                position = position_dodge(0.3),
                colour = "purple")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(limits = c(0,49), breaks = c(0,5,10,15,20,25,30,35,40),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate),
                                         name = "Hatching Success",
                                         breaks = c(0.00,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1)))+
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  
  # statistical groups for EPR
  geom_text(x=0,y=2,label = "c", colour = "white", size = 12)+
  geom_text(x=3,y=2,label = "a", colour = "white", size = 12)+
  geom_text(x=6,y=2,label = "a", colour = "white", size = 12)+
  geom_text(x=10,y=2,label = "a", colour = "white", size = 12)+
  geom_text(x=12,y=2,label = "a", colour = "white", size = 12)+
  geom_text(x=15,y=2,label = "a", colour = "white", size = 12)+
  geom_text(x=25,y=2,label = "b", colour = "white", size = 12)+
  
  # statistical groups for HS
  geom_text(x=0.5,y=43,label = "b", colour = "purple", size = 12)+
  geom_text(x=3,y=42.5,label = "ab", colour = "purple", size = 12)+
  geom_text(x=6,y=36,label = "a", colour = "purple", size = 12)+
  geom_text(x=10,y=38,label = "ab", colour = "purple", size = 12)+
  geom_text(x=12,y=40.5,label = "ab", colour = "purple", size = 12)+
  geom_text(x=15,y=40.5,label = "ab", colour = "purple", size = 12)+
  geom_text(x=25,y=41,label = "ab", colour = "purple", size = 12)


epr.bar.graphs.AA

## AH
epr.bar.graphs.AH <- ggplot(data = eprSumAH, aes(Generation.c, Rate))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8,
           fill = "#009e73")+
  geom_errorbar(aes(ymin=lower.ci, 
                    ymax=upper.ci), 
                size=1.3, 
                width=1.8,
                position = position_dodge())+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), 
       x="Generation")+
  geom_line(aes(y=HF), 
            size=1.3,
            colour="purple")+
  geom_errorbar(aes(ymin=lower.ci2, ymax=upper.ci2), 
                size=1.3,
                width = 1.8,
                position = position_dodge(0.3),
                colour = "purple")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(limits = c(0,45), breaks = c(0,5,10,15,20,25,30,35,40),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate),
                                         name = "Hatching Success",
                                         breaks = c(0.00,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1)))+
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  
  # statistical groups for EPR
  geom_text(x=0,y=2,label = "c", colour = "white", size = 12)+
  geom_text(x=3,y=2,label = "bc", colour = "white", size = 12)+
  geom_text(x=6,y=2,label = "ab", colour = "white", size = 12)+
  geom_text(x=10,y=2,label = "a", colour = "white", size = 12)+
  geom_text(x=12,y=2,label = "ab", colour = "white", size = 12)+
  geom_text(x=15,y=2,label = "a", colour = "white", size = 12)+
  geom_text(x=25,y=2,label = "bc", colour = "white", size = 12)+
  
  # statistical groups for HS
  geom_text(x=0,y=41,label = "ab", colour = "purple", size = 12)+
  geom_text(x=3,y=37.5,label = "ab", colour = "purple", size = 12)+
  geom_text(x=6,y=34,label = "a", colour = "purple", size = 12)+
  geom_text(x=10,y=39,label = "ab", colour = "purple", size = 12)+
  geom_text(x=12,y=37,label = "ab", colour = "purple", size = 12)+
  geom_text(x=15,y=35.5,label = "ab", colour = "purple", size = 12)+
  geom_text(x=25,y=41,label = "b", colour = "purple", size = 12)


epr.bar.graphs.AH


## HA
epr.bar.graphs.HA <- ggplot(data = eprSumHA, aes(Generation.c, Rate))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8,
           fill = "#ffa500")+
  geom_errorbar(aes(ymin=lower.ci, 
                    ymax=upper.ci), 
                size=1.3, 
                width=1.8,
                position = position_dodge())+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), 
       x="Generation")+
  geom_line(aes(y=HF), 
            size=1.3,
            colour="purple")+
  geom_errorbar(aes(ymin=lower.ci2, ymax=upper.ci2), 
                size=1.3,
                width = 1.8,
                position = position_dodge(0.3),
                colour = "purple")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(limits = c(0,45), breaks = c(0,5,10,15,20,25,30,35,40),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate),
                                         name = "Hatching Success",
                                         breaks = c(0.00,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1)))+
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  
  # statistical groups for EPR
  geom_text(x=0,y=2,label = "ab", colour = "black", size = 12)+
  geom_text(x=3,y=2,label = "b", colour = "black", size = 12)+
  geom_text(x=6,y=2,label = "b", colour = "black", size = 12)+
  #geom_text(x=9,y=1.7,label = "ab", colour = "white", size = 12)+
  geom_text(x=12,y=2,label = "ab", colour = "black", size = 12)+
  geom_text(x=15,y=2,label = "a", colour = "black", size = 12)+
  geom_text(x=25,y=2,label = "a", colour = "black", size = 12)+
  
  # statistical groups for HS
  geom_text(x=0,y=38,label = "a", colour = "purple", size = 12)+
  geom_text(x=3,y=42,label = "ab", colour = "purple", size = 12)+
  geom_text(x=6,y=42,label = "ab", colour = "purple", size = 12)+
  #geom_text(x=9,y=42,label = "ab", colour = "purple", size = 12)+
  geom_text(x=12,y=39,label = "ab", colour = "purple", size = 12)+
  geom_text(x=15,y=43,label = "ab", colour = "purple", size = 12)+
  geom_text(x=25,y=44,label = "b", colour = "purple", size = 12)


epr.bar.graphs.HA

## HH
epr.bar.graphs.HH <- ggplot(data = eprSumHH, aes(Generation.c, Rate))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8,
           fill = "#d55e00")+
  geom_errorbar(aes(ymin=lower.ci, 
                    ymax=upper.ci), 
                size=1.3, 
                width=1.8,
                position = position_dodge())+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), 
       x="Generation")+
  geom_line(aes(y=HF), 
            size=1.3,
            colour="purple")+
  geom_errorbar(aes(ymin=lower.ci2, ymax=upper.ci2), 
                size=1.3,
                width = 1.8,
                position = position_dodge(0.3),
                colour = "purple")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(limits = c(0,45), breaks = c(0,5,10,15,20,25,30,35,40),
                     sec.axis = sec_axis(~./max(eprStatsAll$Rate),
                                         name = "Hatching Success",
                                         breaks = c(0.00,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1)))+
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  
  # statistical groups for EPR
  geom_text(x=0,y=2,label = "a", colour = "white", size = 12)+
  geom_text(x=3,y=2,label = "b", colour = "white", size = 12)+
  geom_text(x=6,y=2,label = "ab", colour = "white", size = 12)+
  geom_text(x=9,y=2,label = "b", colour = "white", size = 12)+
  geom_text(x=12,y=2,label = "b", colour = "white", size = 12)+
  geom_text(x=15,y=2,label = "b", colour = "white", size = 12)+
  geom_text(x=25,y=2,label = "ab", colour = "white", size = 12)+
  
  # statistical groups for HS
  geom_text(x=0,y=25,label = "a", colour = "purple", size = 12)+
  geom_text(x=3,y=43,label = "b", colour = "purple", size = 12)+
  geom_text(x=6,y=40,label = "b", colour = "purple", size = 12)+
  geom_text(x=9,y=42.5,label = "b", colour = "purple", size = 12)+
  geom_text(x=12,y=42.5,label = "b", colour = "purple", size = 12)+
  geom_text(x=15,y=42.5,label = "b", colour = "purple", size = 12)+
  geom_text(x=25,y=43.3,label = "b", colour = "purple", size = 12)


epr.bar.graphs.HH




epr.plot.complete <- ggarrange(epr.bar.graphs.AA+rremove("x.title"), epr.bar.graphs.AH+rremove("x.title"),
                               epr.bar.graphs.HA,epr.bar.graphs.HH,
                               ncol = 2, nrow = 2)

epr.plot.complete

ggsave(filename = paste0(epr.directory,"EPR_HS_plot_total.pdf"), plot = epr.plot.complete, height = 210, width = 180, units = "mm")

## EPR and HS boxplots for supplemental material
eprStatsAll$Generation <- as.factor(as.numeric(eprStatsAll$Generation.c))
eprStatsAll$Rate[eprStatsAll$Rate==0] <- NA
eprStatsAll <- eprStatsAll[complete.cases(eprStatsAll),]

eprboxplot <- ggplot(data = EPRtot, aes(x = Generation))+
  geom_boxplot(lwd = 1.1, aes(y = EPRtot, fill = factor(Treatment)),
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = eprStatsAll, aes(Generation, Rate, color=factor(Treatment)),
             size = 3, 
             shape = 18,
             position = position_dodge(width = 0.75)
  ) +
  
  
  scale_fill_manual(values = c("#0c7bdc", #AM
                               "#009e73", #OA
                               "#ffa500", #OW
                               "#d55e00" #OWA
  )
  )+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  )
  )+
  theme_classic()+
  labs(y=expression(atop(Egg~Production~Rate,(female^-1~day^-1))), x="Generation")+
  scale_x_discrete(breaks = c(0,3,6,9,10,12,15,25))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  #ylim(1.402, 1.405)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


eprboxplot

hfStatsAll$Generation <- as.factor(as.numeric(hfStatsAll$Generation2))

hfStatsAll$HF[hfStatsAll$HF==0] <- NA
hfStatsAll <- hfStatsAll[complete.cases(hfStatsAll),]

hfboxplot <- ggplot(data = EPRtot, aes(x = Generation))+
  geom_boxplot(lwd = 1.1, aes(y = HFtot, fill = factor(Treatment)),
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = hfStatsAll, aes(Generation, HF, color=factor(Treatment2)),
             size = 3, 
             shape = 18,
             position = position_dodge(width = 0.75)
  ) +
  
  
  scale_fill_manual(values = c("#0c7bdc", #AM
                               "#009e73", #OA
                               "#ffa500", #OW
                               "#d55e00" #OWA
  )
  )+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  )
  )+
  theme_classic()+
  labs(y="Hatching Frequency", x="Generation")+
  scale_x_discrete(breaks = c(0,3,6,9,10,12,15,25))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


hfboxplot

epr.hf.boxplots <- ggarrange(eprboxplot, hfboxplot, ncol = 1)

ggsave(filename = paste0(epr.directory,"EPR_HF_boxplots.pdf"), plot = epr.hf.boxplots, height = 210, width = 180, units = "mm")


##### Statistics for EPR data #####

# remove for where there is no data
EPRtot$Generation[EPRtot$Generation==10] <-9


EPRtot <- unite(EPRtot,
                Treat.Rep,
                c(Treatment, Rep),
                remove = FALSE)


## create gam models to test if epr and hf change across generations

# s() indicates a smooth function
gam1 <- gam(EPRtot ~ s(Generation.c, by = Treatment, k = 3), data = EPRtot)

epr.gam.stats <- summary(gam1)
epr.gam.stats

tab_model(gam1)

## generalized additive mixed models with replicates as mixed effects
mm <- gamm(EPRtot ~ s(Generation.c, by = Treatment, k = 3), random = list(Treat.Rep=~1), data = EPRtot)
summary(mm$lme)
summary(mm$gam)


## gam for HF
gam2 <- gam(HFtot ~ s(Generation.c, by = Treatment, k = 3), data = EPRtot)
hf.gam.stats <- summary(gam2)
hf.gam.stats

tab_model(gam2)

# plot the gam models

gam.list <- list(gam1, gam2)

lapply(gam.list, function (x) plot_smooth(x, view = "Generation.c", plot_all = "Treatment", rug = FALSE))






# create models for Tukey pairwise comparisons

l <- glmer(EPRtot~Generation.c*Treatment+(1|Treat.Rep), 
           REML = FALSE,
           family = gaussian,
           data = EPRtot)
l
summary(l)
Anova(l)
ranef(l)

plot_model(l, type = 'int')
tab_model(l)



l2 <- glmer(EPRtot~Generation*Treatment+(1|Treat.Rep), 
            REML = FALSE,
            family = gaussian,
            data = EPRtot)
tab_model(l2)

emmip(l2, Treatment ~ Generation)

epr.emm <- emmeans(l2, pairwise ~ Generation | Treatment)

pairs(epr.emm)


l3 <- glmer(HFtot~Generation*Treatment+(1|Treat.Rep), 
            REML = FALSE,
            family = gaussian,
            data = EPRtot)

hf.emm <- emmeans(l3, pairwise ~ Generation | Treatment)



# create an object with a dataframe of the estimated marginal means and the groups for each factor
epr.tukey.groups <- CLD(epr.emm)
hf.tukey.groups <- CLD(hf.emm)


fwrite(epr.tukey.groups, file = paste(epr.directory,"Statistics/EPR_tukey_groups.txt", sep = ""), sep = "\t")
fwrite(hf.tukey.groups, file = paste(epr.directory,"Statistics/HF_tukey_groups.txt", sep = ""), sep = "\t")





######################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################

######################### Lambda calculations ###############################

SurvData <- unite(SurvTot,#the data frame
                  unique, #the name of the new column
                  c(Generation, Treatment, Rep, Beak), #the existing columns you want to combine
                  remove = FALSE)

# change number of individuals to numeric format
SurvData$nx <- as.numeric(SurvData$nx)




SurvData$lx <- as.numeric(SurvData$lx)

## Extract the sex ratio data for scaling fecundity
SurvData$M.Ratio <- as.numeric(SurvData$M.Ratio)
SurvData$F.Ratio <- as.numeric(SurvData$F.Ratio)

SurvDataSex <- SurvData %>%
  filter(F.Ratio >= 0 | M.Ratio >= 0) %>%
  group_by(Generation, Treatment, Rep, Beak) %>%
  summarise(F.Ratio = last(F.Ratio))

SurvDataSex$Generation.c <- as.numeric(SurvDataSex$Generation)

SurvDataSex$Generation <- as.factor(as.numeric(SurvDataSex$Generation))



SurvDataSex.Mean <- SurvDataSex %>%
  group_by(Generation, Treatment, Rep) %>%
  summarise(F.Ratio = mean(F.Ratio, na.rm = TRUE)) %>%
  as.data.frame(SurvDataSex.Mean)

SurvDataSex.Mean

SurvDataSex.Mean$Generation <- as.numeric(SurvDataSex.Mean$Generation)



# Create a vector formatted for day-specific survivorship as it fits along the off-diagonal of a Leslie Matrix
# See Caswell, H. 2001. Matrix Population Models for further details


SurvData1 <- SurvData %>%
  group_by(Treatment, Rep, Beak) %>%
  mutate(lx_diag = if_else(time == 0, 1, as.numeric(lx/lag(lx, default = first(lx))))) %>% ## always start with 100%
  mutate(lx_diag = if_else(lx_diag <= 1.000000, lx_diag, lag(lx_diag, n=1, default =last(lx_diag)))) %>%
  mutate(days = if_else(time == 0, 1, as.numeric(time-lag(time)))) # create a new column for the number of days spent at the respective survivorships


# Survivorship is reflective of prior day. Therefore, there can be no 0 survivorship in this vector.
# If all animals die, then the vector ends and the matrix is truncated at the appropriate time
SurvData1 <- filter(SurvData1, lx_diag > 0) 


# Check if there is any super survivorship. There can be no survivorship > 1.
# No individuals can be lost and reappear
if (any(SurvData1$lx_diag > 1) == TRUE ){
  SurvData1$lx_diag <- if_else(SurvData1$lx_diag > 1, lag(SurvData1$lx_diag), SurvData1$lx_diag)
}

any(SurvData1$lx_diag > 1)

SurvData1 <- as.data.frame(SurvData1)


# elongate the data frame to make it reflect actual days spent over the experiment. This essentially changes the matrix to a daily matrix.

SurvData1 <- SurvData1[rep(seq(nrow(SurvData1)), SurvData1$days),] 


SurvData1$Generation <- as.numeric(SurvData1$Generation)


#################################################################################################################################################################################
##### Dev time data #####

Dev.time <- filter(SurvData, SurvData$Cdev > 0)
Dev.time <- Dev.time[rep(seq(nrow(Dev.time)), Dev.time$Cdev),]

Dev.time.sum <- Dev.time %>%
  group_by(Generation, Treatment, Rep, Beak) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            sd = sd(time, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

Dev.time.sum <- unite(Dev.time.sum,#the data frame
                      unique, #the name of the new column
                      c(Generation, Treatment, Rep, Beak), #the existing columns you want to combine
                      remove = FALSE)


Dev.time.sum$Generation <- as.numeric(Dev.time.sum$Generation)

Dev.time.sum.short <- Dev.time.sum[,-c(1,7:11)]

## Add the dev.time to the survival table
SurvData1 <- inner_join(SurvData1,
                        Dev.time.sum.short,
                        by = c("Generation", "Treatment", "Rep", "Beak"))

SurvData1 <- SurvData1 %>% rename(dev.time = mean) # the new name of the column has to come first

SurvData1 <- unite(SurvData1,#the data frame
                   unique2, #the name of the new column
                   c(Generation, Treatment, Rep), #the existing columns you want to combine
                   remove = FALSE) 


Survival.list <- split(SurvData1, f = SurvData1$unique2) # don't create the list with a list of organizers. Then you create a complete matrix of samples which include some tables with no data



# lists within lists
Survival.list <- lapply(Survival.list, function (x) split(x, f = x$Beak))


#################################################################################################################################################################################
##### EPR data #####

EPR.Data <- unite(EPRtot,#the data frame
                  unique, #the name of the new column
                  c(Generation, Treatment, Rep), #the existing columns you want to combine
                  remove = FALSE)

EPR.data$fecundity <- EPR.data$EPRtot*EPR.data$HFtot



EPR.data$Generation[EPR.data$Generation==10] <- 9

EPR.data$Treatment <- case_when(EPR.data$Treatment == "AA" ~ 1,
                                EPR.data$Treatment == "AH" ~ 2,
                                EPR.data$Treatment == "HA" ~ 3,
                                EPR.data$Treatment == "HH" ~ 4)

EPR.data <- inner_join(EPR.data, #the data frame you want to modify
                       
                       SurvDataSex.Mean, # the data frame that has the column you want to add
                       
                       by = c("Generation", "Treatment", "Rep"))



# Caluclate per capita sex-specific fecundity

EPR.data$sex.spec.fecundity <- EPR.data$fecundity*EPR.data$F.Ratio


which(is.na(EPR.data$sex.spec.fecundity))

# Create a list of lists to use for matrix creations.


### NOTE: lists of survivorship and EPR data MUST have same names and same length of items

EPR.list <- split(EPR.data, f = EPR.data$unique)

# no fecundity data available
Survival.list[["3_1_4"]] <- NULL

names(Survival.list)
names(EPR.list)
names.surv <- names(Survival.list)
names.surv <- as.list(names.surv)
names.surv

#EPR.list[[1]]

# Create a dummy table to add data to

lambda.results <- data.frame(Variables = "dummy", lambda = 0, surv = 0, epr = 0, hf = 0, sex = 0, dev.time = 0)


for (i in 1:length(Survival.list)) { #list of names for each gen and treatment
  
  diag.ls <- Survival.list[[i]]
  
  epr.df <- EPR.list[[i]] # select the appropriate data frame to pull epr values from
  
  for (j in 1: length(diag.ls)) {
    
    diag.df <-  diag.ls[[j]] # temporary list used to index at each point
    
    
    
    diag.v <- diag.df$lx_diag # extract the vector of interest for the diagonal
    
    surv.value <- min(diag.df$lx) # find the corresponding survival value for plotting the fitness landscape
    
    leslie.matrix <- diag(diag.v) # create the matrix without the top row added
    
    dev.time.value <- as.integer(mean(diag.df$dev.time))
    
    zero <- matrix(0, nrow = 1, ncol = (dev.time.value))
    
    fecundity.vector <- epr.df$sex.spec.fecundity # select the vector with sex spec fecundity
    
    epr.vector <- epr.df$EPRtot
    
    hf.vector <- epr.df$HFtot
    
    sex.vector <- epr.df$F.Ratio
    
    epr.count <- as.integer(dim(leslie.matrix)[1]-dev.time.value)
    
    if (epr.count < 1) {
      
      epr.count <- 1
      
    } 
    
    
    for(k in 1:length(epr.vector)) {
      
      fecundity.value <- fecundity.vector[k] # use the sex specific fecundity values in sequence with the same survival matrix
      
      
      if (is.na(fecundity.value) == TRUE) {
        
        fecundity.value <- 0 ## use an arrow when assigning numbers, not a "=="
        
      } 
      
      epr.value <- epr.vector[k]
      
      hf.value <- hf.vector[k]
      
      sex.value <- sex.vector[k]
      
      fecundity.row <- t(c(zero, rep(fecundity.value, epr.count))) # combine the zero row and the epr.value for the matrix and transpose it to make it a row
      
      if (ncol(fecundity.row) > ncol(leslie.matrix)) { # if the dev.time is somehow greater than the survival matrix, then only make the last column reflective of epr
        
        delete <- ncol(fecundity.row)-ncol(leslie.matrix) # find the number of days that the dev time is greater than the survivorship
        
        fecundity.row <- fecundity.row[,-c(1:delete)] # delete those days from the fecundity row to make it the same size as the leslie matrix
        
      }
      
      leslie.matrix1 <- rbind(fecundity.row, leslie.matrix) # add the fecundity row to the matrix
      
      matrix.final <- leslie.matrix1[-nrow(leslie.matrix1),] # eliminate the last row to make it a square
      
      eigen.calcs <- eigen.analysis(matrix.final, zero = FALSE) # calculate the eigen values
      
      lambda.value <- eigen.calcs$lambda1 # extract the dominant eigen values
      
      lambda.row <- data.frame(Variables = names.surv[i], 
                               lambda = lambda.value, 
                               surv = surv.value, 
                               epr = epr.value, 
                               hf = hf.value,
                               sex = sex.value,
                               dev.time = dev.time.value) # create a 1x2 data frame to add to the end of the final data frame
      
      # data frame has to have the same colnames in order to rbind
      
      colnames(lambda.row) <- colnames(lambda.results) # make sure the data frames have the same names
      
      lambda.results <- bind_rows(lambda.results, lambda.row) # append the data frame with new results
      
      
    }
    
    
    
  }
}


lambda.results <- separate(lambda.results, "Variables", into = c("Generation", "Treatment", "Rep"))


# remove the first row as a last step
lambda.results <- lambda.results[-1,] 

#setwd("C:/Users/james/Documents/Grad_school/OA_Project/Fitness/")



lambda.results <- lambda.results %>% 
  group_by(Generation, Treatment, Rep) %>% 
  mutate(lambda.stand = (lambda-mean(lambda))/sd(lambda),
         surv.stand = (surv-mean(surv))/sd(surv),
         epr.stand = (epr-mean(epr))/sd(epr),
         hf.stand = (hf-mean(hf))/sd(hf),
         sex.stand = (sex-mean(sex))/sd(sex),
         dev.stand = (dev.time-mean(dev.time))/sd(dev.time),
         lambda.rel = lambda/mean(lambda),
         surv.rel = surv/mean(surv),
         epr.rel = epr/mean(epr),
         hf.rel = hf/mean(hf),
         sex.rel = sex/mean(sex),
         dev.rel = dev.time/mean(dev.time))


lambda.results[is.na(lambda.results$lambda.rel)] <- 0




##### Selection coefficients #####

## check directional selection gradients
AA <- lm(lambda.rel ~ surv + epr + hf + sex, data = subset(lambda.results, Treatment == 1))

AA.beta <- lm.beta(AA)


summary(AA.0.beta)


AH <- lm(lambda.rel ~ surv + epr + hf + sex, data = subset(lambda.results, Treatment == 2))

AH.beta <- lm.beta(AH)

summary(AH.beta)


HA <- lm(lambda.rel ~ surv + epr + hf + sex, data = subset(lambda.results, Treatment == 3))


HA.beta <- lm.beta(HA)


summary(HA.beta)



HH <- lm(lambda.rel ~ surv + epr + hf + sex, data = subset(lambda.results, Treatment == 4))

HH.beta <- lm.beta(HH)


summary(HH.beta)





lambda.results <- unite(lambda.results,
                        unique,
                        c(Generation, Treatment),
                        remove = F)




lambda.list <- split(lambda.results, f = lambda.results$unique)

lm.list.2 <- lapply(lambda.list, function (x) lm(lambda.rel ~ surv + epr + hf + sex, data = x))


beta.list <- lapply(lm.list.2, function (x) lm.beta(x))


summary(beta.list$`0_1`)

summary(beta.list$`0_2`)

summary(beta.list$`0_3`)

summary(beta.list$`0_4`)


summary(beta.list$`25_1`)

summary(beta.list$`25_2`)

summary(beta.list$`25_3`)

summary(beta.list$`25_4`)

## see if there are differences in selection gradient for each trait at gen's 0 and 25 for each treatment
## no difference in survival, but yes difference in epr and hf for all treatments
Anova(beta.list$`0_1`, beta.list$`25_1`)

Anova(beta.list$`0_2`, beta.list$`25_2`)

Anova(beta.list$`0_3`, beta.list$`25_3`)

Anova(beta.list$`0_4`, beta.list$`25_4`)




library(lme4)

lm.list <- lmList(lambda.rel ~ surv + epr + hf + sex | unique, data = lambda.results)

summary(lm.list$`0_1`)





mult.reg <- lm(lambda.rel~surv*epr*hf*dev.time*sex, data = lambda.results)

summary(mult.reg)
tab_model(mult.reg) # right click on the table to select all, then copy, then paste into word to keep formatting


mult.reg$coefficients

mult.reg.anova <- Anova(mult.reg)
mult.reg.anova$factors <- rownames(mult.reg.anova)

fwrite(mult.reg.anova, file = paste(fit.directory,"Statistics/lambda_mult_regression.txt", sep = ""), sep = "\t")



mult.reg.2 <- lm(lambda~surv+epr+hf+dev.time+sex, data = lambda.results)
coef(mult.reg.2)
tab_model(mult.reg.2)

# relative coefficients for the model. Must use non-standardized data
library(lm.beta)

mult.reg.2.beta <- lm.beta(mult.reg.2)
summary(mult.reg.2.beta)

tab_model(mult.reg.2.beta)

mult.reg.coef <- as.data.frame(coef(mult.reg.2.beta))

mult.reg.coef$factors <- rownames(mult.reg.coef)

fwrite(mult.reg.coef, file = paste(fit.directory,"Statistics/lambda_mult_regression_standardized_coefficients.txt", sep = ""), sep = "\t")

plot_model(mult.reg.2)

## find the coefficients for the models in the landscape figures

names(lambda.list)

lambda.coef <- data.frame(Variables = "dummy", surv = 0, epr = 0, hf = 0)


#i=1
for(i in 1:length(lambda.list)) {
  
  dat <- lambda.list[[i]]
  
  one <- first(dat$unique)
  
  
  
  surv.lm <- lm(lambda.rel ~ surv, data = dat)
  
  surv.coef <- surv.lm$coefficients[2]
  
  
  
  
  epr.lm <- lm(lambda.rel ~ epr, data = dat)
  
  epr.coef <- epr.lm$coefficients[2]
  
  
  
  
  hf.lm <- lm(lambda.rel ~ hf, data = dat)
  
  hf.coef <- hf.lm$coefficients[2]
  
  
  
  
  
  
  
  new.coefs <- data.frame(Variables = one, surv = surv.coef, epr = epr.coef, hf = hf.coef)
  
  lambda.coef <- bind_rows(lambda.coef, new.coefs)
  
  
}


lambda.coef <- separate(lambda.coef, "Variables", into = c("Generation", "Treatment"))


# remove the first row as a last step
lambda.coef <- lambda.coef[-1,] 



## Test to see if models for epr and surv are different from each other between F0 and F25

## AA F0
a <- lambda.results %>% 
  filter(Treatment == 1 & Generation == 0)

a.surv <- lm(lambda.rel~surv, data = a)

a.epr <- lm(lambda.rel~epr, data = a)

summary(a.epr)

summary(a.surv)

## AA F25
b <- lambda.results %>% 
  filter(Treatment == 1 & Generation == 25)

b.surv <- lm(lambda.rel~surv, data = b)

b.epr <- lm(lambda.rel~epr, data = b)

summary(b.surv)

summary(b.epr)


Anova(a.surv, b.surv)
Anova(a.epr, b.epr)


## AH F0
c <- lambda.results %>% 
  filter(Treatment == 2 & Generation == 0)


c.surv <- lm(lambda.rel ~ surv, data = c)

c.epr <- lm(lambda.rel ~ epr, data = c)

summary(c.surv)

summary(c.epr)


## AH F25

d <- lambda.results %>% 
  filter(Treatment == 2 & Generation == 25)


d.surv <- lm(lambda.rel ~ surv, data = d)

d.epr <- lm(lambda.rel ~ epr, data = d)


summary(d.surv)


summary(d.epr)


Anova(c.surv, d.surv)
Anova(c.epr,d.epr)



## HA F0


e <- lambda.results %>% 
  filter(Treatment == 3 & Generation == 0)


e.surv <- lm(lambda.rel ~ surv, data = e)

e.epr <- lm(lambda.rel ~ epr, data = e)


summary(e.surv)

summary(e.epr)

## HA F25

f <- lambda.results %>% 
  filter(Treatment == 3 & Generation == 25)


f.surv <- lm(lambda.rel ~ surv, data = f)

f.epr <- lm(lambda.rel ~ epr, data = f)


summary(f.surv)

summary(f.epr)


Anova(e.surv, f.surv)
Anova(e.epr, f.epr)

## HH F0


g <- lambda.results %>% 
  filter(Treatment == 4 & Generation == 0)


g.surv <- lm(lambda.rel ~ surv, data = g)

g.epr <- lm(lambda.rel ~ epr, data = g)

g.sex <- lm(lambda.rel ~ sex, data = g)

summary(g.surv)

summary(g.epr)

summary(g.sex)

## HH F25


h <- lambda.results %>% 
  filter(Treatment == 4 & Generation == 25)


h.surv <- lm(lambda.rel ~ surv, data = h)

h.epr <- lm(lambda.rel ~ epr, data = h)

h.sex <- lm(lambda.rel ~ sex, data = h)


summary(h.surv)

summary(h.epr)

summary(h.sex)

Anova(g.surv, h.surv)
Anova(g.epr, h.epr)
Anova(g.sex, h.sex)



## multiple regression analyses

f1 <- lm(lambda.rel~surv*factor(Treatment)*factor(Generation), data = lambda.results)
summary(f1)

f1.emm <- emmeans(f1,  ~ Generation | Treatment)
f1.emm




f1.coef <- tidy(coef(f1))
View(f1.coef)


f2 <- lm(lambda.stand~epr*factor(Treatment)*factor(Generation), data = lambda.results)
f2.emm <- emmeans(f2,  ~ Generation | Treatment)

f2.emm

f2.coef <- tidy(coef(f2))
View(f2.coef)




##### PATH analysis ######

# help from: https://rpubs.com/tbihansk/302732

mod <- 'lambda.rel ~ surv + epr + hf + dev.time + sex'

path.fit <- cfa(mod, data = lambda.results)

summary(path.fit)


library(mnormt)
library(lavaan)
library(semPlot)

cisemPaths(path.fit, 'std', layout = "circle", title = T)
semPaths(path.fit, 'std', layout = 'tree')
semPaths(path.fit, 'std', layout = 'spring')
semPaths(path.fit, 'std', layout = 'tree2')
semPaths(path.fit, 'std', layout = 'circle2')


## split by gen and treatment
mod2 <- 'lambda.rel ~ surv + epr + hf + sex'


path.list <- cfaList(mod2, lambda.list)


coef(path.list)



path.list2 <- lapply(lambda.list, function (x) cfa(mod2, x)) ## produces the same results as using SEM

summary(path.list2$`0_1`)
summary(path.list2$`0_2`)
summary(path.list2$`0_3`)
summary(path.list2$`0_4`)
summary(path.list2$`25_1`)
summary(path.list2$`25_2`)
summary(path.list2$`25_3`)
summary(path.list2$`25_4`)

path.list3 <- list(path.list2$`0_1`,
                   path.list2$`0_2`,
                   path.list2$`0_3`,
                   path.list2$`0_4`,
                   path.list2$`25_1`,
                   path.list2$`25_2`,
                   path.list2$`25_3`,
                   path.list2$`25_4`)

sem.figs <- lapply(path.list3, function (x) semPaths(x, 'std', layout = 'circle', title = T))


##### Statistics #####

lambda.results <- unite(lambda.results,#the data frame
                        unique, #the name of the new column
                        c(Generation, Treatment), #the existing columns you want to combine
                        remove = FALSE)

lambda.results <- unite(lambda.results,
                        Treat.Rep,
                        c(Treatment, Rep),
                        remove = FALSE)

lambda.results$Rep.c <- case_when(lambda.results$Treat.Rep == "1_1" ~ 1,
                                  lambda.results$Treat.Rep == "1_2" ~ 2,
                                  lambda.results$Treat.Rep == "1_3" ~ 3,
                                  lambda.results$Treat.Rep == "1_4" ~ 4,
                                  lambda.results$Treat.Rep == "2_1" ~ 5,
                                  lambda.results$Treat.Rep == "2_2" ~ 6,
                                  lambda.results$Treat.Rep == "2_3" ~ 7,
                                  lambda.results$Treat.Rep == "2_4" ~ 8,
                                  lambda.results$Treat.Rep == "3_1" ~ 9,
                                  lambda.results$Treat.Rep == "3_2" ~ 10,
                                  lambda.results$Treat.Rep == "3_3" ~ 11,
                                  lambda.results$Treat.Rep == "3_4" ~ 12,
                                  lambda.results$Treat.Rep == "4_1" ~ 13,
                                  lambda.results$Treat.Rep == "4_2" ~ 14,
                                  lambda.results$Treat.Rep == "4_3" ~ 15,
                                  lambda.results$Treat.Rep == "4_4" ~ 16)





# remove the data for gen 9 and treatment 3 since we do not have data for this generation
lambda.results <- lambda.results %>% 
  filter(unique != "9_3")


lambda.results <- fread(paste(fit.directory,"lambda_results_devtime_surv_epr_hf_sex_standardized_relative.txt", sep = ""))

lambda.results <- as.data.frame(lambda.results)


# create continuous generation vector for anova and plotting
lambda.results$Generation.c <- as.numeric(as.character(lambda.results$Generation)) 
lambda.results$Generation <- as.factor(as.numeric(lambda.results$Generation))
lambda.results$Treatment <- as.factor(lambda.results$Treatment)
lambda.results$Rep.c <- as.numeric(lambda.results$Rep)
lambda.results$Rep <- as.factor(as.numeric(lambda.results$Rep.c))

# check to see what the distribution looks like
hist(lambda.results$lambda, freq = F)

lambda.ls <- split(lambda.results, f = lambda.results$Treatment)

lambda.hist <- lapply(lambda.ls, function (x) hist(x$lambda))

AM.hist <- lambda.hist$`1`
OA.hist <- lambda.hist$`2`
OW.hist <- lambda.hist$`3`
OWA.hist <- lambda.hist$`4`

# save each as 3.5 in x 4.4 in when selecting portrait


# s() indicates a smooth function
gam1 <- gam(lambda ~ s(Generation.c, by = factor(Treatment), k = 3), data = lambda.results)
summary(gam1)
#gam2 <- gam(lambda ~ s(Generation.c, by = Treatment, k = 4), data = lambda.results)
#gam3 <- gam(lambda ~ s(Generation.c, by = Treatment, k = 5), data = lambda.results)
#gam4 <- gam(lambda ~ s(Generation.c, by = Treatment, k = 6), data = lambda.results)
#gam5 <- gam(lambda ~ s(Generation.c, by = Treatment, k = 7), data = lambda.results)

AIC(gam1,gam2,gam3,gam4,gam5)

## try BIC also because that penalizes the more complex models more



plot(gam1, pages = 1)
plot_smooth(gam1, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
gam.stats <- summary(gam1)

gam.stats

# an HTML table of the results
tab_model(gam1)


##### two-part mixed effects model with zero-inflated #####

## use to identify how inflated the zeroes make the data



fit_zigauss <- glmmTMB(lambda~Generation.c*Treatment+(1|Rep.c),
                       data = lambda.results,
                       ziformula = ~.) # this specifies the zero-inflated part of the model

## fixed effects results correspond to when response is >0, and zero-inflation results correspond to when lambda includes 0
fit_zigauss

summary(fit_zigauss)


## create a formatted output of model results with intra-class correlation (ICC)
tab_model(fit_zigauss)


fixef(fit_zigauss)$zi # the fixed-effects results of the zero-inflated model
ranef(fit_zigauss)$zi # the random effects intercepts of the model


## create linear mixed effects models that include zeroes and omit zeroes

## create the predicted values graph
lambda.results <- lambda.results %>% 
  mutate(lambda_zero = if_else(lambda==0, 0, 1))



# model for predicting when lambda is either 0 or >0
gm1 <- glmer(lambda_zero ~ Generation.c*Treatment+(1|Rep.c),
             data = lambda.results, family = binomial)

plot_model(gm1,type='int') + theme_sjplot2()

# model for predicting when lambda is not 0
lambda.nonzero = lambda.results[lambda.results$lambda>0,]

aa = lmer(lambda~Generation.c*Treatment+(1|Rep.c),
          data=lambda.nonzero)


plot_model(aa,type='int') + theme_sjplot2()


################################################################################################################################################
##### Plots of lambda and malthusian parameter #####
#lambda.results$Generation <- as.numeric(lambda.results$Generation)

lambda.mean.c <- lambda.results %>%
  group_by(Generation.c, Treatment) %>%
  summarise(mean = mean(lambda, na.rm = TRUE),
            sd = sd(lambda, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)


lambdaPlotTotal <- ggplot(data = lambda.mean.c, aes(Generation.c, mean, color=factor(Treatment)))+
  geom_line(size=1,
            position = position_dodge(width = 2))+
  geom_point(size = 2, position = position_dodge(width = 2))+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width=0, position = position_dodge(width = 2))+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#0c7bdc", #AM
                                  "#009e73", #OA
                                  "#ffa500", #OW
                                  "#d55e00" #OWA
  ))+
  theme_classic()+
  labs(y=expression(Population~Fitness~(lambda)), x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(limits = c(0.38,1.3),
                     breaks = c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))

lambdaPlotTotal

ggsave(filename = paste0(fit.directory,"Lambda_plot.pdf"), plot = lambdaPlotTotal, height = 101, width = 180, units = "mm")


lambda.mean.c$Generation <- as.factor(as.numeric(lambda.mean.c$Generation.c))
lambdaBoxplot <- ggplot()+
  geom_boxplot(data = lambda.results, aes(x = Generation, y = lambda, fill = factor(Treatment)), lwd = 1.1,
               alpha = 0.3,
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  geom_point(data = lambda.mean.c, aes(Generation, mean, color=factor(Treatment)),
             size = 3, 
             shape = 18,
             position = position_dodge(width = 0.75)
  ) +
  
  
  
  
  scale_fill_manual(values = c("#0c7bdc", #AM
                               "#009e73", #OA
                               "#ffa500", #OW
                               "#d55e00" #OWA
  ))+
  scale_color_manual(values = c("#0c7bdc", #AM
                                "#009e73", #OA
                                "#ffa500", #OW
                                "#d55e00" #OWA
  )
  )+
  theme_classic()+
  labs(y=expression(Population~Fitness~(lambda)), x="Generation")+
  scale_x_discrete(breaks = c(0,3,6,9,12,15,25))+
  #scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12))


lambdaBoxplot

ggsave(filename = paste0(fit.directory,"Lambda_boxplot.pdf"), plot = lambdaBoxplot, height = 101, width = 180, units = "mm")
##### Fitness landscape plots #####



lambda.results.0.full <- filter(lambda.results, Generation == 0)
last.gen <- max(lambda.results$Generation.c)
lambda.results.last.full <- filter(lambda.results, Generation == last.gen)
lambda.results.0.last.full <- rbind(lambda.results.0.full, lambda.results.last.full)



lambda.results.0.last.full$Treatment2 <- case_when(lambda.results.0.last.full$Treatment == 1 ~ "AM",
                                                   lambda.results.0.last.full$Treatment == 2 ~ "OA",
                                                   lambda.results.0.last.full$Treatment == 3 ~ "OW",
                                                   lambda.results.0.last.full$Treatment == 4 ~ "OWA")


# change the order of the factor to be in the order you want. Help from: https://stackoverflow.com/questions/5490638/how-to-change-the-order-of-facet-labels-in-ggplot-custom-facet-wrap-labels
lambda.results.0.last.full <- within(lambda.results.0.last.full, 
                                     Treatment2 <- factor(Treatment2, 
                                                          
                                                          # put the specific levels in the order you want
                                                          
                                                          levels = c("AM", # 1
                                                                     "OA", # 2
                                                                     "OW", # 3
                                                                     "OWA") # 4 
                                                          
                                     ))



### Survival fitness landscape

## Panels are wrapped by Treatment
x.pos <- 0.25

surv.landscape <- ggplot()+
  geom_density(data = SurvTot.agg.0.last, aes(x=lx, 
                                              y=after_stat(count),
                                              group = Generation,
                                              fill = Generation), alpha = 0.7, binwidth = 0.01)+
  
  #geom_freqpoly(data = SurvTot.agg.0.25, aes(x=lx,
  #                                     y=after_stat(count),
  #                                    group = factor(Generation),
  #                                   fill = factor(Generation)), binwidth = 0.10)+
  #geom_point(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 0.25
  #                                       y = mean*3, 
  #                                      color = factor(Generation)), 
  #        alpha = 0.7,
  #       size = 7,
  #      position = position_dodge(width = 0.12))+


#geom_errorbar(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 0.25
#                                          ymin=lower.ci*3,
#                                         ymax=upper.ci*3,
#                                        color = factor(Generation)),
#          position = position_dodge(width = 0.12),
#         size = 1.0,
#        width = 0.07,
#       alpha = 0.7)+

## for a fitness landscape (i.e. how fitness changes with changing phenotype distribution, use this instead of point and errorbar)

geom_smooth(data = lambda.results.0.last.full, aes(x = surv, y = lambda.rel*10, color = factor(Generation), group = factor(Generation)),
            method = "lm"
)+
  
  #geom_point(data = lambda.results.0.last.full, aes(x = surv, y = lambda.stand*10, color = factor(Generation), group = factor(Generation)))+
  
  #geom_vline(xintercept = mean(lambda.results.0.full$surv), color = "black", size = 1.5)+ # add vertical lines where the means are
  
  #geom_vline(xintercept = mean(lambda.results.25.full$surv), color = "gray", size = 1.5) +
  
  
  theme_minimal()+
  
  scale_fill_manual(values = c("black", "red"),
                    
                    labels = c("F0", 
                               bquote(F*.(last.gen))))+ ## reference an object in the legend. Help from: https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
  
  scale_color_manual(guide = "none",
                     values = c("black", "red"))+
  
  
  scale_x_continuous(name = "Survival")+
  
  scale_y_continuous(sec.axis = dup_axis(trans = ~./10, 
                                         name = expression(Relative~Fitness)))+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        strip.text.x = element_text(size = 9),
        legend.position = "none",
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  
  facet_wrap(~Treatment2)


surv.landscape

ggsave(filename = paste0(surv.directory,"Surv_fit_landscape.pdf"), plot = surv.landscape, height = 70, width = 180, units = "mm")
### SAVE AS 6.6 X 16 ############################################################################################


##### EPR fitness landscape #####


EPRtot.0 <- filter(EPRtot, Generation == 0)



EPRtot.last <- filter(EPRtot, Generation == last.gen)

EPRtot.0.last <- rbind(EPRtot.0, EPRtot.last)

EPRtot.0.last <- EPRtot.0.last %>% 
  mutate(Treatment2 <- case_when(Treatment == "AA" ~ "AM",
                                 Treatment == "AH" ~ "OA",
                                 Treatment == "HA" ~ "OW",
                                 Treatment == "HH" ~ "OWA"))

#EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="AA"] <- 1
#EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="AH"] <- 2
#EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="HA"] <- 3
#EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="HH"] <- 4



epr.landscape <- ggplot()+
  geom_density(data = lambda.results.0.last.full, aes(x=epr, 
                                                      y=after_stat(count),
                                                      group = factor(Generation), 
                                                      
                                                      fill = factor(Generation)), alpha = 0.7, binwidth = 0.01)+
  
  #geom_freqpoly(data = EPRtot.0.25, aes(x=EPRtot,
  #                                     y=after_stat(count),
  #                                    group = factor(Generation),
  #                                   fill = factor(Generation)), alpha = 0.7, binwidth = 10)+
  
  ## add the points for where the mean epr value shows up
  
  #geom_point(data = EPRtot.0.25.mean, aes(x = Rate,
  #                                   y = y.pos,
  #                                  color = factor(Generation)),
#         position = position_dodge(width = 15),
#        size = 7) +


#geom_point(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 62.5

#                                       y = mean/scale, # keep this as divided when just using the points and no landscape with scale = 25 for density or 0.6666665 for count

#                                      color = factor(Generation)),

#        position = position_dodge(width = 15),
#       alpha = 0.7,
#      size = 7)+

#geom_errorbar(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 62.5

#                                         ymin=lower.ci/scale,
#                                         ymax=upper.ci/scale, # keep this as divided when just using the points and no landscape with scale = 25 for density or 0.6666665 for count

#                                       color = factor(Generation)),
#          position = position_dodge(width = 15),
#         size = 1.0,
#        width = 8)+

geom_smooth(data = lambda.results.0.last.full, aes(x = epr, y = lambda.rel, color = factor(Generation), group = factor(Generation)),
            method = "lm",# make sure we know that there is only 2 possible data values
            se = TRUE
)+
  
  #geom_point(data = lambda.results.0.last.full, aes(x=epr, y = lambda.rel, color = factor(Generation), group = factor(Generation)))+
  
  
  
  theme_minimal()+
  scale_fill_manual(values = c("black", "red"),
                    
                    labels = c("F0", 
                               bquote(F*.(last.gen))) ## reference an object in the legend. Help from: https://stackoverflow.com/questions/15074127/use-expression-with-a-variable-r
  )+
  
  scale_color_manual(guide = "none",
                     values = c("black", "red"))+
  
  
  scale_x_continuous(name = expression(Egg~Production~Rate~(female^-1~day^-1)))+
  
  scale_y_continuous(labels = function (x) x*10,
                     sec.axis = dup_axis(trans = ~./10, # keep this as multiplied when just using the points and no landscape with scale = 25 for density or 0.6666665 for count
                                         
                                         name = expression(Relative~Fitness)))+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_blank(),
        legend.position = "none",
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.text = element_text(size = 9))+
  facet_wrap(~Treatment2)

epr.landscape

ggsave(filename = paste0(epr.directory,"EPR_fit_landscape.pdf"), plot = epr.landscape, height = 70, width = 180, units = "mm")
##### HF fitness landscape #####


hf.landscape <- ggplot()+
  geom_density(data = lambda.results.0.last.full, aes(x=hf, 
                                                      y=after_stat(count)/25,
                                                      group = factor(Generation), 
                                                      
                                                      fill = factor(Generation)), 
               alpha = 0.7, 
               binwidth = 0.01)+
  
  ## add the points for where the mean epr value shows up
  
  #geom_point(data = EPRtot.0.25.mean, aes(x = Rate,
  #                                   y = y.pos,
  #                                  color = factor(Generation)),
  #         position = position_dodge(width = 15),
  #        size = 7) +
  
  
  #geom_point(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 62.5

#                                       y = mean/scale, # keep this as divided when just using the points and no landscape with scale = 25

#                                      color = factor(Generation)),

#        position = position_dodge(width = 15),
#       alpha = 0.7,
#      size = 7)+

#geom_errorbar(data = lambda.mean.0.25, aes(x=x.pos, # x.pos needs to be 62.5

#                                         ymin=lower.ci/scale,
#                                         ymax=upper.ci/scale, # keep this as divided when just using the points and no landscape with scale = 25

#                                       color = factor(Generation)),
#          position = position_dodge(width = 15),
#         size = 1.0,
#        width = 8)+

geom_smooth(data = lambda.results.0.last.full, aes(x = hf, y = lambda.rel*20, color = factor(Generation), group = factor(Generation)),
            method = "lm"
)+
  
  
  
  theme_minimal()+
  scale_fill_manual(values = c("black", "red"),
                    
                    labels = c("F0", 
                               "F25"),
  )+
  
  scale_color_manual(guide = "none",
                     values = c("black", "red"))+
  
  
  scale_x_continuous(name = "Hatching Frequency")+
  
  
  
  
  scale_y_continuous(limits = c(0,40),
                     
                     name = "count",
                     
                     sec.axis = dup_axis(trans = ~./20, # keep this as multiplied when just using the points and no landscape with scale = 25
                                         
                                         name = "Relative Fitness"))+
  
  
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_blank(),
        legend.position = "none",
        legend.text.align = 0,
        legend.text = element_text(size = 9),
        legend.title = element_blank())+
  facet_wrap(~Treatment2)#+ylim(0,1000)

hf.landscape



### SAVE THE HF FIGURE AS 6.6 x 16 INCHES #########################################################################################################

fit.landscapes <- ggarrange(hf.landscape,epr.landscape,surv.landscape, common.legend = T, ncol = 1)

ggsave(filename = paste0(fit.directory,"Fitness_landscapes_complete.pdf"), plot = fit.landscapes, height = 210, width = 180, units = "mm")
