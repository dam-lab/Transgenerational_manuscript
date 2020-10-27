
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(broom)
library(survminer)
library(survival)
library(survMisc)
library(sjPlot)
library(mgcv)
library(lme4)


wd <- "C:/Users/james/Documents/Grad_school/OA_Project/Survival/SurvDataFiles/"

#wd <- "C:/Users/james/Documents/Grad_school/OA_hudsonica/Survival/SurvDataFiles/"


if (wd == "C:/Users/james/Documents/Grad_school/OA_Project/Survival/SurvDataFiles/") {
  
  
  surv.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Survival/"
  dev.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Development_time/"
  fit.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Fitness/"
  
  
  
   } else if (wd == "C:/Users/james/Documents/Grad_school/OA_hudsonica/Survival/SurvDataFiles/") {
    
    
    surv.directory <- "C:/Users/james/Documents/Grad_school/OA_hudsonica/Survival/"
    dev.directory <- "C:/Users/james/Documents/Grad_school/OA_hudsonica/Development_time/"
    fit.directory <- "C:/Users/james/Documents/Grad_school/OA_hudsonica/Fitness/"
    
  
}


#Take all the lists available and combine them all as you read them in
#received help from: https://stackoverflow.com/questions/38312864/adding-new-column-and-its-value-based-on-file-name?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa

#LF = list.files(pattern = ".*080718.*")
#LF

setwd(wd)


#LF = list.files(pattern = "txt$")
#LF

#SurvTot <- rbindlist(lapply(
#  setNames(LF,LF), #what you want to be "lapplied"
#  fread), #what you want the lapply to do
#  idcol = "source", #what the new column will be named (i.e. the header of the new column)
#  fill = TRUE) #fills missing columns with NA's



#SurvTot <- separate(SurvTot, "source", into = c("Date", 
#                                                "x", 
#                                                "y", 
#                                                "Generation"))
#View(SurvTot)
#fwrite(SurvTot, file = paste(surv.directory, "SurvDataFiles/Survival_data_total.txt", sep = ""), sep = "\t")

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
#SurvTot1 <- filter(SurvTot1, SurvTot1$Generation < 13)

library(car)


SurvTot.agg <- aggregate(SurvTot1$lx, 
                         list(SurvTot1$Treatment, SurvTot1$Rep, SurvTot1$Beak, SurvTot1$Generation), #use a list so you don't have to use "interaction"
                         min)#use the minimum value for each of the individual beakers which should be the final survival

colnames(SurvTot.agg) <- c("Treatment", "Rep", "Beak", "Generation", "lx")

is.even <- function(x) x %% 2 == 0



if (wd == "C:/Users/james/Documents/Grad_school/OA_hudsonica/Survival/SurvDataFiles/") {
  
  SurvTot.agg <- rbind(SurvTot.agg, SurvData.agg)
  
}

SurvTot.agg <- filter(SurvTot.agg, SurvTot.agg$Treatment < 5)

if (wd == "C:/Users/james/Documents/Grad_school/OA_Project/Survival/SurvDataFiles/") {

SurvTot.agg$Temp <- if_else(SurvTot.agg$Treatment < 3, 18, 22)

SurvTot.agg$pH <- if_else(is.even(SurvTot.agg$Treatment), 7.5, 8.2)

} else if (wd == "C:/Users/james/Documents/Grad_school/OA_hudsonica/Survival/SurvDataFiles/") {
  
  SurvTot.agg$Temp <- if_else(SurvTot.agg$Treatment < 3, 13, 15)
  
  SurvTot.agg$pH <- if_else(is.even(SurvTot.agg$Treatment), 7.8, 8.2)
  
  
}

SurvTot.agg$Generation.c <- as.numeric(SurvTot.agg$Generation)

SurvTot.agg$Generation <- as.factor(as.numeric(SurvTot.agg$Generation))

SurvTot.agg$Treatment <- as.factor(as.numeric(SurvTot.agg$Treatment))


SurvTot.agg <- unite(SurvTot.agg,
                     Treat.Rep,
                     c(Treatment, Rep),
                     remove = FALSE)

#SurvTot.agg$line <- if_else(SurvTot.agg$Treatment == 1 | SurvTot.agg$Treatment == 6, 
#                                          "AA", "HH")


SurvTot.mean <- SurvTot.agg %>%
  group_by(Treatment, Generation.c) %>%
  summarise(mean = mean(lx, na.rm = TRUE),
            sd = sd(lx, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)



SurvTot.mean[is.na(SurvTot.mean)] <- 0

survPlotTotal <- ggplot(data = SurvTot.mean, aes(Generation.c, mean, color=factor(Treatment)))+
  geom_line(size=1,
            position = position_dodge(width = 2))+
  geom_point(size = 2, position = position_dodge(width = 2))+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width = 0,
                position = position_dodge(width = 2))+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Survival", x="Generation")+
  scale_x_continuous(breaks = c(0,2,4,11))+
  #scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

survPlotTotal


## create density plots for survival for each treatment
## help from: https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html



library(ggridges)

#lapply(SurvTot.agg.ls, function (x) 
ggplot(data = SurvTot.agg, aes(x=lx, y=Generation, fill = Generation)) +
  geom_density_ridges() +
  theme_classic()+
  facet_wrap(~Treatment)
#facet_wrap(~Generation)





SurvTot.agg.0 <- filter(SurvTot.agg, Generation == 0)

last.gen <- max(SurvTot.agg$Generation.c)

SurvTot.agg.last <- filter(SurvTot.agg, Generation == last.gen)
SurvTot.agg.0.last <- rbind(SurvTot.agg.0, SurvTot.agg.last)

SurvTot.agg.0.last$Treatment2 <- case_when(SurvTot.agg.0.last$Treatment == 1 ~ "Ambient",
                                 SurvTot.agg.0.last$Treatment == 2 ~ "Acidification",
                                 SurvTot.agg.0.last$Treatment == 3 ~ "Warming",
                                 SurvTot.agg.0.last$Treatment == 4 ~ "Greenhouse")

# change the order of the factor to be in the order you want. Help from: https://stackoverflow.com/questions/5490638/how-to-change-the-order-of-facet-labels-in-ggplot-custom-facet-wrap-labels
SurvTot.agg.0.last <- within(SurvTot.agg.0.last, 
                             Treatment2 <- factor(Treatment2, 
                                                  
                                                  # put the specific levels in the order you want
                                                  
                                                  levels = c("Ambient", # 1
                                                             "Acidification", # 2
                                                             "Warming", # 3
                                                             "Greenhouse") # 4 
                                                  
                                                  ))


SurvTot.agg.0.last %>% 
  group_by(Treatment, Generation) %>% 
  summarise(n=n())

###### Fitness landscape plots #####

lambda.results <- fread(paste(fit.directory,"lambda_results_devtime_surv_epr_hf_sex_standardized_relative.txt", sep = ""))
lambda.results.0.full <- filter(lambda.results, Generation == 0)
lambda.results.last.full <- filter(lambda.results, Generation == last.gen)
lambda.results.0.last.full <- rbind(lambda.results.0.full, lambda.results.last.full)




lambda.mean <- lambda.results %>%
  group_by(Generation, Treatment) %>%
  summarise(mean = mean(lambda, na.rm = TRUE),
            sd = sd(lambda, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)


lambda.mean.0 <- filter(lambda.mean, Generation == 0)
lambda.mean.last <- filter(lambda.mean, Generation == last.gen)
lambda.mean.0.last <- rbind(lambda.mean.0, lambda.mean.last)

lambda.mean.0.last$Treatment.scale <- lambda.mean.0.last$Treatment/4
lambda.mean.0.last$mean.scale <- lambda.mean.0.last$mean*3
lambda.mean.0.last$lower.ci.scale <- lambda.mean.0.last$lower.ci*3
lambda.mean.0.last$upper.ci.scale <- lambda.mean.0.last$upper.ci*3

# create a variable to order panels by
#f <- factor(Treatment2, levels = c("Ambient", "Acidification", "Warming", "Greenhouse"))


## Panels are wrapped by Treatment
x.pos <- 0.25

ggplot()+
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
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_blank(),
        legend.position = "bottom",
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.text = element_text(size = 20))+
  
  facet_wrap(~Treatment2)


### SAVE AS 6.6 X 16 ############################################################################################



## Panels wrapped by generation = 2 panels
ggplot()+
  
  geom_density(data = SurvTot.agg.0.last, aes(x=lx, fill = Treatment), alpha = 0.3, binwidth = 0.01)+
  
  geom_point(data = lambda.mean.0.last, aes(x=Treatment/4, y = mean.scale, color = factor(Treatment)), 
             alpha = 0.7,
             size = 7)+
  
  geom_errorbar(data = lambda.mean.0.last, aes(x=Treatment/4,
                                               ymin=lower.ci.scale,
                                               ymax=upper.ci.scale,
                                               color = factor(Treatment)),
                size = 1.0,
                width = 0.07,
                alpha = 0.7)+
  
  #geom_line(data = lambda.mean.0.25, aes(x=Treatment.scale, y = mean.scale, color = factor(Treatment), group = Treatment))+
  
  theme_minimal()+
  scale_fill_manual(values = c("blue", "forestgreen", "orange", "red"),
                    labels = c("Control", 
                               expression("Effect of C"*O[2]), 
                               "Effect of Temperature", 
                               expression("Effect of Temperature and C"*O[2])))+
  scale_color_manual(guide = "none",
                     values = c("blue", "forestgreen", "orange", "red"))+
  scale_x_continuous(name = "Survivorship")+
  scale_y_continuous(sec.axis = dup_axis(trans = ~./3, 
                                         name = expression(Mean~Absolute~Fitness~(generation^-1))))+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_blank(),
        legend.position = "bottom",
        legend.text.align = 0)+
  facet_wrap(~Generation)





##### Statistics for Survival #####

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




###############################################################################################################################################################


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
  geom_point(size = 2, position = position_dodge(width = 2))+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), size=0.3, width=0, position = position_dodge(width = 2))+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Total Development Time (days)", x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

CdevPlotTotal





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




gam.sexratio <- gam(F.Ratio ~ s(Generation.c, by = factor(Treatment), k = 3), data = SurvTotSex)
summary(gam.sexratio)


plot(gam.sexratio, pages = 1)
plot_smooth(gam.sexratio, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)
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





##### Naupliar development time #####

SurvTot.Ndev <- filter(SurvTot, SurvTot$Ndev > 0)

SurvTot.Ndev <- SurvTot.Ndev[rep(seq(nrow(SurvTot.Ndev)), SurvTot.Ndev$Ndev),]

Ndev.sum <- SurvTot.Ndev %>%
  group_by(Treatment, Generation) %>%
  summarise(mean = mean(time, na.rm = TRUE),
            sd = sd(time, na.rm = TRUE),
            n.count = n()) %>%
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)



Ndev.sum$Generation.c <- as.numeric(Ndev.sum$Generation)

NdevPlotTotal <- ggplot(data = Ndev.sum, aes(Generation.c, mean, color=factor(Treatment)))+
  geom_line(size=1.3)+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), size=1.3, width=0.3)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Naupliar Development Time (days)", x="Generation")+
  scale_x_continuous(breaks = c(0,2,4))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

NdevPlotTotal



SurvTot.Ndev$time <- as.numeric(SurvTot.Ndev$time)

SurvTot.Ndev$Treatment <- as.factor(SurvTot.Ndev$Treatment)

SurvTot.Ndev$Generation <- as.factor(as.numeric(SurvTot.Ndev$Generation))

pairwise.t.test(SurvTot.Ndev$time, SurvTot.Ndev$Generation:SurvTot.Ndev$Treatment, p.adjust.method = "none")

Ndev.pairwise <- tidy(pairwise.t.test(SurvTot.Ndev$time, SurvTot.Ndev$Generation:SurvTot.Ndev$Treatment, p.adjust.method = "none"))

Ndev.pairwise <- filter(Ndev.pairwise, Ndev.pairwise$p.value > 0.05)

SurvTot.Ndev$Generation[SurvTot.Ndev$Generation == 10] <- 9
#Ndev.anova <- SurvTot.Ndev %>%
#  group_by(Treatment) %>%
#  do(tidy(aov(time ~ Temp*pH*Generation, data = .)))

library(emmeans)
#options(contrasts = c("contr.sum","contr.poly"))
Ndev.lm.lsm <- lm(time ~ Generation+Treatment, data = SurvTot.Ndev)
#anova(Ndev.aov)


#reference grid
#Ndev.rg <- ref_grid(Ndev.lm)
#summary(Ndev.rg)

Ndev.lsm <- lsmeans(Ndev.lm.lsm, "Generation", by = "Treatment",
                    at = list(Treatment = c(1,2,3,4), 
                              Generation = c("0", "12", "15", "3", "6", "9"))) #find the least squares means for all the generations, using temp or pH just gives one mean

Ndev.lsm.sum <- summary(Ndev.lsm) #a summary data frame of the ls means for each treatment by generation


Ndev.lm <- lm(time ~ (Generation+Temp+pH)^2, data = SurvTot.Ndev)
Ndev.anova <- Anova(Ndev.lm, contrasts=list(topic=contr.sum, sys=contr.sum), 
                    type = 3)#make it a type III test which considers all factors significant at the beginning instead of in sequential order (i.e. effect of temp, then pH, then generation)


Ndev.anova


#Cdev.lm.lsm <- lm(time ~ Generation.c*Treatment, data = SurvTot.Cdev)
#Cdev.anova <- Anova(Cdev.lm.lsm)
#Cdev.anova$factors <- rownames(Cdev.anova)
#fwrite(Cdev.anova, file = "Dev_time_anova.txt", sep = "\t")
#summary(Cdev.lm.lsm)
#hist(Cdev.lm.lsm$residuals)
#foo <- plot_model(Cdev.lm.lsm, 'int') # use this to visualize the models for development time for each treatment
#foo


#anova(Cdev.aov)


#reference grid
#Cdev.rg <- ref_grid(Cdev.lm)
#summary(Cdev.rg)

#Cdev.lsm <- lsmeans(Cdev.lm.lsm, "Generation", by = "Treatment",
#                    at = list(Treatment = c(1,2,3,4), 
#                              Generation = c("0", "12", "15", "3", "6", "9"))) #find the least squares means for all the generations, using temp or pH just gives one mean

#Cdev.lsm.sum <- summary(Cdev.lsm) #a summary data frame of the ls means for each treatment by generation


#Cdev.lm <- lm(time ~ (Generation+Temp+pH)^2, data = SurvTot.Cdev)
#Cdev.anova <- Anova(Cdev.lm, contrasts=list(topic=contr.sum, sys=contr.sum), 
#                    type = 3)#make it a type III test which considers all factors significant at the beginning instead of in sequential order (i.e. effect of temp, then pH, then generation)


#Cdev.anova


#Cdev.anova <- SurvTot.Cdev %>%
#  group_by(Treatment) %>%
#  do(tidy(aov(time~Generation*Temp*pH, data = .)))






#create a function that will take the partial contents of one column to add conditional contents of another column
#taken from https://stackoverflow.com/questions/19747384/create-new-column-in-dataframe-based-on-partial-string-matching-other-column?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
#ff = function(x, patterns, replacements = patterns, fill = NA, ...)
#{
#  stopifnot(length(patterns) == length(replacements))#stop if the length of the things you're searching for is not the same as the thing you want to replace it with

#  ans = rep_len(as.character(fill), length(x))    
#  empty = seq_along(x)
#  
#  for(i in seq_along(patterns)) {
#    greps = grepl(patterns[[i]], x[empty], ...)
#    ans[empty[greps]] = replacements[[i]]  
#    empty = empty[!greps]
#  }

#  return(ans)
#}





#SurvTot$Generation <- ff(SurvTot$source, #what the first column you're searching through is 
#                         c("F0", "F3", "F6", "F9", "F12", "F15", "F25"), #what string you're looking for in that first column 
#                         c("0", "3", "6", "9", "12", "15", "25"), #what you want the output column to include in order of the first column string
#                         ignore.case = TRUE)


#SurvTot[2302:2876,12] <- 10






######################################################################################################################################
#####Old data analysis################################################################################################################
######################################################################################################################################

setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/F0")
survDataF0 <- read.table(file = "032618_SurvDataF0.txt", header = TRUE, sep = "\t")

setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/F3")
survDataF3 <- read.table(file = "032218_SurvData_F3.txt", header = TRUE, sep = "\t")

setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/F6")
survDataF6 <- read.table(file = "032318_SurvData_F6.txt", header = TRUE, sep = "\t")

setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/F9")
survDataF9 <- read.table(file = "032518_SurvData_F9.txt", header = TRUE, sep = "\t")

setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/F12")
survDataF12 <- read.table(file = "032618_SurvData_F12.txt", header = TRUE, sep = "\t")

setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/F15")
survDataF15 <- read.table(file = "032218_SurvDataF15.txt", header = TRUE, sep = "\t")

setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/F25")
survDataF25 <- read.table(file = "SurvDataF25.txt", header = TRUE, sep = "\t")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



library(dplyr)
#get rid of all times where survival is 0 or negative
survDataF0 <- filter(survDataF0, survDataF0$lx > 0)
survDataF3 <- filter(survDataF3, survDataF3$lx > 0)
survDataF6 <- filter(survDataF6, survDataF6$lx > 0)
survDataF9 <- filter(survDataF9, survDataF9$lx > 0)
survDataF12 <- filter(survDataF12, survDataF12$lx > 0)
survDataF15 <- filter(survDataF15, survDataF15$lx > 0)
survDataF25 <- filter(survDataF25, survDataF25$lx >0)

survDataF0Ndev <- filter(survDataF0, survDataF0$Ndev > 0)
survDataF0Cdev <- filter(survDataF0, survDataF0$Cdev > 0)
survDataF3Ndev <- filter(survDataF3, survDataF3$Ndev > 0)
survDataF3Cdev <- filter(survDataF3, survDataF3$Cdev > 0)
survDataF6Ndev <- filter(survDataF6, survDataF6$Ndev > 0)
survDataF6Cdev <- filter(survDataF6, survDataF6$Cdev > 0)
survDataF9Ndev <- filter(survDataF9, survDataF9$Ndev > 0)
survDataF9Cdev <- filter(survDataF9, survDataF9$Cdev > 0)
survDataF12Ndev <- filter(survDataF12, survDataF12$Ndev > 0)
survDataF12Cdev <- filter(survDataF12, survDataF12$Cdev > 0)
survDataF15Ndev <- filter(survDataF15, survDataF15$Ndev > 0)
survDataF15Cdev <- filter(survDataF15, survDataF15$Cdev > 0)
survDataF25Ndev <- filter(survDataF25, survDataF25$Ndev >0)
survDataF25Cdev <- filter(survDataF25, survDataF25$Cdev >0)
detach("package:dplyr", unload = TRUE)


survDataF25.stats <- summarySE(survDataF25, measurevar = "lx", groupvars = c("Treatment", "time"))

library(ggplot2)
survPlotF25 <- ggplot(survDataF25.stats, aes(time, lx, color=factor(Treatment)))+
  labs(x="Time (Days)", y="Survival")+
  scale_x_continuous(breaks = seq(0,18,2))+
  scale_y_continuous(limits = c(0,1.5))+
  geom_smooth(position = position_dodge(), method = loess, fill=NA)+
  geom_errorbar(aes(ymin=lx-ci, ymax=lx+ci), size=1, width=2.0, position = position_dodge())+
  scale_color_manual(values = c("orange", "red"))+ 
  theme_classic()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20)
  )

survPlotF25

surv.anova.F25 <- aov(lx~Temp+pH+Temp:pH, data = survDataF25)
summary(surv.anova.F25)

pairwise.t.test(survDataF25$lx, survDataF25$Treatment, p.adjust.method = "none")

#####F0 development#####
survDataF0.Stats<-summarySE(data = survDataF0, measurevar = "lx", groupvars = c("Treatment"))
survDataF0.Stats[is.na(survDataF0.Stats)] <- 0
survDataF0.Stats$Generation <- c(0,0,0,0)




survDataF0Ndev <- survDataF0Ndev[rep(seq(nrow(survDataF0Ndev)), survDataF0Ndev$Ndev),]

#remove unnecessary columns
survDataF0Ndev <- survDataF0Ndev[-grep('Cdev', colnames(survDataF0Ndev))]
survDataF0Ndev <- survDataF0Ndev[-grep('Ndev', colnames(survDataF0Ndev))]
survDataF0Ndev <- survDataF0Ndev[-grep('nx', colnames(survDataF0Ndev))]
survDataF0Ndev <- survDataF0Ndev[-grep('lx', colnames(survDataF0Ndev))]
survDataF0Ndev$Generation <- c(0)

F0Ndev.stats <- summarySE(data = survDataF0Ndev, measurevar = "time", groupvars = "Treatment")

F0Ndev.stats$Temp <- c(18, 18, 22, 22)
F0Ndev.stats$pH <- c(8.2,7.5,8.2,7.5)
F0Ndev.stats$Generation <- c(0,0,0,0)

#create object for x axis label for development time plots
xlab <- "Temperature (°C)"

F0NdevPlot <- ggplot(F0Ndev.stats, aes(Temp, time, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=time-ci, ymax=time+ci), size=1, width=1)+
  guides(size = guide_legend(override.aes = list(size=20)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Naupliar Development Time (days)", x=xlab)+
  scale_x_continuous(breaks = unique(F0Ndev.stats$Temp))


F0NdevPlot



survDataF0Cdev <- survDataF0Cdev[rep(seq(nrow(survDataF0Cdev)), survDataF0Cdev$Cdev),]

#remove unnecessary columns
survDataF0Cdev <- survDataF0Cdev[-grep('Ndev', colnames(survDataF0Cdev))]
survDataF0Cdev <- survDataF0Cdev[-grep('nx', colnames(survDataF0Cdev))]
survDataF0Cdev <- survDataF0Cdev[-grep('lx', colnames(survDataF0Cdev))]
survDataF0Cdev <- survDataF0Cdev[-grep('Cdev', colnames(survDataF0Cdev))]
survDataF0Cdev$Generation <- c(0)


#create summary table
F0Cdev.stats <- summarySE(data = survDataF0Cdev, measurevar = "time", groupvars = "Treatment")

#find the actual development time for the copepodite stages
F0Cdev.stats$time <- F0Cdev.stats$time-F0Ndev.stats$time

AA <- F0Ndev.stats[1,3]
AH <- F0Ndev.stats[2,3]
HA <- F0Ndev.stats[3,3]
HH <- F0Ndev.stats[4,3]

library(dplyr)
#nested ifelse statement to make the actual time during copepodite development time
#help from: https://www.listendata.com/2017/03/if-else-in-r.html
survDataF0Cdev$time=if_else(survDataF0Cdev$Treatment==1, survDataF0Cdev$time-AA,
                             if_else(survDataF0Cdev$Treatment==2, survDataF0Cdev$time-AH,
                                     if_else(survDataF0Cdev$Treatment==3, survDataF0Cdev$time-HA,
                                             if_else(survDataF0Cdev$Treatment==4, survDataF0Cdev$time-HH, 0))))
detach("package:dplyr", unload = TRUE)

F0Cdev.stats <- summarySE(data = survDataF0Cdev, measurevar = "time", groupvars = "Treatment")

#add temp and pH to summary table
F0Cdev.stats$Temp <- c(18, 18, 22, 22)
F0Cdev.stats$pH <- c(8.2, 7.5, 8.2, 7.5)
F0Cdev.stats$Generation <- c(0,0,0,0)

F0CdevPlot <- ggplot(F0Cdev.stats, aes(Temp, time, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=time-ci, ymax=time+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Copepodite Development Time (days)", x=xlab)+
  scale_x_continuous(breaks = unique(F0Cdev.stats$Temp))


F0CdevPlot


#####F3 development#####
survDataF3.Stats<-summarySE(data = survDataF3, measurevar = "lx", groupvars = c("Treatment"))
survDataF3.Stats[is.na(survDataF3.Stats)] <- 0
survDataF3.Stats$Generation <- c(3,3,3,3)

survDataF3Ndev <- survDataF3Ndev[rep(seq(nrow(survDataF3Ndev)), survDataF3Ndev$Ndev),]

#remove unnecessary columns
survDataF3Ndev <- survDataF3Ndev[-grep('Cdev', colnames(survDataF3Ndev))]
survDataF3Ndev <- survDataF3Ndev[-grep('Ndev', colnames(survDataF3Ndev))]
survDataF3Ndev <- survDataF3Ndev[-grep('nx', colnames(survDataF3Ndev))]
survDataF3Ndev <- survDataF3Ndev[-grep('lx', colnames(survDataF3Ndev))]
survDataF3Ndev$Generation <- c(3)

F3Ndev.stats <- summarySE(data = survDataF3Ndev, measurevar = "time", groupvars = "Treatment")

F3Ndev.stats$Temp <- c(18, 18, 22, 22)
F3Ndev.stats$pH <- c(8.2,7.5,8.2,7.5)
F3Ndev.stats$Generation <- c(3,3,3,3)





survDataF3Cdev <- survDataF3Cdev[rep(seq(nrow(survDataF3Cdev)), survDataF3Cdev$Cdev),]

#remove unnecessary columns
survDataF3Cdev <- survDataF3Cdev[-grep('Ndev', colnames(survDataF3Cdev))]
survDataF3Cdev <- survDataF3Cdev[-grep('nx', colnames(survDataF3Cdev))]
survDataF3Cdev <- survDataF3Cdev[-grep('lx', colnames(survDataF3Cdev))]
survDataF3Cdev <- survDataF3Cdev[-grep('Cdev', colnames(survDataF3Cdev))]
survDataF3Cdev$Generation <- c(3)


#create summary table
F3Cdev.stats <- summarySE(data = survDataF3Cdev, measurevar = "time", groupvars = "Treatment")

#find the actual development time for the copepodite stages
F3Cdev.stats$time <- F3Cdev.stats$time-F3Ndev.stats$time

AA <- F3Ndev.stats[1,3]
AH <- F3Ndev.stats[2,3]
HA <- F3Ndev.stats[3,3]
HH <- F3Ndev.stats[4,3]

library(dplyr)
#nested ifelse statement to make the actual time during copepodite development time
#help from: https://www.listendata.com/2017/03/if-else-in-r.html
survDataF3Cdev$time=if_else(survDataF3Cdev$Treatment==1, survDataF3Cdev$time-AA,
                             if_else(survDataF3Cdev$Treatment==2, survDataF3Cdev$time-AH,
                                     if_else(survDataF3Cdev$Treatment==3, survDataF3Cdev$time-HA,
                                             if_else(survDataF3Cdev$Treatment==4, survDataF3Cdev$time-HH, 0))))
detach("package:dplyr", unload = TRUE)

F3Cdev.stats <- summarySE(data = survDataF3Cdev, measurevar = "time", groupvars = "Treatment")

#add temp and pH to summary table
F3Cdev.stats$Temp <- c(18, 18, 22, 22)
F3Cdev.stats$pH <- c(8.2, 7.5, 8.2, 7.5)
F3Cdev.stats$Generation <- c(3,3,3,3)


#####F6 development#####
survDataF6.Stats<-summarySE(data = survDataF6, measurevar = "lx", groupvars = c("Treatment"))
survDataF6.Stats[is.na(survDataF6.Stats)] <- 0
survDataF6.Stats$Generation <- c(6,6,6,6)



survDataF6Ndev <- survDataF6Ndev[rep(seq(nrow(survDataF6Ndev)), survDataF6Ndev$Ndev),]

#remove unnecessary columns
survDataF6Ndev <- survDataF6Ndev[-grep('Cdev', colnames(survDataF6Ndev))]
survDataF6Ndev <- survDataF6Ndev[-grep('Ndev', colnames(survDataF6Ndev))]
survDataF6Ndev <- survDataF6Ndev[-grep('nx', colnames(survDataF6Ndev))]
survDataF6Ndev <- survDataF6Ndev[-grep('lx', colnames(survDataF6Ndev))]
survDataF6Ndev$Generation <- c(6)

F6Ndev.stats <- summarySE(data = survDataF6Ndev, measurevar = "time", groupvars = "Treatment")

F6Ndev.stats$Temp <- c(18, 18, 22, 22)
F6Ndev.stats$pH <- c(8.2,7.5,8.2,7.5)
F6Ndev.stats$Generation <- c(6,6,6,6)





survDataF6Cdev <- survDataF6Cdev[rep(seq(nrow(survDataF6Cdev)), survDataF6Cdev$Cdev),]

#remove unnecessary columns
survDataF6Cdev <- survDataF6Cdev[-grep('Ndev', colnames(survDataF6Cdev))]
survDataF6Cdev <- survDataF6Cdev[-grep('nx', colnames(survDataF6Cdev))]
survDataF6Cdev <- survDataF6Cdev[-grep('lx', colnames(survDataF6Cdev))]
survDataF6Cdev <- survDataF6Cdev[-grep('Cdev', colnames(survDataF6Cdev))]
survDataF6Cdev$Generation <- c(6)


#create summary table
F6Cdev.stats <- summarySE(data = survDataF6Cdev, measurevar = "time", groupvars = "Treatment")

#find the actual development time for the copepodite stages
F6Cdev.stats$time <- F6Cdev.stats$time-F6Ndev.stats$time

AA <- F6Ndev.stats[1,3]
AH <- F6Ndev.stats[2,3]
HA <- F6Ndev.stats[3,3]
HH <- F6Ndev.stats[4,3]

library(dplyr)
#nested ifelse statement to make the actual time during copepodite development time
#help from: https://www.listendata.com/2017/03/if-else-in-r.html
survDataF6Cdev$time=if_else(survDataF6Cdev$Treatment==1, survDataF6Cdev$time-AA,
                             if_else(survDataF6Cdev$Treatment==2, survDataF6Cdev$time-AH,
                                     if_else(survDataF6Cdev$Treatment==3, survDataF6Cdev$time-HA,
                                             if_else(survDataF6Cdev$Treatment==4, survDataF6Cdev$time-HH, 0))))
detach("package:dplyr", unload = TRUE)

F6Cdev.stats <- summarySE(data = survDataF6Cdev, measurevar = "time", groupvars = "Treatment")

#add temp and pH to summary table
F6Cdev.stats$Temp <- c(18, 18, 22, 22)
F6Cdev.stats$pH <- c(8.2, 7.5, 8.2, 7.5)
F6Cdev.stats$Generation <- c(6,6,6,6)

#####F9 development#####
survDataF9.Stats<-summarySE(data = survDataF9, measurevar = "lx", groupvars = c("Treatment"))
survDataF9.Stats[is.na(survDataF9.Stats)] <- 0
survDataF9.Stats$Generation <- c(10,10,10,9)




survDataF9Ndev <- survDataF9Ndev[rep(seq(nrow(survDataF9Ndev)), survDataF9Ndev$Ndev),]

#remove unnecessary columns
survDataF9Ndev <- survDataF9Ndev[-grep('Cdev', colnames(survDataF9Ndev))]
survDataF9Ndev <- survDataF9Ndev[-grep('Ndev', colnames(survDataF9Ndev))]
survDataF9Ndev <- survDataF9Ndev[-grep('nx', colnames(survDataF9Ndev))]
survDataF9Ndev <- survDataF9Ndev[-grep('lx', colnames(survDataF9Ndev))]
survDataF9Ndev$Generation <- c(9)

F9Ndev.stats <- summarySE(data = survDataF9Ndev, measurevar = "time", groupvars = "Treatment")

F9Ndev.stats$Temp <- c(18, 18, 22, 22)
F9Ndev.stats$pH <- c(8.2,7.5,8.2,7.5)
F9Ndev.stats$Generation <- c(10,10,10,9)





survDataF9Cdev <- survDataF9Cdev[rep(seq(nrow(survDataF9Cdev)), survDataF9Cdev$Cdev),]

#remove unnecessary columns
survDataF9Cdev <- survDataF9Cdev[-grep('Ndev', colnames(survDataF9Cdev))]
survDataF9Cdev <- survDataF9Cdev[-grep('nx', colnames(survDataF9Cdev))]
survDataF9Cdev <- survDataF9Cdev[-grep('lx', colnames(survDataF9Cdev))]
survDataF9Cdev <- survDataF9Cdev[-grep('Cdev', colnames(survDataF9Cdev))]
survDataF9Cdev$Generation <- c(9)


#create summary table
F9Cdev.stats <- summarySE(data = survDataF9Cdev, measurevar = "time", groupvars = "Treatment")

#find the actual development time for the copepodite stages
F9Cdev.stats$time <- F9Cdev.stats$time-F9Ndev.stats$time

AA <- F9Ndev.stats[1,3]
AH <- F9Ndev.stats[2,3]
HA <- F9Ndev.stats[3,3]
HH <- F9Ndev.stats[4,3]

library(dplyr)
#nested ifelse statement to make the actual time during copepodite development time
#help from: https://www.listendata.com/2017/03/if-else-in-r.html
survDataF9Cdev$time=if_else(survDataF9Cdev$Treatment==1, survDataF9Cdev$time-AA,
                             if_else(survDataF9Cdev$Treatment==2, survDataF9Cdev$time-AH,
                                     if_else(survDataF9Cdev$Treatment==3, survDataF9Cdev$time-HA,
                                             if_else(survDataF9Cdev$Treatment==4, survDataF9Cdev$time-HH, 0))))
detach("package:dplyr", unload = TRUE)

F9Cdev.stats <- summarySE(data = survDataF9Cdev, measurevar = "time", groupvars = "Treatment")

#add temp and pH to summary table
F9Cdev.stats$Temp <- c(18, 18, 22, 22)
F9Cdev.stats$pH <- c(8.2, 7.5, 8.2, 7.5)
F9Cdev.stats$Generation <- c(10,10,10,9)


#####F12 development#####
survDataF12.Stats<-summarySE(data = survDataF12, measurevar = "lx", groupvars = c("Treatment"))
survDataF12.Stats[is.na(survDataF12.Stats)] <- 0
survDataF12.Stats$Generation <- c(12,12,12)

survDataF12Ndev <- survDataF12Ndev[rep(seq(nrow(survDataF12Ndev)), survDataF12Ndev$Ndev),]

#remove unnecessary columns
survDataF12Ndev <- survDataF12Ndev[-grep('Cdev', colnames(survDataF12Ndev))]
survDataF12Ndev <- survDataF12Ndev[-grep('Ndev', colnames(survDataF12Ndev))]
survDataF12Ndev <- survDataF12Ndev[-grep('nx', colnames(survDataF12Ndev))]
survDataF12Ndev <- survDataF12Ndev[-grep('lx', colnames(survDataF12Ndev))]
survDataF12Ndev$Generation <- c(12)

F12Ndev.stats <- summarySE(data = survDataF12Ndev, measurevar = "time", groupvars = "Treatment")

F12Ndev.stats$Temp <- c(18, 22, 22)
F12Ndev.stats$pH <- c(8.2,8.2,7.5)
F12Ndev.stats$Generation <- c(12,12,12)





survDataF12Cdev <- survDataF12Cdev[rep(seq(nrow(survDataF12Cdev)), survDataF12Cdev$Cdev),]

#remove unnecessary columns
survDataF12Cdev <- survDataF12Cdev[-grep('Ndev', colnames(survDataF12Cdev))]
survDataF12Cdev <- survDataF12Cdev[-grep('nx', colnames(survDataF12Cdev))]
survDataF12Cdev <- survDataF12Cdev[-grep('lx', colnames(survDataF12Cdev))]
survDataF12Cdev <- survDataF12Cdev[-grep('Cdev', colnames(survDataF12Cdev))]
survDataF12Cdev$Generation <- c(12)


#create summary table
F12Cdev.stats <- summarySE(data = survDataF12Cdev, measurevar = "time", groupvars = "Treatment")

#find the actual development time for the copepodite stages
F12Cdev.stats$time <- F12Cdev.stats$time-F12Ndev.stats$time

AA <- F12Ndev.stats[1,3]
HA <- F12Ndev.stats[2,3]
HH <- F12Ndev.stats[3,3]


library(dplyr)
#nested ifelse statement to make the actual time during copepodite development time
#help from: https://www.listendata.com/2017/03/if-else-in-r.html
survDataF12Cdev$time=if_else(survDataF12Cdev$Treatment==1, survDataF12Cdev$time-AA,
                             if_else(survDataF12Cdev$Treatment==2, survDataF12Cdev$time-AH,
                                     if_else(survDataF12Cdev$Treatment==3, survDataF12Cdev$time-HA,
                                             if_else(survDataF12Cdev$Treatment==4, survDataF12Cdev$time-HH, 0))))
detach("package:dplyr", unload = TRUE)

F12Cdev.stats <- summarySE(data = survDataF12Cdev, measurevar = "time", groupvars = "Treatment")

#add temp and pH to summary table
F12Cdev.stats$Temp <- c(18, 22, 22)
F12Cdev.stats$pH <- c(8.2, 8.2, 7.5)
F12Cdev.stats$Generation <- c(12,12,12)




#####F15 development######
survDataF15.Stats<-summarySE(data = survDataF15, measurevar = "lx", groupvars = c("Treatment"))
survDataF15.Stats[is.na(survDataF15.Stats)] <- 0
survDataF15.Stats$Generation <- c(15,15,15,15)

survDataF15Ndev <- survDataF15Ndev[rep(seq(nrow(survDataF15Ndev)), survDataF15Ndev$Ndev),]

#remove unnecessary columns
survDataF15Ndev <- survDataF15Ndev[-grep('Cdev', colnames(survDataF15Ndev))]
survDataF15Ndev <- survDataF15Ndev[-grep('Ndev', colnames(survDataF15Ndev))]
survDataF15Ndev <- survDataF15Ndev[-grep('nx', colnames(survDataF15Ndev))]
survDataF15Ndev <- survDataF15Ndev[-grep('lx', colnames(survDataF15Ndev))]
survDataF15Ndev$Generation <- c(15)

F15Ndev.stats <- summarySE(data = survDataF15Ndev, measurevar = "time", groupvars = "Treatment")

F15Ndev.stats$Temp <- c(18, 18, 22, 22)
F15Ndev.stats$pH <- c(8.2,7.5,8.2,7.5)
F15Ndev.stats$Generation <- c(15,15,15,15)





survDataF15Cdev <- survDataF15Cdev[rep(seq(nrow(survDataF15Cdev)), survDataF15Cdev$Cdev),]

#remove unnecessary columns
survDataF15Cdev <- survDataF15Cdev[-grep('Ndev', colnames(survDataF15Cdev))]
survDataF15Cdev <- survDataF15Cdev[-grep('nx', colnames(survDataF15Cdev))]
survDataF15Cdev <- survDataF15Cdev[-grep('lx', colnames(survDataF15Cdev))]
survDataF15Cdev <- survDataF15Cdev[-grep('Cdev', colnames(survDataF15Cdev))]
survDataF15Cdev$Generation <- c(15)


#create summary table
F15Cdev.stats <- summarySE(data = survDataF15Cdev, measurevar = "time", groupvars = "Treatment")

#find the actual development time for the copepodite stages
F15Cdev.stats$time <- F15Cdev.stats$time-F15Ndev.stats$time

AA <- F15Ndev.stats[1,3]
AH <- F15Ndev.stats[2,3]
HA <- F15Ndev.stats[3,3]
HH <- F15Ndev.stats[4,3]

library(dplyr)
#nested ifelse statement to make the actual time during copepodite development time
#help from: https://www.listendata.com/2017/03/if-else-in-r.html
survDataF15Cdev$time=if_else(survDataF15Cdev$Treatment==1, survDataF15Cdev$time-AA,
                              if_else(survDataF15Cdev$Treatment==2, survDataF15Cdev$time-AH,
                                      if_else(survDataF15Cdev$Treatment==3, survDataF15Cdev$time-HA,
                                              if_else(survDataF15Cdev$Treatment==4, survDataF15Cdev$time-HH, 0))))
detach("package:dplyr", unload = TRUE)

F15Cdev.stats <- summarySE(data = survDataF15Cdev, measurevar = "time", groupvars = "Treatment")

#add temp and pH to summary table
F15Cdev.stats$Temp <- c(18, 18, 22, 22)
F15Cdev.stats$pH <- c(8.2, 7.5, 8.2, 7.5)
F15Cdev.stats$Generation <- c(15,15,15,15)

#####Graphs#####

survData.Total <- rbind(survDataF0.Stats, 
                        survDataF3.Stats, 
                        survDataF6.Stats, 
                        survDataF9.Stats,
                        survDataF12.Stats,
                        survDataF15.Stats)
library(dplyr)
#nested ifelse statement to make the actual time during copepodite development time
#help from: https://www.listendata.com/2017/03/if-else-in-r.html
survData.Total$Temp=if_else(survData.Total$Treatment==1, 18,
                             if_else(survData.Total$Treatment==2, 18,
                                     if_else(survData.Total$Treatment==3, 22,
                                             if_else(survData.Total$Treatment==4, 22, 0))))

survData.Total$pH=if_else(survData.Total$Treatment==1, 8.2,
                            if_else(survData.Total$Treatment==2, 7.5,
                                    if_else(survData.Total$Treatment==3, 8.2,
                                            if_else(survData.Total$Treatment==4, 7.5, 0))))


detach("package:dplyr", unload = TRUE)


survPlotTotal <- ggplot(data = survData.Total, aes(Generation, lx, color=factor(Treatment)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=lx-ci, ymax=lx+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Survival", x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15))+
  scale_y_continuous(limits = c(0.75,1), breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

survPlotTotal


surv.anova <- aov(lx~Temp*pH*Generation, data = survData.Total)
surv.anova2 <- update(surv.anova, . ~ . - Temp:pH:Generation)
anova(surv.anova, surv.anova2)



summary(surv.anova)


Ndev.total <- rbind(F0Ndev.stats, 
                    F3Ndev.stats, 
                    F6Ndev.stats, 
                    F9Ndev.stats,
                    F12Ndev.stats, 
                    F15Ndev.stats)

NdevPlotTotal <- ggplot(data = Ndev.total, aes(Generation, time, color=factor(Treatment)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=time-ci, ymax=time+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Naupliar Development Time (days)", x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20), 
              axis.text.x = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              axis.text.y = element_text(size = 20))

NdevPlotTotal

Cdev.total <- rbind(F0Cdev.stats, 
                    F3Cdev.stats, 
                    F6Cdev.stats,
                    F9Cdev.stats,
                    F12Cdev.stats, 
                    F15Cdev.stats)
CdevPlotTotal <- ggplot(data = Cdev.total, aes(Generation, time, color=factor(Treatment)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=time-ci, ymax=time+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Copepodite Development Time (days)", x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

CdevPlotTotal


Ndev.complete <- rbind(survDataF0Ndev, 
                        survDataF3Ndev, 
                        survDataF6Ndev, 
                        survDataF9Ndev, 
                        survDataF12Ndev,
                        survDataF15Ndev)

Ndev.anova <- aov(time~Temp*pH*Generation, data = Ndev.complete)
#Ndev.anova2 <- update(Ndev.anova, . ~ . - Temp:pH:Generation)
#anova(Ndev.anova, Ndev.anova2)

summary(Ndev.anova)


Cdev.complete <- rbind(survDataF0Cdev, 
                    survDataF3Cdev, 
                    survDataF6Cdev, 
                    survDataF9Cdev, 
                    survDataF12Cdev,
                    survDataF15Cdev)


Cdev.anova <- aov(time~Temp*pH*Generation, data = Cdev.complete)
summary(Cdev.anova)
