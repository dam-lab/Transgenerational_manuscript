
wd <- "C:/Users/james/Documents/Grad_school/OA_Project/EPR/Data_frames/"

#wd <- "C:/Users/james/Documents/Grad_school/OA_hudsonica/EPR/Data_frames/"



if (wd == "C:/Users/james/Documents/Grad_school/OA_Project/EPR/Data_frames/") {
  
  
  fit.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Fitness/"
  epr.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/EPR/"

  
} else if (wd == "C:/Users/james/Documents/Grad_school/OA_hudsonica/EPR/Data_frames/") {
  
  
  fit.directory <- "C:/Users/james/Documents/Grad_school/OA_hudsonica/Fitness/"
  epr.directory <- "C:/Users/james/Documents/Grad_school/OA_hudsonica/EPR/"
  
  
}






library(ggplot2)

library(data.table)

library(tidyr)

library(broom)

library(dplyr)

setwd(paste(epr.directory,"Data_frames/", sep = ""))

#Take all the lists available and combine them all as you read them in
#received help from: https://stackoverflow.com/questions/38312864/adding-new-column-and-its-value-based-on-file-name?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
#LF = list.files(pattern = "EPR_F_")
#LF

#EPRtot <- rbindlist(lapply(
#  setNames(LF,LF), #what you want to be "lapplied"
#  fread), #what you want the lapply to do
#  idcol = "source", #what the new column will be named (i.e. the header of the new column)
#  fill = TRUE) #fills missing columns with NA's

#fwrite(EPRtot, file = paste(epr.directory,"Data_frames/EPR_HF_data_total.txt", sep = "\t"))


EPRtot <- fread(paste(epr.directory,"Data_frames/EPR_HF_data_total.txt", sep = ""))


#remove the last row because it adds an extra row in the summary of NaNs
last.row <- nrow(EPRtot)
EPRtot <- EPRtot[-last.row,]


if (EPRtot$pH == 8.13) {

  EPRtot$pH[EPRtot$pH == '8.13'] <- 8.20

  }


eprStatsAll <- EPRtot %>% 
  group_by(Treatment, Generation) %>% 
  summarise(Rate = mean(EPRtot, na.rm = TRUE),
            sd = sd(EPRtot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = Rate - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = Rate + qt(1 - (0.05 / 2), n.count-1)*se)

if (wd == "C:/Users/james/Documents/Grad_school/OA_hudsonica/EPR/Data_frames/") {
new.df <- data.frame(Treatment = c("AA", "HH"),
                     Generation = c(11, 11),
                     Rate = c(23.8889, 9.0000),
                     sd = c(7.868019, 5.829219),
                     n.count = c(12, 12),
                     se = c(2.271301, 1.682751),
                     lower.ci = c(18.889788, 5.296291),
                     upper.ci = c(28.88799, 12.70371))

eprStatsAll <- rbind(eprStatsAll, new.df)
}



hfStatsAll <- EPRtot %>% 
  group_by(Treatment, Generation) %>% 
  summarise(HF = mean(Hftot, na.rm = TRUE),
            sd = sd(Hftot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = HF - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = HF + qt(1 - (0.05 / 2), n.count-1)*se)

if (wd == "C:/Users/james/Documents/Grad_school/OA_hudsonica/EPR/Data_frames/") {
new.df <- data.frame(Treatment = c("AA", "HH"),
                     Generation = c(11, 11),
                     HF = c(0.9120532, 0.7157698),
                     sd = c(0.1399012, 0.3545944),
                     n.count = c(12, 12),
                     se = c(0.04038601, 0.10236260),
                     lower.ci = c(0.8231642, 0.4904713),
                     upper.ci = c(1.0009422, 0.9410684))


hfStatsAll <- rbind(hfStatsAll, new.df)
}



eprStatsAll <- eprStatsAll[-1,]
hfStatsAll <- hfStatsAll[-1,]


#####combine all plots with HF on the same plot#####

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

epr.bar.graphs <- ggplot(data = eprSumComplete, aes(Generation, Rate, fill = factor(Treatment)))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8)+
  geom_errorbar(aes(ymin=lower.ci, 
                    ymax=upper.ci), 
                size=1.3,
                width = 1.0, 
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,2,4,11))+
  geom_line(aes(y=HF), 
            size=1.3,
            colour="purple")+
  geom_errorbar(aes(ymin=lower.ci2, ymax=upper.ci2), 
                size=1.3,
                width = 1,
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Frequency"))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

epr.bar.graphs

epr.bar.graphs+
  facet_wrap(~Treatment)+
  theme(
    strip.text.x = element_blank()#remove white strips with labelling
  )


#### separate bar graphs for each treatment

## AA
epr.bar.graphs.AA <- ggplot(data = eprSumAA, aes(Generation, Rate, fill = factor(Treatment)))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8)+
  geom_errorbar(aes(ymin=Rate-ci, 
                    ymax=Rate+ci), 
                size=1.3, 
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("blue"))+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  geom_line(aes(y=HF), 
            size=1.3,
            colour="purple")+
  geom_errorbar(aes(ymin=HF-ci2, ymax=HF+ci2), 
                size=1.3, 
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Frequency"))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


epr.bar.graphs.AA

## AH
epr.bar.graphs.AH <- ggplot(data = eprSumAH, aes(Generation, Rate, fill = factor(Treatment)))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8)+
  geom_errorbar(aes(ymin=Rate-ci, 
                    ymax=Rate+ci), 
                size=1.3, 
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("forestgreen"))+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  geom_line(aes(y=HF), 
            size=1.3,
            colour="purple")+
  geom_errorbar(aes(ymin=HF-ci2, ymax=HF+ci2), 
                size=1.3, 
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Frequency"))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


epr.bar.graphs.AH


## HA
epr.bar.graphs.HA <- ggplot(data = eprSumHA, aes(Generation, Rate, fill = factor(Treatment)))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8)+
  geom_errorbar(aes(ymin=Rate-ci, 
                    ymax=Rate+ci), 
                size=1.3, 
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("orange"))+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  geom_line(aes(y=HF), 
            size=1.3,
            colour="purple")+
  geom_errorbar(aes(ymin=HF-ci2, ymax=HF+ci2), 
                size=1.3, 
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Frequency"))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


epr.bar.graphs.HA

## HH
epr.bar.graphs.HH <- ggplot(data = eprSumHH, aes(Generation, Rate, fill = factor(Treatment)))+
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           size =4,
           width = 1.8)+
  geom_errorbar(aes(ymin=Rate-ci, 
                    ymax=Rate+ci), 
                size=1.3, 
                position = position_dodge())+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("red"))+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), 
       x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  geom_line(aes(y=HF), 
            size=1.3,
            colour="purple")+
  geom_errorbar(aes(ymin=HF-ci2, ymax=HF+ci2), 
                size=1.3, 
                position = position_dodge(0.3),
                colour = "purple")+
  scale_y_continuous(sec.axis = sec_axis(~./max(eprStatsAll$Rate), 
                                         name = "Hatching Frequency"))+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


epr.bar.graphs.HH




#eprStatsAll$Rate[eprStatsAll$Rate==0] <- NA
#hfStatsAll$HF[hfStatsAll$HF==0] <- NA

eprPlotTotal <- ggplot(data = eprStatsAll, aes(Generation, EPRtot, color=factor(Treatment)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=EPRtot-ci, ymax=EPRtot+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), x="Generation")+
  scale_x_continuous(breaks = c(0,2,4))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

eprPlotTotal

hfPlotTotal <- ggplot(data = hfStatsAll, aes(Generation2, HF, color=factor(Treatment2)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci2, ymax=HF+ci2), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Hatching Frequency", x="Generation")+
  scale_x_continuous(breaks = c(0,2,4))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


hfPlotTotal 



##### Fitness landscapes for epr #####


lambda.results <- fread(paste(fit.directory,"lambda_results_devtime_surv_epr_hf_sex_standardized_relative.txt", sep = ""))

EPRtot.0 <- filter(EPRtot, Generation == 0)

last.gen <- max(lambda.results$Generation)

EPRtot.last <- filter(EPRtot, Generation == last.gen)

EPRtot.0.last <- rbind(EPRtot.0, EPRtot.last)


EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="AA"] <- 1
EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="AH"] <- 2
EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="HA"] <- 3
EPRtot.0.last$Treatment[EPRtot.0.last$Treatment=="HH"] <- 4




EPRtot.0.last.HH <- filter(EPRtot.0.last, EPRtot.0.last$Treatment == "HH")

# test if frequency distributions are different between gens 0 and 25
chisq.test(EPRtot.0.last.HH$EPRtot, EPRtot.0.last.HH$Generation) # not significant :(
chisq.test(EPRtot.0.last.HH$HFtot, EPRtot.0.last.HH$Generation) # significant!




lambda.results.0.full <- filter(lambda.results, Generation == 0)
lambda.results.last.full <- filter(lambda.results, Generation == last.gen)
lambda.results.0.last.full <- rbind(lambda.results.0.full, lambda.results.last.full)

lambda.results.0.last.full$Treatment2 <- case_when(lambda.results.0.last.full$Treatment == 1 ~ "Ambient",
                                      lambda.results.0.last.full$Treatment == 2 ~ "Acidification",
                                      lambda.results.0.last.full$Treatment == 3 ~ "Warming",
                                      lambda.results.0.last.full$Treatment == 4 ~ "Greenhouse")

# change the order of the factor to be in the order you want. Help from: https://stackoverflow.com/questions/5490638/how-to-change-the-order-of-facet-labels-in-ggplot-custom-facet-wrap-labels
lambda.results.0.last.full <- within(lambda.results.0.last.full, 
                        Treatment2 <- factor(Treatment2, 
                                             
                                             # put the specific levels in the order you want
                                             
                                             levels = c("Ambient", # 1
                                                        "Acidification", # 2
                                                        "Warming", # 3
                                                        "Greenhouse") # 4 
                                             
                        ))

scale <- 0.0166666666666666666666666666666666666667



## Panels are wrapped by Treatment
x.pos <- 62.5
y.pos <- 0.02


EPRtot.0.last.mean <- EPRtot.0.last %>% 
  group_by(Treatment, Generation) %>% 
  summarise(Rate = mean(EPRtot, na.rm = TRUE),
            sd = sd(EPRtot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = Rate - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = Rate + qt(1 - (0.05 / 2), n.count-1)*se)

## calculate pearson correlation coefficient

lambda.results.0.last.full <- unite(lambda.results.0.last.full,
                Gen.treat,
                c(Generation, Treatment),
                remove = FALSE)


gp = group_by(lambda.results.0.last.full, Gen.treat)

summarise(gp, cor(epr, lambda))


ggplot()+
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
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_blank(),
        legend.position = "bottom",
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.text = element_text(size = 20))+
  facet_wrap(~Treatment2)



## try it with Hatching frequency
HF.0.last.mean <- EPRtot.0.last %>% 
  group_by(Treatment, Generation) %>% 
  summarise(Freq = mean(HFtot, na.rm = TRUE),
            sd = sd(HFtot, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = Freq - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = Freq + qt(1 - (0.05 / 2), n.count-1)*se)



ggplot()+
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
  
  
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_blank(),
        legend.position = "bottom",
        legend.text.align = 0,
        legend.text = element_text(size = 20),
        legend.title = element_blank())+
  facet_wrap(~Treatment2)#+ylim(0,1000)

### SAVE THE HF FIGURE AS 6.6 x 16 INCHES #########################################################################################################

ggplot()+
  geom_smooth(data = lambda.results.0.last.full, aes(x = hf, y = lambda.rel, color = factor(Generation), group = factor(Generation)),
              method = "lm")+
  facet_wrap(~Treatment)

## Panels are wrapped by Generation
ggplot()+
  geom_density(data = EPRtot.0.last, aes(x=EPRtot, 
                                         #group = Generation, 
                                         fill = Treatment), alpha = 0.3, binwidth = 0.01)+
  geom_point(data = lambda.mean.0.last, aes(x=Treatment*scale, y = mean/scale, color = factor(Treatment)), 
             #alpha = 0.7,
             size = 7)+
  geom_errorbar(data = lambda.mean.0.last, aes(x=Treatment*scale,
                                               ymin=lower.ci/scale,
                                               ymax=upper.ci/scale,
                                               color = factor(Treatment)),
                size = 1.0,
                width = 8,
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
  
  
  scale_x_continuous(name = expression(Egg~Production~Rate~(female^-1~day^-1)))+
  
  scale_y_continuous(sec.axis = dup_axis(trans = ~.*scale, 
                                         name = expression(Mean~Absolute~Fitness~(generation^-1))))+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text.x = element_text(size = 24),
        strip.text.y = element_blank(),
        legend.position = "bottom",
        legend.text.align = 0)+
  facet_wrap(~Generation)




##### F0 plots #####

Temp <- c(13,13,15,15)
pH <- c(8.2,7.8,8.2,7.8)

eprStatsF0 <- cbind(eprStatsAll, Temp, pH)
hfStatsF0 <- cbind(hfStatsAll, Temp, pH)

eprANOVAF0 <- aov(EPRtot~Temp+pH+Temp:pH, data = EPRtot)
F0_epr_Anova <- tidy(eprANOVAF0)

hfANOVAF0 <- aov(Hftot~Temp+pH+Temp:pH, data = EPRtot)
F0_hf_Anova <- tidy(hfANOVAF0)

fwrite(F0_epr_Anova, file = paste(epr.directory,"F0_EPR_anova.txt", sep = ""), sep = "\t")
fwrite(F0_hf_Anova, file = paste(epr.directory,"F0_HF_anova.txt", sep = ""), sep = "\t")

xlab <- "Temperature (°C)"

eprPlotF0 <- ggplot(eprStatsF0, aes(Temp, EPRtot, color=factor(pH)))+
  geom_point(size=1.3)+
  geom_line(size=0.5)+
  geom_errorbar(aes(ymin=EPRtot-ci, ymax=EPRtot+ci), size=0.5, width=0.3)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), x=xlab)+
  scale_x_continuous(breaks = unique(eprStatsF0$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F0 EPR')

eprPlotF0

hfPlotF0 <- ggplot(hfStatsF0, aes(Temp, Hftot, color=factor(pH)))+
  geom_point(size=1.3)+
  geom_line(size=0.5)+
  geom_errorbar(aes(ymin=Hftot-ci, ymax=Hftot+ci), size=0.5, width=0.3)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Hatching Frequency", x=xlab)+
  scale_x_continuous(breaks = unique(hfStatsF0$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F0 Hatching Frequency')

hfPlotF0
###############################################################################################################################################################
##### Stats for EPR and HF #####

## Notes: gams can't be used as a three-way interaction. 
## To do a pairwise comparison, you need to create linear models to feed into emmeans.


library(dplyr)
library(broom)
library(car)
library(sjPlot)

EPRtot <- fread(paste(epr.directory,"Data_frames/EPR_HF_data_total.txt", sep = ""))

#create linear models that are tested against each other grouped by generation and analyzes EPR/HF by temp, pH, and temp*pH
EPRtot <- EPRtot[!is.na(EPRtot),]
EPRtot <- EPRtot[!is.na(HFtot),]
EPRtot$Generation[EPRtot$Generation==10] <-9


#EPRtot <- na.omit(EPRtot)
EPRtot$Generation.c <- as.numeric(EPRtot$Generation)
EPRtot$Generation <- as.factor(as.numeric(EPRtot$Generation))
EPRtot$Treatment <- as.factor(EPRtot$Treatment)


EPRtot <- unite(EPRtot,
                Treat.Rep,
                c(Treatment, Rep),
                remove = FALSE)

EPRtot <- unite(EPRtot,
                Gen.treat,
                c(Generation, Treatment),
                remove = FALSE)


# remove for where there is no data
EPRtot <- subset(EPRtot, Gen.treat != "9_HA")

## create gam models to test if epr and hf change across generations

library(mgcv) # s() indicates a smooth function
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
library(itsadug)

gam.list <- list(gam1, gam2)

lapply(gam.list, function (x) plot_smooth(x, view = "Generation.c", plot_all = "Treatment", rug = FALSE))



#plot(gam1, pages = 1)
#plot_smooth(gam1, view="Generation.c", plot_all="Treatment", rug=FALSE)
#plot_smooth(gam1, view="Generation.c", plot_all="Beak", rug=FALSE)

#plot_smooth(mm$gam, view="Generation.c", plot_all="Treatment", rug=FALSE)


#plot(gam2, pages = 1)
#plot_smooth(gam2, view="Generation.c", plot_all="Treatment", rug=FALSE)

#anova.gam(gam2)






# create models for Tukey pairwise comparisons

library(lme4)
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






# pairwise comparisons with grouping variables
epr.pairwise.all <- tidy(pairwise.t.test(EPRtot$EPRtot, EPRtot$Generation:EPRtot$Treatment, p.adj = "none"))

epr.pairwise.all <- filter(epr.pairwise.all, epr.pairwise.all$p.value < 0.05)

hf.pairwise.all <- tidy(pairwise.t.test(EPRtot$HFtot, EPRtot$Generation:EPRtot$Treatment, p.adj = "none"))

hf.pairwise.all <- filter(hf.pairwise.all, hf.pairwise.all$p.value < 0.05)

fwrite(epr.pairwise.all, file = paste(epr.directory,"Statistics/EPR.pairwise.all.txt", sep = ""), sep = "\t")

fwrite(hf.pairwise.all, file = paste(epr.directory,"Statistics/hf.pairwise.all.txt", sep = ""), sep = "\t")




EPR.lm2 <- lm(EPRtot ~ Generation*Treatment, data = EPRtot)
emm_EPR <- emmeans(EPR.lm2, pairwise ~ Generation | Treatment)

p <- pairs(emm_EPR)
p
epr.groups <- CLD(emm_EPR)


HF.lm2 <- lm(HFtot ~ Generation*Treatment, data = EPRtot)

emm_HF <- emmeans(HF.lm2, pairwise ~ Generation | Treatment)


p2 <- pairs(emm_HF)
p2

hf.groups <- CLD(emm_HF)

fwrite(epr.groups, file = paste(epr.directory,"Statistics/EPR_groups.txt", sep = ""), sep = "\t")
fwrite(hf.groups, file = paste(epr.directory,"Statistics/HF_groups.txt", sep = ""), sep = "\t")



EPRtot$Temp <- as.factor(as.numeric(EPRtot$Temp))
EPRtot$pH <- as.factor(as.numeric(EPRtot$pH))


library(car)
## 3 way anova

epr.3.way <- aov(EPRtot~Generation.c*Temp*pH, data = EPRtot)
summary(epr.3.way)
Three.way.anova.epr <- Anova(epr.3.way)
Three.way.anova.epr$factors <- rownames(Three.way.anova.epr)
fwrite(Three.way.anova.epr, file = paste(epr.directory, "Statistics/Three_way_anova_epr.txt", sep = ""), sep = "\t")

hf.3.way <- aov(HFtot~Generation*Temp*pH, data = EPRtot)
summary(hf.3.way)
Three.way.anova.hf <- Anova(hf.3.way)
Three.way.anova.hf$factors <- rownames(Three.way.anova.hf)
fwrite(Three.way.anova.hf, file = paste(epr.directory, "Statistics/Three_way_anova_hf.txt", sep = ""), sep = "\t")



##########################################################################################################################################################################################################
##### Outdated stats methods #####

# test to see if random effects are significant
#library(lmerTest)

#rand(l)



#l2 <- lm(EPRtot~Generation.c*Treat.Rep, data = EPRtot)

#summary(l2)

#AIC(l, l2)

#drop1(l, test = "Chisq")

#l3 <- lm(EPRtot~Generation*Treat.Rep, data = EPRtot)
#pairs(emmeans(l, pairwise ~ Treatment | Generation))

#hist(l2$residuals)




# create linear models for pairwise interactions 
# default pairwise comparison is Tukey HSD for emmeans

#EPR.lm <- lm(EPRtot~Generation.c*Treatment, data = EPRtot)
#summary(EPR.lm)

#AIC(gam1, EPR.lm)


#HF.lm <- lm(HFtot ~ Generation.c*Treatment, data = EPRtot)
#summary(HF.lm)
#hf.anova <- Anova(HF.lm)
#hf.anova$factors <- rownames(hf.anova$factors)
#plot_model(HF.lm, 'int')
#fwrite(hf.anova, file = "HF_anova.txt", sep = "\t")

#epr.anova <- Anova(EPR.lm)
#epr.anova$factors <- rownames(epr.anova)
#plot_model(EPR.lm, 'int')
#fwrite(epr.anova, file = "EPR_anova.txt", sep = "\t")


# create a model with three way interactions to test for individual effects of temp and pH on EPR and HF



#Epr.broom <- EPRtot %>%
#  group_by(Treatment) %>%
#  do(Anova(aov(EPRtot ~ Generation, data = .)))




# levenes test to test for differences in variance
# values less than 0.05 indicate that the variances are not random and require a non-parametric anova
# values greater than 0.05 indicate random variance and homogeneity
#Epr.levene <- EPRtot %>%
#  group_by(Treatment) %>%
#  do(leveneTest(EPRtot~Generation, data = .))

# use a non-paramentric test for the HH treatment only because Levene's test is >0.05
#HH.kw.test <- kruskal.test(Rate~Generation, data = EPRtot, subset = Treatment=="HH")

#EPRtot$Generation <- as.factor(EPRtot$Generation)
#EPRtot$Treatment <- as.factor(EPRtot$Treatment)

#epr.gen.lm <- lm(Rate ~ Generation*Treatment, data = EPRtot)
#leveneTest(Rate~Generation*Treatment, data = EPRtot)

#epr.gen.anova <- Anova(epr.gen.lm)

#Epr.tukey <- EPRtot %>%
#  group_by(Treatment) %>%
#  do(TukeyHSD(aov(Rate ~ as.factor(Generation), data = .)))

#epr.pH.temp <- EPRtot %>%
#  group_by(Generation) %>%
#  do(tidy(aov(Rate ~ Temp*pH, data = .)))

#Epr.broom

#HF.broom <- EPRtot %>%
#  group_by(Treatment) %>%
#  do(tidy(aov(HF ~ Generation, data = .)))


#HF.tukey <- EPRtot %>%
#  group_by(Treatment) %>%
#  do(tidy(TukeyHSD(aov(HF ~ as.factor(Generation), data = .))))

#HF.levene <- EPRtot %>%
#  group_by(Treatment) %>%
#  do(leveneTest(HF~Generation, data = .))

#AH.kw.HF.test <- kruskal.test(HF~Generation, data = EPRtot, subset = Treatment=="AH")

fwrite(Epr.tukey, file = "C:/users/james/Documents/Grad_school/OA_Project/EPR/epr.tukey.txt")
fwrite(HF.tukey, file = "C:/users/james/Documents/Grad_school/OA_Project/EPR/hf.tukey.txt")


hf.pH.temp <- EPRtot %>%
  group_by(Generation) %>%
  do(tidy(aov(HF ~ Temp*pH, data = .)))

HF.broom
epr.gen <- aov(Rate ~ as.factor(Generation), data = EPRtot)
tidy(epr.gen)


gen.posthoc <- TukeyHSD(x=epr.gen, "Generation")
epr.posthoc <- tidy(gen.posthoc)


fwrite(epr.pH.temp, file = "C:/users/james/Documents/Grad_school/OA_Project/EPR/EPR.pH.temp.stats.txt")
fwrite(hf.pH.temp, file = "C:/users/james/Documents/Grad_school/OA_Project/EPR/hf.pH.temp.stats.txt")


fwrite(Epr.broom, file = "C:/users/james/Documents/Grad_school/OA_Project/EPR/EPR.stats.txt")
fwrite(HF.broom, file = "C:/users/james/Documents/Grad_school/OA_Project/EPR/hf.stats.txt")


EPRtot$Generation <- as.factor(EPRtot$Generation)
EPRtot$Treatment <- as.factor(EPRtot$Treatment)





###### Attempts with gam models #####


gam3 <- gam(EPRtot ~ Generation.c*Temp*pH, data = EPRtot)
gam3 <- gam(EPRtot ~ s(Generation.c, by = Temp, k = 7), data = EPRtot) # k = 7 seems to be the largest that is able to use

gam.check(gam3)
summary(gam3)
anova.gam(gam3)
plot_smooth(gam3, view="Generation.c", plot_all="Temp", rug=FALSE)
#plot_model(gam3, 'int')


gam4 <- gam(EPRtot ~ s(Generation.c, by = pH, k = 7), data = EPRtot)
gam.check(gam4)


summary(gam4)
anova.gam(gam4)
plot_smooth(gam4,view="Generation.c", plot_all="pH", rug=FALSE)

EPRtot$Temp.n <- as.numeric(EPRtot$Temp)
EPRtot$pH.n <- as.numeric(EPRtot$pH)
m <- loess(EPRtot~Generation.c*Temp.n*pH.n, data = EPRtot)
m
summary(m)
predict(m)
loess.smooth(EPRtot$EPRtot, EPRtot$Generation.c)

library(interactions)
m2 <- lm(EPRtot ~ Generation.c * Temp * pH, data = EPRtot)
ss <- sim_slopes(m2, pred = Generation.c, modx = Temp, mod2 = pH, johnson_neyman = FALSE)
plot(ss)


########################################################################################################################################
#####Old data analysis method#####
#######################################################################################################################################
setwd("C:/users/james/Documents/Grad_school/OA_Project/EPR/")
eprDataAll <- read.table(file = "022818_EPR_all.txt", header = TRUE, sep = "\t")

setwd("C:/users/james/Documents/Grad_school/OA_Project/EPR/F0")
eprDataF0 <- read.table(file = "022818_EPR_F0.txt", header = TRUE, sep = "\t")

setwd("C:/users/james/Documents/Grad_school/OA_Project/EPR/F3")
eprDataF3 <- read.table(file = "022818_EPR_F3.txt", header = TRUE, sep = "\t")

setwd("C:/users/james/Documents/Grad_school/OA_Project/EPR/F6")
eprDataF6 <- read.table(file = "022818_EPR_F6.txt", header = TRUE, sep = "\t")

setwd("C:/users/james/Documents/Grad_school/OA_Project/EPR/F9")
eprDataF9 <- read.table(file = "022818_EPR_F9.txt", header = TRUE, sep = "\t")

setwd("C:/users/james/Documents/Grad_school/OA_Project/EPR/F12")
eprDataF12 <- read.table(file = "022818_EPR_F12.txt", header = TRUE, sep = "\t")

setwd("C:/users/james/Documents/Grad_school/OA_Project/EPR/F15")
eprDataF15 <- read.table(file = "022818_EPR_F15.txt", header = TRUE, sep = "\t")

setwd("C:/users/james/Documents/Grad_school/OA_Project/EPR/F25")
eprDataF25 <- read.table(file = "072518_EPR_F25.txt", header = TRUE, sep = "\t")

##### old summaryse #####
#change the HA F25 generation column to 25 for some reason
#EPRtot[949:982,2] <- 25
detach("package:dplyr", unload = TRUE)
#have to set na.rm to TRUE because of the N/A values
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=TRUE) {
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
##### Continue #####

eprStatsAll <- summarySE(data = eprDataAll, measurevar = "Rate", groupvars = "Treatment")
hfStatsAll <- summarySE(data = eprDataAll, measurevar = "HF", groupvars = "Treatment")




eprStatsAll$Temp <- c(18, 18, 22, 22)
eprStatsAll$pH <- c(8.2,7.5,8.2,7.5)

xlab <- "Temperature (°C)"

eprPlotAll <- ggplot(eprStatsAll, aes(Temp, Rate, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Rate-ci, ymax=Rate+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(day^-1)), x=xlab)+
  scale_x_continuous(breaks = unique(eprStatsAll$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())


eprPlotAll

hfStatsAll$Temp <- c(18, 18, 22, 22)
hfStatsAll$pH <- c(8.2,7.5,8.2,7.5)

hfPlot <- ggplot(hfStatsAll, aes(Temp, HF, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci, ymax=HF+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Hatching Frequency", x=xlab)+
  scale_x_continuous(breaks = unique(hfStatsAll$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())


hfPlot

#eprANOVA <- aov(Rate~Temp+pH+Temp:pH, data = eprData)
#summary(eprANOVA)

#hfANOVA <- aov(HF~Temp+pH+Temp:pH, data = eprData)
#summary(hfANOVA)


library(gridExtra)
library(ggpubr)

Temp <- c(18,18,22,22)
pH <- c(8.2,7.5,8.2,7.5)
#####F0#####

eprStatsF0 <- summarySE(data = eprDataF0, measurevar = "Rate", groupvars = "Treatment")
hfStatsF0 <- summarySE(data = eprDataF0, measurevar = "HF", groupvars = "Treatment")
eprStatsF0$Generation <- c(0,0,0,0)
hfStatsF0$Generation <- c(0,0,0,0)
eprStatsF0 <- cbind(eprStatsF0, Temp, pH)
hfStatsF0 <- cbind(hfStatsF0, Temp, pH)

eprANOVAF0 <- aov(Rate~Temp+pH+Temp:pH, data = eprDataF0)
summary(eprANOVAF0)

#other ways of comparing statistical data
#from: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
compare_means(Rate~Treatment, data = eprDataF0)

hfANOVAF0 <- aov(HF~Temp+pH+Temp:pH, data = eprDataF0)
summary(hfANOVAF0)

#eprF0table <- cbind(eprANOVAF0)

eprPlotF0 <- ggplot(eprStatsF0, aes(Temp, Rate, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Rate-ci, ymax=Rate+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), x=xlab)+
  scale_x_continuous(breaks = unique(eprStatsF0$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F0 EPR')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

eprPlotF0

hfPlotF0 <- ggplot(hfStatsF0, aes(Temp, HF, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci, ymax=HF+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Hatching Frequency", x=xlab)+
  scale_x_continuous(breaks = unique(hfStatsF0$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F0 Hatching Frequency')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

hfPlotF0



#####F3#####
eprStatsF3 <- summarySE(data = eprDataF3, measurevar = "Rate", groupvars = "Treatment")
hfStatsF3 <- summarySE(data = eprDataF3, measurevar = "HF", groupvars = "Treatment")
eprStatsF3$Generation <- c(3,3,3,3)
hfStatsF3$Generation <- c(3,3,3,3)
eprStatsF3 <- cbind(eprStatsF3, Temp, pH)
hfStatsF3 <- cbind(hfStatsF3, Temp, pH)

eprPlotF3 <- ggplot(eprStatsF3, aes(Temp, Rate, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Rate-ci, ymax=Rate+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(day^-1)), x=xlab)+
  scale_x_continuous(breaks = unique(eprStatsF3$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F3 EPR')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

eprPlotF3

hfPlotF3 <- ggplot(hfStatsF3, aes(Temp, HF, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci, ymax=HF+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Hatching Frequency", x=xlab)+
  scale_x_continuous(breaks = unique(hfStatsF3$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F3 Hatching Frequency')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

hfPlotF3

#####F6#####
eprStatsF6 <- summarySE(data = eprDataF6, measurevar = "Rate", groupvars = "Treatment")
hfStatsF6 <- summarySE(data = eprDataF6, measurevar = "HF", groupvars = "Treatment")
eprStatsF6$Generation <- c(6,6,6,6)
hfStatsF6$Generation <- c(6,6,6,6)
eprStatsF6 <- cbind(eprStatsF6, Temp, pH)
hfStatsF6 <- cbind(hfStatsF6, Temp, pH)

eprPlotF6 <- ggplot(eprStatsF6, aes(Temp, Rate, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Rate-ci, ymax=Rate+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(day^-1)), x=xlab)+
  scale_x_continuous(breaks = unique(eprStatsF6$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F6 EPR')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

eprPlotF6

hfPlotF6 <- ggplot(hfStatsF6, aes(Temp, HF, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci, ymax=HF+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Hatching Frequency", x=xlab)+
  scale_x_continuous(breaks = unique(hfStatsF6$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F6 Hatching Frequency')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

hfPlotF6

#####F9#####
eprStatsF9 <- summarySE(data = eprDataF9, measurevar = "Rate", groupvars = "Treatment")
hfStatsF9 <- summarySE(data = eprDataF9, measurevar = "HF", groupvars = "Treatment")
eprStatsF9$Generation <- c(10,10,10,9)
hfStatsF9$Generation <- c(10,10,10,9)
eprStatsF9 <- cbind(eprStatsF9, Temp, pH)
hfStatsF9 <- cbind(hfStatsF9, Temp, pH)

eprPlotF9 <- ggplot(eprStatsF9, aes(Temp, Rate, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Rate-ci, ymax=Rate+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(day^-1)), x=xlab)+
  scale_x_continuous(breaks = unique(eprStatsF9$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F9 EPR')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


eprPlotF9

hfPlotF9 <- ggplot(hfStatsF9, aes(Temp, HF, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci, ymax=HF+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Hatching Frequency", x=xlab)+
  scale_x_continuous(breaks = unique(hfStatsF9$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F9 Hatching Frequency')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

hfPlotF9

#####F12#####
eprStatsF12 <- summarySE(data = eprDataF12, measurevar = "Rate", groupvars = "Treatment")
hfStatsF12 <- summarySE(data = eprDataF12, measurevar = "HF", groupvars = "Treatment")
eprStatsF12$Generation <- c(12,12,12,12)
hfStatsF12$Generation <- c(12,12,12,12)
eprStatsF12 <- cbind(eprStatsF12, Temp, pH)
hfStatsF12 <- cbind(hfStatsF12, Temp, pH)

eprPlotF12 <- ggplot(eprStatsF12, aes(Temp, Rate, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Rate-ci, ymax=Rate+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(day^-1)), x=xlab)+
  scale_x_continuous(breaks = unique(eprStatsF12$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F12 EPR')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

eprPlotF12

hfPlotF12 <- ggplot(hfStatsF12, aes(Temp, HF, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci, ymax=HF+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Hatching Frequency", x=xlab)+
  scale_x_continuous(breaks = unique(hfStatsF12$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F12 Hatching Frequency')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

hfPlotF12

#####F15#####
eprStatsF15 <- summarySE(data = eprDataF15, measurevar = "Rate", groupvars = "Treatment")
hfStatsF15 <- summarySE(data = eprDataF15, measurevar = "HF", groupvars = "Treatment")
eprStatsF15$Generation <- c(15,15,15,15)
hfStatsF15$Generation <- c(15,15,15,15)
eprStatsF15 <- cbind(eprStatsF15, Temp, pH)
hfStatsF15 <- cbind(hfStatsF15, Temp, pH)

eprPlotF15 <- ggplot(eprStatsF15, aes(Temp, Rate, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Rate-ci, ymax=Rate+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12),
        legend.text = element_text(size = 20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), x=xlab)+
  scale_x_continuous(breaks = unique(eprStatsF15$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F15 EPR')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

eprPlotF15

hfPlotF15 <- ggplot(hfStatsF15, aes(Temp, HF, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci, ymax=HF+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Hatching Frequency", x=xlab)+
  scale_x_continuous(breaks = unique(hfStatsF15$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F15 Hatching Frequency')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))
  

hfPlotF15


#####F25#####

eprStatsF25 <- summarySE(data = eprDataF25, measurevar = "Rate", groupvars = "Treatment")
hfStatsF25 <- summarySE(data = eprDataF25, measurevar = "HF", groupvars = "Treatment")
eprStatsF25$Generation <- 25
hfStatsF25$Generation <- 25
eprStatsF25$Temp <- c(18,18,22,22)
eprStatsF25$pH <- c(8.2, 7.5,8.2,7.5)
hfStatsF25$Temp <- c(18,18,22,22)
hfStatsF25$pH <- c(8.2, 7.5,8.2,7.5)
#eprStatsF25 <- cbind(eprStatsF25, Temp, pH)
#hfStatsF25 <- cbind(hfStatsF25, Temp, pH)

######All the data combined into plot over time#####
eprStatsTotal <- rbind(eprStatsF0, eprStatsF3, eprStatsF6, eprStatsF9, eprStatsF12, eprStatsF15, eprStatsF25)
hfStatsTotal <- rbind(hfStatsF0, hfStatsF3, hfStatsF6, hfStatsF9, hfStatsF12, hfStatsF15, hfStatsF25)

#get rid of 0 values so that it just runs through the data but still uses generation 0
eprStatsTotal[,1:6][eprStatsTotal[, 1:6]==0] <- NA
hfStatsTotal[,1:6][hfStatsTotal[, 1:6]==0] <- NA

#f25 reaction norm plots

eprPlotF25 <- ggplot(eprStatsF25, aes(Temp, Rate, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Rate-ci, ymax=Rate+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), x=xlab)+
  scale_x_continuous(breaks = unique(eprStatsF0$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F25 EPR')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

eprPlotF25

hfPlotF25 <- ggplot(hfStatsF25, aes(Temp, HF, color=factor(pH)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci, ymax=HF+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_discrete(name="pH")+
  theme_classic()+
  labs(y="Hatching Frequency", x=xlab)+
  scale_x_continuous(breaks = unique(hfStatsF0$Temp))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  ggtitle('F25 Hatching Frequency')+
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

hfPlotF25

eprPlotTotal <- ggplot(data = eprStatsTotal, aes(Generation, Rate, color=factor(Treatment)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Rate-ci, ymax=Rate+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y=expression(Egg~Production~Rate~(female^-1~day^-1)), x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))
  
eprPlotTotal


hfPlotTotal <- hfPlotTotal <- ggplot(data = hfStatsTotal, aes(Generation, HF, color=factor(Treatment)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=HF-ci, ymax=HF+ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Hatching Frequency", x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


hfPlotTotal  
  
epr.anova <- aov(Rate~Temp*pH*Generation, data = eprDataAll)
summary(epr.anova)

hf.anova <- aov(HF~Temp*pH*Generation, data = eprDataAll)
summary(hf.anova)


eprDataF0 <- eprDataF0[-grep('X', colnames(eprDataF0))]
eprF0F15 <- rbind(eprDataF0, eprDataF15)
eprF0F15$Generation <- factor(eprF0F15$Generation)
epr.pairwise <- pairwise.t.test(eprF0F15$Rate, eprF0F15$Generation:eprF0F15$Treatment)

epr.pairwise


F0F15anova <- aov(Rate~Temp*pH, data = eprF0F15)
summary(F0F15anova)

F0F15hatching.anova <- aov(HF~Temp*pH, data = eprF0F15)
summary(F0F15hatching.anova)

pairwise.t.test(eprF0F15$HF, eprF0F15$Generation:eprF0F15$Treatment)

#fileList <- list(eprStatsF0, eprStatsF3, eprStatsF6, eprStatsF9, eprStatsF12, eprStatsF15, 
                # hfStatsF0, hfStatsF3, hfStatsF6, hfStatsF9, hfStatsF12, hfStatsF15)

#library(dplyr)
#Map(cbind, fileList, Temp = temp, pH = pH)

#for(i in fileList) {
 # do.call(cbind, (i, Temp, pH))
#}

#final <- do.call(cbind, fileList)

#lapply(fileList, function(x) cbind(x, Temp, pH))
