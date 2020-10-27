setwd("C:/users/james/Documents/Grad_school/OA_Project")

library(ggplot2)

library(data.table)

library(dplyr)

library(car)

tonsa.physical.data <- fread(file = "OA_Physical_data.txt")




Atonsa.pH.mean <- tonsa.physical.data %>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(pH, na.rm = TRUE),
            sd = sd(pH, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)


Atonsa.temp.mean <- tonsa.physical.data %>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(Temperature, na.rm = TRUE),
            sd = sd(Temperature, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)


Atonsa.pCO2.mean <- tonsa.physical.data %>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(pCO2, na.rm = T),
            sd = sd(pCO2, na.rm = T),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)

#Atonsa.pH.mean <- summarySE(data = tonsa.physical.data, measurevar = c("pH"), groupvars = "Treatment")
#Atonsa.temp.mean <- summarySE(data = tonsa.physical.data, measurevar = c("Temperature"), groupvars = "Treatment")
#Atonsa.pco2.mean <- summarySE(data = tonsa.physical.data, measurevar = "pCO2", groupvars = "Treatment")


fwrite(x = Atonsa.temp.mean, file = "A_tonsa_temp_mean.txt", sep = "\t")
fwrite(x = Atonsa.pH.mean, file = "A_tonsa_pH_mean.txt", sep = "\t")
fwrite(x = Atonsa.pco2.mean, file = "A_tonsa_pCO2_mean.txt", sep = "\t")


tonsa.physical.data$Treatment <- as.factor(tonsa.physical.data$Treatment)
tonsa.physical.data$Date <- as.Date(tonsa.physical.data$Date, "%m/%d/%Y")
tonsa.physical.data$pCO2 <- as.numeric(tonsa.physical.data$pCO2)

tonsa.physical.data <- tonsa.physical.data[!grep("AAHH", tonsa.physical.data$Treatment),]
tonsa.physical.data <- tonsa.physical.data[!grep("HHAA", tonsa.physical.data$Treatment),]



co2.low <- tonsa.physical.data %>% 
  filter(pH > 8)

co2.low.comp <- lm(pH ~ Treatment, data = co2.low)
summary(co2.low.comp)
  
co2.high <- tonsa.physical.data %>% 
  filter(pH < 8)


co2.model <- glm(co2 ~ Treatment, data = tonsa.physical.data)



temp.low <- tonsa.physical.data %>% 
  filter(Temperature < 20)

temp.low.comp <- lm(Temperature ~ Treatment, data = temp.low)
summary(temp.low.comp)
Anova(temp.low.comp)


temp.high <- tonsa.physical.data %>% 
  filter(Temperature > 20)

temp.high.comp <- lm(Temperature ~ Treatment, data = temp.high)
summary(temp.high.comp)
Anova(temp.high.comp)



Atonsa.temp.plot <- ggplot(data = tonsa.physical.data, aes(Date, Temperature, color = factor(Treatment)))+
  geom_point(size = 1)+
  scale_y_continuous(name = "Temperature (°C)",
                     breaks = c(16,17,18,19,20,21,22,23), limits = c(16,23))+
  scale_x_date(date_labels = "%b/%Y")+
  scale_color_manual(name = NULL,
                     values = c("blue", "forestgreen", "orange", "red"),
                     labels = c("Control",
                                expression("High C"*O[2]),
                                "High Temp",
                                expression("High Temp High C"*O[2])))+
  theme_light()+
  theme(legend.position = "bottom",
        legend.text.align = 0)
  
  
  
Atonsa.temp.plot



Atonsa.pH.plot <- ggplot(data = tonsa.physical.data, aes(Date, pH, color = factor(Treatment)))+
  geom_point(size = 1)+
  scale_y_continuous(name = "pH")+
  scale_x_date(date_labels = "%b/%Y")+
  scale_color_manual(name = NULL,
                     values = c("blue", "forestgreen", "orange", "red"),
                     labels = c("Control",
                                expression("High C"*O[2]),
                                "High Temp",
                                expression("High Temp High C"*O[2])))+
  theme_light()+
  theme(legend.position = "bottom",
        legend.text.align = 0)
Atonsa.pH.plot


Atonsa.pCO2.plot <- ggplot(data = tonsa.physical.data, aes(Date, pCO2, color = factor(Treatment)))+
  geom_point(size = 1)+
  scale_y_continuous(name = "pCO2")+
  scale_x_date(date_labels = "%b/%Y")+
  scale_color_manual(name = NULL,
                     values = c("blue", "forestgreen", "orange", "red"),
                     labels = c("Control",
                                expression("High C"*O[2]),
                                "High Temp",
                                expression("High Temp High C"*O[2])))+
  theme_light()+
  theme(legend.position = "bottom",
        legend.text.align = 0)

Atonsa.pCO2.plot



Atonsa.pco2 <- fread("A_tonsa_alkalinity_data.txt")

pco2.mean <- Atonsa.pco2 %>% 
  group_by(Temp, pH) %>% 
  summarise(mean = mean(pCO2, na.rm = TRUE),
            sd = sd(pCO2, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)


TA.mean <- Atonsa.pco2 %>% 
  group_by(Exp, Temp, pH) %>% 
  summarise(mean = mean(TA, na.rm = TRUE),
            sd = sd(TA, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

sal.mean<- Atonsa.pco2 %>% 
  group_by(Exp, Temp, pH) %>% 
  summarise(mean = mean(sal, na.rm = TRUE),
            sd = sd(sal, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

fco2.mean <- Atonsa.pco2 %>% 
  group_by(Exp, Temp, pH) %>% 
  summarise(mean = mean(fCO2, na.rm = TRUE),
            sd = sd(fCO2, na.rm = TRUE),
            n.count = n()) %>% 
  mutate(se = sd/sqrt(n.count),
         lower.ci = mean - qt(1 - (0.05 / 2), n.count-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n.count-1)*se)

fwrite(pco2.mean, file = "Atonsa_pCO2_mean.txt", sep = "\t")
fwrite(TA.mean, file = "Atonsa_TA_mean.txt", sep = "\t")
fwrite(fco2.mean, file = "Atonsa_fCO2_mean.txt", sep = "\t")

  
co2.data <- fread("A_tonsa_alkalinity_data.txt")
co2.data.low <- co2.data %>% 
  filter(pH == 8.2)

co2.model.low <- lm(pCO2 ~ Treatment, data = co2.data.low)
summary(co2.model.low)
Anova(co2.model.low)


co2.data.high <- co2.data %>% 
  filter(pH == 7.5)

co2.model.high <- lm(pCO2 ~ Treatment, data = co2.data.high)
summary(co2.model.high)
Anova(co2.model.high)


##### Cost experiments #####


cost.dir <- "C:/Users/james/Documents/Grad_school/OA_Project/Cost_experiments/"

tonsa.cost.physical.data <- fread(file = paste(cost.dir,"Physical_data.txt",sep = ""))




Atonsa.pH.mean <- tonsa.cost.physical.data %>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(pH, na.rm = TRUE),
            sd = sd(pH, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)


Atonsa.temp.mean <- tonsa.cost.physical.data %>% 
  group_by(Treatment) %>% 
  summarise(mean = mean(Temperature, na.rm = TRUE),
            sd = sd(Temperature, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)


fwrite(Atonsa.pH.mean, file = paste(cost.dir,"Cost_pH_mean.txt", sep = ""), sep = "\t")
fwrite(Atonsa.temp.mean, file = paste(cost.dir,"Cost_temp_mean.txt", sep = ""), sep = "\t")




co2.low <- tonsa.cost.physical.data %>% 
  filter(pH > 8)

co2.low.comp <- lm(pH ~ Treatment, data = co2.low)
summary(co2.low.comp)
Anova(co2.low.comp)

co2.high <- tonsa.cost.physical.data %>% 
  filter(pH < 8)


co2.high.comp <- lm(pH ~ Treatment, data = co2.high)
summary(co2.high.comp)
Anova(co2.high.comp)


temp.low <- tonsa.cost.physical.data %>% 
  filter(Temperature < 20)

temp.low.comp <- lm(Temperature ~ Treatment, data = temp.low)
summary(temp.low.comp)
Anova(temp.low.comp)


temp.high <- tonsa.cost.physical.data %>% 
  filter(Temperature > 20)

temp.high.comp <- lm(Temperature ~ Treatment, data = temp.high)
summary(temp.high.comp)
Anova(temp.high.comp)



