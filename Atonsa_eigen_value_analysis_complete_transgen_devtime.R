#setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/SurvDataFiles/")

surv.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Survival/"
fit.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/Fitness/"
epr.directory <- "C:/Users/james/Documents/Grad_school/OA_Project/EPR/"

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


#LF = list.files(pattern = paste(surv.directory,"SurvDataFiles/*.080718.*", sep = ""))
#LF

#SurvData <- rbindlist(lapply(
#  setNames(LF,LF), #what you want to be "lapplied"
#  fread), #what you want the lapply to do
#  idcol = "source", #what the new column will be named (i.e. the header of the new column)
#  fill = TRUE) #fills missing columns with NA's

#SurvData <- separate(SurvData, "source", into = c("Date", "x", "Generation"))

SurvData <- fread(paste(surv.directory, "SurvDataFiles/Survival_data_total.txt", sep = ""))
# Create uniqe sorting column to group each technical and biological replicate

SurvData <- unite(SurvData,#the data frame
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


#setwd("C:/Users/james/Documents/Grad_school/OA_Project/EPR/Data_frames/")

EPR.data <- fread(paste(epr.directory,"Data_frames/EPR_HF_data_total.txt", sep = ""))


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

EPR.data <- unite(EPR.data,#the data frame
                  unique, #the name of the new column
                  c(Generation, Treatment, Rep), #the existing columns you want to combine
                  remove = FALSE) 

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





fwrite(lambda.results, file = paste(fit.directory,"lambda_results_devtime_surv_epr_hf_sex_standardized_relative.txt", sep = ""), sep = "\t")
lambda.results <- fread(file = paste(fit.directory,"lambda_results_devtime_surv_epr_hf_sex_standardized_relative.txt", sep = ""))
lambda.results[is.na(lambda.results$lambda.rel)] <- 0
#lambda.results.GH <- filter(lambda.results, Treatment == 4)
#lambda.results.GH.0 <- filter(lambda.results.GH, Generation == 0)
#hist(lambda.results.GH.0$lambda)


# calculate malthusian parameter
lambda.results$Malthusian <- log10(lambda.results$lambda)



library(lm.beta)


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

HH9 <- lambda.results %>% 
  filter(Treatment == 4 & Generation == 9)


summary(lm(lambda.rel ~ epr, data = HH9))
summary(lm(lambda.rel ~ surv, data = HH9))
summary(lm(lambda.rel ~ hf, data = HH9))


# remove the first row as a last step
lambda.coef <- lambda.coef[-1,] 


fwrite(lambda.coef, file = paste(fit.directory,"Statistics/lambda_regression_coefficients.txt", sep = ""), sep = "\t")

#summary(lm(lambda.rel~epr, data = lambda.results.0.full.HH))
#lm(lambda.rel~epr, data = lambda.results.last.full.HH)



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

fwrite(f1.coef, file = paste(surv.directory,"Statistics/survival_regression_coefficients.txt", sep = ""), sep = "\t")

f2 <- lm(lambda.stand~epr*factor(Treatment)*factor(Generation), data = lambda.results)
f2.emm <- emmeans(f2,  ~ Generation | Treatment)

f2.emm

f2.coef <- tidy(coef(f2))
View(f2.coef)

fwrite(f2.coef, file = paste(epr.directory,"Statistics/epr_regression_coefficients.txt", sep = ""), sep = "\t")



##### PATH analysis ######

# help from: https://rpubs.com/tbihansk/302732

mod <- 'lambda.rel ~ surv + epr + hf + dev.time + sex'

library(mnormt)
library(lavaan)
path.fit <- cfa(mod, data = lambda.results)

summary(path.fit)


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


#Anova(path.list2$`0_4`, path.list2$`25_4`)

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


lambda.results <- fread(paste(fit.directory,"lambda_results_devtime.txt", sep = ""))

lambda.results <- as.data.frame(lambda.results)


# create continuous generation vector for anova and plotting
lambda.results$Generation.c <- as.numeric(as.character(lambda.results$Generation)) 


lambda.results$Generation <- as.factor(as.numeric(lambda.results$Generation))
lambda.results$Treatment <- as.factor(lambda.results$Treatment)
lambda.results$Rep.c <- as.numeric(lambda.results$Rep)
lambda.results$Rep.c <- as.factor(as.numeric(lambda.results$Rep.c))

# check to see what the distribution looks like
hist(lambda.results$lambda, freq = F)

lambda.ls <- split(lambda.results, f = lambda.results$Treatment)

lambda.hist <- lapply(lambda.ls, function (x) hist(x$lambda))




# use a continuous generation model for the anova

malthusian.model <- lm(lambda~Generation.c*Treatment, data=lambda.results)


m <- Anova(malthusian.model)
m

m$factors <- rownames(m)
fwrite(m, file = paste(fit.directory,"Statistics/lambda_devtime_Anova.txt", sep = ""), sep = "\t")

shapiro.test(malthusian.model$residuals)

hist(malthusian.model$residuals)


#leveneTest(malthusian.model)

foo <- plot_model(malthusian.model, 'int')

foo

# Pairwise comparison of malthusian parameters

# create the model to use for Tukey pairwise comparisons -- generation must be a factor, not a continuous variable
malthusian.model.2 <- lm(lambda~Generation*Treatment, data=lambda.results)

emmip(malthusian.model.2, Treatment ~ Generation)

emm_1 <- emmeans(malthusian.model.2, pairwise ~ Generation | Treatment)

p <- pairs(emm_1)
p

# model for writing the tukey results
malthusian.model.3 <- aov(lambda~unique, data = lambda.results)
p2 <- TukeyHSD(malthusian.model.3, "unique")
p2.tukey <- tidy(p2)

p3 <- tidy(pairwise.t.test(lambda.results$lambda, lambda.results$Generation:lambda.results$Treatment, p.adjust.method = "none"))



lambda.results$Temp <- case_when(lambda.results$Treatment == "1" ~ 18,
                                 lambda.results$Treatment == "2" ~ 18,
                                 lambda.results$Treatment == "3" ~ 22,
                                 lambda.results$Treatment == "4" ~ 22)


lambda.results$pH <- case_when(lambda.results$Treatment == "1" ~ 8.2,
                                 lambda.results$Treatment == "2" ~ 7.5,
                                 lambda.results$Treatment == "3" ~ 8.2,
                                 lambda.results$Treatment == "4" ~ 7.5)


malthusian.model.5 <- aov(lambda ~ Generation.c*as.factor(Temp)*as.factor(pH), data = lambda.results)

m5 <- Anova(malthusian.model.5)
m5
m5$factors <- rownames(m5)



fwrite(lambda.results, file = paste(fit.directory,"lambda_results_devtime.txt", sep = ""), sep = "\t")
fwrite(p2.tukey, file = paste(fit.directory,"Statistics/lambda_devtime_Tukey.txt", sep = ""), sep = "\t")
fwrite(m5, file = paste(fit.directory,"Statistics/lambda_3way_anova.txt", sep = ""), sep = "\t")
fwrite(p3, file =paste(fit.directory,"Statistics/lambda_pairwise.txt", sep = ""), sep = "\t")




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

library(glmmTMB)

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
library(lme4)
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




# test if the replicates affect the lambda
#malthusian.model.4 <- aov(lambda~Generation.c*Treat.Rep, data = lambda.results)
#m4 <- Anova(malthusian.model.4)
#m4



## how to do pairwise comparisons with zero-inflated models
## help from: https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html

## and https://rcompanion.org/handbook/J_01.html


## can't use a poisson model for non-integer data https://github.com/atahk/pscl/issues/5
#library(pscl)
#lambda.results$Treatment.c <- as.numeric(lambda.results$Treatment)
#lambda.results$lambda.i <- as.integer(lambda.results$lambda)

#m.zi <- zeroinfl(lambda ~ Generation.c*Treatment,
#                 data = lambda.results,
#                 offset = log(lambda),
#

#dist = "poisson")
#summary(m.zi)

## try to do a pairwise comparison with these results

#fit_zigauss.2 <- glmmTMB(lambda~Generation*Treatment + (1|Rep.c),
#                         data = lambda.results,
#                         ziformula = ~.)

#emmip(fit_zigauss.2, Treatment ~ Generation)

#emm_2 <- emmeans(fit_zigauss.2, pairwise ~ Generation | Treatment)

#p2 <- pairs(emm_2)
#p2


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
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="Population Fitness (??)", x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  #ylim(1.402, 1.405)+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

lambdaPlotTotal


#setwd("C:/Users/james/Documents/Grad_school/OA_Project/Fitness/")
fwrite(lambda.mean.c, file = paste(fit.directory,"lambda_mean_devtime.txt", sep = ""), sep = "\t")
#lambda.results$Malthusian <- log10(lambda.results$lambda)

# plot to see how data is distributed among treatments and reps
ggplot(data = lambda.results, aes(x = Generation.c,
                                  y = lambda,
                                  color = factor(Rep)))+
  geom_smooth(method = "lm", se = F)+
  
  geom_jitter(width = 0.7, alpha = 0.2)+
  
  geom_point()+
  
  facet_wrap(~Treatment)

# plot how many zeroes change across the experiment and how many positive lambda values change
lambda.results.zero <- lambda.results %>% 
  group_by(Generation.c, Treatment, Rep) %>% 
  filter(lambda == 0) %>% 
  summarise(zero_n = n())


ggplot(data = lambda.results.zero, aes(x = Generation.c,
                                       y = zero_n,
                                       color = factor(Rep)))+
  geom_smooth(method = "lm", se = F)+
  
  #geom_jitter(width = 0.7, alpha = 0.2)+
  
  geom_point()+
  
  facet_wrap(~Treatment)
  
## plot to see how many non-zero values for each replicate
lambda.results.pos <- lambda.results %>% 
  group_by(Generation.c, Treatment, Rep) %>% 
  filter(lambda > 0) %>% 
  summarise(zero_n = n())

ggplot(data = lambda.results.pos, aes(x = Generation.c,
                                       y = zero_n,
                                       color = factor(Rep)))+
  geom_smooth(method = "lm", se = F)+
  
  #geom_jitter(width = 0.7, alpha = 0.2)+
  
  geom_point()+
  
  facet_wrap(~Treatment)


Malthusian.mean <- lambda.results %>%
  group_by(Generation, Treatment) %>%
  summarise(mean = mean(Malthusian, na.rm = TRUE),
            sd = sd(Malthusian, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)





## plots of traits vs. lambda

setwd("C:/Users/james/Documents/Grad_school/OA_Project/Fitness/")

lambda.results.traits <- fread(paste(fit.directory,"lambda_results_devtime_surv_epr.txt", sep = ""))

epr.lambda.plot <- ggplot(data = lambda.results.traits, aes(x = epr,
                                                            y = lambda
                                                            )) +
  geom_smooth(method = "lm", se = F) +
  
  #geom_point()+
  
  facet_wrap(~Treatment)
  
epr.lambda.plot


surv.lambda.plot <- ggplot(data = lambda.results.traits, aes(x = surv,
                                                             y = lambda,
                                                             color = factor(Generation))) +
  geom_smooth( se = F) +
  
  ylim(0,2)+
  #geom_point()+
  
  facet_wrap(~Treatment)

surv.lambda.plot






lambdaBoxplot <- ggplot(data = lambda.results, aes(Generation, lambda))+
  geom_boxplot(lwd = 1.1, aes(fill = factor(Treatment)),
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y="lambda", x="Generation")+
  scale_x_discrete(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  #ylim(1.402, 1.405)+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


lambdaBoxplot






MalthusianPlotTotal <- ggplot(data = Malthusian.mean, aes(Generation, mean, color=factor(Treatment)))+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), size=1, width=1)+
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y=expression("Malthusian Parameter "*(generation^-1)), x="Generation")+
  scale_x_continuous(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  #ylim(1.402, 1.405)+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))

MalthusianPlotTotal







MalthusianBoxplot <- ggplot(data = lambda.results, aes(Generation, Malthusian))+
  geom_boxplot(lwd = 1.1, aes(fill = factor(Treatment)),
               outlier.size = 1.5)+ # can't make a continuous x axis with boxplots
  
  theme(legend.title = element_text(colour = "black", size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("blue", "forestgreen", "orange", "red"))+
  theme_classic()+
  labs(y=expression("Malthusian Parameter "*(generation^-1)), x="Generation")+
  scale_x_discrete(breaks = c(0,3,6,9,12,15,25))+
  scale_y_continuous(breaks = waiver(), minor_breaks = waiver())+
  #ylim(1.402, 1.405)+
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))


MalthusianBoxplot




################################################################################################################################################
##### Old analysis methods #####################






#malthusian.pairwise <- pairwise.t.test(lambda.results$Malthusian, lambda.results$unique)

#malthusian.pairwise


#malthusian.pairwise.table <- lambda.results %>%
#  group_by(Generation, Treatment) %>%
#  do(tidy(pairwise.t.test(lambda.results$Malthusian, lambda.results$unique)))

#malthusian.pairwise.table.adj <- filter(malthusian.pairwise.table, p.value < 0.05)




# Two-way ANOVA
malthusian.model <- aov(Malthusian~Generation*Treatment, data=lambda.results)

summary(malthusian.model)


malthusian.model.table <- lambda.results %>%
  group_by(Generation, Treatment) %>%
  do(tidy(aov(Malthusian~Generation*Treatment, data=lambda.results)))



## find the r2 values
Models <- lambda.results %>%
  group_by(Treatment) %>%
  do(model = lm(Malthusian~Generation, data=.))


model.results <- Models %>% glance(model)




lambda.lists <- split(lambda.results, lambda.results$Treatment) #create different dataframes based on the unique feature
lambda.lists


lambda.model.list <- lapply(lambda.lists,
                           function(x) lm(Malthusian~as.numeric(Generation), data = x))



lambda.slope.list <- lapply(lambda.model.list,
                            function(x) coef(x)["as.numeric(Generation)"])


fwrite(lambda.results, file = "Lambda_results_extended_method_dev_time.txt", sep = "\t")
fwrite(malthusian.pairwise.table, file = "Malthusian_pairwise_table_extended.txt", sep = "\t")
fwrite(malthusian.pairwise.table.adj, file = "Malthusian_pairwise_table_extended_adj.txt", sep = "\t")
fwrite(malthusian.model.table, file = "Malthusian_anova_table_extended.txt", sep = "\t")


###### Malthusian Point Plot #####



setwd("C:/Users/james/Documents/Grad_school/OA_Project/Survival/Cost_experiment/")
lambda.results <- fread("Malthusian_results_all_reps.txt")






lambda.list.mean <- lambda.results %>%
  group_by(Treatment, Food) %>%
  summarise(mean = mean(Malthusian, na.rm = TRUE),
            sd = sd(Malthusian, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)