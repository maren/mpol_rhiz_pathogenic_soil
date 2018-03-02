# FLM018: Rhizobia confer protection against pathogenic effect of soil microbiome in M. polymorpha

# ====== Project Description ========
# This is a project initiated by Katie Wozniak where she grew 12 genotypes of Medicago polymorpha in the presence of soil inoculates and/or rhizobia. The initial project (Trapping experiment) used soil from 8 different locations in either Florida or Portugal. The second part of the experiment, Katie used only two MP genotypes and location matching soil at low and high concentrations. In both experiments, WSM was the rhizobia used for inoculations.

# ===== Packages =====
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(car)
library(multcomp)
library(data.table)
library(lme4)
library(lmerTest)

# ==== Data prepping =====
setwd("~/Documents/Projects/Rhizobia_Pathogenic_Soil/results/processed_data/")

# I'm first trying to untangle all of the data and put it back in a logical order. I know that the Soil_RhizMort dataset has all of the treatments, I will use that to populate with all tissue data. Anything missing will be supplemented by checking the scans.

SRMort <- read.csv("SoilRhizMort.csv")

# First will start with Root data
SRroot <- read.csv("RootData.csv")

# Will match up data using plant pot id
SRMort$RootWeight <- SRroot[match(SRMort$Pot, SRroot$Pot), "Weight"]

# Do the same thing for the shoot and nodule data
SRShoot <- read.csv("ShootData.csv")
SRNod <- read.csv("NoduleData.csv")

SRMort$ShootWeight <- SRShoot[match(SRMort$Pot, SRShoot$Pot), "Weight"]
SRMort$NodWeight <- SRNod[match(SRMort$Pot, SRNod$Pot), "NoduleWeight"]
SRMort$NodNum <- SRNod[match(SRMort$Pot, SRNod$Pot), "NoduleNum"]
SRMort$PodNum <- SRNod[match(SRMort$Pot, SRNod$Pot), "PodNum"]

SRMort <- within(SRMort, Shoot_Root <- ShootWeight / RootWeight)
# Saving this dataset:
write.table(SRMort, file = "Soil_Rhiz_All_Data20Feb2018.csv", row = F, sep = ",")


# Now doing the same thing for the Trapping experiment
TrapMort <- read.csv("TrappingMortalityFeb12018.csv")
TR <- read.csv("TrapRoot.csv")
TS <- read.csv("TrapShoot.csv")
TN <- read.csv("TrapNod.csv")

TrapMort$ShootWeight <- TS[match(TrapMort$Pot.., TS$Pot_no), "final_weight..g."]
TrapMort$RootWeight <- TR[match(TrapMort$Pot.., TR$Pot_no), "final_weight..g."]
TrapMort$NodWeight <- TN[match(TrapMort$Pot.., TN$Pot_no), "final_weight..g."]
TrapMort$NodNum <- TN[match(TrapMort$Pot.., TN$Pot_no), "nod_num"]

# Adding Range Data and Shoot: Root Ratio
TrapMort <- within(TrapMort, {
  Range <- if_else(Location %in% c("Australia" , "California" , "Chile" , "Florida"), "Invaded", "Native")
  Shoot_Root <- ShootWeight / RootWeight
})

write.table(TrapMort, "Trapping_Alll_Data20Feb2018.csv", row = F, sep = ",")

# ====== Beginning actual data analysis and figures for results sections =====

# Trapping Exp: Mortality ====
MortMod <- glmer(as.factor(Dead_0) ~ Treat * Range + (1 | Genotype), family = binomial(link = "logit"), data = TrapMort)
lmerTest::anova(MortMod, test = "Chisq") # Doesn't work
car::Anova(MortMod) # works


MortMod2 <- glmer(as.factor(Dead_0) ~ Treatment * Range + (1 | Genotype), family = binomial(link = "logit"), data = TrapMort)
car::Anova(MortMod2) # This shows significance but it runs out of DFs...

MM <- glm(as.factor(Dead_0) ~ Treatment * Range, family = binomial(link = "logit"), data = TrapMort)
anova(MM, test = "Chisq")

MortMod3 <- glmer(as.factor(Dead_0) ~ Treatment  + (1 | Genotype), family = binomial(link = "logit"), data = TrapMort)

car::Anova(MortMod3) # works

ggplot(TrapMort, aes(x = Treatment, y = Dead_0, colour = Treat)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Survival Rate", title = "Trapping Experiment Mortality 21 Feb 2018")
ggsave(file = "SurvivalRateTrapExp21Feb2018.pdf")

summary(glht(MortMod3, mcp(Treatment = "Tukey"))) # NS between Soil treatments, only soil and buffer/wsm

ddply(TrapMort, "Treatment", summarise, Mean = mean(Dead_0))
ddply(TrapMort, "Treat", summarise, Mean = mean(Dead_0))

# Trapping Exp: Biomass Data ====

# Manually fixed some errors in the csv file
TrapMort <- read.csv("Trapping_Alll_Data20Feb2018.csv")
# Only the buffer and WSM plants were harvested
TrapBio <- subset(TrapMort, Treat != "Soil")

# Shoot Data
ShootMod <- lmer(ShootWeight ~ Treatment * Range + (1| Genotype), data = TrapBio)
anova(ShootMod)
ggplot(TrapBio, aes(x = Treatment, y = ShootWeight, colour = Range)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Shoot Weight (mg)", title = "Trapping Experiment Shoot Weight 21 Feb 2018") 
ggsave(file = "TrappingShootWeight21Feb2018.pdf")

# Root Data
RootMod <- lmer(RootWeight ~ Treatment * Range + (1| Genotype), data = TrapBio)
anova(RootMod)
ggplot(TrapBio, aes(x = Treatment, y = RootWeight, colour = Range)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Root Weight (mg)", title = "Trapping Experiment Root Weight 21 Feb 2018") 
ggsave(file = "TrappingRootWeight21Feb2018.pdf")

# Nodule Data
#   Number
NNMod <- lmer(NodNum ~ Range + (1| Genotype), data = TrapBio)
anova(NNMod)
ggplot(TrapBio, aes(x = Range, y = NodNum)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Nodule Number", title = "Trapping Experiment Nodule Number 21 Feb 2018")  + coord_cartesian(ylim = c(0,30))
ggsave(file = "TrappingNoduleNumber21Feb2018.pdf") 

#   Weight
NWMod <- lmer(NodWeight ~ Range + (1| Genotype), data = TrapBio)
anova(NWMod)
ggplot(TrapBio, aes(x = Range, y = NodWeight)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Nodule Weight (mg)", title = "Trapping Experiment Nodule Weight 21 Feb 2018")  + coord_cartesian(ylim = c(0,.0035))
ggsave(file = "TrappingWeight21Feb2018.pdf")


# Soil x Rhiz Exp: Mortality ====
SRMort$Treatment <- as.factor(SRMort$Treatment)
SRMortMod <- glm(DeadOrAlive ~ Treatment , family = binomial(link = "logit"), data = SRMort)
anova(c, test = "Chisq")
summary(glht(SRMortMod, mcp(Treatment = "Tukey"))) # NS between Soil treatments, only soil and buffer/wsm


SRMort$Survival <- ifelse(SRMort$DeadOrAlive == "Alive", 1,0)
ggplot(SRMort, aes(x = Treatment, y = Survival, colour = Inoculate, shape = SoilLocation)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Survival Rate", title = "SxR Experiment Mortality 21 Feb 2018") + facet_grid(. ~ SoilConc, scales = "free_x")

ggplot(SRMort, aes(x = SoilConc, y = Survival, colour = Inoculate)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Soil Concentration", y = "Survival Rate", title = "SxR Experiment Mortality 21 Feb 2018") #+ facet_grid(. ~ Match, scales = "free_x")
ggsave(file = "SurvivalRateSRExperiment21Feb2018.pdf")


ggplot(SRMort, aes(x = Soil, y = Survival)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Soil Concentration", y = "Survival Rate", title = "SxR Experiment Mortality 21 Feb 2018") + facet_grid(. ~ Range, scales = "free_x")

ggplot(SRMort, aes(x = Soil, y = Survival, colour = SoilLocation)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Soil Concentration", y = "Survival Rate", title = "SxR Experiment Mortality 21 Feb 2018") + facet_grid(. ~ Range, scales = "free_x")

SRMortMod2 <- glm(DeadOrAlive ~ SoilConc * Inoculate , family = binomial(link = "logit"), data = SRMort)
anova(SRMortMod2, test = "Chisq")

SRMort <- within(SRMort, Soil <- as.factor(paste(SoilConc, Inoculate)))

SRMortMod3 <- glm(DeadOrAlive ~ Soil * Range * SoilLocation, family = binomial(link = "logit"), data = SRMort)
anova(SRMortMod3, test = "Chisq")

SRMort <- within(SRMort, SoilRSL_Factor <- as.factor(paste(Soil, Range, SoilLocation)))
levels(SRMort$SoilRSL_Factor)

SRSLContLevels <- as.data.frame(levels(SRMort$SoilRSL_Factor))
write.table(SRSLContLevels, file = "ContrastLevelsSoil_Range_SLMort23Feb2018.csv", row = F, sep = ",")

# Lets do some contrasts!
SRMortMod4 <- glm(DeadOrAlive ~ SoilRSL_Factor, family = binomial(link = "logit"), data = SRMort)
anova(SRMortMod4, test = "Chisq")

tt = lsmeans(SRMortMod3, specs = ~ SoilLocation | Soil:Range)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(dd, by = "Range")

tt = lsmeans(SRMortMod3, specs = ~ Range | Soil:SoilLocation)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(dd, by = "Soil")

SoilContLevels <- as.data.frame(unique(SRMort$Soil))
write.table(SoilContLevels, file = "ContrastLevelsSoilMort21Feb2018.csv", row = F, sep = ",")

SoilMortCont <- read.csv("ContrastLevelsSoilMort21Feb2018.csv", header = F)
H.BvW <- SoilMortCont$V1
L.BvW <- SoilMortCont$V2
N.BvW <- SoilMortCont$V3

mat.temp <- rbind(constant = 1/6, H.BvW, L.BvW, N.BvW)
mat <- ginv(mat.temp)
mat <- mat[,-1]
H.BvW <- mat[,1]
L.BvW <- mat[,2]
N.BvW <- mat[,3]

MyContMat <- rbind("High: Buffer vs WSM" = H.BvW, "Low: Buffer vs WSM" = L.BvW, "None: Buffer vs WSM" = N.BvW)
colnames(MyContMat) = levels(SRMort$Soil)
Mod_GLHT <- glht(SRMortMod3, linfct = mcp(Soil = MyContMat))
summary(Mod_GLHT, test = adjusted("fdr"))
summary(Mod_GLHT)

# SR Biomass Data =====
#   Shoot Data
ShootModSR <- lm(ShootWeight ~ Treatment * Range, data = SRMort)
anova(ShootModSR) #  Treatment and Range are S (p < 0.001; p = 0.003)

ShootModSR <- lm(ShootWeight ~ Soil * Range, data = SRMort)
anova(ShootModSR)

ShootModSR <- lm(ShootWeight ~ Soil , data = SRMort)
anova(ShootModSR)

summary(glht(ShootModSR, mcp(Soil = "Tukey")))

ggplot(SRMort, aes(x = Soil, y = ShootWeight, colour = Range)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Shoot Weight (mg)", title = "S x R Experiment Shoot Weight 21 Feb 2018") 

ggplot(SRMort, aes(x = Soil, y = ShootWeight)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Shoot Weight (mg)", title = "S x R Experiment Shoot Weight 21 Feb 2018") + facet_grid(. ~ Range, scales = "free_x")
ggsave(file = "SR_ShootWeight21Feb2018.pdf")

SRMort_N <- subset(SRMort, Range == " Native")
ShootModSRs <- lm(ShootWeight ~ Soil , data = SRMort_N)
anova(ShootModSRs)
summary(glht(ShootModSRs, mcp(Soil = "Tukey")))


SRMort_I <- subset(SRMort, Range == "Invasive")
ShootModSRs <- lm(ShootWeight ~ Soil , data = SRMort_I)
Anova(ShootModSRs, type = "3")
summary(glht(ShootModSRs, mcp(Soil = "Tukey")))

# SR Biomass Data =====
#   Root Data
RootModSR <- lm(RootWeight ~ Treatment * Range, data = SRMort)
anova(RootModSR) #  All S (p < 0.001)

RootModSR <- lm(RootWeight ~ Soil * Range, data = SRMort)
anova(RootModSR)

RootModSR <- lm(RootWeight ~ Soil * Match, data = SRMort)
anova(RootModSR)

ggplot(SRMort, aes(x = Soil, y = RootWeight, colour = Range)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Root Weight (mg)", title = "S x R Experiment root Weight 21 Feb 2018") 

ggplot(SRMort, aes(x = Soil, y = RootWeight)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Root Weight (mg)", title = "S x R Experiment Root Weight 21 Feb 2018") + facet_grid(. ~ Range, scales = "free_x")


RootModSRs <- lm(RootWeight ~ Soil , data = SRMort_N)
anova(RootModSRs)
summary(glht(RootModSRs, mcp(Soil = "Tukey")))


RootModSRs <- lm(RootWeight ~ Soil , data = SRMort_I)
Anova(ShootModSRs, type = "3")
summary(glht(ShootModSRs, mcp(Soil = "Tukey")))

#   Nodule Number Data
NNModSR <- lm(NodNum ~ Treatment * Range, data = SRMort)
anova(NNModSR) #  All S (p < 0.001; p < 0.01)

NNModSR <- lm(NodNum ~ Soil * Range, data = SRMort)
anova(NNModSR) #  All S (p < 0.001; p < 0.01)

ggplot(SRMort, aes(x = Soil, y = NodNum, colour = Range)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Nodule Number", title = "S x R Experiment Nodule Number 21 Feb 2018") 

ggplot(SRMort, aes(x = Soil, y = NodNum)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Nodule Number", title = "S x R Experiment Nodule Number 21 Feb 2018") + facet_grid(. ~ Range, scales = "free_x")

NNModSRs <- lm(NodNum ~ Soil , data = SRMort_N)
anova(NNModSRs)
summary(glht(NNModSRs, mcp(Soil = "Tukey")))

NNModSRs <- lm(subset(NodNum, S ~ Soil , data = SRMort_I))
Anova(NNModSRs, type = "3")
summary(glht(NNModSRs, mcp(Soil = "Tukey")))





NWModSR <- lm(NodWeight ~ Soil * Range * SoilLocation, data = SRMort)
Anova(NWModSR, type = "2") 



mod <- lm(ShootWeight ~ Soil * Range, data = SRMort)
tt = lsmeans(mod, specs = ~ Range | Soil)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(dd, by = "Soil")

mod <- glm(DeadOrAlive ~ SoilConc * Inoculate, family = binomial(link = "logit"), data = SRMort)
tt = lsmeans(mod, specs = ~ SoilConc | Inoculate)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(dd, by = "Inoculate")


 SRMort_no_none <- subset(SRMort, SoilConc != "none")
mod <- glm(DeadOrAlive ~ SoilConc * Inoculate, family = binomial(link = "logit"), data = SRMort_no_none)

anova(mod, test = "Chisq")

mod <- glm(DeadOrAlive ~ Soil, family = binomial(link = "logit"), data = SRMort_no_none)
summary(glht(mod, mcp(Soil = "Tukey")))

tt = lsmeans(mod, specs = ~ Soil)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(dd)



# 26 Feb 2018
# testing the impact of soil on rhizobia
# SRMort_no_none has soil only removed

# Mortality
mod <- glm(DeadOrAlive ~ Soil, family = binomial(link = "logit"), data = SRMort_no_none)
anova(mod, test = "Chisq")
summary(glht(mod, mcp(Soil = "Tukey")))

mod <- glm(DeadOrAlive ~ SoilConc * Inoculate, family = binomial(link = "logit"), data = SRMort_no_none)
anova(mod, test = "Chisq")

#SHoot - will continue with this b/c of line 320
mod <- lm(ShootWeight ~ SoilConc * Inoculate, data = SRMort_no_none)
anova(mod)
ggplot(SRMort_no_none, aes(x = Inoculate, y = ShootWeight, colour = SoilConc)) + stat_summary(fun.data = "mean_se")

mod <- lm(ShootWeight ~ Soil, data = SRMort_no_none)
anova(mod)
ggplot(SRMort_no_none, aes(x = Soil, y = ShootWeight)) + stat_summary(fun.data = "mean_se")

# All high buffer plants died!

SRMortSub2 <- subset(SRMort_no_none, Inoculate == "WSM")
SRMortSub3 <- subset(SRMort, Inoculate == "WSM")

#Mortality
mod <- glm(DeadOrAlive ~ SoilConc * SoilLocation * Range, family = binomial(link = "logit"), data = SRMortSub2)
anova(mod, test = "Chisq")

tt = lsmeans(mod, specs = ~ SoilLocation | SoilConc:Range)
dd = pairs(tt)
summary(dd, by = NULL)

tt = lsmeans(mod, specs = ~ SoilConc | SoilLocation:Range)
dd = pairs(tt)
summary(dd, by = NULL)

tt = lsmeans(mod, specs = ~ Range | SoilLocation:SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)



ggplot(SRMortSub2, aes(x = SoilConc, y = Survival, colour = SoilConc)) + stat_summary(fun.data = "mean_se") + coord_cartesian(ylim = c(0,1)) +theme_bw()

mod <- glm(DeadOrAlive ~ SoilConc * SoilLocation * Range, family = binomial(link = "logit"), data = SRMortSub3)
anova(mod, test = "Chisq")

ggplot(SRMortSub3, aes(x = SoilConc, y = Survival, colour = SoilConc)) + stat_summary(fun.data = "mean_se") + coord_cartesian(ylim = c(0,1)) +theme_bw()
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = lsmeans(mod, specs = ~ SoilLocation | SoilConc:Range)
dd = pairs(tt)
summary(dd, by = NULL)

tt = lsmeans(mod, specs = ~ SoilConc | SoilLocation:Range)
dd = pairs(tt)
summary(dd, by = NULL)

tt = lsmeans(mod, specs = ~ Range | SoilLocation:SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)



# I you just compare high vs low WSM, there is a significant difference for survival. If no soil is included, then it becomes marginal due to ns difference from none treatment

# Shoot weight

mod <- lm(ShootWeight ~ SoilConc, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = ShootWeight, colour = SoilConc)) + stat_summary(fun.data = "mean_se")


mod <- lm(ShootWeight ~ SoilConc, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = ShootWeight, colour = SoilConc)) + stat_summary(fun.data = "mean_se")
summary(glht(mod, mcp(SoilConc = "Tukey")))

# No sig difference for either subset for shoot weight


# Root weight

mod <- lm(RootWeight ~ SoilConc, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = RootWeight, colour = SoilConc)) + stat_summary(fun.data = "mean_se")


mod <- lm(RootWeight ~ SoilConc, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = RootWeight, colour = SoilConc)) + stat_summary(fun.data = "mean_se")
summary(glht(mod, mcp(SoilConc = "Tukey")))

# No sig diff for either subset, but there is a trend for increased root mass with decreasing soil


# Root:Shoot
SRMort <- within(SRMort, Root_Shoot <- 1 / Shoot_Root)

mod <- lm(Root_Shoot ~ SoilConc, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = Root_Shoot, colour = SoilConc)) + stat_summary(fun.data = "mean_se")


mod <- lm(Root_Shoot ~ SoilConc, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = Root_Shoot, colour = SoilConc)) + stat_summary(fun.data = "mean_se")
summary(glht(mod, mcp(SoilConc = "Tukey")))

# no sig diff for first subset, but second is significant ( p = 0.043). There is no diff between H-L, sig between N-L, marg between N-H.

# Nod number

mod <- lm(NodNum ~ SoilConc, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = NodNum, colour = SoilConc)) + stat_summary(fun.data = "mean_se")


mod <- lm(NodNum ~ SoilConc, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = NodNum, colour = SoilConc)) + stat_summary(fun.data = "mean_se")
summary(glht(mod, mcp(SoilConc = "Tukey")))

# No sig diff for either subset, but there is a trend for increased nod number with decreasing soil

# Nod Weight

mod <- lm(NodWeight ~ SoilConc, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = NodWeight, colour = SoilConc)) + stat_summary(fun.data = "mean_se")


mod <- lm(NodWeight ~ SoilConc, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = NodWeight, colour = SoilConc)) + stat_summary(fun.data = "mean_se")
summary(glht(mod, mcp(SoilConc = "Tukey")))

# No sig diff for either subset, but there is a trend for increased nod weight with decreasing soil


# Now running models to see if there is a difference between genotypes
#Mortality
mod <- glm(DeadOrAlive ~ SoilConc * Genotype, family = binomial(link = "logit"), data = SRMortSub2)
anova(mod, test = "Chisq")

ggplot(SRMortSub2, aes(x = SoilConc, y = Survival, colour = Genotype)) + stat_summary(fun.data = "mean_se") + coord_cartesian(ylim = c(0,1)) +theme_bw() + facet_wrap(~Genotype)

mod <- glm(DeadOrAlive ~ SoilConc * Genotype, family = binomial(link = "logit"), data = SRMortSub3)
anova(mod, test = "Chisq")

ggplot(SRMortSub3, aes(x = SoilConc, y = Survival, colour = Genotype)) + stat_summary(fun.data = "mean_se") + coord_cartesian(ylim = c(0,1)) +theme_bw() + facet_wrap(~ Genotype)


# Genotype NS for mortality for either

# Shoot weight

mod <- lm(ShootWeight ~ SoilConc* Genotype, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = ShootWeight, colour = Genotype)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)


mod <- lm(ShootWeight ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = ShootWeight, colour = Genotype)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)

tt = emmeans(mod,  ~ SoilConc|Genotype)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(tt, simple = "Genotype")

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(tt, simple = "Genotype")

# Using emmeans, there is a difference between ST. Aug and PI493292 at low WSM levels.




# Root weight

mod <- lm(RootWeight ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = RootWeight, colour = Genotype)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)


mod <- lm(RootWeight ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = RootWeight, colour = Genotype)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)

tt = emmeans(mod,  ~ SoilConc|Genotype)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(tt, simple = "Genotype")

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)
# There is a sig diff in first subset between soil conc, genotype and interaction. In the second, the interaction becomes marginal. Then ran emmeans and found that in  Native Geno: H-L and H-N are sig. None in Invasive. But the two genotypes vary sig within low and none for root weight.


# Root:Shoot

mod <- lm(Root_Shoot ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = Root_Shoot, colour = Genotype)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)


mod <- lm(Root_Shoot ~ SoilConc *Genotype, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = Root_Shoot, colour = Genotype)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)

tt = emmeans(mod,  ~ SoilConc|Genotype)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(tt, simple = "Genotype")

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

# In the first subset, only genotype is s ( p = 0.025). In the second, soil and genotype but not interaction are s (p = 0.022 for both). Using emmeans, in native sig dif between H-N and L-N (p = 0.0216, 0.0208) non for invasive. Between genotypes, only none is sig ( p= 0.0167)

# Nod number

mod <- lm(NodNum ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = NodNum, Genotype = SoilConc)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)


mod <- lm(NodNum ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = NodNum, colour = Genotype)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)

tt = emmeans(mod,  ~ SoilConc|Genotype)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(tt, simple = "Genotype")

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

# Genotype and interaction are sig in both subsets. Within Native, H-L is sig (0.0278). Between genotypes, low is sig ( p < 0.001)

# Nod Weight

mod <- lm(NodWeight ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)
ggplot(SRMortSub2, aes(x = SoilConc, y = NodWeight, colour = Genotype)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)


mod <- lm(NodWeight ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
ggplot(SRMortSub3, aes(x = SoilConc, y = NodWeight, colour = Genotype)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype)

tt = emmeans(mod,  ~ SoilConc|Genotype)
dd = pairs(tt)
summary(dd, by = NULL)
pairs(tt, simple = "Genotype")

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)
# First subset, genotype is sig ( p < 0.001). 2nd subset, soil and genotype are sig (p = 0.049, p < 0.001). Within genotype, native H-N is sig (p = 0.0224). between genotypes, sig for low and none ( p = 0.002, p = 0.0081)

# Going back to trapping data to check for match x no match for mortality

# Create match variable
TrapMort <- within(TrapMort, MatchNoMatch <- ifelse(Range == "Native" & grepl("PT", Treatment), "Match", ifelse(Range == "Invaded" & grepl("FL", Treatment), "Match", "NoMatch")))

Trap

ggplot(TrapMortSub, aes(x = Genotype, y = Dead_0, colour = MatchNoMatch)) + stat_summary(fun.data = "mean_se") + coord_cartesian(ylim = c(0,1)) +theme_bw() + facet_wrap(~ Range, scales = "free_x")

mod <- glm(as.factor(Dead_0) ~ MatchNoMatch * Range, family = binomial(link = "logit"), data = TrapMortSub)
anova(mod, test = "Chisq")

tt = emmeans(mod,  ~ Range| MatchNoMatch)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ MatchNoMatch | Range)
dd = pairs(tt)
summary(dd, by = NULL)

# No Significance

# March 2, 2018 Analysis for comparison chart
mod <- glm(DeadOrAlive ~ SoilConc * Genotype, family = binomial(link = "logit"), data = SRMortSub2)
anova(mod, test = "Chisq")

mod <- glm(DeadOrAlive ~ SoilConc * Genotype, family = binomial(link = "logit"), data = SRMortSub3)
anova(mod, test = "Chisq")
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL)

mod <- lm(ShootWeight ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)

mod <- lm(ShootWeight ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL)

mod <- lm(RootWeight ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)

mod <- lm(RootWeight ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL)    

mod <- lm(RootWeight ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)

mod <- lm(RootWeight ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL) 

mod <- lm(Root_Shoot ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)

mod <- lm(Root_Shoot ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL) 

mod <- lm(NodNum ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)

mod <- lm(NodNum ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL)

mod <- lm(NodWeight ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)

mod <- lm(NodWeight ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL)

SRMortSub2 <- within(SRMortSub2, {
  NNperRoot <- NodNum / RootWeight
  NWperRoot <- NodWeight / RootWeight
  ShootperNW <- ShootWeight / NodWeight
  })

SRMortSub3 <- within(SRMortSub3, {
  NNperRoot <- NodNum / RootWeight
  NWperRoot <- NodWeight / RootWeight
  ShootperNW <- ShootWeight / NodWeight
})

mod <- lm(NNperRoot ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)

mod <- lm(NNperRoot ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL)

mod <- lm(NWperRoot ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)

mod <- lm(NWperRoot ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL)

mod <- lm(ShootperNW ~ SoilConc * Genotype, data = SRMortSub2)
anova(mod)

mod <- lm(ShootperNW ~ SoilConc * Genotype, data = SRMortSub3)
anova(mod)
summary(glht(mod, mcp(SoilConc = "Tukey")))

tt = emmeans(mod,  ~ Genotype| SoilConc)
dd = pairs(tt)
summary(dd, by = NULL)

tt = emmeans(mod,  ~ SoilConc | Genotype)
dd = pairs(tt)
summary(dd, by = NULL)

# PLots 
P1 <- ggplot(SRMortSub2, aes(x = SoilConc, y = Survival)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Mortality")

P2 <- ggplot(SRMortSub3, aes(x = SoilConc, y = Survival)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Mortality")

P3 <- ggplot(SRMortSub2, aes(x = SoilConc, y = ShootWeight)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Shoot Weight")

P4 <- ggplot(SRMortSub3, aes(x = SoilConc, y = ShootWeight)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Shoot Weight")

P5 <- ggplot(SRMortSub2, aes(x = SoilConc, y = RootWeight)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Root Weight")

P6 <- ggplot(SRMortSub3, aes(x = SoilConc, y = RootWeight)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Root Weight")

P7 <- ggplot(SRMortSub2, aes(x = SoilConc, y = Root_Shoot)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Root:Shoot Ratio")

P8 <- ggplot(SRMortSub3, aes(x = SoilConc, y = Root_Shoot)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Root:Shoot Ratio")

P9 <- ggplot(SRMortSub2, aes(x = SoilConc, y = NodNum)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Nodule Number")

P10 <- ggplot(SRMortSub3, aes(x = SoilConc, y = NodNum)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Nodule Number")

P11 <- ggplot(SRMortSub2, aes(x = SoilConc, y = NodWeight)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Nodule Weight")

P12 <- ggplot(SRMortSub3, aes(x = SoilConc, y = NodWeight)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Nodule Weight")

P13 <- ggplot(SRMortSub2, aes(x = SoilConc, y = NNperRoot)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Nodule Number per Root Biomass")

P14 <- ggplot(SRMortSub3, aes(x = SoilConc, y = NNperRoot)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Nodule Number per Root Biomass")

P15 <- ggplot(SRMortSub2, aes(x = SoilConc, y = NWperRoot)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Nodule Weight per Root Biomass")

P16 <- ggplot(SRMortSub2, aes(x = SoilConc, y = NWperRoot)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Nodule Weight per Root Biomass")

P17 <- ggplot(SRMortSub2, aes(x = SoilConc, y = ShootperNW)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Shoot Weight per Nodule Weight ")

P18 <- ggplot(SRMortSub3, aes(x = SoilConc, y = ShootperNW)) + stat_summary(fun.data = "mean_se") + facet_wrap(~ Genotype) + theme_bw() + ggtitle("Shoot Weight per Nodule Weight ")

multi.page <- ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18, nrow = 2, ncol = 2)
multi.page[1]
ggexport(multi.page, filename = "multi.page.ggplot2.pdf")
