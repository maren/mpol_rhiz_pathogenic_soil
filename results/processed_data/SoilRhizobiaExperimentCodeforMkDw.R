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

SRMort <- read.csv("Soil_Rhiz_All_Data20Feb2018.csv")




# Now doing the same thing for the Trapping experiment
TrapMort <- read.csv("Trapping_Alll_Data20Feb2018.csv")


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
#ggsave(file = "TrappingShootWeight21Feb2018.pdf")

# Root Data
RootMod <- lmer(RootWeight ~ Treatment * Range + (1| Genotype), data = TrapBio)
anova(RootMod)
ggplot(TrapBio, aes(x = Treatment, y = RootWeight, colour = Range)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Root Weight (mg)", title = "Trapping Experiment Root Weight 21 Feb 2018") 
#ggsave(file = "TrappingRootWeight21Feb2018.pdf")

# Nodule Data
#   Number
NNMod <- lmer(NodNum ~ Range + (1| Genotype), data = TrapBio)
anova(NNMod)
ggplot(TrapBio, aes(x = Range, y = NodNum)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Nodule Number", title = "Trapping Experiment Nodule Number 21 Feb 2018")  + coord_cartesian(ylim = c(0,30))
#ggsave(file = "TrappingNoduleNumber21Feb2018.pdf") 

#   Weight
NWMod <- lmer(NodWeight ~ Range + (1| Genotype), data = TrapBio)
anova(NWMod)
ggplot(TrapBio, aes(x = Range, y = NodWeight)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Nodule Weight (mg)", title = "Trapping Experiment Nodule Weight 21 Feb 2018")  + coord_cartesian(ylim = c(0,.0035))
#ggsave(file = "TrappingWeight21Feb2018.pdf")


# Soil x Rhiz Exp: Mortality ====
SRMort$Treatment <- as.factor(SRMort$Treatment)
SRMortMod <- glm(DeadOrAlive ~ Treatment , family = binomial(link = "logit"), data = SRMort)
anova(SRMortMod, test = "Chisq")
summary(glht(SRMortMod, mcp(Treatment = "Tukey"))) # NS between Soil treatments, only soil and buffer/wsm


SRMort$Survival <- ifelse(SRMort$DeadOrAlive == "Alive", 1,0)
ggplot(SRMort, aes(x = Treatment, y = Survival, colour = Inoculate, shape = SoilLocation)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Survival Rate", title = "SxR Experiment Mortality 21 Feb 2018") + facet_grid(. ~ SoilConc, scales = "free_x")

ggplot(SRMort, aes(x = SoilConc, y = Survival, colour = Inoculate)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Soil Concentration", y = "Survival Rate", title = "SxR Experiment Mortality 21 Feb 2018") #+ facet_grid(. ~ Match, scales = "free_x")
#ggsave(file = "SurvivalRateSRExperiment21Feb2018.pdf")


ggplot(SRMort, aes(x = Soil, y = Survival)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Soil Concentration", y = "Survival Rate", title = "SxR Experiment Mortality 21 Feb 2018") + facet_grid(. ~ Range, scales = "free_x")

ggplot(SRMort, aes(x = Soil, y = Survival, shape = Match)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Soil Concentration", y = "Survival Rate", title = "SxR Experiment Mortality 21 Feb 2018") + facet_grid(. ~ Range, scales = "free_x")

SRMortMod2 <- glm(DeadOrAlive ~ SoilConc * Inoculate , family = binomial(link = "logit"), data = SRMort)
anova(SRMortMod2, test = "Chisq")

SRMort <- within(SRMort, Soil <- as.factor(paste(SoilConc, Inoculate)))

SRMortMod3 <- glm(DeadOrAlive ~ Soil * Match , family = binomial(link = "logit"), data = SRMort)
anova(SRMortMod3, test = "Chisq")


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


# SR Biomass Data =====
#   Shoot Data
ShootModSR <- lm(ShootWeight ~ Treatment * Range, data = SRMort)
anova(ShootModSR) #  Treatment and Range are S (p < 0.001; p = 0.003)

ShootModSR <- lm(ShootWeight ~ Soil * Range, data = SRMort)
anova(ShootModSR)

ShootModSR <- lm(ShootWeight ~ Soil * Match, data = SRMort)
anova(ShootModSR)

ggplot(SRMort, aes(x = Soil, y = ShootWeight, colour = Range)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Shoot Weight (mg)", title = "S x R Experiment Shoot Weight 21 Feb 2018") 

ggplot(SRMort, aes(x = Soil, y = ShootWeight)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Shoot Weight (mg)", title = "S x R Experiment Shoot Weight 21 Feb 2018") + facet_grid(. ~ Range, scales = "free_x")
#ggsave(file = "SR_ShootWeight21Feb2018.pdf")

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

ggplot(SRMort, aes(x = Soil, y = RootWeight, colour = Range)) + stat_summary(fun.data = "mean_se") + theme_classic() + labs(x = "Treatment", y = "Shoot Weight (mg)", title = "S x R Experiment Shoot Weight 21 Feb 2018") 

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

NNModSRs <- lm(NodNum ~ Soil , data = SRMort_I)
Anova(NNModSRs, type = "3")
summary(glht(NNModSRs, mcp(Soil = "Tukey")))
