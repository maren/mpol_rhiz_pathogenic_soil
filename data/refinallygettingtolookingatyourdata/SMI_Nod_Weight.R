read.csv("SMIOpt_Ratios_KW.csv")
dat <- read.csv("SMIOpt_Ratios_KW.csv")

#effect of WSM (y/n) on nod no.
boxplot(dat$Nod_no ~dat$WSM., data=dat, xlab="Presence of WSM", ylab="Nod no.", main="Effect of WSM on Nod No.")
NW.lm <- lm(dat$Nod_no ~dat$WSM., data=dat)
Anova(NW.lm)
WSM <- dat$WSM.
class(WSM)
WSM1 <- as.character(WSM)

#WSM and geno effect on nod no
NWG.lm <- lm(dat$Nod_no ~dat$WSM+Genotype, data=dat)
Anova(NWG.lm)

#plotting num of nods according to treatment
boxplot(dat$Nod_no ~dat$Treatment, data=dat, ylab="Nod number", col="forestgreen", main="Nod no. vs. Treatment", las=2)

Treatment <- dat$Treatment
nod_treat <- lm(dat$Nod_no ~Treatment, data=dat)
effect("Treatment", nod_treat) 
glht(nod_treat, mcp(Treatment = "Tukey"))
cld(glht(nod_treat, mcp(Treatment = "Tukey")))

TN.lm<- lm(dat$Nod_no ~dat$Treatment+Genotype, data=dat)
summary(TN.lm)
Anova(TN.lm)

#plotting num of nods according to genotype location *PT is signif*
boxplot(dat$Nod_no ~dat$Location, data=dat, xlab="Genotype", ylab="Nod number", main="Nod no. vs. Genotype")
NL.lm <- lm(dat$Nod_no ~dat$Location, data=dat)
summary(NL.lm)
Anova(NL.lm)

#subsetting out -WSMs
dat.pos <- subset(dat, dat$WSM.=="1")

#Geno effect on root biomass +WSM only
RWG.lm <- lm(dat.pos$Root_wt ~dat.pos$Genotype, data=dat.pos)
Anova(RWG.lm)

#Geno effect on shoot weight?
SWG.lm <- lm(dat.pos$Shoot_wt ~dat.pos$Genotype, data=dat.pos)
Anova(SWG.lm)

#Genotype effect on nod biomass
boxplot(dat.pos$Nod_wt ~dat.pos$Genotype, xlab="Genotype", ylab="Nod Biomass (g)", main="Genotype vs. Nod Biomass")
gn.lm <- lm(dat.pos$Nod_wt ~dat.pos$Genotype, data=dat.pos)
Anova(gn.lm) 

gnt.lm <- lm(dat.pos$Nod_wt ~dat.pos$Genotype+dat.pos$Treatment)
Anova(gnt.lm) #signific effect of geno*** and treatment*

treatment <- dat.pos$Treatment
gnt.lm <- lm(dat.pos$Nod_wt ~treatment, data=dat.pos)
effect("treatment", gnt.lm) 

#Testing if there are differences between treatments (affecting nod wt)
glht(gnt.lm, mcp(treatment = "Tukey"))
cld(glht(gnt.lm, mcp(treatment = "Tukey"))) # **All performing the same***
anova(gnt.lm) 
boxplot(dat.pos$Nod_wt ~dat.pos$Treatment, data=dat.pos, las=2, xlab="Treatment", ylab="Nod weight (g)", main="Nodule weight vs. Treatment")

#Pairwise comparison differences (accord to zero)
pair.root.treat<- glht(gnt.lm, linfct = mcp(treatment = "Tukey"))
summary(pair.root.treat)
--

#Geno effect on root biomass
boxplot(dat.pos$Root_wt ~dat.pos$Genotype, xlab="Genotype", ylab="Root biomass (g)", main="Genotype vs. Root Biomass")
gr.lm <- lm(dat.pos$Root_wt ~dat.pos$Genotype+dat.pos$Treatment, data=dat.pos)
Anova(gr.lm) #geno and treat having significant effect on root biomass

Treatment1 <- dat.pos$Treatment
GR.lm <- lm(dat.pos$Root_wt ~Treatment1, data=dat.pos)
effect("Treatment1", GR.lm) 

#Testing if there are differences between two (diff than zero)
glht(GR.lm, mcp(Treatment1 = "Tukey"))
cld(glht(GR.lm, mcp(Treatment1 = "Tukey")))
anova(GR.lm) #Treat has effect on leaf ct

#Pairwise comparison differences (accord to zero)
pair.root.treat<- glht(GR.lm, linfct = mcp(Treatment1 = "Tukey"))
summary(pair.root.treat)


---
Treatment <- dat.pos$Treatment
GT.lm <- lm(dat.pos$Shoot_wt ~Treatment, data=dat.pos)
effect("Treatment", GT.lm) 

#Testing if there are differences between two (diff than zero)
glht(GT.lm, mcp(Treatment = "Tukey"))
cld(glht(GT.lm, mcp(Treatment = "Tukey")))
anova(GT.lm) 

#Pairwise comparison differences (accord to zero)
pair.shoot.treat<- glht(GT.lm, linfct = mcp(Treatment = "Tukey"))
summary(pair.shoot.treat)
-----
#Treatment effect on nod number
treat <- dat.pos$Treatment
NT.lm <- lm(dat.pos$Nod_no ~treat, data=dat.pos)
effect("treat", NT.lm)
glht(NT.lm, mcp(treat= "Tukey"))
cld(glht(NT.lm, mcp(treat = "Tukey")))

