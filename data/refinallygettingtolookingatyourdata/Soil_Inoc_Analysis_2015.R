read.csv('TrappingData_19May15_KW-MB - Sheet1.csv',as.is=T)->dat
head(dat); dim(dat)

summary(as.factor(dat$Treatment))
# correct neg -> Neg
dat[dat$Treatment=="neg","Treatment"]<-"Neg"
summary(as.factor(dat$Treatment))

# add N vs I for plant, soil
# add column for soil inoc

soil.inoc<-dat$Treatment!="WSM"&dat$Treatment!="Neg"
head(soil.inoc)
rhiz.inoc<-rep(NA,length(dat[,1]))
rhiz.inoc[dat$Treatment=="WSM"]<-"WSM"
rhiz.inoc[dat$Treatment=="Neg"]<-"Neg"

soil.type<-soil.inoc
soil.type[dat$Treatment=="#1"|dat$Treatment=="#2"|dat$Treatment=="#3"|dat$Treatment=="#4"]<-"Inv.s"
soil.type[dat$Treatment=="#5"|dat$Treatment=="#6"|dat$Treatment=="#7"|dat$Treatment=="#8"]<-"Nat.s"
soil.type[soil.type==F]<-NA
head(soil.type)
summary(as.factor(soil.type))

summary(as.factor(dat$Location))
plant.type<-dat$Location
plant.type[dat$Location=="Australia"|dat$Location=="California"|dat$Location=="Chile"|dat$Location=="Florida"]<-"Inv.p"
plant.type[dat$Location=="France "|dat$Location=="Italy"|dat$Location=="Portugal"|dat$Location=="Spain"|dat$Location=="Malta "]<-"Nat.p"
summary(as.factor(plant.type))

hist(dat$Leaf.Count..13.May.15[dat$Leaf.Count..13.May.15>0],breaks=20)

boxplot(dat$Leaf.Count..13.May.15~dat$Treatment)

last.fungus<-dat$Fungus..0.no..1.yes.3
last.fungus[is.na(dat$Fungus..0.no..1.yes.3)]<-0
summary(as.factor(last.fungus))
aggregate(last.fungus,by=list(dat$Treatment),sum)

ones<-rep(1,length(dat[,1]))
alive<-dat$Dead.0
alive[is.na(dat$Dead.0)]<-1
summary(as.factor(alive))
aggregate(alive,by=list(dat$Location),sum)->alive.loc
aggregate(ones,by=list(dat$Location),sum)->alive.loc.num
cbind(alive.loc[,1],alive.loc[,2]/alive.loc.num[,2])

aggregate(alive,by=list(dat$Treatment),sum)->alive.trt
aggregate(ones,by=list(dat$Treatment),sum)->alive.trt.num
cbind(alive.trt[,1],alive.trt[,2]/alive.trt.num[,2])

aggregate(alive,by=list(soil.type,plant.type),sum)->alive.comb
aggregate(ones,by=list(soil.type,plant.type),sum)->alive.comb.num
cbind(alive.comb,alive.comb.num[,3],alive.comb[,3]/alive.comb.num[,3])

library(car)
Anova(lm(alive~soil.type*plant.type)->lm1)
Anova(lm(dat$Leaf.Count..13.May.15~ soil.type*dat$Treatment*plant.type*dat$Genotype)->lm1)

aggregate(alive,by=list(dat$Treatment,dat$Genotype),sum)

get.counts.tot<-function(phen,my.groups){
	aggregate(phen,by=my.groups,sum)->by.trt
	aggregate(ones,by=my.groups,sum)->by.trt.num
	ll<-length(my.groups)
cbind(by.trt,by.trt.num[,ll+1],by.trt[,ll+1]/by.trt.num[,ll+1])
}

get.counts.tot(alive,list(soil.inoc))->soil.inoc.table
fisher.test(soil.inoc.table[,2:3])

get.counts.tot(alive,list(soil.type,plant.type))->NvsI.table

glm.fit<-glm(alive~soil.inoc,family='binomial')
summary(glm.fit)

glm.fit<-glm(alive~soil.typedat$Treatment+plant.type,family='binomial')
summary(glm.fit)

library(lme4)
leaf.5<-dat$Leaf.Count..13.May.15
leaf.2<-dat$Leaf.Count..22.April.15.
leaf.3<-dat$X4.28.2015
leaf.4<-dat$Leaf.Count..05.May.15
fungus.2<-dat$Fungus..0.no..1.yes
fungus.3<-dat$Fungus..0.no..1.yes.1
fungus.4<-dat$Fungus..0.no..1.yes.2
fungus.5<-dat$Fungus..0.no..1.yes.3
fungus.2[is.na(fungus.2)&!is.na(leaf.2)]<-0
fungus.3[is.na(fungus.3)&!is.na(leaf.2)]<-0
fungus.4[is.na(fungus.4)&!is.na(leaf.2)]<-0
fungus.5[is.na(fungus.5)&!is.na(leaf.2)]<-0
data.frame(alive,leaf.2,leaf.3,leaf.4,leaf.5,fungus.2,fungus.3,fungus.4,fungus.5, soil.type,dat$Treatment,plant.type,dat$Genotype,rhiz.inoc,soil.inoc)->alive.data

lmer(alive~plant.type*soil.type +(1|dat.Treatment) +(1|dat.Genotype) ,data=alive.data)->lmer1
summary(lmer1)

summary(lmer(leaf.2~plant.type*soil.type +(1|dat.Treatment) +(1|dat.Genotype) ,data=alive.data))
summary(lmer(leaf.3~plant.type*soil.type +(1|dat.Treatment) +(1|dat.Genotype) ,data=alive.data))
summary(lmer(leaf.4~plant.type*soil.type +(1|dat.Treatment) +(1|dat.Genotype) ,data=alive.data))
summary(lmer(leaf.5~plant.type*soil.type +(1|dat.Treatment) +(1|dat.Genotype) ,data=alive.data))

summary(lmer(fungus.2~plant.type*soil.type +(1|dat.Treatment) +(1|dat.Genotype) ,data=alive.data))
summary(lmer(fungus.3~plant.type*soil.type +(1|dat.Treatment) +(1|dat.Genotype) ,data=alive.data))
summary(lmer(fungus.4~plant.type*soil.type +(1|dat.Treatment) +(1|dat.Genotype) ,data=alive.data))
summary(lmer(fungus.5~plant.type*soil.type +(1|dat.Treatment) +(1|dat.Genotype) ,data=alive.data))

summary(lmer(leaf.5~plant.type*rhiz.inoc+(1|dat.Genotype),data=alive.data))
summary(lmer(leaf.4~plant.type*rhiz.inoc+(1|dat.Genotype),data=alive.data))
summary(lmer(leaf.3~plant.type*rhiz.inoc+(1|dat.Genotype),data=alive.data))
summary(lmer(leaf.2~plant.type*rhiz.inoc+(1|dat.Genotype),data=alive.data))

boxplot(leaf.5~rhiz.inoc*dat.Genotype,data=alive.data,las=2)
aggregate(leaf.5,by=list(rhiz.inoc,dat$Genotype),data=alive.data,mean,na.rm=T)->mean.leaf
mean.leaf[mean.leaf[,1]=="Neg",]->mean.leaf.neg
mean.leaf[mean.leaf[,1]=="WSM",]->mean.leaf.pos
cbind(mean.leaf.pos[,2],mean.leaf.pos[,3]/mean.leaf.neg[,3])->rhiz.effect


## not working
geno.inv<-merge(unique(dat$Genotype),alive.data[,c('plant.type','dat.Genotype')],by.x=1,by.y=2,all=F)
merge(rhiz.effect,geno.inv,by.x=1,by.y=1,all=F)->rhiz.effect.2
my.col<-rep("gray",length(rhiz.effect.2[,1]))
my.col[rhiz.effect.2[,3]=="Nat.p"]<-"black"
barplot(as.numeric(as.character(rhiz.effect.2[,2])),names.arg=rhiz.effect.2[,1],col=my.col,las=2)
