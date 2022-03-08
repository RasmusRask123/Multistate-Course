library(survival)
library(Epi)
library(popEpi)

data(cancer)
data(lung)
lung$sex <- factor(lung$sex,levels=1:2,labels=c("M","W"))
lung$time <- lung$time/(365.25/12)
head(lung)
km <- survfit(Surv(time,status==2)~1,data=lung)
km
summary(km)
plot(km)
km_sex <- survfit(Surv(time,status==2)~sex,data=lung)
plot(km_sex,col=1:2)
plot(km_sex, col = c("blue", "red"), lwd = 1, conf.int = TRUE)
lines(km_sex, col = c("blue", "red"), lwd = 3)
with(lung,tapply(age,sex,mean))
survdiff(Surv(time,status==2)~sex,data=lung)

c0 <- coxph(Surv(time,status==2)~sex,data=lung)
c1 <- coxph(Surv(time,status==2)~sex+age,data=lung)
c2 <- coxph(Surv(time,status==2)~sex+I(age/10),data=lung)
ci.exp(c0)
ci.exp(c1)
ci.exp(c2)
summary(c1)


p1 <- glm(cbind(status==2,time)~sex+age,
          family=poisreg,
          data=lung)
px <- glm(status == 2 ~ sex + age + offset(log(time)),
          family = poisson,
          data = lung)
ci.exp(p1)
ci.exp(c2)
ci.exp(px)

L1 <-  Lexis(exit = list(tfl = time),
              exit.status = factor(status,
                                     levels = 1:2,
                                     labels = c("Alive","Dead")),
              data = lung)
L1_cox <- coxph(Surv(tfl,tfl + lex.dur,lex.Xst == "Dead") ~ sex + age,
             data = L1)
boxes(L1, boxpos = TRUE)

cL <- coxph.Lexis(L1, tfl ~ sex + age)
pL <- glm.Lexis(L1, ~ sex + age)

cL
pL

S1 <- splitMulti(L1,tfl=0:36)
summary(L1)
summary(S1)
wh <- names(L1)[1:10] # names of variables in some order
subset(L1, lex.id == 10)[,wh]
subset(S1,lex.id==10)[,wh]
ps <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ Ns(tfl, knots = seq(0, 36, 12)) + sex + age,
           family = poisreg,
           data = S1)
pc <- glm(cbind(lex.Xst == "Dead", lex.dur) ~ sex + age,
           family = poisreg,
           data = L1)

ci.exp(ps)
summary(ps)

round(cbind(ci.exp(cl),
            ci.exp(ps, subset = c("sex","age")),
            ci.exp(pc, subset = c("sex","age"))), 3)

plot(ps)
data(DMlate)
# str(DMlate)
set.seed(1952)
DMlate <- DMlate[sample(1:nrow(DMlate), 2000),]
str(DMlate)
Ldm <- Lexis(entry = list(per = dodm,
                          age = dodm - dobth,
                          tfd = 0),
              exit = list(per = dox),
              exit.status = factor(!is.na(dodth),
                                     labels = c("DM","Dead")),
              data = DMlate)
summary(Ldm)

Cdm <- cutLexis(Ldm,
                 cut = Ldm$dooad,
                 timescale = "per",
                 new.state = "OAD")
summary(Cdm)
Adm <- subset(Cdm, lex.Cst == "DM")
summary(Adm)
boxes(Adm,
      boxpos = TRUE,
      scale.R = 100,
      show.BE = T)
boxes(Cdm,
      boxpos = TRUE,
      scale.R = 100,
      show.BE = T)

m3 <- survfit(Surv(tfd,
                    tfd + lex.dur,
                    lex.Xst) ~ 1,
               id = lex.id,
               data = Adm)

head(cbind(time = m3$time, m3$pstate))
tail(cbind(time = m3$time, m3$pstate))
names(m3)
m3$transitions
summary(Adm)
m3$pstate
par( mfrow=c(1,2) )
 matplot(m3$time, m3$pstate,
           type="s", lty=1, lwd=4,
           col=c("ForestGreen","red","black"),
           xlim=c(0,15), xaxs="i",
           ylim=c(0,1), yaxs="i" )
stackedCIF(m3, lwd=3, xlim=c(0,15), xaxs="i", yaxs="i" )
text( rep(12,3), c(0.9,0.3,0.6), levels(Cdm) )
box()

 
 
Sdm <- splitMulti(Adm, tfd = seq(0,20,0.1) )
summary(Adm) 
round(cbind(
with(subset(Sdm, lex.Xst == "OAD" ), quantile(tfd + lex.dur, 0:10/10)),
with(subset(Sdm, lex.Xst == "Dead"), quantile(tfd + lex.dur, 0:10/10))),3)
okn <- c(0,0.5,3,6)
dkn <- c(0,2.0,5,9)
OAD.glm <- glm.Lexis(Sdm, ~ Ns(tfd, knots = okn), from = "DM", to = "OAD" )
Dead.glm <- glm.Lexis(Sdm, ~ Ns(tfd, knots = dkn), from = "DM", to = "Dead")

int <- 0.01
nd <- data.frame(tfd = seq(0, 15, int))
l.glm <- ci.pred( OAD.glm, nd)
m.glm <- ci.pred(Dead.glm, nd)
matshade(nd$tfd,
          cbind(l.glm, m.glm) * 100,
          plot = TRUE,
          log = "y", ylim = c(2, 20),
          col = rep(c("red","black"), 2), lwd = 3)

cR <- ci.Crisk(mods = list(OAD = OAD.glm,
                           Dead = Dead.glm), nd = nd)

matshade(as.numeric(dimnames(cR$Crisk)[[1]]),
         cbind(cR$Crisk[,1,],
               cR$Crisk[,2,],
               cR$Crisk[,3,]), plot = TRUE,
         lwd = 2, col = c("limegreen","red","black"))

mat2pol(cR$Crisk[,3:1,1], col = c("forestgreen","red","black")[3:1])
matshade(as.numeric(dimnames(cR$Srisk)[[1]]),
          cbind(cR$Srisk[,1,],
                cR$Srisk[,2,]), plot = TRUE,
         lwd = 2, col = c("black","red"),
         ylim = 0:1, yaxs = "i")


data(steno2)
steno2 <- cal.yr(steno2)
steno2 <- transform(steno2,
                       doEnd = pmin(doEnd, doDth, na.rm = TRUE))
str(steno2)
L2 <- Lexis(entry = list(per = doBase,
                          age = doBase - doBth,
                          tfi = 0),
             exit = list(per = doEnd),
             exit.status = factor(deathCVD + !is.na(doDth),
                                    labels=c("Mic","D(oth)","D(CVD)")),
             id = id,
             data = steno2)
summary(L2,t=TRUE)
boxes(L2, boxpos = TRUE, show.BE = TRUE)

data(st2alb)
cut2 <- cal.yr(st2alb)
names(cut2)
names(cut2) <- c("lex.id","cut","new.state")
str(cut2)
addmargins(table(table(cut2$lex.id)))
cut2$lex.id %>% table %>% table %>% addmargins
cut2$cut <- as.numeric(cut2$cut)
L2$per <- as.numeric(L2$per)
L3 <- rcutLexis(L2, cut2,timescale = "per")
summary(L3)
boxes(L3, boxpos = TRUE, cex = 0.8)

(jump <- subset(L3, (lex.Cst == "Norm" & lex.Xst == "Mac") |
                (lex.Xst == "Norm" & lex.Cst == "Mac"))[,c("lex.id", "per", "lex.dur","lex.Cst", "lex.Xst")])

set.seed(1952)
xcut <- select(transform(jump, cut = per + lex.dur * runif(per, 0.1, 0.9),
                            new.state = "Mic"),
                  c(lex.id, cut, new.state))
L4 <- rcutLexis(L3, xcut)
L4 <- Relevel(L4, c("Norm","Mic","Mac","D(CVD)","D(oth)"))
summary(L4)
clr <- c("forestgreen","orange","red","blue",gray(0.3))
boxes(L4, boxpos = list(x = c(20,20,20,80,80),
                           y = c(10,50,90,75,25)),
         show.BE = "nz",
         scale.R = 100,
         digits.R = 2,
         cex = 0.9,
         pos.arr = 0.3,
         col.bg = clr,
         col.border = clr,
         col.txt = c("white","black")[c(1,2,1,1,1)])
S4 <- splitMulti(L4, tfi = seq(0, 25, 1/2))
ma <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                   Ns(age, knots = seq(50, 80, 10)) +
                   lex.Cst)
ma_other <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                  Ns(age, knots = seq(50, 80, 10)) +
                  lex.Cst,to= "D(oth)")
ma_CVD <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                  Ns(age, knots = seq(50, 80, 10)) +
                  lex.Cst,
                 to= "D(CVD)")

m1 <- glm(cbind(lex.Xst %in% c("D(oth)", "D(CVD)")
                 & lex.Cst != lex.Xst,
                 lex.dur)
           ~ Ns(tfi, knots = seq( 0, 20, 5)) +
             Ns(age, knots = seq(50, 80, 10)) +
             lex.Cst,
           family = poisreg,
           data = subset(S4, lex.Cst %in% c("Norm","Mic","Mac")))
Wald(ma,subset="lex.Cst")

round(ci.exp(ma), 2)
round(ci.exp(ma_other), 2)
round(ci.exp(ma_CVD), 2)
expand.grid(tfi = c(NA, seq(0, 20, 5)),
             ain = c(45, 55, 65))
prf <- transform(expand.grid(tfi = c(NA, seq(0, 20, 0.5)),
                              ain = c(45, 55, 65))[-1,],
                  age = ain + tfi,
                  lex.Cst = "Mic")
prf[ 1: 5,]
matshade(prf$age, cbind(ci.pred(ma_other, prf),
                         ci.pred(ma_CVD, prf)) * 100,
          lwd = 3, col = c("black","blue"),
          log = "y", ylim = c(0.02,20), plot = TRUE)

matshade(prf$age, pmin(pmax(
   cbind(ci.pred(ma_other, prf),
           ci.pred(ma_CVD, prf)) * 100, 0.01), 40),
   lwd = 3, col = c("black","blue"),
   log = "y", ylim = c(0.02,20), plot = TRUE)
par(mfrow=c(1,3))
for(st in c("Norm","Mic","Mac")){
     matshade(prf$age, pmin(pmax(
       cbind(ci.pred(ma_other, transform(prf, lex.Cst = st)),
               ci.pred(ma_CVD, transform(prf, lex.Cst = st))) * 100,
       0.05), 60),
       lwd = 3, col = c("black","blue"),
       log = "y", ylim = c(0.1,50), plot = TRUE)
     text(60, 50, st, adj = 0)
}

mix <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                    Ns(age, knots = seq(50, 80, 10)) +
                    lex.Cst * allo,
                  to = "D(oth)")
mix_cvd <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                   Ns(age, knots = seq(50, 80, 10)) +
                   lex.Cst * allo,
                 to = "D(CVD)")
mox <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                   Ns(age, knots = seq(50, 80, 10)) +
                    lex.Cst / allo,
                  to = "D(oth)")
mCx <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                    Ns(age, knots = seq(50, 80, 10)) +
                    lex.Cst / allo,
                  to = "D(CVD)")

mix <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                                     lex.Cst * allo,
                 to = "D(oth)")
mix_cvd <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                       lex.Cst * allo,
                     to = "D(CVD)")
mox <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                    lex.Cst / allo,
                 to = "D(oth)")
mCx <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                   lex.Cst / allo,
                 to = "D(CVD)")



c(deviance(mox), deviance(mix))
round(ci.exp(mix), 3)
round(ci.exp(mox), 3)
round(ci.exp(mix_cvd), 3)
round(ci.exp(mCx), 3)

det <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                    Ns(age, knots = seq(50, 80, 10)) +
                    lex.Cst / allo,
                  from = c("Norm","Mic"),
                  to = c("Mic","Mac"))
imp <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                    Ns(age, knots = seq(50, 80, 10)) +
                    lex.Cst / allo,
                  from = c("Mac","Mic"),
                  to = c("Mic","Norm"))

det <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                                   lex.Cst / allo,
                 from = c("Norm","Mic"),
                 to = c("Mic","Mac"))
imp <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                    lex.Cst / allo,
                 from = c("Mac","Mic"),
                 to = c("Mic","Norm"))


Mima <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                   Ns(age, knots = seq(50, 80, 10)) +
                   allo,
                 from = "Mic",
                 to = "Mac")

nomi <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                    Ns(age, knots = seq(50, 80, 10)) +
                    allo,
                  from = "Norm",
                  to = "Mic")
mima <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                    Ns(age, knots = seq(50, 80, 10)) +
                    allo,
                  from = "Mic",
                  to = "Mac")
mami <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                    Ns(age, knots = seq(50, 80, 10)) +
                    allo,
                  from = "Mac",
                  to = "Mic")
round(ci.exp(det), 3)
round(ci.exp(imp), 3)
round(ci.exp(Mimp), 3)
round(ci.exp(nimp), 3) 
Tr <- list(Norm = list("Mic" = det,
                                               "D(oth)" = mox,
                                               "D(CVD)" = mCx),
                                   Mic = list("Mac" = det,
                                                "Norm" = imp,
                                                "D(oth)" = mox,
                                               "D(CVD)" = mCx),
                                   Mac = list("Mic" = imp,
                                                "D(oth)" = mox,
                                                "D(CVD)" = mCx))
lapply(Tr, names)
ini <- L2[,c("per", "age", "tfi")]
 ini <- rbind(transform(ini, lex.Cst = factor("Mic"), allo = factor("Int")),
              transform(ini, lex.Cst = factor("Mic"), allo = factor("Conv")))
str(ini)
set.seed(1952)
system.time(
   Sorg <- simLexis(Tr = Tr, # models for each transition
                      init = ini, # cohort of straters
                      N = 10, # how many copies of each person in ini
                      t.range = 20, # how long should we simulate before censoring
                      n.int = 200)
)
summary(Sorg)
Sorg <- Relevel(Sorg, c("Norm", "Mic", "Mac", "D(CVD)", "D(oth)"))
summary(Sorg, t = T)
addmargins(table(table(Sorg$lex.id)))
system.time(
   Nst <- nState(Sorg,
                   at = seq(0, 20, 0.2),
                   from = 0,
                   time.scale = "tfi"))
Nint <- nState(subset(Sorg, allo == "Int"),
                at = seq(0, 20, 0.1),
                from = 0,
                time.scale = "tfi")
Nconv<- nState(subset(Sorg, allo == "Conv"),
                  at = seq(0, 20, 0.1),
                  from = 0,
                  time.scale = "tfi")
cbind(
   head(Nint), NA,
   head(Nconv))
Pint <- pState(Nint )
 Pconv <- pState(Nconv)
 clr <- c("forestgreen", "orange", "red", "blue", gray(0.4))
  par(mfrow = c(1,2), mar=c(3,3,2,2))
  plot(Pint , col = clr, xlim = c(0, 20))
  plot(Pconv, col = clr, xlim = c(20, 0))

  clr <- c("forestgreen", "orange", "red", "blue", gray(0.4))
   par(mfrow = c(1,2), mar=c(3,3,2,2))
   plot(Pint, col = clr, xlim = c(0, 20))
   # the survival curve
     lines(as.numeric(rownames(Pint)), Pint[,"Mac"], lwd = 4, col = "white")
   lines(as.numeric(rownames(Pint)), Pint[,"Mac"], lwd = 1, col = "black")
   text(rownames(Pint)[150],
          Pint[150,] - diff(c(0, Pint[150,]))/2,
          colnames(Pint), col = "white", cex = 0.8)
   plot(Pconv, col = clr, xlim = c(20, 0))
   # the survival curve
     lines(as.numeric(rownames(Pconv)), Pconv[,"Mac"], lwd = 4, col = "white")
   lines(as.numeric(rownames(Pconv)), Pconv[,"Mac"], lwd = 1, col = "black")
   text(rownames(Pconv)[150],
         Pconv[150,] - diff(c(0, Pconv[150,]))/2,
          colnames(Pint), col = "white", cex = 0.8)
   mtext(c("Intensive care","Conventional care"),
           side = 3, at = c(1,3)/4, outer = TRUE, line = -2) 


   ini <- S4[1:10,c("lex.id", "per", "age", "tfi", "lex.Cst", "allo")]
   ini[,"lex.id"] <- 1:10
   ini[,"per"] <- 1993 # not used but it is a time scale in S4
   ini[,"age"] <-
      ini[,"ain"] <- rep(seq(45,65,5), 2)
    ini[,"tfi"] <- 0
    ini[,"lex.Cst"] <- factor("Mic",
                                levels = c("Norm","Mic","Mac","D(CVD)","D(oth)"))
    ini[,"allo"] <- factor(rep(c("Int","Conv"), each = 5))
    ini
    
    system.time(
      Sdef <- simLexis(Tr = Tr,
                         + init = ini,
                         + N = 100,
                         + t.range = 20,
                         + n.int = 200))

    
    
    data(st2clin)
    st2clin <- cal.yr(st2clin)
    
    
    
    
    
    
    ## R code from vignette source 'ms-steno2.rnw'
    
    ###################################################
    ### code chunk number 1: ms-steno2.rnw:2-5
    ###################################################
    options( width=90,
             SweaveHooks=list( fig=function()
               par(mar=c(3,3,1,1),mgp=c(3,1,0)/1.6,las=1,bty="n") ) )
    
    
    ###################################################
    ### code chunk number 2: ms-steno2.rnw:14-23
    ###################################################
    library(survival)
    library(Epi)
    library(popEpi)
    # popEpi::splitMulti returns a data.frame rather than a data.table
    options("popEpi.datatable" = FALSE)
    library(tidyverse)
    setwd("c:/bendix/teach/AdvCoh/courses/Aalborg.2022/pracs")
    getwd()
    clear()
    
    
    ###################################################
    ### code chunk number 3: ms-steno2.rnw:28-35
    ###################################################
    nround <-
      function(df, dec = 2)
      {
        wh.num <- sapply(df, is.numeric)
        df[,wh.num] <- round(df[,wh.num], dec)
        print(df)
      }
    
    
    ###################################################
    ### code chunk number 4: ms-steno2.rnw:47-52
    ###################################################
    data(steno2)
    steno2 <- cal.yr(steno2)
    steno2 <- transform(steno2,
                        doEnd = pmin(doEnd, doDth, na.rm = TRUE))
    str(steno2)
    
    
    ###################################################
    ### code chunk number 5: ms-steno2.rnw:61-71
    ###################################################
    L2 <- Lexis(entry = list(per = doBase,
                             age = doBase - doBth,
                             tfi = 0),
                exit = list(per = doEnd),
                exit.status = factor(deathCVD + !is.na(doDth),
                                     labels=c("Mic","D(oth)","D(CVD)")),
                id = id,
                data = steno2)
    summary(L2, t = TRUE)
    boxes(L2, boxpos = TRUE, show.BE = TRUE)
    
    
    ###################################################
    ### code chunk number 6: ms-steno2.rnw:105-111
    ###################################################
    data(st2alb)
    cut2 <- cal.yr(st2alb)
    names(cut2)
    names(cut2) <- c("lex.id","cut","new.state")
    str(cut2)
    head(cut2)
    
    
    ###################################################
    ### code chunk number 7: ms-steno2.rnw:115-117
    ###################################################
    addmargins(table(table(cut2$lex.id)))
    cut2$lex.id %>% table %>% table %>% addmargins
    
    
    ###################################################
    ### code chunk number 8: boxL3
    ###################################################
    cut2$cut <- as.numeric(cut2$cut)
    L2$per <- as.numeric(L2$per)
    L3 <- rcutLexis(L2, cut2)
    summary(L3)
    boxes(L3, boxpos = TRUE, cex = 0.8)
    
    
    ###################################################
    ### code chunk number 9: ms-steno2.rnw:149-153
    ###################################################
    (jump <-
       subset(L3, (lex.Cst == "Norm" & lex.Xst == "Mac") |
                (lex.Xst == "Norm" & lex.Cst == "Mac"))[,
                                                        c("lex.id", "per", "lex.dur","lex.Cst", "lex.Xst")])
    
    
    ###################################################
    ### code chunk number 10: ms-steno2.rnw:161-167
    ###################################################
    set.seed(1952)
    xcut <- select(transform(jump,
                             cut = per + lex.dur * runif(per, 0.1, 0.9),
                             new.state = "Mic"),
                   c(lex.id, cut, new.state))
    xcut
    
    
    ###################################################
    ### code chunk number 11: ms-steno2.rnw:174-177
    ###################################################
    L4 <- rcutLexis(L3, xcut)
    L4 <- Relevel(L4, c("Norm","Mic","Mac","D(CVD)","D(oth)"))
    summary(L4)
    
    
    ###################################################
    ### code chunk number 12: b4
    ###################################################
    clr <- c("forestgreen","orange","red","blue",gray(0.3))
    boxes(L4, boxpos = list(x = c(20,20,20,80,80),
                            y = c(10,50,90,75,25)),
          show.BE = "nz",
          scale.R = 100,
          digits.R = 2,
          cex = 0.9,
          pos.arr = 0.3,
          col.bg = clr,
          col.border = clr,
          col.txt = c("white","black")[c(1,2,1,1,1)])
    
    
    ###################################################
    ### code chunk number 13: ms-steno2.rnw:219-222
    ###################################################
    S4 <- splitMulti(L4, tfi = seq(0, 25, 1/2))
    summary(L4)
    summary(S4)
    
    
    ###################################################
    ### code chunk number 14: ms-steno2.rnw:231-234
    ###################################################
    ma <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                      Ns(age, knots = seq(50, 80, 10)) +
                      lex.Cst)
    
    
    ###################################################
    ### code chunk number 15: ms-steno2.rnw:252-260 (eval = FALSE)
    ###################################################
    ## m1 <- glm(cbind(lex.Xst %in% c("D(oth)", "D(CVD)")
    ##                 & lex.Cst != lex.Xst,
    ##                 lex.dur)
    ##           ~ Ns(tfi, knots = seq( 0, 20, 5)) +
    ##             Ns(age, knots = seq(50, 80, 10)) +
    ##             lex.Cst,
    ##           family = poisreg,
    ##             data = subset(S4, lex.Cst %in% c("Norm","Mic","Mac")))
    
    
    ###################################################
    ### code chunk number 16: ms-steno2.rnw:263-270 (eval = FALSE)
    ###################################################
    ## m2 <- glm((lex.Xst %in% c("D(oth)", "D(CVD)") & lex.Cst != lex.Xst)
    ##           ~ Ns(tfi, knots = seq( 0, 20, 5)) +
    ##             Ns(age, knots = seq(50, 80, 10)) +
    ##             lex.Cst,
    ##           offset = log(lex.dur),
    ##           family = poisson,
    ##             data = subset(S4, lex.Cst %in% c("Norm","Mic","Mac")))
    
    
    ###################################################
    ### code chunk number 17: ms-steno2.rnw:282-283
    ###################################################
    round(ci.exp(ma), 2)
    
    
    ###################################################
    ### code chunk number 18: ms-steno2.rnw:292-293
    ###################################################
    Wald(ma, subset = "lex.Cst")
    
    
    ###################################################
    ### code chunk number 19: ms-steno2.rnw:299-309
    ###################################################
    mo <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                      Ns(age, knots = seq(50, 80, 10)) +
                      lex.Cst,
                    to = "D(oth)")
    round(ci.exp(mo), 3)
    mC <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                      Ns(age, knots = seq(50, 80, 10)) +
                      lex.Cst,
                    to = "D(CVD)")
    round(ci.exp(mC), 3)
    
    
    ###################################################
    ### code chunk number 20: ms-steno2.rnw:315-317
    ###################################################
    Wald(mo, subset = "Cst")
    Wald(mC, subset = "Cst")
    
    
    ###################################################
    ### code chunk number 21: ms-steno2.rnw:325-326
    ###################################################
    summary(L2$age)
    
    
    ###################################################
    ### code chunk number 22: ms-steno2.rnw:334-336
    ###################################################
    expand.grid(tfi = c(NA, seq(0, 20, 5)),
                ain = c(45, 55, 65))
    
    
    ###################################################
    ### code chunk number 23: mort1
    ###################################################
    prf <- transform(expand.grid(tfi = c(NA, seq(0, 20, 0.5)),
                                 ain = c(45, 55, 65))[-1,],
                     age = ain + tfi,
                     lex.Cst = "Mic")
    prf[ 1: 5,]
    prf[40:44,]
    matshade(prf$age, cbind(ci.pred(mo, prf),
                            ci.pred(mC, prf)) * 100,
             lwd = 3, col = c("black","blue"),
             log = "y", ylim = c(0.02,20), plot = TRUE)
    
    
    ###################################################
    ### code chunk number 24: ms-steno2.rnw:362-367
    ###################################################
    matshade(prf$age, pmin(pmax(
      cbind(ci.pred(mo, prf),
            ci.pred(mC, prf)) * 100, 0.01), 40),
      lwd = 3, col = c("black","blue"),
      log = "y", ylim = c(0.02,20), plot = TRUE)
    
    
    ###################################################
    ### code chunk number 25: mort3
    ###################################################
    par(mfrow=c(1,3))
    for(st in c("Norm","Mic","Mac"))
    {
      matshade(prf$age, pmin(pmax(
        cbind(ci.pred(mo, transform(prf, lex.Cst = st)),
              ci.pred(mC, transform(prf, lex.Cst = st))) * 100,
        0.05), 60),
        lwd = 3, col = c("black","blue"),
        log = "y", ylim = c(0.1,50), plot = TRUE)
      text(60, 50, st, adj = 0)
    }
    
    
    ###################################################
    ### code chunk number 26: ms-steno2.rnw:424-429
    ###################################################
    mix <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                       Ns(age, knots = seq(50, 80, 10)) +
                       lex.Cst * allo,
                     to = "D(oth)")
    round(ci.exp(mix), 3)
    
    
    ###################################################
    ### code chunk number 27: ms-steno2.rnw:435-441
    ###################################################
    mox <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                       Ns(age, knots = seq(50, 80, 10)) +
                       lex.Cst / allo,
                     to = "D(oth)")
    round(ci.exp(mox), 3)
    c(deviance(mox), deviance(mix))
    
    
    ###################################################
    ### code chunk number 28: ms-steno2.rnw:449-454
    ###################################################
    mCx <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                       Ns(age, knots = seq(50, 80, 10)) +
                       lex.Cst / allo,
                     to = "D(CVD)")
    round(ci.exp(mCx), 3)
    
    
    ###################################################
    ### code chunk number 29: ms-steno2.rnw:474-488
    ###################################################
    det <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                       Ns(age, knots = seq(50, 80, 10)) +
                       lex.Cst / allo,
                     from = c("Norm","Mic"),
                     to = c("Mic","Mac"))
    imp <- glm.Lexis(S4, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                       Ns(age, knots = seq(50, 80, 10)) +
                       lex.Cst / allo,
                     from = c("Mac","Mic"),
                     to = c("Mic","Norm"))
    round(ci.exp(det), 3)
    round(ci.exp(imp), 3)
    round(ci.exp(det, subset="al"), 1)
    round(ci.exp(imp, subset="al"), 1)
    
    
    ###################################################
    ### code chunk number 30: ms-steno2.rnw:524-535
    ###################################################
    Tr <- list(Norm = list("Mic" = det,
                           "D(oth)" = mox,
                           "D(CVD)" = mCx),
               Mic = list("Mac" = det,
                          "Norm" = imp,
                          "D(oth)" = mox,
                          "D(CVD)" = mCx),
               Mac = list("Mic" = imp,
                          "D(oth)" = mox,
                          "D(CVD)" = mCx))
    lapply(Tr, names)
    
    
    ###################################################
    ### code chunk number 31: ms-steno2.rnw:581-585
    ###################################################
    ini <- L2[,c("per", "age", "tfi")]
    ini <- rbind(transform(ini, lex.Cst = factor("Mic"), allo = factor("Int")),
                 transform(ini, lex.Cst = factor("Mic"), allo = factor("Conv")))
    str(ini)
    
    
    ###################################################
    ### code chunk number 32: ms-steno2.rnw:595-602
    ###################################################
    set.seed(1952)
    system.time(
      Sorg <- simLexis(Tr = Tr,  # models for each transition
                       init = ini, # cohort of straters
                       N = 10,  # how many copies of each person in ini
                       t.range = 20,  # how long should we simulate before censoring
                       n.int = 200))# how many intervals for evaluating rates
    
    
    ###################################################
    ### code chunk number 33: ms-steno2.rnw:606-609
    ###################################################
    Sorg <- Relevel(Sorg, c("Norm", "Mic", "Mac", "D(CVD)", "D(oth)"))
    summary(Sorg, t = T)
    nround(subset(Sorg, lex.id %in% 28:32), 2)
    
    
    ###################################################
    ### code chunk number 34: ms-steno2.rnw:614-615
    ###################################################
    addmargins(table(table(Sorg$lex.id)))
    
    
    ###################################################
    ### code chunk number 35: ms-steno2.rnw:622-629
    ###################################################
    system.time(
      Nst <- nState(Sorg,
                    at = seq(0, 20, 0.2),
                    from = 0,
                    time.scale = "tfi"))
    str(Nst)
    head(Nst)
    
    
    ###################################################
    ### code chunk number 36: ms-steno2.rnw:634-645
    ###################################################
    Nint <- nState(subset(Sorg, allo == "Int"),
                   at = seq(0, 20, 0.1),
                   from = 0,
                   time.scale = "tfi")
    Nconv<- nState(subset(Sorg, allo == "Conv"),
                   at = seq(0, 20, 0.1),
                   from = 0,
                   time.scale = "tfi")
    cbind(
      head(Nint), NA,
      head(Nconv))
    
    
    ###################################################
    ### code chunk number 37: ms-steno2.rnw:653-657
    ###################################################
    Pint  <- pState(Nint )
    Pconv <- pState(Nconv)
    str(Pint)
    head(Pint)
    
    
    ###################################################
    ### code chunk number 38: ms-steno2.rnw:665-669
    ###################################################
    clr <- c("forestgreen", "orange", "red", "blue", gray(0.4))
    par(mfrow = c(1,2), mar=c(3,3,2,2))
    plot(Pint , col = clr, xlim = c(0, 20))
    plot(Pconv, col = clr, xlim = c(20, 0))
    
    
    ###################################################
    ### code chunk number 39: pStates
    ###################################################
    clr <- c("forestgreen", "orange", "red", "blue", gray(0.4))
    par(mfrow = c(1,2), mar=c(3,3,2,2))
    
    plot(Pint, col = clr, xlim = c(0, 20))
    # the survival curve
    lines(as.numeric(rownames(Pint)), Pint[,"Mac"], lwd = 4, col = "white")
    lines(as.numeric(rownames(Pint)), Pint[,"Mac"], lwd = 1, col = "black")
    text(rownames(Pint)[150],
         Pint[150,] - diff(c(0, Pint[150,]))/2,
         colnames(Pint), col = "white", cex = 0.8)
    
    plot(Pconv, col = clr, xlim = c(20, 0))
    # the survival curve
    lines(as.numeric(rownames(Pconv)), Pconv[,"Mac"], lwd = 4, col = "white")
    lines(as.numeric(rownames(Pconv)), Pconv[,"Mac"], lwd = 1, col = "black")
    text(rownames(Pconv)[150],
         Pconv[150,] - diff(c(0, Pconv[150,]))/2,
         colnames(Pint), col = "white", cex = 0.8)
    
    mtext(c("Intensive care","Conventional care"),
          side = 3, at = c(1,3)/4, outer = TRUE, line = -2)
    
    
    ###################################################
    ### code chunk number 40: ms-steno2.rnw:729-740
    ###################################################
    ini <- S4[1:10,c("lex.id", "per", "age", "tfi", "lex.Cst", "allo")]
    ini[,"lex.id"]  <- 1:10
    ini[,"per"]     <- 1993 # not used but it is a time scale in S4
    ini[,"age"]     <-
      ini[,"ain"]     <- rep(seq(45,65,5), 2)
    ini[,"tfi"]     <- 0
    ini[,"lex.Cst"] <- factor("Mic",
                              levels = c("Norm","Mic","Mac","D(CVD)","D(oth)"))
    ini[,"allo"]    <- factor(rep(c("Int","Conv"), each = 5))
    ini
    str(ini)
    
    
    ###################################################
    ### code chunk number 41: ms-steno2.rnw:752-760
    ###################################################
    system.time(
      Sdef <- simLexis(Tr = Tr,
                       init = ini,
                       N = 100,
                       t.range = 20,
                       n.int = 200))
    summary(Sdef)
    nround(head(Sdef))
    
    
    ###################################################
    ### code chunk number 42: ms-steno2.rnw:768-774
    ###################################################
    P45i <- nState(subset(Sdef, ain == 45 & allo == "Int"),
                   at = seq(0, 20, 0.1),
                   from = 0,
                   time.scale = "tfi")
    head(P45i)
    head(pState(P45i))
    
    
    ###################################################
    ### code chunk number 43: ms-steno2.rnw:778-818
    ###################################################
    P45c <- pState(nState(subset(Sdef, ain == 45 & allo == "Conv"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    P45i <- pState(nState(subset(Sdef, ain == 45 & allo == "Int"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    P50c <- pState(nState(subset(Sdef, ain == 55 & allo == "Conv"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    P50i <- pState(nState(subset(Sdef, ain == 55 & allo == "Int"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    P55c <- pState(nState(subset(Sdef, ain == 55 & allo == "Conv"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    P55i <- pState(nState(subset(Sdef, ain == 55 & allo == "Int"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    P60c <- pState(nState(subset(Sdef, ain == 55 & allo == "Conv"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    P60i <- pState(nState(subset(Sdef, ain == 55 & allo == "Int"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    P65c <- pState(nState(subset(Sdef, ain == 65 & allo == "Conv"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    P65i <- pState(nState(subset(Sdef, ain == 65 & allo == "Int"),
                          at = seq(0, 20, 0.1),
                          from = 0,
                          time.scale = "tfi"))
    
    
    ###################################################
    ### code chunk number 44: panel5
    ###################################################
    par(mfrow = c(5,2), mar = c(1,1,0,0),
        oma = c(3,3,1,0), mgp=c(3,1,0)/1.6)
    
    plot(P45i, col = clr, xlim = c(0,20))
    plot(P45c, col = clr, xlim = c(20,0))
    
    plot(P50i, col = clr, xlim = c(0,20))
    plot(P50c, col = clr, xlim = c(20,0))
    
    plot(P55i, col = clr, xlim = c(0,20))
    plot(P55c, col = clr, xlim = c(20,0))
    
    plot(P60i, col = clr, xlim = c(0,20))
    plot(P60c, col = clr, xlim = c(20,0))
    
    plot(P65i, col = clr, xlim = c(0,20))
    plot(P65c, col = clr, xlim = c(20,0))
    
    mtext(c("Int","Conv"), side = 3, at = c(1,3)/4, outer = TRUE, line = 0)
    mtext(paste(seq(45,65,5)), side = 2, at = (5:1*2-1)/10,
          outer = TRUE, line = 0)
    
    
    ###################################################
    ### code chunk number 45: ms-steno2.rnw:856-862
    ###################################################
    (ain <- seq(45, 65, 5))
    (all <- levels(S4$allo))
    pdef <- NArray(c(list(ain = ain,
                          allo = all),
                     dimnames(P45i)))
    str(pdef)
    
    
    ###################################################
    ### code chunk number 46: ms-steno2.rnw:869-878
    ###################################################
    for(aa in ain)
      for(gg in all)
        pdef[paste(aa), gg, ,] <-
      nState(subset(Sdef, ain == aa & allo == gg),
             at = as.numeric(dimnames(pdef)[["when"]]),
             from = 0,
             time.scale = "tfi")
    pdef <- sweep(pdef, 1:2, pdef[,,1,"Mic"], "/")
    str(pdef)
    
    
    ###################################################
    ### code chunk number 47: panel3
    ###################################################
    ain <- seq(45, 65, 10)
    par(mfrow = c(length(ain),2),
        mar = c(3,3,1,1),
        oma = c(0,2,1,0),
        mgp = c(3,1,0) / 1.6)
    for(aa in ain)
    {
      mat2pol(pdef[paste(aa),"Int" ,,], col = clr, xlim = c(0,20))
      mat2pol(pdef[paste(aa),"Conv",,], col = clr, xlim = c(20,0))
    }
    mtext(c("Int","Conv"), side = 3, at = c(1,3)/4, outer = TRUE, line = 0)
    mtext(ain, side = 2, at = (length(ain):1 * 2 - 1) / (length(ain) * 2),
          outer = TRUE, line = 0)
    
    
    ###################################################
    ### code chunk number 48: ms-steno2.rnw:927-932
    ###################################################
    AaJ <- survfit(Surv(tfi, tfi + lex.dur, lex.Xst) ~ 1,
                   id = lex.id,
                   data = L4)
    AaJ$transitions
    summary(L4)
    
    
    ###################################################
    ### code chunk number 49: ms-steno2.rnw:936-940 (eval = FALSE)
    ###################################################
    ## survfit(Surv(tfi, tfi + lex.dur, lex.Xst) ~ 1,
    ##         id = lex.id,
    ##     istate = lex.Cst,
    ##       data = L4)
    
    
    ###################################################
    ### code chunk number 50: ms-steno2.rnw:951-959
    ###################################################
    R4 <- sortLexis(L4)
    last <- rev(!duplicated(rev(R4$lex.id)))
    R4$lex.Xst <- ifelse(last & R4$lex.Cst == R4$lex.Xst,
                         "cens",
                         as.character(R4$lex.Xst))
    R4 <- Relevel(factorize(R4), "cens")
    summary(L4)
    summary(R4)
    
    
    ###################################################
    ### code chunk number 51: ms-steno2.rnw:966-970
    ###################################################
    AaJest <- survfit(Surv(tfi, tfi + lex.dur, lex.Xst) ~ 1,
                      id = lex.id,
                      istate = lex.Cst,
                      data = R4)
    
    
    ###################################################
    ### code chunk number 52: ms-steno2.rnw:974-977
    ###################################################
    AaJest$transitions[,c(6,1:5)]
    summary(R4)
    summary(L4)
    
    
    ###################################################
    ### code chunk number 53: ms-steno2.rnw:981-982
    ###################################################
    source("http://bendixcarstensen.com/AdvCoh/courses/Aalborg-2022/R/AaJ.Lexis.R")
    
    
    ###################################################
    ### code chunk number 54: ms-steno2.rnw:986-989
    ###################################################
    AaJepi <- AaJ(L4, timeScale = "tfi")
    AaJepi
    AaJest
    
    
    ###################################################
    ### code chunk number 55: ms-steno2.rnw:995-1000
    ###################################################
    names(AaJepi)
    AaJepi$states
    head(AaJepi$pstate)
    head(AaJepi$lower)
    head(AaJepi$upper)
    
    
    ###################################################
    ### code chunk number 56: AaJ
    ###################################################
    par(mfrow = c(1,1))
    mat2pol(AaJepi$pstate, perm = c(2,1,3,5,4), x = AaJepi$time,
            col = clr)
    lines(AaJepi$time, apply(AaJepi$pstate[,1:3], 1, sum), lwd = 5)
    
    
    ###################################################
    ### code chunk number 57: ms-steno2.rnw:1013-1021
    ###################################################
    AaJallo <- survfit(Surv(tfi, tfi + lex.dur, lex.Xst) ~ allo,
                       id = lex.id,
                       istate = lex.Cst,
                       data = R4)
    AaJallo                 
    AaJallo <- AaJ(L4, ~ allo, timeScale = "tfi")
    AaJallo                 
    names(AaJallo)
    
    
    ###################################################
    ### code chunk number 58: ms-steno2.rnw:1027-1031
    ###################################################
    AaJallo$states
    AaJallo$strata
    wh <- rep(substr(names(AaJallo$strata), 6, 9), AaJallo$strata)
    table(wh)
    
    
    ###################################################
    ### code chunk number 59: AaJstates
    ###################################################
    par(mfrow = c(1,2), mar=c(3,3,2,2))
    mat2pol(AaJallo$pstate[wh=="Int",],
            perm = c(2,1,3:5),
            x = AaJallo$time[wh=="Int"],
            col = clr, xlim = c(0,21), xaxs = "i", yaxs = "i")
    lines(AaJallo$time[wh=="Int"],
          apply(AaJallo$pstate[,1:3], 1, sum)[wh=="Int"], lwd = 4)
    
    mat2pol(AaJallo$pstate[wh=="Conv",],
            perm = c(2,1,3:5),
            x = AaJallo$time[wh=="Conv"],
            col = clr, xlim = c(21,0), xaxs = "i", yaxs = "i")
    lines(AaJallo$time[wh=="Conv"],
          apply(AaJallo$pstate[,1:3], 1, sum)[wh=="Conv"], lwd = 4)
    
    mtext(c("Int","Conv"), side = 3, at = c(1,3)/4, outer = TRUE, line = -2)
    
    
    ###################################################
    ### code chunk number 60: ms-steno2.rnw:1096-1101
    ###################################################
    tLive <- xtabs(lex.dur ~ ain + allo + lex.Cst, data = Sdef) /
      nid(Sdef) * 10 
    str(mtLive <- addmargins(tLive, 3))
    round(ftable(mtLive          , col.vars = c(3,2)), 1)
    round(ftable(mtLive[,,-(4:5)], col.vars = c(3,2)), 1) 
    
    
    ###################################################
    ### code chunk number 61: ms-steno2.rnw:1105-1106
    ###################################################
    round((mtLive[,"Int",-(4:5)] - mtLive[,"Conv",-(4:5)]), 1)
    
    
    ###################################################
    ### code chunk number 62: ms-steno2.rnw:1116-1122
    ###################################################
    tDead <- xtabs((20 - tfi - lex.dur) ~ ain + allo + lex.Xst, 
                   data = subset(Sdef, lex.Xst %in% c("D(oth)", "D(CVD)"))) /
      nid(Sdef) * 10 
    str(mtDead <- addmargins(tDead[,,4:5], 3))
    round(ftable(mtDead          , col.vars = c(3,2)), 1)
    round(ftable(mtDead[,,-(4:5)], col.vars = c(3,2)), 1) 
    
    
    ###################################################
    ### code chunk number 63: ms-steno2.rnw:1135-1142
    ###################################################
    data(st2clin)
    str(st2clin)
    st2clin <- cal.yr(st2clin)
    names(st2clin)
    names(st2clin)[1:2] <- c("lex.id","per")
    summary(st2clin)
    addmargins(table(table(st2clin$lex.id)))
    
    
    ###################################################
    ### code chunk number 64: ms-steno2.rnw:1148-1157
    ###################################################
    S5 <- addCov.Lexis(S4, st2clin, "per")
    tt <- table(st2clin$lex.id)
    (who <- names(tt[tt == 3])[1])
    subset(st2clin, lex.id == who)
    nround(subset(S5,
                  lex.id == who,
                  select = c(lex.id,per,tfi,tfc,exnam,a1c,chol,crea)))
    timeScales(S5)
    timeSince(S5)
    
    
    ###################################################
    ### code chunk number 65: ms-steno2.rnw:1176-1202
    ###################################################
    detc <- glm.Lexis(S5, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                        Ns(age, knots = seq(50, 80, 10)) +
                        lex.Cst / allo +
                        a1c + chol + log2(crea),
                      from = c("Norm","Mic"),
                      to = c("Mic","Mac"))
    impc <- glm.Lexis(S5, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                        Ns(age, knots = seq(50, 80, 10)) +
                        lex.Cst / allo +
                        a1c + chol + log2(crea),
                      from = c("Mic","Mac"),
                      to = c("Norm","Mic")
                      )
    round(ci.exp(detc), 3)
    round(ci.exp(impc), 3)
    moc <- glm.Lexis(S5, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                       Ns(age, knots = seq(50, 80, 10)) +
                       lex.Cst / allo +
                       a1c + chol + log2(crea),
                     to = "D(oth)")
    mCc <- glm.Lexis(S5, ~ Ns(tfi, knots = seq( 0, 20, 5)) +
                       Ns(age, knots = seq(50, 80, 10)) +
                       lex.Cst / allo +
                       a1c + chol + log2(crea),
                     to = "D(CVD)")
    round(ci.exp(moc), 3)
    round(ci.exp(mCc), 3)
    
    
    ###################################################
    ### code chunk number 66: ms-steno2.rnw:1248-1253
    ###################################################
    St4 <- stack(S4)
    c(nrow(S4), nrow(St4))
    table(S4$lex.Cst)
    table(St4$lex.Tr, St4$lex.Cst)
    ftable(St4$lex.Tr, St4$lex.Xst, St4$lex.Fail, col.vars = 2:3)
    
    
    ###################################################
    ### code chunk number 67: ms-steno2.rnw:1259-1261
    ###################################################
    nround(subset(S4 , lex.id == 102)[,1:8], 1)
    nround(subset(St4, lex.id == 102)[,1:9], 1)
    
    
    ###################################################
    ### code chunk number 68: ms-steno2.rnw:1270-1271
    ###################################################
    cbind(with(subset(St4, grepl("D", lex.Tr)), table(lex.Tr)))
    
    
    ###################################################
    ### code chunk number 69: ms-steno2.rnw:1275-1283
    ###################################################
    stD <- glm(cbind(lex.Fail, lex.dur)
               ~  Ns(tfi, knots = seq( 0, 20,  5)) * lex.Tr +
                 Ns(age, knots = seq(50, 80, 10)) * lex.Tr +
                 lex.Tr / allo + sex,
               family = poisreg,
               offset = log(lex.dur),
               data = subset(St4, grepl("D", lex.Tr)))
    round(ci.exp(stD)[,1,drop=F],3)
    