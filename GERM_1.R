G_R_1 <- read.csv("C:/Users/uqmsabam/OneDrive - The University of Queensland/RESEARCH/M-2/R-FILES/GERMINATION/G_R_1.csv")
View(G_R_1)
library(MASS)
library(Matrix)
library(lme4)
library(tidyverse)
library(ggplot2)
library(tibble)
library(tidyr)
library(readr)
library(nls2)
library(proto)


G_R_1 <- read.csv("C:/Users/uqmsabam/OneDrive - The University of Queensland/RESEARCH/M-2/R-FILES/GERMINATION/G_R_1.csv")
comment(G_R_1)
levels(G_R_1$GTP)
str(test <- subset(G_R_1, GTP =="QPL1001"))
test2 <- subset(test,G_R>0)
ggplot(test, aes(x=TEM, y=G_R)) + geom_point()
summary(lm1 <- lm(G_R ~ TEM, data = test))

d <- 30
coefficients(lm(G_R ~ TEM, data = test, subset = TEM < d))
coefficients(lm(G_R ~ TEM, data = test, subset = TEM >= d))

nls1 <- nls(G_R ~ (b1 + b2 * TEM)*(TEM<d) + 
             (b1 + (b2 - b3)*d + b3*TEM)* (TEM>=d),
           start=list(b1=-.34,b2=0.04, b3=-0.11, d=30), data=test)
summary(nls1)

nls2 <- nls(G_R ~ (b1 + b2 * TEM)*(TEM<d) + 
              (b1 + (b2 - b3)*d + b3*TEM)* (TEM>=d),
            start=list(b1=-.34,b2=0.04, b3=-0.11, d=30), data=test2)
summary(nls2)
coeffnls2 <- coefficients(nls2)
(minnls2 <-  -coeffnls2[1]/coeffnls2[2])
(maxnls2 <- -(coeffnls2[1]-coeffnls2[3]*coeffnls2[4])/(coeffnls2[2]+coeffnls2[3]))
(ttsubnls2 <- 1/coeffnls2[2])
(ttsupnls2 <- -1/(coeffnls2[2]+coeffnls2[3]))

x.seq <- seq(5,50, by=0.1)
p1 <- predict(nls1, newdata=list(TEM=x.seq), interval = "confidance")
plot(G_R ~ TEM, data=test, ylim = c(0, 1.3))
lines(x = x.seq, y = p1)

p2 <- predict(nls2, newdata=list(TEM=x.seq), interval = "confidance")
plot(G_R ~ TEM, data=test, ylim = c(0, 1.3))
lines(x = x.seq, y = p2)

df1 <- data.frame(x = x.seq, y = p1)
ggplot(test, aes(x=TEM, y=G_R)) + geom_point() + geom_path(aes(x = x.seq, y = p1), data = df1)


#Adapting code to document
nls2 <- nls(G_R ~ (b1 + b2 * TEM)*(TEM<d) + 
              (b1 - b3*d + (b2+b3)*TEM)* (TEM>=d),
            start=list(b1=-.34,b2=0.04, b3=-0.11, d=30), data=test2)
summary(nls2)
coeffnls2 <- coefficients(nls2)
(minnls2 <-  -coeffnls2[1]/coeffnls2[2])
(maxnls2 <- -(coeffnls2[1]-coeffnls2[3]*coeffnls2[4])/(coeffnls2[2]+coeffnls2[3]))
(ttsubnls2 <- 1/coeffnls2[2])
(ttsupnls2 <- -1/(coeffnls2[2]+coeffnls2[3]))

#multiple gtp

Gtest <- subset(G_R_1,G_R>0)

str(test <- subset(G_R_1, GTP =="QPL1001"))
test2 <- subset(test,G_R>0)
ggplot(test, aes(x=TEM, y=G_R)) + geom_point()
summary(lm1 <- lm(G_R ~ TEM, data = test))
str(test3 <- subset(G_R_1, GTP =="Quest"))
test4 <- subset(test,G_R>0)
d <- 32
coefficients(lm(G_R ~ TEM, data = test3, subset = TEM < d))
coefficients(lm(G_R ~ TEM, data = test3, subset = TEM >= d))

startcoeff <- data.frame(b1=c(-.34,-0.43),b2=c(0.04,0.04),b3=c(-0.11,-0.12),d=c(30,32))

library(nlsLoop)
nlsCrazy <- nlsLoop(G_R ~ (b1 + b2 * TEM)*(TEM<d) + 
                       (b1 - b3*d + (b2+b3)*TEM)* (TEM>=d),id_col = "GTP",
                    tries=500,
                    supp_errors = "Y",
                    data=Gtest)
nlsCrazy
