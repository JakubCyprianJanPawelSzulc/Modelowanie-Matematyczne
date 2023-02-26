library(moments)
library(fitdistrplus)
library(ggplot2)
library(ggExtra)
library(mnormt)
library(MASS)
library(QRM)
library(evir)

dane <- read.csv("C:/Modelowanie/car_d.csv")
dane_zamkniecie <- dane$Zamkniecie
dane2 <-read.csv("C:/Modelowanie/dbc_d.csv")
dane_zamkniecie2 <- dane2$Zamkniecie

#1
plot(dane_zamkniecie)
hist(dane_zamkniecie, prob=TRUE)

#2
srednia <- mean(dane_zamkniecie)
odchylenie <- sd(dane_zamkniecie)
skosnosc <- skewness(dane_zamkniecie)
kurtoza <- kurtosis(dane_zamkniecie)

#Skośność dodatnia, czyli ma długi ogon z prawej strony.
#Kurtoza dodatnia oznacza, że w danych jest więcej skrajnych wartości odstających niż w rozkładzie normalnym.

#3

#rozkład normalny

fit <- fitdist(dane_zamkniecie,"norm")
fit
plot(fit)

#rozkład wykładniczy

fit_1 <- fitdist(dane_zamkniecie,"exp")
fit_1
plot(fit_1)

#rozkład log-norm

fit_2 <- fitdist(dane_zamkniecie,"lnorm")
fit_2
plot(fit_2)

#4

gofstat(list(fit, fit_1, fit_2), fitnames=c('normal', 'exponential', 'lognormal'))

#wybieram lognormal bo najmniejsze wartosci

denscomp(list(fit, fit_1, fit_2), legendtext = c('normal', 'exponential', 'lognormal'))
qqcomp(list(fit, fit_1, fit_2), legendtext = c('normal', 'exponential', 'lognormal'))
cdfcomp(list(fit, fit_1, fit_2), legendtext = c('normal', 'exponential', 'lognormal'))
ppcomp(list(fit, fit_1, fit_2), legendtext = c('normal', 'exponential', 'lognormal'))

#5

N <- 10000
n <- length(dane_zamkniecie); n


Dln <- c()

for (i in 1:N) { 
  
  Yln <- rlnorm(n,fit_2$estimate[1],fit_2$estimate[2])
  
  Dln[i] <-  ks.test(Yln,plnorm, fit_2$estimate[1],fit_2$estimate[2],exact=TRUE)$statistic
}

#2. Obliczamy dn, czyli wartosc statystyki Dn, dla danych Loss i rozkładu F0 wybranego w punkcie A.
dn_ln <-  ks.test(dane_zamkniecie,plnorm,fit_2$estimate[[1]],fit_2$estimate[[2]],exact=TRUE)$statistic
dn_ln

#wyniki z punktow 1.2 na histogramie
par(mfrow=c(1,1))
hist(Dln,prob=T)
points(dn_ln,0,pch=19,col=2)

#Odleglosc dystrybuanty empirycznej (wartosc statystyki Dn) dla Loss, oraz dystrybuanty F0
# jest istotnie większa od odleglosci obserwowanych dla probek tej samej licznosci z rozkladu F0.

#3. Obliczamy p-value.
p_value_ln <- length(Dln[Dln>dn_ln])/N; p_value_ln

#4. Przyjmujemy poziom istotnosci alpha=0.05
alpha <- 0.05
p_value_ln <= alpha

#Wartosc p-value jest mniejsza od przyjetego poziomu istotnosci, 
#zatem hipoteze o rownosci dystrybuant (F=F0, gdzie F poszukiwany rozklad) odrzucamy.




#6
#diff


#A
#1
Intercars <- diff(log(dane_zamkniecie))
Debica <- diff(log(dane_zamkniecie2))

df <- data.frame(Intercars=Intercars,Debica=Debica)
p <-  ggplot(df, aes(x=Intercars, y=Debica)) + geom_point()
ggMarginal(p, type="histogram")

#to jest super

#2
mu <- colMeans(df); mu #wektor srednich

Sigma <- cov(df); Sigma #macierz kowariancji, estymator nieobciążony

Sigma[2]

mi1 <- mu[1];mi1
mi2 <- mu[2];mi2


a1=sqrt(0.0008403938); a1
a2=sqrt(0.0001667264); a2

p=0.0001193947/(a1*a2); p #Wspolczynnik korelacji

P <- cor(df); P #macierz korelacji

#3
#N (srednia1, srednia2, odchylenie_R, odchylenie_R2, wspolczynnik_korelacji)

s1 <- sd(Intercars)
s2 <- sd(Debica)

x <- seq(-3*s1, 3*s1, 0.005) 
y <- seq(-3*s2, 3*s2, 0.005)

#gestosc rozkladu normalnego o sredniej mu i macierzy kowariancji S

f <- function(x, y) dmnorm(cbind(x, y), mu, Sigma)  
z <- outer(x, y, f) 

persp(x,y,z,theta = 30,phi = 30,col="lightblue")


persp(x, y, z, theta = -30, phi = 25, 
      shade = 0.75, col = "lightblue", expand = 0.5, r = 2, 
      ltheta = 25, ticktype = "detailed")

#B
#1
n <- nrow(df); n
set.seed(100)
Z <- MASS::mvrnorm(n,mu=mu,Sigma=Sigma);Z

par(mfrow=c(1,2))
plot(df, xlim=c(-0.15,0.15),ylim=c(-0.10,0.10))
plot(Z,xlim=c(-0.15,0.15),ylim=c(-0.10,0.10))

#2
Z <- rmnorm(n,mu,Sigma)
dM <- mahalanobis(df,mu,Sigma)
dM_probki <- mahalanobis(Z,mu,Sigma)

hist(dM)

n <- dim(df)[1]; n
alpha <- ppoints(n)
q_emp <- quantile(dM,alpha)
q_teo <- qchisq(alpha,df=2)

plot(q_emp,q_teo,pch=19)
abline(a=0,b=1,col=2)

q_emp_probki <- quantile(dM_probki, alpha)
q_teo_probki <- qchisq(alpha, df=2)

plot(q_emp_probki, q_teo_probki, pch=19)
abline(a=0, b=1,col=2)

ks.test(dM, 'pchisq', 2, exact=TRUE)
ks.test(dM, 'pchisq', 2)

#####################################################
#Regresja Liniowa

#Intercars
Intercars <- diff(log(dane_zamkniecie))
alpha <- 0.05
wo_Intercars <- mean(Intercars)
sd_Intercars <- sd(Intercars)
x <- qnorm(1-alpha/2, 0, 1)
lower_CI_Intercars <- wo_Intercars - (x * (sd_Intercars/sqrt(length(Intercars))))
upper_CI_Intercars <- wo_Intercars + (x * (sd_Intercars/sqrt(length(Intercars))))
length_CI_Intercars <- upper_CI_Intercars - lower_CI_Intercars
c(lower=lower_CI_Intercars, upper=upper_CI_Intercars, length=length_CI_Intercars)


#Debica
Debica <- diff(log(dane_zamkniecie2))
alpha <- 0.05
wo_Debica <- mean(Debica)
sd_Debica <- sd(Debica)
x <- qnorm(1-alpha/2, 0, 1)
lower_CI_Debica <- wo_Debica - (x * (sd_Debica/sqrt(length(Debica))))
upper_CI_Debica <- wo_Debica + (x * (sd_Debica/sqrt(length(Debica))))
length_CI_Debica <- upper_CI_Debica - lower_CI_Debica
c(lower=lower_CI_Debica, upper=upper_CI_Debica, length=length_CI_Debica)



library(ggplot2)
library(boot)  


df <- data.frame(Intercars, Debica)

#obejrzyjmy dane na wykresie
qplot(Intercars,Debica, data = df,
      main = "") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(colour = "blue", size = 1.5) 


#estymatory wspolczynnikow
beta1 <- cov(Intercars,Debica)/var(Intercars)
beta0 <- mean(Debica)-mean(Intercars)*beta1
beta1; beta0

#linia regresji na  wykresie
qplot(Intercars, Debica, data = df,
      ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(colour = "blue", size = 1.5) +
  geom_abline(intercept = beta0, slope = beta1, color="red",size=1)



data <- lm(Debica~Intercars,data=df)
data

summary(data)
confint(data)
m_p <- plot(data); m_p

sum <- summary(data)
sum

#reszty (residuals) 
reszty <- data$residuals; reszty

#histogram i qq-ploty
hist(reszty)
qqnorm(reszty)
qqline(reszty,col=2)

m <- mean(reszty); m
s <- sd(reszty); s

ks.test(reszty,'pnorm',m,s)

#test Shapiro-Wilka
shapiro.test(reszty)

#p-value=0.1381, na poziomie 5% nie ma podstaw
#do odrzucenia hipotezy o normalnosci rozkladu reszt

#RSE - blad standardowy reszt
RSE <- sqrt(sum(reszty^2)/(length(Intercars)-2))
RSE

m <- mean(Intercars)
beta0+beta1*m

summary(data)$r.squared 


##Predykcja

df <- data.frame(Intercars=Intercars,Debica=Debica)

data <- lm(Debica~Intercars-1,data=df)
data 

sum2 <- summary(data)
sum2 

#Kiedy m to średnia z danych
m <- mean(Intercars)

beta1_model2 <- data$coefficients

beta1_model2*m


#Przyklad 1.g (Predykcja i przedzialy ufnosci dla predykcji)
#=========
nowe_dane <- data.frame(Intercars=m)

predict(data, nowe_dane, interval="confidence") #model 2




