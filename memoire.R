library(readxl)
# CORRECTION : chemin en dur remplacé par un chemin relatif au dossier de travail.
# Placez meanstemp.xlsx dans le meme dossier que ce script, ou modifiez le chemin ci-dessous.
meanstemp <- read_excel("meanstemp.xlsx")
library(styler)
library(ggplot2)
library(forecast)
library(ggfortify)
library(tidyverse)
library(timeSeries)
library(xts)
library(ggpmisc)
library(zoo)
# CORRECTION : library(roll) supprimée -> package non utilisé
# (rollapply() provient bien du package zoo, chargé ci-dessus).
library(GA) #R-package for genetic algorithm (on s'en servira plus tard pour la détection de ruptures)
library(microbenchmark)

#meanstemp est le tibble tel qu'importé du MetOffice
meanstemp <- meanstemp[-363,]   #Suppression de l'année 2021, les données ne sont pas encore mises à jour
#Comme ce sont des chaines de caractere, je les transforme en numeric, au cas où on en aura besoin apres
meanstemp[,"JAN"] <- as.numeric(unlist(meanstemp[,"JAN"]))
meanstemp[,"FEB"] <- as.numeric(unlist(meanstemp[,"FEB"]))
meanstemp[,"MAR"] <- as.numeric(unlist(meanstemp[,"MAR"]))
meanstemp[,"APR"] <- as.numeric(unlist(meanstemp[,"APR"]))
meanstemp[,"MAY"] <- as.numeric(unlist(meanstemp[,"MAY"]))
meanstemp[,"JUN"] <- as.numeric(unlist(meanstemp[,"JUN"]))
meanstemp[,"JUL"] <- as.numeric(unlist(meanstemp[,"JUL"]))
meanstemp[,"AUG"] <- as.numeric(unlist(meanstemp[,"AUG"]))
meanstemp[,"SEP"] <- as.numeric(unlist(meanstemp[,"SEP"]))
meanstemp[,"OCT"] <- as.numeric(unlist(meanstemp[,"OCT"]))
meanstemp[,"NOV"] <- as.numeric(unlist(meanstemp[,"NOV"]))
meanstemp[,"DEC"] <- as.numeric(unlist(meanstemp[,"DEC"]))
meanstemp[,"YEAR"] <- as.numeric(unlist(meanstemp[,"YEAR"]))
meanstempts <- ts(meanstemp["YEAR"], start = 1, frequency = 1)
annual <- meanstemp %>% select(Annee, YEAR)  #data.frame
annualts <- as.ts(annual)  

toto <- as.data.frame(meanstemp) #transformation du tb en data.frame

plot(meanstemp[,"Annee", drop = TRUE], meanstemp[,"YEAR", drop = TRUE], type="l")
ggplot(data = as_tibble(toto), aes(x = Annee, y = YEAR)) + geom_line() 
pacf(toto$JAN) #Autocorrélation Janvier toutes années    #ceci n'a aucun sens
pacf(toto$JUL) #Autocorrélation Juillet toutes années    #ceci aussi 

acf(toto$YEAR) #Autocorrélation moyennes annuelles

ggAcf(toto, lag.max = NULL, plot = TRUE) #Plot autocorrélation
autoplot(annualts, facets = TRUE)
autoplot(annualts, facets = FALSE)
#ts.plot(meanstempts) 
k <- 50
moyennemobile <- rollapply(annualts[,2], width = 2*k+1, mean, align = "center")
dates<-annualts[,1]
mb<-data.frame(dates[(1:length(moyennemobile))+k],moyennemobile) #Il faut commencer à k+1
plot(mb[,1],mb[,2], type="l", col = "red")
# mb <- ma(annualts, order = 100, centre = FALSE)

ecarttypemobile <- rollapply(annualts[,2],  width = 2*k+1, sd, align="center")
binf<-moyennemobile-(2/sqrt(2*k+1))*ecarttypemobile
bsup<-moyennemobile+(2/sqrt(2*k+1))*ecarttypemobile
ic<-data.frame(dates[(1:length(ecarttypemobile))+k],binf,bsup)

ggplot(ic, aes(ic[,1])) +                    # basic graphical object
  geom_line(aes(y=ic[,2]), colour="red") +  # first layer
  geom_line(aes(y=ic[,3]), colour="green")+
  geom_line(aes(y=moyennemobile))

par(mfrow=c(1,1))
plot(meanstemp[,"Annee", drop = TRUE], meanstemp[,"YEAR", drop = TRUE], type="l")
lines(mb[,1],mb[,2], type="l", col = "red")

#Problème d'echelle pour la superposition les graphes

lines(ic[,1],ic[,2],type="l",col="blue")
lines(ic[,1],ic[,3],type="l",col="green")


E<-summary(lm(meanstempts~time(meanstempts), data=meanstempts))

# f fonction de regression avec les coeff estimés
# CORRECTION : "a" et "b" étaient inversés (l'intercept était multiplié par x,
# et la pente ajoutée comme constante), ce qui donnait des valeurs aberrantes
# (ex. ~3169°C en fin de série au lieu de ~9.8°C). Renommage explicite pour éviter
# toute confusion, et pente correctement multipliée par x.
f <- function(x,E) {
  intercept <- E$coefficients[1,1]
  pente     <- E$coefficients[2,1]
  return(pente*x + intercept)}
autoplot(meanstempts)
plot(meanstempts)
# CORRECTION : l'appel original f(seq(1,362,0.01)) ne passait pas E, ce qui
# provoquait une erreur ("argument E is missing, with no default").
lines(seq(1,362,0.01),f(seq(1,362,0.01),E),col="red")


meanstempts.trend.corr <- meanstempts-predict(lm(meanstempts~time(meanstempts), data=meanstempts))
summary(lm(meanstempts.trend.corr ~time(meanstempts.trend.corr ), data=meanstempts.trend.corr ))
autoplot(meanstempts.trend.corr, colour="blue")



### Programmation dynamique pour la détection des points de rupture dans la moyenne :
###ou dans la tendance

##Fonction de coût pour les ruptures en moyennes :

cost <- function(y,s,t) { #y est un vecteur, s et t des indices
  if (s<t) {
    yst <- y[s:t] 
    return(sum((yst-mean(yst))^2))
  }
  else {return(0)}   #Triangulaire supérieure
}

##Matrice de coût pour les ruptures en moyennes : 
cout<- function(y){ #y est un vecteur
  J <- length(y)
  Q0 <- matrix(numeric(J*J),ncol = J, nrow = J)
  for(i in 1:J) {
    for(j in i:J){
      Q0[i,j] <- cost(y,i,j)
    }
  }
  return(Q0)
}

cost(annual[[2]],1,329)+cost(annual[[2]],330,362)


##Estimation des coefficients de regression linéraire pour les ruptures dans la tendance :


estim_a_b<-function(y,s,t){# y est un data frame (pour avoir les estimations en fonctions des années), s et t des indices 
  data <- data.frame(y = y[[2]][s:t], x = y[[1]][s:t]) #y[[2]] sont dans notre cas par exemple les moyennes de temperature, y[[1]] les années
  E<-summary(lm(y~x, data = data)) 
  if (s+1<t){#(lorsque t >leght(X) fonction renvoie meme cout que pour t =lenght(X))
    b<-E$coefficients[1,1] 
    a<-E$coefficients[2,1]
  }
  else {
    b<-E$coefficients[1,1]
    a<-0
  }
  return(c(a,b))
}



##Fonction de coût pour les ruptures dans la tendance: 

cost_t<-function(y,s,t){ #y est un data frame, s et t des indices
  if(s<t){
    res <- estim_a_b(y,s,t) 
    a<-res[1]
    b<-res[2]
    M<-y[[2]][s:t]
    I<-y[[1]][s:t]
    return(sum((M-(a*I+b))^2))}
  else{return(0)}
}    

#Matrice de coût pour les ruptures dans la tendance :
cout_t<- function(y){ #y est un data frame
  J <- length(y[[2]])
  Q0 <- matrix(numeric(J*J),ncol = J, nrow = J)
  for(i in 1:J) {
    for(j in i:J){
      Q0[i,j] <- cost_t(y,i,j)
    }
  }
  return(Q0)
}


#Fonction auxiliaire pour le sapply :

fdynamic1<-function(j,k,Q,y,Q0){ #Cette fonction ne sert qu'à remplacer le for de prog_dyn ci dessus par un sapply 
  return(min(Q[(k-1),1:(j-1)]+ Q0[(2:j),j]))
}


#Fonction qui retourne toutes les configurations à k points de ruptures 
#et le coût de parcours de la série avec k ruptures pour tous les k dans 1:K

prog_dyn <- function(y,K,Q0) { # y est un vecteur, K le nombre maximal de ruptures à considérer, Q0 une matrice de cout 
  t <- sapply(c(1:K), numeric)
  J <- length(y)
  Q <- matrix(numeric((K+1)*J), nrow = K+1 , ncol = J)
  Q[1,] <- Q0[1,]
  for(k in 2:(K+1)) {
    Q[k,-1]<-t(sapply(c(2:J),fdynamic1,k,Q,y,Q0)) #Les couts dans une configuration à k points de rupture
    t[[k-1]][1] <- which.min(Q[(k-1),1:(J-1)]+ Q0[(2:J),J]) #Dernier point de rupture
    if (k>2){
      for (i in 2:(k-1)) { #Les autres points de ruptures par récurrence
        t[[k-1]][i] <- which.min(Q[(k-i),1:(t[[k-1]][i-1]-1)]+ Q0[(2:t[[k-1]][i-1]),t[[k-1]][i-1]])
      }
      
      
    }
  }
  return(list(configuration_optimale=t,couts=Q[-1,J])) #L'output est une liste contenant toutes les configation à k ruptures
  #k in 1:K, et les couts de parcours de la série correspondants
}


#Estimation des moyennes sur chaque segment :

estim_mean <- function (y,k,s) {#y est un vecteur, k un réel (ce serait celui qui minimisera le bic) et Q lé résultat de prog_dyn
  J <- length(y)
  L <- numeric(k+1) #un vecteur de 0 à k+1 composantes, on le remplit composante par composante par les moyennes sur chaque segment
  if (k>1) {     #autrement le for ne marche pas, donc on fait le cas k=1 à part
    L[1] <- mean(y[1:s[1]])      #on estime les moyennes sur chaque segment
    for(i in 1:(k-1)) {
      L[i+1] <- mean(y[(s[i]+1):s[i+1]])
    }
    L[k+1] <- mean(y[(s[k]+1):J])
  }
  else { #Cas k=1 i.e une seule rupture
    L[1] <- mean(y[1:s[1]])
    L[2] <- mean(y[(s[1]+1):J])
  }
  return(L)
}

#Estimation de la tendance sur chaque segment :

#La fonction ci dessous renvoie une matrice (p+1)*2 où sont stockés les estimations de a et b
#sur chaque segment d'une configuration à p ruptures
estim_trend <- function(y,p,s) {#y est un data frame, s un vecteurs de ruptures, p le nombre de ruptures
  J <- length(y[[2]])
  M <- matrix(numeric(2*(p+1)), nrow = p+1, ncol = 2) #Matrice de taille (p+1)*2 de 0, qu'on remplira par les coef
  #On mettra les coefficients a sur la premiere colonne et les intercepts (b) sur la seconde
  res1 <- estim_a_b(y,1,s[1])
  M[1,1] <- res1[1]
  M[1,2] <- res1[2]
  if (p>1) { #Le cas p=1 sera traité à part
    for(i in 1:(p-1)) {
      resi <- estim_a_b(y,(s[i]+1),s[i+1]) 
      M[(i+1),1] <- resi[1]
      M[(i+1),2] <- resi[2]
    }
    resf <- estim_a_b(y,s[p]+1,J)
    M[(p+1),1] <- resf[1]
    M[(p+1),2] <- resf[2]

  }
  else { #Cas p=1 ie une seule rupture
    res2 <- estim_a_b(y,(s[1]+1),J)
    M[2,1] <- res2[1]
    M[2,2] <- res2[2]
  }
  return(M)
}




#Pénalisation pour choisir la mailleure configuration par le BIC : 

bicc <- function(k,y,Q,W) { #y est un data frame, k un nombre de ruptures donné, Q un resultat de prog_dyn 
  #et W prend "mean" ou "trend" dependemment de si on fait l'etude des ruptures en moyenne ou en tendance 
  J <- length(y[[2]])
  if (W=="mean"){
    return(Q[[2]][k]+(k+1)*log(J))  #Q[[2]][k] est le cout de parcours de la série en presence de k ruptures estimées
  }
  if (W=="trend") {
    return(Q[[2]][k]+2*(k+1)*log(J)) #*2 car il y a 2 paramètre qu'on estime dans le modèle (a et b)
  }
}



#Estimation de la variance des erreurs :

estim_var<-function(y,k,s,W){ #y est un data frame, k un nombre de ruptures donné, Q un resultat de prog_dyn 
  #et W prend "mean" ou "trend" dependemment de si on fait l'etude des ruptures en moyenne ou en tendance 
  J <- length(y[[2]])
  if (W=="mean") { #Dans le cas de ruptures en moyenne
    m <- estim_mean(y[[2]],k,s) #vecteur des moyenne sur les segments de la configuration s à k ruptures
    G <- numeric(k+1) #vecteur de 0
    if (k>1) { #On traitera le cas k=1 à part 
      G[1] <-sum((y[[2]][1:s[1]]-m[1])^2)
      for (i in 1:(k-1)) {
        G[i+1] <- sum((y[[2]][(s[i]+1):s[i+1]]-m[i+1])^2)
      }
      G[k+1] <-sum((y[[2]][(s[k]+1):J]-m[k+1])^2)
    }
    else { #Cas k=1 
      G[1] <-sum((y[[1]][1:s[1]]-m[1])^2)
      G[2] <-sum((y[[1]][(s[1]+1):J]-m[2])^2)
    }
    return((1/J)*sum(G))
  }
  
  if (W=="trend") {
    coef <- estim_trend(y,k,s) #Matrice de coefficients de regression linéaire sur chaque segment d'une configuration
    #à k ruptures, les coef a son sur la premiere colonne et les b sur la seconde (c'est une matrice (k+1)*2)
    G <- numeric(k+1) #vecteur de 0
    if (k>1) { #On traitera le cas k=1 à part
      G[1] <-sum((y[[2]][1:s[1]]-coef[1,1]*y[[1]][(1:s[1])]-coef[1,2])^2)
      for (i in 1:(k-1)) {
        G[i+1] <- sum((y[[2]][(s[i]+1):s[i+1]]-coef[(i+1),1]*y[[1]][((s[i]+1):s[i+1])]-coef[(i+1),2])^2)
      }
      G[k+1] <-sum((y[[2]][(s[k]+1):J]-coef[(k+1),1]*y[[1]][((s[k]+1):J)]-coef[(k+1),2])^2)
      
    }
    else { #Cas k=1 ie dans le cas d'une seule rupture
      G[1] <- sum((y[[2]][1:s[1]]-coef[1,1]*y[[1]][(1:s[1])]-coef[1,2])^2)
      G[2] <- sum((y[[2]][(s[1]+1):J]-coef[2,1]*y[[1]][((s[1]+1):J)]-coef[2,2])^2)
    }
    return((1/J)*sum(G))
    
  }
  
}



#Correction du BIC par la variance :

biccorr <- function(k,y,Q,W) { #y est un data frame, k un nombre de ruptures donné, Q un resultat de prog_dyn 
  #et W prend "mean" ou "trend" dependemment de si on fait l'etude des ruptures en moyenne ou en tendance 
  J <- length(y[[2]])
  # CORRECTION : Q[[1]][k] (simple crochet) renvoie une SOUS-LISTE de longueur 1
  # au lieu du vecteur de ruptures lui-meme -> plantage en aval ("NA/NaN argument"
  # dans estim_var/estim_trend/estim_a_b). Q[[1]][[k]] (double crochet) renvoie
  # le vecteur numerique attendu.
  s <- rev(Q[[1]][[k]])
  if (W=="mean") {
    return(J*log(estim_var(y,k,s,W))+(k+1)*log(J)) 
  }
  if(W=="trend") {
    return(J*log(estim_var(y,k,s,W))+2*(k+1)*log(J))
  }
  
}


#Le plot :
#Le graphe fonctionne assez bien pour le moment même si on a beaucoup trop de ruptures
graphe <- function(y,p,rupt,W){ #y est un data frame, p un nombre de ruptures donné, Q un resultat de prog_dyn 
  #et W prend "mean" ou "trend" dependemment de si on fait l'etude des ruptures en moyenne ou en tendance 
  s <- y[[1]][1]+rupt-1 #ici rev(Q[[1]][[p]] est le vecteur des points de ruptures dans un ordre croissant
  #auquel on ajoute y[[1]][1]-1 pour avoir les points de ruptures correspondants à notre série (par exemple des années entre 1659 et 2020)
 print(s)
   if (W=="mean") {
    m <- estim_mean(y[[2]],p,rupt)
    print(m)
    plot(y[[1]],y[[2]], type='l')
    lines(y[[1]][1]:s[1],rep(m[1],(s[1])-y[[1]][1]+1),col= "darkgreen", lwd=4)
    abline(v=s[1],lty='dotted',col="darkgreen",lwd=2)
    for(i in 2:p) {
      lines((s[i-1]:s[i]),rep(m[i],(s[i]-s[i-1])+1),col="darkgreen", lwd=4)
      abline(v=s[i],lty='dotted',col="darkgreen",lwd=2)
    }
    
    lines(s[p]:(length(y[[1]])+y[[1]][1]), rep(m[p+1],y[[1]][1]+length(y[[1]])-s[p]+1), col="darkgreen",lwd=4)
  }
  if (W=="trend") {
    coef <- estim_trend(y,p,rupt)
    plot(y[[1]],y[[2]], type='l')
    lines(x <- (y[[1]][1]:s[1]),(coef[1,1]*x+coef[1,2]),col= "turquoise4", lwd=3)
    abline(v=s[1],lty='dotted',col="turquoise4",lwd=2)
    if (p>1) {
    for(i in 2:p) {
      lines(x <- (s[i-1]:s[i]),(coef[i,1]*x+coef[i,2]),col="turquoise4", lwd=3)
      abline(v=s[i],lty='dotted',col="turquoise4",lwd=2)
    }
    }
    
    lines(x <- s[p]:(length(y[[1]])+y[[1]][1]), coef[(p+1),1]*x+coef[(p+1),2], col="turquoise4",lwd=3)
    
    
  }
}






###Application à des données simulées :
##Moyenne :
X <- data.frame(1:400,c(rnorm(100,0,1),rnorm(100,5,1),rnorm(100,10,1),rnorm(100,2,1)))
J<-length(X)
coutX <- cout(X[,2])
Xkmax <- 10
progX <- prog_dyn(X[,2],Xkmax,coutX)
r <- which.min(sapply(c(1:Xkmax),bicc,X,progX,"mean"))

graphe(X,r,rev(progX[[1]][[r]]),"mean")
r

Z <- data.frame(1:400,c(runif(100,0,1),runif(100,3,4),runif(100,1,2),runif(100,0,1)))
J1 <- length(Z)
coutZ <- cout(Z[,2])
Zkmax <- 6
progZ <- prog_dyn(Z[,2],Zkmax,coutZ)
r1 <- which.min(sapply(c(1:Zkmax),bicc,Z,progZ,"mean"))

graphe(Z,r1,rev(progZ[[1]][[r1]]),"mean")
r1

##Tendance : 
X <- 1:400
Y <- data.frame(X,c(2*X[1:100]+rnorm(100,0,10), 3*X[101:300]+2+rnorm(200,0,10),4*X[301:400]+10+rnorm(100,0,10)))


plot(Y,type="l")

coutY <- cout_t(Y)   #Cout_t prend un peu de temps, dans les environs de 3 minutes
Ykmax <- 10
progY <- prog_dyn(Y[,2],Ykmax,coutY)
r <- which.min(sapply(c(1:Ykmax),biccorr,Y,progY,"trend"))
graphe(Y,r,rev(progY[[1]][[r]]),"trend")


### Application à notre série CET :

annual <- meanstemp %>% select(Annee, YEAR) 

annual <- as.data.frame(annual)
Kmax <- 10

##Moyenne :

P0 <- cout(annual[[2]])
P <- prog_dyn(annual[[2]],Kmax,P0)
ro <- which.min(sapply(c(1:Kmax),bicc,annual,P,"mean"))
graphe(annual,ro,rev(P[[1]][[ro]]),"mean")


##Tendance :

Q <- cout_t(annual)    #Ceci prendra un peu de temps, quelques minutes, 3 ou 4 environ
P1 <- prog_dyn(annual[[2]],Kmax,Q)
ro1 <- which.min(sapply(c(1:Kmax),biccorr,annual,P1,"trend"))
graphe(annual,ro1,rev(P1[[1]][[ro1]]),"trend")


##Autocorrélation des erreurs :

#Vecteur des erreurs :

erreur <- function(y,r,Q,W) { #data.frame
  J <- length(y[[2]])
  conf <- rev(Q[[1]][[r]])
  if(W=="mean") {
    m <- estim_mean(y[[2]],r,conf)
    v <- y[[2]][1:conf[1]]-m[1]
    
    i <- 2
    while(i<r+1) {
      v <- c(v,y[[2]][(conf[i-1]+1):conf[i]]-m[i])
      i <- i+1
    }
    v <- c(v,y[[2]][(conf[r]+1):J]-m[r+1])
    return(v)
  }
  
  if(W=="trend") {
    coef <- estim_trend(y,r,conf)
    v <- y[[2]][1:conf[1]]-coef[1,1]*y[[1]][1:conf[1]]-coef[1,2]
    i <- 2
    while(i<r+1) {
      v <- c(v,y[[2]][(conf[i-1]+1):conf[i]]-coef[i,1]*y[[1]][(conf[i-1]+1):conf[i]]-coef[i,2])
      i <- i+1
    }
    v <- c(v,y[[2]][(conf[r]+1):J]-coef[(r+1),1]*y[[1]][(conf[r]+1):J]-coef[(r+1),2])
    return(v)
  }
}

erreur(annual,ro1,P1,"trend")

erreur(annual,ro,P,"mean")


#Autocorrélogramme des erreurs :

#En moyenne :
pacf(erreur(annual,ro,P,"mean"))

#En tendance : 
pacf(erreur(annual,ro1,P,"trend"))



###Détection des ruptures dans la moyenne via l'algo génétique :



#Fonction auxiliaire pour le sapply dans costga()
# CORRECTION : la ligne "J <- length(v)" référençait une variable "v" qui n'est
# ni un argument de rec() ni une variable globale définie -> selon la portée
# lexicale de R, ceci provoquait une erreur ("object 'v' not found") dès le
# premier appel de costga()/costga2() dans une session propre. Cette ligne
# était de toute façon un reliquat mort (J n'était jamais utilisé) : supprimée.
rec <- function(i, s, y) {
  return(cost(y, (s[i]+1), s[i+1]))
}

costga <- function(v,y,minseglen) {
  k <- sum(v)   #Nb de ruptures
  J <- length(y)
  s <- which(v==1)*minseglen  #Vecteur d'emplacements des points de ruptures
  s <- c(1, s, J)
  P <- sapply(1:(k+1), rec, s=s, y=y)
  return(-(sum(P)+(k+1)*log(J)))  
  
}



#Fonction de cout 2 (BIC) : 
costga2 <- function(v,y,minseglen) {
  k <- sum(v)   #Nb de ruptures
  J <- length(y)
  s <- which(v==1)*minseglen  #Vecteur d'emplacements des points de ruptures
  s <- c(1, s, J)
  P <- sapply(1:(k+1), rec, s=s, y=y)
  return(-J*log(sum(P)/J)-(k+1)*log(J))
}

X <- data.frame(1:400,c(rnorm(100,0,1),rnorm(100,5,1),rnorm(100,10,1),rnorm(100,2,1)))

minseglen = 5

microbenchmark(resultX2 <- ga(type = "binary",fitness =costga,X[[2]],minseglen, nBits = length(X[[2]])%/%minseglen,
                              population = gabin_Population, maxiter = 300, popSize = 500, pmutation = 0.2, pcrossover = 0.9), times=1)

sum(resultX2@solution)
which(resultX2@solution==1)*minseglen

graphe(X,sum(resultX2@solution),which(resultX2@solution==1)*minseglen,"mean")



temp <- annual[[2]]

minseglen <- 10

microbenchmark(result <- ga(type = "binary",fitness =costga2,temp,minseglen, nBits = length(temp)%/%minseglen, population = gabin_Population, maxiter = 400 , popSize = 250 ), times=1)

sum(result@solution)
which(result@solution==1)*minseglen+1659-1


nrupt <- sum(result@solution)
rupt <- which(result@solution==1)*minseglen

graphe(annual,nrupt,rupt,"mean")

plot(result)



##    Les ruptures dans la tendance     ##

#Fonction auxiliaire pour costga_t et costga_t2
rec_t <- function(i, s, y) { #i un indice, s un vecteur d'emplacements de ruptures, y un data.frame
  return(cost_t(y, (s[i]+1), s[i+1]))
}

#Fonction de cout à maximiser par ga(), correspond à (-BIC)
costga_t <- function(v,y,minseglen) { #v un vecteur binaire, y un data.frame, 
  #minseglen le nombre de points minimal entre 2 ruptures successives
  k <- sum(v)   #Nb de ruptures
  J <- length(y[[2]])
  s <- which(v==1)*minseglen  #Vecteur d'emplacements des ruptures
  s <- c(1, s, J)
  P <- sapply(1:(k+1), rec_t, s=s, y=y)
  return(-J*log(sum(P)/J)-2*(k+1)*log(J))
  
}

costga_t2 <- function(v,y,minseglen) {
  k <- sum(v)   #Nb de ruptures
  J <- length(y[[2]])
  s <- which(v==1)*minseglen  #Vecteur d'emplacements des points de ruptures
  s <- c(1, s, J)
  P <- sapply(1:(k+1), rec_t, s=s, y=y)
  return(-sum(P)-2*(k+1)*log(J))
  
}



X <- 1:400
Y <- data.frame(X,c(2*X[1:100]+rnorm(100,0,10), 3*X[101:300]+2+rnorm(200,0,10),4*X[301:400]+10+rnorm(100,0,10)))
plot(Y,type="l")


minseglen <- 10
resultY <- ga(type = "binary",fitness =costga_t,Y,minseglen, nBits = length(Y[[2]])%/%minseglen, population = gabin_Population, maxiter = 300, popSize = 150)

sum(resultY@solution)

which(resultY@solution==1)*minseglen

graphe(Y,sum(resultY@solution),which(resultY@solution==1)*minseglen,"trend")


temp <- annual[[2]]

minseglen <- 10

microbenchmark(result_t <- ga(type = "binary",fitness =costga_t,annual,minseglen, nBits = length(annual[[2]])%/%minseglen, population = gabin_Population, maxiter = 300 , popSize = 250 ), times=1)

sum(result_t@solution)
which(result_t@solution==1)*minseglen+1659-1


nrupt <- sum(result_t@solution)
rupt <- which(result_t@solution==1)*minseglen

graphe(annual,nrupt,rupt,"trend")

plot(result@solution)




##Modèle AR : 
#Pour la moyenne : 

######Rupture en moyenne########

Xtilde<-function(y,s){
  k <- length(s) #Nb de ruptures
  J <- length(y)
  m<-estim_mean(y,k,s)
  z<-list()
  z<-append(z,c(rep(m[1],s[1])))
  for (i in 2:k){
    z<-append(z,c(rep(m[i],s[i]-s[i-1]))) }
  z<-append(z,c(rep(m[k],J-s[k])))
  xtilde<-y-unlist(z)
  return (data.frame(t=seq(0:(length(xtilde)-1)),xtilde=xtilde))
}



# estimation phi:
estim_phi<-function (Xtilde){
  l<-dim(Xtilde)[1]
  xt<-Xtilde[-1,2]
  xt1<-Xtilde[-l,2]
  P<-(sum(xt*xt1)/l)
  Q<-(sum(xt1^2)/l)
  return(P/Q)
}

#Fonction de cout non pénalisée :
cost_ar <- function(y,s,t,phi){
  if (s<t) {
    return(sum((y[s:t]-mean(y[s:t])-phi*(y[(s-1):(t-1)]-mean(y[s:t])))^2))
  }
  else {return(0)}   #Triangulaire supérieure
}

#Fonction auxiliaire pour costga_ar et costga_ar2
rec_ar <- function(i,s,y,phi){
  return(cost_ar(y, s[i], s[i+1],phi))
  
}

costga_ar <- function(v,y,minseglen) {
  k <- sum(v)   #Nb de ruptures
  J <- length(y)
  s <- which(v==1)*minseglen  #Vecteur d'emplacements des points de ruptures
  s <- c(2,s,J)
  phi <- estim_phi(Xtilde(y,s))
  P <- sapply(1:(k+1), rec_ar, s=s, y=y,phi=phi)
  return(-(J*(sum(P)/J)+(k+1)*log(J)))  
}

costga_ar2 <- function(v,y,minseglen) {
  k <- sum(v)   #Nb de ruptures
  J <- length(y)
  s <- which(v==1)*minseglen  #Vecteur d'emplacements des points de ruptures
  s <- c(2,s,J)
  phi <- estim_phi(Xtilde(y,s))
  P <- sapply(1:(k+1), rec_ar, s=s, y=y,phi=phi)
  return(-(sum(P)+(k+1)*log(J)))  
}


minseglen <- 5

result_ar <- ga(type = "binary",fitness =costga_ar,temp,minseglen, nBits = length(temp)%/%minseglen, population = gabin_Population, maxiter = 400 , popSize = 250 )



sum(result_ar@solution)
which(result_ar@solution==1)*minseglen+1659-1


nrupt_ar <- sum(result_ar@solution)
rupt_ar <- which(result_ar@solution==1)*minseglen

graphe(annual,nrupt_ar,rupt_ar,"mean")


library(tseries)

epsilon1 <- arima.sim(model = list(ar=c(0.5)), n=100, rand.gen = rnorm, sd=1)
epsilon2 <- arima.sim(model = list(ar=c(0.5)),n=200, rand.gen = rnorm, sd=1)
epsilon3 <- arima.sim(model = list(ar=c(0.5)),n=100, rand.gen = rnorm, sd=1)

XAR <- data.frame(1:400,c(epsilon1,10+epsilon2,2+epsilon3))
plot(XAR, type='l')


minseglen <- 10

result_XAR <- ga(type = "binary",fitness =costga_ar,XAR[[2]],minseglen, nBits = length(XAR[[2]])%/%minseglen, population = gabin_Population, maxiter = 300 , popSize = 250 )



sum(result_XAR@solution)
which(result_XAR@solution==1)*minseglen
graphe(XAR,sum(result_XAR@solution),which(result_XAR@solution==1)*minseglen,"mean")


##    Rupture dans la tendance AR

Xtilde_t <- function(y,s){
  J <- length(y[[2]])
  k <- length(s)  #Nombre de ruptures
  coef <- estim_trend(y,k,s)
  G <- numeric(J)
  if(k>1){
    G[1:s[1]] <- y[[2]][1:s[1]]-coef[1,1]*y[[1]][(1:s[1])]-coef[1,2]
    for (i in 2:k) {
      if(s[i]<J){
        G[(s[i-1]+1):s[i]] <- y[[2]][(s[i-1]+1):s[i]]-coef[i,1]*y[[1]][(s[i-1]+1):s[i]] -coef[i,2]
      }
    }
    if(s[k]<J){
      G[(s[k]+1):J] <- y[[2]][(s[k]+1):J]-coef[(k+1),1]*y[[1]][(s[k]+1):J]-coef[(k+1),2]
    }
  }
  else {
    G[1:s[1]] <- y[[2]][1:s[1]]-coef[1,1]*y[[1]][(1:s[1])]-coef[1,2]
    G[(s[1]+1):J] <- y[[2]][(s[1]+1):J]-coef[2,1]*y[[1]][((s[1]+1):J)]-coef[2,2]
  }
  return(data.frame(y[[1]],G))
}



cost_ar_t <- function(y,r,t,phi){
  if(r<t){
    res <- estim_a_b(y,r,t) 
    a<-res[1]
    b<-res[2]
    M<-y[[2]][r:t]
    I<-y[[1]][r:t]
    N <- y[[2]][(r-1):(t-1)]
    return(sum((M-(a*I+b)-phi*(N-(a*I+b)))^2))}
  else{return(0)}
  
}


rec_ar_t <- function(i,s,y,phi){
  return(cost_ar_t(y, (s[i]+1), s[i+1],phi))
  
}

costga_ar_t <- function(v,y,minseglen) {
  k <- sum(v)
  J <- length(y[[2]])
  s <- which(v==1)*minseglen  #Vecteur d'emplacements des points de ruptures
  s <- c(2,s,J)
  phi <- estim_phi(Xtilde_t(y,s))
  P <- sapply(1:(k+1), rec_ar_t, s=s, y=y,phi=phi)
  return(-(J*log(sum(P)/J)+2*(k+1)*log(J)))  
}


minseglen <- 2

result_ar_t <- ga(type = "binary",fitness =costga_ar_t,annual,minseglen, nBits = length(annual[[2]])%/%minseglen, population = gabin_Population, maxiter = 300 , popSize = 250 )

nrupt_ar_t <- sum(result_ar_t@solution)
rupt_ar_t <- which(result_ar_t@solution==1)*minseglen


graphe(annual,nrupt_ar_t,rupt_ar_t,"trend")


