library(readxl)
library(rlist)
library(sp)
library(tmap)
library(rgeos)
library(rgdal)
library(dplyr) 
library(Gifi)
library(fitdistrplus)
library(raster)
library(ggplot2)

pesos1<-as.data.frame(cbind(sectorvalQuito$SectorN,sectorvalQuito$Area/sectorvalQuito$NumeroViviendas))
V3<-sectorvalQuito$NumeroViviendas/sectorvalQuito$NumeroHogares
pesos1<-as.data.frame(cbind(pesos1,V3)) ###pesos

#pesos<-as.data.frame(cbind(sectorvalQuito$SectorN,sectorvalQuito$NumeroViviendas))
pesos1<-pesos1[order(pesos1$V1), ]
#pesos1<-pesos1$V2
pesos<-pesos2*pesos1$V3

dendrite <- function(D)##funcion del arbol de minima expansion
{
  library("ape")
  D <- as.matrix(D)
  M <- mst(D)
  lw=nrow(M)
  den <- matrix(0,lw-1,3)
  k <- 0
  for (i in 1:(lw-1))
    for (j in (i+1):lw)
    {
      if (M[i,j]==1)
      {
        k <- k+1 
        den[k,1] <- i
        den[k,2] <- j
        den[k,3] <- D[i,j]
      }
    }
  par <- matrix(0,5,1)
  par[1,1] <- mean(den[,3])#media de los valores de i a j
  par[2,1] <- sd(den[,3])#desvest
  par[3,1] <- par[1,1]+par[2,1]#media+desvest
  par[4,1] <- par[1,1]+2*par[2,1]#media+2desvest
  par[5,1] <- par[1,1]+3*par[2,1]#media+3desvest
  wynik <- list(den=den,par=par)
  return(wynik)
}
library(fda)
D <- base2 #lista de matrices


p <- ncol(D[[1]])#numero de variables
T<-nrow(D[[1]]) #numero de observaciones por region
N <- length(D) ## numero re regiones

F <- matrix(0,nrow=N,ncol=N) ## calculo de la norma de frobenius
for (i in 1:N)
  for (j in 1:N)
  {
    m <- t(as.matrix(D[[i]]-D[[j]]))%*%(as.matrix(D[[i]]-D[[j]]))
    F[i,j] <- sqrt(sum(diag(m)))
  }

K1 <- matrix(0,nrow=p,ncol=p)
K1 <- exp(-F^2/(sum(F^2)/(N*(N-1)))) ##calculo del kernel de gauss

I <- diag(N)
J <- matrix(1,nrow=N,ncol=N)
J <- (1/N)*J
P <- I-J
K1 <- P%*%K1%*%P ##calculo del kernel centrado

u1 <- eigen(K1) ##valores y vectores propios
uK1 <- u1$values
Abi<-5 ##numero de componentes para tener el 70% de representatividad
m1 <- matrix(0,nrow=N,ncol=Abi) ##aqui se almacenan las proyecciones

for (k in 1:Abi) 
{
  V<- u1$vectors[,k]
  for (i in 1:N)
  {
    suma <- 0
    for (j in 1:N)
    {
      
      suma <- suma+V[j]*K1[i,j]
    }
    m1[i,k] <- suma
  }
}#calculo de las proyecciones

A<-NULL
for(i in 1:191){A[i]<-sum(m1[i,])} #suma de las componentes para el indicador
A1<-A-min(A)
A1<-as.vector(A1/max(A1))#escalar el indicador a la escala 0-1
A1<-1-A1 ## calculo del indicador


D <- dist(m1)
wynik <- dendrite(D)
te <- uK1/sum(uK1)*100
te <- round(te,2)
plot(m1[,1],m1[,2],type='n',xlab=bquote(paste(U[1]," (")~.(te[1])*"% )"),ylab=bquote(paste(U[2]," (")~.(te[2])*"% )"),pch=19,col="blue",main="")
dendr <- wynik
ded <- dendr$den
dep <- dendr$par
w <- nrow(ded)
lev <- 3
for (i in 1:w)
{
  if (ded[i,3]<dep[lev,1])
  {
    punktx <- c(m1[ded[i,1],1],m1[ded[i,2],1])
    punkty <- c(m1[ded[i,1],2],m1[ded[i,2],2])
    segments(m1[ded[i,1],1],m1[ded[i,1],2],m1[ded[i,2],1],m1[ded[i,2],2])
    points(punktx, punkty, pch = 21, col = "black", bg = "white", cex = 3)
  }
  else
  {
    punktx <- c(m1[ded[i,1],1],m1[ded[i,2],1])
    punkty <- c(m1[ded[i,1],2],m1[ded[i,2],2])
    segments(m1[ded[i,1],1],m1[ded[i,1],2],m1[ded[i,2],1],m1[ded[i,2],2],lty=2)
    points(punktx, punkty, pch = 21, col = "black", bg = "white", cex = 3)
  }
}
text(m1[,1],m1[,2],cex=0.8) #grafico del arbol de minima expansion


####CASO CON PESOS GEOGRAFICOS


Abi<-5
DK=as.numeric(pesos) #VECTOR DE PESOS
SDK <- sum(DK)
#DK<-(DK-mean(DK))/sd(DK)
DK <- DK/SDK #ESCALAR PESOS
#DK<-1-DK

D <- base2
p <- ncol(D[[1]])#numero de variables
T<-nrow(D[[1]]);
N <- length(D)
X<-D

for (i in 1:N)
{
  X[[i]] <- X[[i]]*DK[i]
}


F <- matrix(0,nrow=N,ncol=N)
for (i in 1:N)
  for (j in 1:N)
  {
    m <- t(as.matrix(X[[i]]-X[[j]]))%*%(as.matrix(X[[i]]-X[[j]]))
    F[i,j] <- sqrt(sum(diag(m)))
  }

K1 <- matrix(0,nrow=N,ncol=N)
K1 <- exp(-F^2/(sum(F^2)/(N*(N-1))))

I <- diag(N)
J <- matrix(1,nrow=N,ncol=N)
J <- (1/N)*J
P <- I-J
K1 <- P%*%K1%*%P

u1 <- eigen(K1)
uK1 <- u1$values

m <- matrix(0,nrow=N,ncol=Abi)

for (k in 1:Abi) 
{
  V<- u1$vectors[,k]
  for (i in 1:N)
  {
    suma <- 0
    for (j in 1:N)
    {
      
      suma <- suma+V[j]*K1[i,j]
    }
    m[i,k] <- suma
  }
}


A2<-NULL
Aux<-NULL
for(i in 1:length(secvalQ)){A2[i]<-sum(m[i,])}
A3<-A2-min(A2)
A3<-as.vector(A3/max(A3))
A3<-1-A3

D <- dist(m)
wynik <- dendrite(D)
te <- uK1/sum(uK1)*100
te <- round(te,2)
plot(m[,1],m[,2],type='n',xlab=bquote(paste(U[1]," (")~.(te[1])*"% )"),ylab=bquote(paste(U[2]," (")~.(te[2])*"% )"),pch=19,col="blue",main="")
dendr <- wynik
ded <- dendr$den
dep <- dendr$par
w <- nrow(ded)
lev <- 3
puntos1<-list()
puntos2<-list()
for (i in 1:w)
{
  if (ded[i,3]<dep[lev,1])
  {
    punktx <- c(m[ded[i,1],1],m[ded[i,2],1])##GRAFICA LOS DOS PUNTOS AL MISMO TIEMPO
    punkty <- c(m[ded[i,1],2],m[ded[i,2],2])
    segments(m[ded[i,1],1],m[ded[i,1],2],m[ded[i,2],1],m[ded[i,2],2])
    points(punktx, punkty, pch = 21, col = "black", bg = "white", cex = 3)
    
  }
  else
  {       ##proyeccion punto 1 proyeccion punto 2
    punktx <- c(m[ded[i,1],1],m[ded[i,2],1])
    punkty <- c(m[ded[i,1],2],m[ded[i,2],2])
    segments(m[ded[i,1],1],m[ded[i,1],2],m[ded[i,2],1],m[ded[i,2],2],lty=2)
    points(punktx, punkty, pch = 21, col = "black", bg = "white", cex = 3)
  }
}
text(m[,1],m[,2],cex=0.8)


sectorIndex2<-sectorvalQuito
sectorIndex2<-sectorIndex2[order(sectorIndex2$SectorN), ]

sectorIndex2$Caton<-A3

save(sectorIndex2, file = "SectorW.Rdata")

## Adjunto SectorW. 





