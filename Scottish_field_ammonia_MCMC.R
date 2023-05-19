library(sp)
library(pscl)
set.seed(345)

#reading in data
setwd("~/JRRSC_paper_code_files")

ammonia_set=read.csv("ammonia_set.csv")

ammonia=ammonia_set$ammonia
ammonia_x=ammonia_set$ammonia_x
ammonia_y=ammonia_set$ammonia_y
ammonia1=log(ammonia)

load("~/cowfid/corner_field_OS.Rda")
corner_field_OS=rbind(corner_field_OS,corner_field_OS[1,])

#transforming
Y=5+log(ammonia)

#scaling
ammonia_x=ammonia_x/10000
ammonia_y=ammonia_y/10000

corner_field_x=corner_field_OS[,1]/10000
corner_field_y=corner_field_OS[,2]/10000

#building a grid
n=length(Y)


minY1=min(ammonia_y)
minX1=min(ammonia_x)

maxY1=max(ammonia_y)
maxX1=max(ammonia_x)

xdist1=(maxX1-minX1)/12
ydist1=(maxY1-minY1)/12

minY=minY1-ydist1
minX=minX1-xdist1

maxX=maxX1+xdist1
maxY=maxY1+ydist1

xleng=40
yleng=40



xlocs=seq(from=minX,to=maxX,length=xleng)
ylocs=seq(from=minY, to=maxY,length=yleng)
square_grid_points=expand.grid(xlocs,ylocs)


in_field=point.in.polygon(square_grid_points[,1],square_grid_points[,2],corner_field_x,corner_field_y)
field_grid=square_grid_points[which(in_field==1),]


n=length(Y)
N=dim(field_grid)[1]

allpoints=field_grid
irreg_locs=cbind(ammonia_x,ammonia_y)
tot_locs=rbind(as.matrix(allpoints),as.matrix(irreg_locs))
alldistmat=as.matrix(dist(tot_locs))
distmat_RR=alldistmat[1:N,1:N]
distmat_IR=alldistmat[1:N,(N+1):(N+n)]
distmat_II=alldistmat[(N+1):(N+n),(N+1):(N+n)]

#Utility
# alpha=2
# beta=50
# delta=3

maxdist=max(distmat_RR)

#defining a utility function
log_utility=function(Z,D,alpha,beta,M){

  
  desdist= M[D,D]
  diag(desdist)=10
  minmean= mean(apply(desdist,1,min))
  maxmean=mean(apply(M[,D],1,min))
  spread=-maxmean+minmean
  
  alpha*sum(Z[D])+ 20000*beta*spread
}


#iterations etc.
num_iter=25000
Kreps=1500
burnin=500
Design_storeZ=matrix(0, nrow=n,ncol=Kreps)
Design_storeZ[,1]=sample(1:N,n)
Prop_Design_storeZ=matrix(0,nrow=n,ncol=Kreps)
Prop_Design_storeZ[,1]=sample(1:N,n)
design_accept=0

####Actual prediction

#storage and starting values for MCMC
ssq_vec=rep(0,num_iter) #sigma
ssq_vec[1]=10
logphi_vec=rep(0,num_iter)#phi
logphi_vec[1]=-1.6
mp_vec=rep(0,num_iter)
mp_vec[1]=1.5
Cormat=exp(-alldistmat/exp(logphi_vec[1]))
cor_RR=exp(-distmat_RR/exp(logphi_vec[1]))
cor_II=exp(-distmat_II/exp(logphi_vec[1]))
cor_IR=exp(-distmat_IR/exp(logphi_vec[1]))
Z_irreg_store=matrix(0,n,num_iter) #Z_I
Z_irreg_store[,1]=Y
Z_reg_store=matrix(0,N,num_iter)#Z_R
Z_reg_store[,1]=mean(Y)+as.vector( (cor_IR)%*%solve(cor_II)%*%(Y-mean(Y)))

alpha_vec=rep(0,num_iter)
alphamean=0
alphavar=20
alpha_step=1.2
betamean=0
betavar=20
beta_step=1
beta_vec=rep(0,num_iter)
beta_vec[1]=5
tsq_vec=rep(0,num_iter)
tsq_vec[1]=0.05

#Hyperparams


a=4
b=0.1

f=1
g=1


h=2.25
k=0.0625

phi_var=0.2

#
gradlogU=function(Z,d,alpha,DesignSample){
  
  nreps=dim(DesignSample)[2]

  chosen_actual=tabulate(d,nbins=N)
  chosen_current=tabulate(DesignSample,nbins=N)/nreps
  
  alpha*(chosen_actual-chosen_current)}

#starting correlation matrices
cor_IR=exp(-distmat_IR/exp(logphi_vec[1]))
cor_II=exp(-distmat_II/exp(logphi_vec[1]))
cor_RR=exp(-distmat_RR/exp(logphi_vec[1]))
ch_II=chol(cor_II)
inv_II=chol2inv(ch_II)
ch_RR=chol(cor_RR)
inv_RR=chol2inv(ch_RR)
R_bar=cor_RR-(cor_IR)%*%inv_II%*%t(cor_IR)
ch_R_bar=chol(R_bar)
#A_RR is the inverse of R bar
A_RR=chol2inv(ch_R_bar)
A_IR=-A_RR%*%cor_IR%*%inv_II
A_II=inv_II+inv_II%*%t(cor_IR)%*%A_RR%*%cor_IR%*%inv_II
tot_cor=exp(-alldistmat/exp(logphi_vec[1]))
ch_tot_cor=chol(tot_cor)
inv_tot_cor=chol2inv(ch_tot_cor)
inv_tot_cor=0.5*(inv_tot_cor+t(inv_tot_cor))
inv_RR=chol2inv(chol(cor_RR))
inv_II=chol2inv(chol(cor_II))
phi_accept=0
alpha_vec[1]=0

#a starting design sample
Kreps=1500
design_accept=0
Design_storeZ=matrix(0,n,Kreps)
Design_storeZ[,1]=sample(1:N,n,replace=TRUE)
utilvec=rep(0,Kreps)
utilvec[1]=log_utility(Z_reg_store[,1],Design_storeZ[,1],alpha_vec[1],beta_vec[1],distmat_RR)
for(l in 2:Kreps){
  change_sites=sample(1:n,1)
  newsite=sample(1:N,1)
  changedsite=(Design_storeZ[,l-1])[change_sites]
  Design_prop=c(Design_storeZ[-change_sites,l-1],newsite)
  
 
  proputil=log_utility(Z_reg_store[,1],Design_prop,alpha_vec[1],beta_vec[1],distmat_RR)
  MH_rat1=proputil-utilvec[l-1]
  if(log(runif(1))<MH_rat1){Design_storeZ[,l]=Design_prop
  utilvec[l]=proputil
  design_accept=design_accept+1
  }else{
    utilvec[l]=utilvec[l-1]
    Design_storeZ[,l]=Design_storeZ[,l-1]}
}

mp_vec=rep(0, num_iter)
mp_vec[1]=1.5
mpmean=2
mp_accept=0
mpvar=1

ZRaccept=0

design_vec=rep(n,0)
for(i in 1:n){
  cur_mon=irreg_locs[i,]
  distvec=rep(0,N)
  for(j in 1:N){distvec[j]= (cur_mon[1]-allpoints[j,1])^2 +(cur_mon[2]-allpoints[j,2])^2}
  design_vec[i]=which(distvec==min(distvec))
}

desdist= distmat_RR[design_vec,design_vec]
diag(desdist)=100000000 #trick so this will never be min val 
minmean= mean(apply(desdist,1,min))
maxmean=mean(apply(distmat_RR[,design_vec],1,min))
realspread=-maxmean+minmean


#MCMC sampling loop
for(i in 2:num_iter){
  
  #ZI
  Sig1= (1/ssq_vec[i-1])*( cor_II-t(cor_IR)%*% inv_RR%*%(cor_IR)) # AII
  Sig2=diag(1/tsq_vec[i-1],n)
  sig3=chol2inv(chol(Sig1+Sig2))
  Mu1=mp_vec[i-1]+t(cor_IR)%*%inv_RR%*%(Z_reg_store[,i-1]-mp_vec[i-1])
  Mu2=Y
  Mu3=sig3%*%Sig1%*%Mu1+sig3%*%Sig2%*%Mu2
  
  
  Z_irreg_store[,i]=Mu3+chol(sig3)%*%(rnorm(n,0,1))
  
  
  #ZR
  
  mu_RR=mp_vec[i-1]+cor_IR%*%solve(cor_II)%*%(Z_irreg_store[,i]-mp_vec[i-1])
  cholcor_II=chol(cor_II)
  sig_RR=cor_RR-cor_IR%*%chol2inv(cholcor_II)%*%t(cor_IR) #rbar
  sig_RR=0.5*(sig_RR+t(sig_RR))
  cholsigRR=chol(sig_RR)
  A_RR=chol2inv(cholsigRR)
  
  dir=gradlogU(Z_reg_store[,i-1],design_vec,alpha_vec[i-1],Design_storeZ)- (1/ssq_vec[i-1])*A_RR%*%((Z_reg_store[,i-1])-(mu_RR))
  
  propnorm=rnorm(N,0,1)
  tp=0.1
  propZ=Z_reg_store[,i-1]+tp*propnorm+(0.5*tp^2)*dir
  
  
  
  
  
  
  #generating designs, first with propZ
  Design_storeZP=Design_storeZ
  
  #generating designs, first with propZ
  utilvecp=rep(0,Kreps)
  Design_storeZP=Design_storeZ
  
  Design_storeZP[,1]=Design_storeZ[,Kreps]
  utilvecp[1]=log_utility(propZ,Design_storeZP[,1],alpha_vec[i-1],beta_vec[i-1],distmat_RR)
  design_accept=0
  for(l in 2:Kreps){
    
    change_sites=sample(1:n,1)
    newsite=sample(1:N,1)
    changedsite=(Design_storeZP[,l-1])[change_sites]
    Design_prop=c(Design_storeZP[-change_sites,l-1],newsite)
    #
    proputil=log_utility(propZ,Design_prop,alpha_vec[i-1],beta_vec[i-1],distmat_RR)
    MH_rat1z=proputil-utilvecp[l-1]
    if(log(runif(1))<MH_rat1z){Design_storeZP[,l]=Design_prop
    
    utilvecp[l]=proputil
    design_accept=design_accept+1}else{
      Design_storeZP[,l]=Design_storeZP[,l-1]
      utilvecp[l]=utilvecp[l-1]
    } }
  
  
  
  
  
  #accept/reject
  
  
  utilrat=log_utility(propZ,design_vec,alpha_vec[i-1],beta_vec[i-1],distmat_RR)-log_utility(Z_reg_store[,i-1],design_vec,alpha_vec[i-1],beta_vec[i-1],distmat_RR)
  
  #Krat
  

  KR_1=0
  
  KB=Kreps-burnin
  KR1_vec=rep(0,KB)
  PZZ=alpha_vec[i-1]*(Z_reg_store[,i-1]-propZ)
  for(t in 1:KB){
    CD=Design_storeZP[,t+burnin]
    KR1_vec[t]=sum(PZZ[CD])}
  
  logKR1=-log(KB)+max(KR1_vec)+log(sum(exp(KR1_vec-max(KR1_vec))))
  
  
  
  #now the transition probabilities and the MVN part of the acceptance
  
  triangle_y=gradlogU(propZ,design_vec,alpha_vec[i-1], Design_storeZP)- (1/ssq_vec[i-1])*A_RR%*%((propZ)-mu_RR)
  
  log_qyx=-0.5*N*log(2*pi*tp^2)-(1/(2*tp^2))*sum(( Z_reg_store[,i-1]- propZ -(0.5*tp^2)*triangle_y   )^2)  
  
  log_qxy=-0.5*N*log(2*pi*tp^2)-(1/(2*tp^2))*sum(( -Z_reg_store[,i-1]+ propZ -(0.5*tp^2)*dir   )^2)  
  
  ## The MVN part of the likelihood ratio
  logmvnLR_top=-(0.5/ssq_vec[i-1])*t( propZ-mu_RR   )%*%A_RR%*%(propZ-mu_RR)
  logmvnLR_bottom=-(0.5/ssq_vec[i-1])*t( Z_reg_store[,i-1]-mu_RR   )%*%A_RR%*%(Z_reg_store[,i-1]-mu_RR)
  logmvnLR=logmvnLR_top-logmvnLR_bottom
  #putting them all together
  accept=log_qyx-log_qxy+logmvnLR+ utilrat+logKR1
  
  
  
  if( log(runif(1))<  accept  ){
    Z_reg_store[,i]=propZ
    
    Design_storeZ=Design_storeZP
    
    ZRaccept=ZRaccept+1
    utilvec=utilvecp
    
  }else{Z_reg_store[,i]=Z_reg_store[,i-1]}
  
  
  #tsq
  new_h=h+0.5*n
  new_k=k+0.5*crossprod((Z_irreg_store[,i]-Y),(Z_irreg_store[,i]-Y))
  tsq_vec[i]=(1/rgamma(1,new_h,new_k))
  
  
  
  
  #alpha
  propA=alpha_vec[i-1] +rnorm(1,0,alpha_step)
  anewaold=propA-alpha_vec[i-1]
  curz=Z_reg_store[,i]
  utilrat=anewaold*sum(curz[design_vec])
  a_pri=-(1/(2*alphavar))*( (propA-alphamean)^2-(alpha_vec[i-1]-alphamean)^2 )
  # #krat
  
  KR1=0
  KR1_vec=rep(0,KB)
  
  for(t in 1:KB){
    CD=Design_storeZ[,t+burnin]
    KR1_vec[t]=anewaold*sum(curz[CD])}
  
  logKR1=-log(KB)+max(KR1_vec)+log(sum(exp((KR1_vec-max(KR1_vec))   )))
  

  #
  alpha_accept=utilrat+a_pri-logKR1
  #
  if(log(runif(1))<alpha_accept){
    alpha_vec[i]=propA
    #we need a new design for dir
    Design_storeZ[,1]=Design_storeZ[,Kreps]
    utilvec[1]=log_utility(Z_reg_store[,i],Design_storeZ[,1],alpha_vec[i],beta_vec[i-1],distmat_RR)
    for(l in 2:Kreps){
      change_sites=sample(1:n,1)
      newsite=sample(1:N,1)
      changedsite=(Design_storeZ[,l-1])[change_sites]
      Design_prop=c(Design_storeZ[-change_sites,l-1],newsite)
      
      proputil=log_utility(Z_reg_store[,i],Design_prop,alpha_vec[i],beta_vec[i-1],distmat_RR)
      MH_rat1a=proputil-utilvec[l-1]
      if(log(runif(1))<MH_rat1a){Design_storeZ[,l]=Design_prop
      utilvec[l]=proputil}else{
        Design_storeZ[,l]=Design_storeZ[,l-1]
        utilvec[l]=utilvec[l-1]}
    }
    #  
    #  
  }else{alpha_vec[i]=alpha_vec[i-1]}
  
  

  # #beta
  propB=beta_vec[i-1] +rnorm(1,0,beta_step)
  bnewbold=propB-beta_vec[i-1]
  utilrat=20000*(bnewbold)*realspread
  
  b_pri=-(1/(2*betavar))*( (propB-betamean)^2-(beta_vec[i-1]-betamean)^2 )
  
  
  KR1=0
  KR1_vec=rep(0,KB)
  for(t in 1:KB){
    CD=Design_storeZ[,t+burnin]
    desdist= distmat_RR[CD,CD]
    diag(desdist)=10000000000000
    minmean= mean(apply(desdist,1,min))
    maxmean=mean(apply(distmat_RR[,CD],1,min))
    CDspread=-maxmean+minmean
    #KR1_vec[t]=log_utility(Z_reg_store[,i],CD,alpha_vec[i],propB,distmat_RR)- utilvec[t]
    KR1_vec[t]=20000*bnewbold*CDspread
  }
  
  
  logKR1=-log(KB)+max(KR1_vec)+log(sum(exp((KR1_vec-max(KR1_vec)))) )
  
  
  
  
  # KR2=KR2+exp(anewaold*sum(sortZ*tabulate(Design_storeZS[,t],nbins=N))) }
  #
  # Krat=(1/(Kreps))*(KR1+KR2)
  Krat=KR1/Kreps
  #
  beta_accept=utilrat+b_pri-logKR1
  #
  if(log(runif(1))<beta_accept){
    beta_vec[i]=propB
    #we need a new design for dir
    Design_storeZ[,1]=Design_storeZ[,Kreps]
    utilvec[1]=log_utility(Z_reg_store[,i],Design_storeZ[,1],alpha_vec[i],beta_vec[i],distmat_RR)
    
    for(l in 2:Kreps){
      change_sites=sample(1:n,1)
      newsite=sample(1:N,1)
      changedsite=(Design_storeZ[,l-1])[change_sites]
      Design_prop=c(Design_storeZ[-change_sites,l-1],newsite)
      
      proputil=log_utility(Z_reg_store[,i],Design_prop,alpha_vec[i],beta_vec[i],distmat_RR)
      MH_rat1b=proputil-utilvec[l-1]
      if(log(runif(1))<MH_rat1b){Design_storeZ[,l]=Design_prop
      utilvec[l]=proputil}else{
        Design_storeZ[,l]=Design_storeZ[,l-1]
        utilvec[l]=utilvec[l-1]}
    }
    #
    #
  }else{beta_vec[i]=beta_vec[i-1]}
  
  #beta_vec[i]=0
  
  
  
  
  #ssq

  
  
  R_C=f+0.5*(N)
  bigZR=c(Z_reg_store[,i])
  R_D=g+0.5*t(bigZR-mp_vec[i-1])%*%inv_RR%*%(bigZR-mp_vec[i-1])
  

  ssq_vec[i]=1/rgamma(1,R_C,R_D)
  
  prop_logphi=logphi_vec[i-1]+rnorm(1,0,sqrt(phi_var))
  
  
  #what these would make the correlation matrices


    prop_cor_RR=exp(-distmat_RR/exp(prop_logphi))
  prop_ch_RR=chol(prop_cor_RR)
  prop_inv_RR=chol2inv(prop_ch_RR)
  
  Zr=Z_reg_store[,i]

  oldloglike=-sum(log(diag(ch_RR))) - (0.5/ssq_vec[i])*( t(Zr-mp_vec[i-1])%*%inv_RR%*%(Zr-mp_vec[i-1])  )- b*exp(-logphi_vec[i-1])-a*logphi_vec[i-1]
  newloglike=-sum(log(diag(prop_ch_RR)))-(0.5/ssq_vec[i])*(t(Zr-mp_vec[i-1])%*%prop_inv_RR%*%(Zr-mp_vec[i-1]))-b*exp(-prop_logphi)-a*prop_logphi


  
  MH_rat2=newloglike-oldloglike
  
  if(log(runif(1))<MH_rat2){ 
    logphi_vec[i]=prop_logphi
  cor_IR=exp(-distmat_IR/exp(prop_logphi))
  cor_II=exp(-distmat_II/exp(prop_logphi))
  cor_RR=exp(-distmat_RR/exp(prop_logphi))
  
  ch_II=chol(cor_II)
  inv_II=chol2inv(ch_II)
  ch_RR=chol(cor_RR)
  inv_RR=chol2inv(ch_RR)
  R_bar=cor_RR-(cor_IR)%*%inv_II%*%t(cor_IR)
  R_bar=(R_bar+t(R_bar))*.5
  ch_R_bar=chol(R_bar)
  #A_RR is the inverse of R bar
  A_RR=chol2inv(ch_R_bar)
  A_IR=-A_RR%*%(cor_IR)%*%inv_II
  A_II=inv_II+inv_II%*%t(cor_IR)%*%A_RR%*%(cor_IR)%*%inv_II
  

  
  Sig1=(1/ssq_vec[i])*(cor_II - t(cor_IR)%*%inv_RR%*%cor_IR)
  
  phi_accept=phi_accept+1
  
  }else{logphi_vec[i]=logphi_vec[i-1]}
  
  
  
  
  #predicting MP
  propMP=mp_vec[i-1]+rnorm(1,0,1.5)
  ZMB=Z_reg_store[,i]-propMP
  loglikepropMP=-0.5*(t(ZMB)%*%inv_RR%*%ZMB) - (1/(2*mpvar))*(propMP-mpmean)^2
  
  
  #repeating for current value
  ZMB=Z_reg_store[,i]-mp_vec[i-1]
  loglikeMP=-0.5*(t(ZMB)%*%inv_RR%*%ZMB)  - (1/(2*mpvar))*(mp_vec[i-1]-mpmean)^2
  
  MP_MH=loglikepropMP-loglikeMP
  
  if(log(runif(1))< MP_MH ){mp_vec[i]=propMP
  mp_accept=mp_accept+1}else{
    mp_vec[i]=mp_vec[i-1]
  } 
  
  
  print(i) #counter

  
}

#plotting the results
PredZ=rowMeans(Z_reg_store[,1:i])
predZR_range=range(PredZ)

emptymat=matrix(0,nrow=yleng,ncol=xleng)
N2=xleng*yleng
nearest_pt=rep(0,N)
for(k in 1:N){
  current_pt=allpoints[k,]
  distvec=rep(0,N2)
  for(j in 1:N2){distvec[j]=(current_pt[1]-square_grid_points[j,1])^2+(current_pt[2]-square_grid_points[j,2])^2
  }
  distvec=as.numeric(distvec)
  nearest_pt[k]=which(distvec==min(distvec))
  print(k)
}

squareZ=rep(0,N2)
for(k in 1:N){squareZ[nearest_pt[k]]=PredZ[k]}
colr=heat.colors(128)
xlocs=seq(from=minX,to=maxX,length=xleng)
ylocs=seq(from=minY, to=maxY,length=yleng)
zrange=range(PredZ)
layout(matrix(1:2,1), widths = c(9,2))
image(xlocs,ylocs,matrix(squareZ,xleng,yleng),col=colr,zlim=zrange, xlab="Longitude", ylab="Latitude")
points(ammonia_x,ammonia_y)
image(c(0,1), seq(predZR_range[1],predZR_range[2],length=length(colr)), matrix(seq_along(colr), 1),
      axes = FALSE, xlab = "", ylab = "", col=colr)
box()
axis(4)

par(mfrow=c(1,1))





