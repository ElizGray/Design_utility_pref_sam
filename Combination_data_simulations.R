#trying it lots of times with the M one  THIS WORKS FOR five iters

set.seed(234)
library(pscl)
repeats=50
n=30
sse_vec_full=rep(0,repeats)
sse_vec_mult=rep(0,repeats)
alpha_vec_full=rep(0,repeats)
alpha_vec_mult=rep(0,repeats)
beta_vec_full=rep(0,repeats)

sigma_vec_full=rep(0,repeats)
sigma_vec_mult=rep(0,repeats)

logphi_vec_full=rep(0,repeats)
logphi_vec_mult=rep(0,repeats)


tau_vec_full=rep(0,repeats)
tau_vec_mult=rep(0,repeats)

mp_vec_full=rep(0,repeats)
mp_vec_mult=rep(0,repeats)


xleng=15
yleng=15

N=xleng*yleng

predZ_store_full=matrix(0, N, repeats)

predZ_store_mult=matrix(0, N, repeats)

Z_store=matrix(0, N, repeats)


num_iter=1000
alpha_store_full=matrix(0, num_iter,repeats)
beta_store_full=matrix(0, num_iter, repeats)
alpha_store_mult=matrix(0, num_iter,repeats)
beta_store_mult=matrix(0, num_iter, repeats)

sigma_store_full=matrix(0, num_iter,repeats)
sigma_store_mult=matrix(0, num_iter,repeats)

logphi_store_full=matrix(0, num_iter,repeats)
logphi_store_mult=matrix(0, num_iter,repeats)

tau_store_full=matrix(0, num_iter,repeats)
tau_store_mult=matrix(0, num_iter,repeats)

mp_store_full=matrix(0, num_iter,repeats)
mp_store_mult=matrix(0, num_iter,repeats)




Z_array_full=array(0,dim=c(N,num_iter,repeats))
Z_array_mult=array(0, dim=c(N,num_iter,repeats))

ZI_array_full=array(0,dim=c(n,num_iter,repeats))
ZI_array_mult=array(0, dim=c(n,num_iter,repeats))

alpha_sd_full=rep(0,repeats)
beta_sd_full=rep(0,repeats)
sigma_sd_full=rep(0,repeats)
logphi_sd_full=rep(0, repeats)
tau_sd_full=rep(0,repeats)
mp_sd_full=rep(0,repeats)

alpha_sd_mult=rep(0,repeats)
beta_sd_mult=rep(0,repeats)
sigma_sd_mult=rep(0,repeats)
logphi_sd_mult=rep(0, repeats)
tau_sd_mult=rep(0,repeats)
mp_sd_mult=rep(0,repeats)




for(p in 1:repeats){
  
  xleng=15
  yleng=15
  
  N=xleng*yleng
  
  logphi=1
  ssq=4
  
  alpha=4.5
  beta=2.75
  
  xlocs=seq(from=0, to=1, length=xleng)
  ylocs=seq(from=0, to=1, length=yleng)
  
  allpoints=expand.grid(xlocs,ylocs)
  
  distmat=as.matrix(dist(allpoints))
  
  Cormat=exp(-distmat/exp(logphi))
  
  big_sig=ssq*Cormat
  
  chol_sig=chol(big_sig)
  mp=5
  
  Z=mp+rnorm(N,0,1)%*%chol_sig
  
  cols=heat.colors(128)
  image(xlocs,ylocs,matrix(Z,xleng,yleng), col=cols)
  
  maxdist=max(distmat)
  log_utility=function(Z,D,alpha,beta,M){
   
    tabd=tabulate(D,nbins=N)
    facbit=sum(log(factorial(tabd)))

 
    
    desdist= M[D,D]
    diag(desdist)=10
    minmean= mean(apply(desdist,1,min))
    maxmean=mean(apply(M[,D],1,min))
    spread=-maxmean+minmean
    
    -facbit+alpha*(sum(Z[D]))+ 500*beta*(spread)
    
  }
  
  
  log_utility2=function(Z,D,alpha,beta,M){

    
    desdist= M[D,D]
    diag(desdist)=10
    minmean= mean(apply(desdist,1,min))
    maxmean=mean(apply(M[,D],1,min))
    spread=-maxmean+minmean
    alpha*(sum(Z[D]))+500*beta*spread}
  
  
  
  
  
  
  
  
  n=30
  Kreps=5000
  design_accept=0
  Design_storeZ=matrix(0, nrow=n,ncol=Kreps)
  
  Design_storeZ[,1]=sample(1:N,n,replace=TRUE)
  for(l in 2:Kreps){
    change_sites=sample(1:n,1)
    newsite=sample(1:N,1)
    changedsite=(Design_storeZ[,l-1])[change_sites]
    Design_prop=c(Design_storeZ[-change_sites,l-1],newsite)
    
    dup_in_new=length(which(Design_prop==newsite))
    dup_in_old=length(which(Design_storeZ[,l-1]==changedsite))
    scalefac=dup_in_new/dup_in_old
    
    MH_rat1=log(scalefac)+log_utility(Z,Design_prop,alpha,beta,distmat)-log_utility(Z,Design_storeZ[,l-1],alpha,beta,distmat)
    if(log(runif(1))<MH_rat1){Design_storeZ[,l]=Design_prop
    design_accept=design_accept+1
    }else{
      Design_storeZ[,l]=Design_storeZ[,l-1]}
    
   # if(l%%100==0){print(l/Kreps)}
  }
  
  design_vec=Design_storeZ[,Kreps]
  
  

  

  
  
  
  corner_points_sel=cbind(allpoints[design_vec,1]-0.5/xleng,allpoints[design_vec,2]-0.5/xleng)
  irreg_locs=as.matrix(corner_points_sel+cbind(runif(n,0,0.5/xleng),runif(n,0,0.5/xleng)))
  
  points(irreg_locs)
  
  allpoints=as.matrix(allpoints)
  distmat_II=as.matrix(dist(irreg_locs))
  totdist=as.matrix(dist(rbind(allpoints,irreg_locs)))
  distmat_IR=totdist[(N+1):(n+N),1:N]
  
  cor_RR=exp(-distmat/exp(logphi))
  cor_II=exp(-distmat_II/exp(logphi))
  cor_IR=exp(-distmat_IR/exp(logphi))
  
  cor_RR=0.5*(cor_RR+t(cor_RR))
  ch_RR=chol(cor_RR)
  inv_RR=chol2inv(ch_RR)
  ZI=mp+as.vector( (cor_IR)%*%inv_RR%*%t(Z-mp))
  tsq=0.01
  Y=ZI+rnorm(n,0,tsq)
  
  #so now we have data of Y and irreg. locs
  
  
  ##############################################
  #FITTING FROM HERE
  
  
  
  n=length(Y)
  
  
  
  
  minY=0
  minX=0
  
  maxX=1
  maxY=1

  xlocs=seq(from=minX,to=maxX,length=xleng)
  ylocs=seq(from=minY, to=maxY,length=yleng)
  square_grid_points=expand.grid(xlocs,ylocs)
  
  
  
  
  n=length(Y)
N=xleng*yleng
  
  allpoints=square_grid_points
  
  tot_locs=rbind(as.matrix(allpoints),as.matrix(irreg_locs))
  alldistmat=as.matrix(dist(tot_locs))
  distmat_RR=alldistmat[1:N,1:N]
  distmat_IR=alldistmat[1:N,(N+1):(N+n)]
  distmat_II=alldistmat[(N+1):(N+n),(N+1):(N+n)]
  
  #Utility
  
  
  log_utility=function(Z,D,alpha,beta,M){
  
    
    tabd=tabulate(D,nbins=N)
    facbit=sum(log(factorial(tabd)))
    
    desdist= M[D,D]
    diag(desdist)=10
   minmean= mean(apply(desdist,1,min))
   maxmean=mean(apply(M[,D],1,min))
    spread=-maxmean+minmean

    -facbit+ alpha*(sum(Z[D]))+500*beta*spread
  }
  
  
  num_iter=1000
  Kreps=1500
  Design_storeZ=matrix(0, nrow=n,ncol=Kreps)
  
  Design_storeZ[,1]=sample(1:N,n)
  Prop_Design_storeZ=matrix(0,nrow=n,ncol=Kreps)
  Prop_Design_storeZ[,1]=sample(1:N,n)
  design_accept=0
  
  ####Actual prediction
  ssq_vec=rep(0,num_iter) #sigma
  ssq_vec[1]=1
  logphi_vec=rep(0,num_iter)#phi
  logphi_vec[1]=5
  mp_vec=rep(0,num_iter)
  mp_vec[1]=mean(Y)
  Cormat=exp(-alldistmat/exp(logphi_vec[1]))
  cor_RR=exp(-distmat_RR/exp(logphi_vec[1]))
  cor_II=exp(-distmat_II/exp(logphi_vec[1]))
  cor_IR=exp(-distmat_IR/exp(logphi_vec[1]))
  Z_irreg_store=matrix(0,n,num_iter) #Z_I
  Z_irreg_store[,1]=Y
  Z_reg_store=matrix(0,N,num_iter)#Z_R
  Z_reg_store[,1]=mean(Y)+as.vector( (cor_IR)%*%solve(cor_II)%*%(Y-mean(Y)))+rnorm(N,0,0.5)
  
  alpha_vec=rep(0,num_iter)
  alpha_vec[1]=3
  alphamean=0
  alphavar=1
  alpha_step=0.25
  betamean=0
  betavar=1
  beta_step=0.5
  beta_vec=rep(0,num_iter)
  beta_vec[1]=1
  tsq_vec=rep(0,num_iter)
  tsq_vec[1]=0.001
  

  
  
  
  
  
  a=10
  b=500
  
 
  f=60
  g=100
  
  
  Mt=0.01
  Vt=0.05
  
  h=(Mt^2)/Vt +2
  k=Mt*(h-1)
  
  
  
  alphamean=0
  alphavar=1
  alpha_step=0.5
  mp1_step=2
  mp1var=20
  mp1mean=mean(Y)
  
  mp2_step=0.05
  mp2_var=10
  mp2_mean=-5
  
  mp3_step=0.5
  mp3_var=5
  mp3_mean=5
  
  
  
  
  phi_var=0.3
  gradlogU=function(Z,d,alpha,DesignSample){
    
    nreps=dim(DesignSample)[2]
    chosen_actual=tabulate(d,nbins=N)
    chosen_current=tabulate(DesignSample,nbins=N)/nreps
    alpha*(chosen_actual-chosen_current)
    
    
  }
  
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
  
  
  
  
  alpha_vec[1]=3
  #a starting design sample
  Kreps=350
  design_accept=0
  Design_storeZ=matrix(0,n,Kreps)
  Design_storeZ[,1]=sample(1:N,n,replace=TRUE)
  utilvec=rep(0,Kreps)
  utilvec[1]=log_utility2(Z_reg_store[,1],Design_storeZ[,1],alpha_vec[1],beta_vec[1],distmat_RR)
  for(l in 2:Kreps){
    
    change_sites=sample(1:n,1)
    newsite=sample(1:N,1)
    changedsite=(Design_storeZ[,l-1])[change_sites]
    Design_prop=c(Design_storeZ[-change_sites,l-1],newsite)
    
    
    newutil=log_utility2(Z_reg_store[,1],Design_prop,alpha_vec[1],beta_vec[1],distmat_RR)
    MH_rat1=newutil-utilvec[l-1]
    if(log(runif(1))<MH_rat1){Design_storeZ[,l]=Design_prop
    design_accept=design_accept+1
    utilvec[l]=newutil
    
    }else{
      Design_storeZ[,l]=Design_storeZ[,l-1]
      utilvec[l]=utilvec[l-1]}
  }
 
  #mean parameter mp
  mp_vec=rep(0, num_iter)
  mp_vec[1]=mean(Y)
  mpmean=mean(Y)
  mp_accept=0
  mpvar=5
  #The loop.
  ZRaccept=0
  
  design_vec=rep(n,0)
  for(i in 1:n){
    cur_mon=irreg_locs[i,]
    distvec=rep(0,N)
    for(j in 1:N){distvec[j]= (cur_mon[1]-allpoints[j,1])^2 +(cur_mon[2]-allpoints[j,2])^2}
    design_vec[i]=which(distvec==min(distvec))
  }
  
  
  
  #realspread 
  desdist= distmat_RR[design_vec,design_vec]
  diag(desdist)=10
  minmean= mean(apply(desdist,1,min))
  maxmean=mean(apply(distmat_RR[,design_vec],1,min))
  realspread=-maxmean+minmean

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
    
    mu_RR=mp_vec[i-1]+cor_IR%*%inv_II%*%(Z_irreg_store[,i]-mp_vec[i-1])

    
    dir=gradlogU(Z_reg_store[,i-1],design_vec,alpha_vec[i-1],Design_storeZ)- (1/ssq_vec[i-1])*A_RR%*%((Z_reg_store[,i-1])-(mu_RR))
    
    propnorm=rnorm(N,0,1)
    tp=0.15
    propZ=Z_reg_store[,i-1]+tp*propnorm+(0.5*tp^2)*dir
    
    
    
    
    
    
    #generating designs, first with propZ
    Design_storeZP=Design_storeZ
    
    Design_storeZP[,1]=Design_storeZ[,Kreps]
    design_accept=0
    utilvecp=rep(0,Kreps)
    utilvecp[1]=log_utility2(propZ,Design_storeZP[,1],alpha_vec[i-1],beta_vec[i-1],distmat_RR)
    for(l in 2:Kreps){
      
      change_sites=sample(1:n,1)
      newsite=sample(1:N,1)
      changedsite=(Design_storeZP[,l-1])[change_sites]
      Design_prop=c(Design_storeZP[-change_sites,l-1],newsite)
      #
      
      
      newutil=log_utility2(propZ,Design_prop,alpha_vec[i-1],beta_vec[i-1],distmat_RR)
      MH_rat1z=newutil-utilvecp[l-1]
      if(log(runif(1))<MH_rat1z){Design_storeZP[,l]=Design_prop
      utilvecp[l]=newutil
      
      design_accept=design_accept+1}else{
        Design_storeZP[,l]=Design_storeZP[,l-1]
        utilvecp[l]=utilvecp[l-1]}
      
      
    }
    
    
    #accept/reject
    
    
    ### The utility part of the likelihood ratio
   
    utilrat=alpha_vec[i-1]*(sum(propZ[design_vec])- sum(Z_reg_store[design_vec,i-1]))
    
    #Krat
 
    KR1_vec=rep(0,Kreps)
    
    for(t in 1:Kreps){KR1_vec[t]=-utilvecp[t]+log_utility2(Z_reg_store[,i-1],Design_storeZP[,t],alpha_vec[i-1],beta_vec[i-1],distmat_RR)}
    
    logKR1=-log(Kreps)+max(KR1_vec)+log(sum(exp(KR1_vec-max(KR1_vec))))
    
    
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
    
    utilrat=(anewaold)*sum(Z_reg_store[design_vec,i])
    
    a_pri=-(1/(2*alphavar))*( (propA-alphamean)^2-(alpha_vec[i-1]-alphamean)^2 )
    # #krat
    
    #KR1=0
    KR1_vec=rep(0,Kreps)
    for(t in 1:Kreps){
      CD=Design_storeZ[,t]
      KR1_vec[t]=anewaold*sum(Z_reg_store[CD,i])}
    
    logKR1=-log(Kreps)+max(KR1_vec)+log(sum(exp((KR1_vec-max(KR1_vec))   )))
    

    
    #
    alpha_accept=utilrat+a_pri-logKR1
    #
    if(log(runif(1))<alpha_accept){
      alpha_vec[i]=propA
      #we need a new design for dir
      Design_storeZ[,1]=Design_storeZ[,Kreps]
      utilvec[1]=log_utility2(Z_reg_store[,1],Design_prop,alpha_vec[i],beta_vec[i-1],distmat_RR)
      for(l in 2:Kreps){
        change_sites=sample(1:n,1)
        newsite=sample(1:N,1)
        changedsite=(Design_storeZ[,l-1])[change_sites]
        Design_prop=c(Design_storeZ[-change_sites,l-1],newsite)
        newutil=log_utility2(Z_reg_store[,i],Design_prop,alpha_vec[i],beta_vec[i-1],distmat_RR)
        MH_rat1a=newutil-utilvec[l-1]
        if(log(runif(1))<MH_rat1a){Design_storeZ[,l]=Design_prop
        utilvec[l]=newutil
        }else{
          Design_storeZ[,l]=Design_storeZ[,l-1]
          utilvec[l]=utilvec[l-1]}
      }
      #  
      #  
    }else{alpha_vec[i]=alpha_vec[i-1]}
    
    
    
    # #beta
    propB=beta_vec[i-1] +rnorm(1,0,beta_step)
    bnewbold=propB-beta_vec[i-1]
    
    utilrat=(500*bnewbold)*realspread
    
    b_pri=-(1/(2*betavar))*( (propB-betamean)^2-(beta_vec[i-1]-betamean)^2 )
    
    KR1_vec=rep(0,Kreps)
    
    for(t in 1:Kreps){
      CD=Design_storeZ[,t]
      desdist= distmat_RR[CD,CD]
      diag(desdist)=10
      minmean= mean(apply(desdist,1,min))
      maxmean=mean(apply(distmat_RR[,CD],1,min))
      CDspread=-maxmean+minmean
      KR1_vec[t]=(500*bnewbold)*CDspread }
    
    
    logKR1=-log(Kreps)+max(KR1_vec)+log(sum(exp((KR1_vec-max(KR1_vec)))) )
    
    
    
    
   
    #
    beta_accept=utilrat+b_pri-logKR1
    #
    if(log(runif(1))<beta_accept){
      beta_vec[i]=propB
      #we need a new design for dir
      Design_storeZ[,1]=Design_storeZ[,Kreps]
      utilvec[1]=log_utility2(Z_reg_store[,1],Design_prop,alpha_vec[i],beta_vec[i],distmat_RR)
      for(l in 2:Kreps){
        change_sites=sample(1:n,1)
        newsite=sample(1:N,1)
        changedsite=(Design_storeZ[,l-1])[change_sites]
        Design_prop=c(Design_storeZ[-change_sites,l-1],newsite)
        
        newutil=log_utility2(Z_reg_store[,i],Design_prop,alpha_vec[i],beta_vec[i],distmat_RR)
        
        MH_rat1b=newutil-utilvec[l-1]
        if(log(runif(1))<MH_rat1b){Design_storeZ[,l]=Design_prop
        utilvec[i]=newutil
        }else{
          Design_storeZ[,l]=Design_storeZ[,l-1]
          utilvec[l]=utilvec[l-1]
        }
      }
      #
      #
    }else{beta_vec[i]=beta_vec[i-1]}
    
    #beta_vec[i]=0
    
    
    
    
    #ssq
    

    
    R_C=f+0.5*N
    ZM=as.vector(Z_reg_store[,i]-mp_vec[i-1])
    R_D=g+0.5*t(ZM)%*%inv_RR%*%ZM
    
    # ssq_vec[i]=min(c(rigamma(1,I_C,I_D),rigamma(1,R_C,R_D),rigamma(1,T_C,T_D)))
    ssq_vec[i]=rigamma(1,alpha=R_C,beta=R_D)
    #ssq_vec[i]=ssq
    #logphi
    phi_var=1.5
    prop_logphi=logphi_vec[i-1]+rnorm(1,0,(phi_var))
    
    
    #what these would make the correlation matrices
    
    prop_cor_RR=exp(-distmat_RR/exp(prop_logphi))
    prop_ch_RR=chol(prop_cor_RR)
    prop_inv_RR=chol2inv(prop_ch_RR)
    
    Zr=Z_reg_store[,i]
    Zi=Z_irreg_store[,i]
    
    bigZ=c(Zr,Zi)
    
    
    oldloglike=-sum(log(diag(ch_RR))) - (0.5/ssq_vec[i])*( t(ZM)%*%inv_RR%*%(ZM)  )- b*exp(-logphi_vec[i-1])-a*logphi_vec[i-1]
    newloglike=-sum(log(diag(prop_ch_RR)))-(0.5/ssq_vec[i])*(t(ZM)%*%prop_inv_RR%*%(ZM))-b*exp(-prop_logphi)-a*prop_logphi
    
    
    MH_rat2=newloglike-oldloglike
    
    if(log(runif(1))<MH_rat2){ logphi_vec[i]=prop_logphi
    
    
    cor_IR=exp(-distmat_IR/exp(prop_logphi))
    cor_II=exp(-distmat_II/exp(prop_logphi))
    #cor_RR=exp(-distmat_RR/exp(prop_logphi))
    cor_RR=prop_cor_RR
    ch_II=chol(cor_II)
    inv_II=chol2inv(ch_II)
    # ch_RR=chol(cor_RR)
    ch_RR=prop_ch_RR
    inv_RR=prop_inv_RR
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
    propMP=mp_vec[i-1]+rnorm(1,0,0.5)
    
    mubarZ=propMP+cor_IR%*%inv_II%*%(Z_irreg_store[,i]-propMP) #cholcor already?
    
    
    ZMB=Z_reg_store[,i]-mubarZ
    XMB=Z_irreg_store[,i]-propMP
    loglikepropMP=-0.5*(t(ZMB)%*%A_RR%*%ZMB) -0.5*(t(XMB)%*%Sig1%*%XMB) - (1/(2*mpvar))*(propMP-mpmean)^2
    
    #repeating for current value
    mubarZ=mp_vec[i-1]+cor_IR%*%inv_II%*%(Z_irreg_store[,i]-mp_vec[i-1]) #cholcor already?
    ZMB=Z_reg_store[,i]-mubarZ
    XMB=Z_irreg_store[,i]-mp_vec[i-1]
    
    loglikeMP=-0.5*(t(ZMB)%*%A_RR%*%ZMB) -0.5*(t(XMB)%*%Sig1%*%XMB) - (1/(2*mpvar))*(mp_vec[i-1]-mpmean)^2
    
    MP_MH=loglikepropMP-loglikeMP
    
    if(log(runif(1))< MP_MH ){mp_vec[i]=propMP
    mp_accept=mp_accept+1}else{
      mp_vec[i]=mp_vec[i-1]
    } 
    
    
   if(i%%50==0){print(i)}
    
    
  }
  
  PredZ=rowMeans(Z_reg_store)
  
  predZ_store_full[,p]=PredZ
  Z_store[,p]=Z
  
  alpha_store_full[,p]=alpha_vec
  beta_store_full[,p]=beta_vec
  sigma_store_full[,p]=ssq_vec
  logphi_store_full[,p]=logphi_vec
  tau_store_full[,p]=tsq_vec
  mp_store_full[,p]=mp_vec
  
  
  sse_vec_full[p]=sum((PredZ-Z)^2)
  alpha_vec_full[p]=mean(alpha_vec)
  beta_vec_full[p]=mean(beta_vec)
  sigma_vec_full[p]=mean(ssq_vec)
  logphi_vec_full[p]=mean(logphi_vec)
  tau_vec_full[p]=mean(tsq_vec)
  mp_vec_full[p]=mean(mp_vec)
  alpha_sd_full[p]=sd(alpha_vec)
  beta_sd_full[p]=sd(beta_vec)
  sigma_sd_full[p]=sd(ssq_vec)
  logphi_sd_full[p]=sd(logphi_vec)
  tau_sd_full[p]=sd(tsq_vec)
  mp_sd_full[p]=sd(mp_vec)
  
  
  
  Z_array_full[,,p]=Z_reg_store
  ZI_array_full[,,p]=Z_irreg_store
  
  
  
  print(p)
  print("full model done")
  print(sse_vec_full[p])
  
  #multinomial one
  
  multlogU=function(D,alpha,z,beta,Mdist){
    d=tabulate(D,nbins=N)
    # nfac=log(factorial(sum(D)))
    logUa <- alpha*sum(z[D])
    logfac <- log(sum(factorial(d)))
    
    AZ=alpha*z
    maxAZ=max(AZ)
    pvec=exp(AZ-maxAZ)
    
    K=n*(log(sum(pvec))+maxAZ)
    # logfac <- sum(lgamma(D+1))
    logUa  - logfac -K
  }
  
  
  gradlogmultU=function(Z,D,alpha){
    d=tabulate(D,nbins=N)
    n=sum(d)
    AZ=alpha*Z
    maxAZ=max(AZ)
    pvec=exp(AZ-maxAZ)
    sump=sum(pvec)
    (alpha*d)-(n*alpha*pvec/sump)
  }
  #so now we have data of Y and irreg. locs
  
  
  ##############################################
  #FITTING FROM HERE
  
  
  

  Cormat=exp(-alldistmat/exp(logphi_vec[1]))
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
  
  
  
  #At the moment this is just getting the largest ones
  #mean parameter mp
  mp_vec=rep(0, num_iter)
  mp_vec[1]=mean(Y)
  mpmean=mean(Y)
  mp_accept=0
  mpvar=5
  #The loop.
  ZRaccept=0
  

  
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
    
    mu_RR=mp_vec[i-1]+cor_IR%*%inv_II%*%(Z_irreg_store[,i]-mp_vec[i-1])
  
    dir=gradlogmultU(Z_reg_store[,i-1],design_vec,alpha_vec[i-1])- (1/ssq_vec[i-1])*A_RR%*%((Z_reg_store[,i-1])-(mu_RR))
    
    propnorm=rnorm(N,0,1)
    tp=0.15
    propZ=Z_reg_store[,i-1]+tp*propnorm+(0.5*tp^2)*dir
    
    
    ### The utility part of the likelihood ratio
    utilrat=multlogU(design_vec,alpha_vec[i-1],propZ,beta_vec[i-1],distmat_RR)-multlogU(design_vec,alpha_vec[i-1],Z_reg_store[,i-1],beta_vec[i-1],distmat_RR)
    
    
    
    #now the transition probabilities and the MVN part of the acceptance
    
    triangle_y=gradlogmultU(propZ,design_vec,alpha_vec[i-1])- (1/ssq_vec[i-1])*A_RR%*%((propZ)-mu_RR)
    
    log_qyx=-0.5*N*log(2*pi*tp^2)-(1/(2*tp^2))*sum(( Z_reg_store[,i-1]- propZ -(0.5*tp^2)*triangle_y   )^2)  
    
    log_qxy=-0.5*N*log(2*pi*tp^2)-(1/(2*tp^2))*sum(( -Z_reg_store[,i-1]+ propZ -(0.5*tp^2)*dir   )^2)  
    
    ## The MVN part of the likelihood ratio
    logmvnLR_top=-(0.5/ssq_vec[i-1])*t( propZ-mu_RR   )%*%A_RR%*%(propZ-mu_RR)
    logmvnLR_bottom=-(0.5/ssq_vec[i-1])*t( Z_reg_store[,i-1]-mu_RR   )%*%A_RR%*%(Z_reg_store[,i-1]-mu_RR)
    logmvnLR=logmvnLR_top-logmvnLR_bottom
    #putting them all together
    accept=log_qyx-log_qxy+logmvnLR+ utilrat
    
    
    
    if( log(runif(1))<  accept  ){
      Z_reg_store[,i]=propZ
      
      ZRaccept=ZRaccept+1
    }else{Z_reg_store[,i]=Z_reg_store[,i-1]}
    
    
    #tsq
    new_h=h+0.5*n
    new_k=k+0.5*crossprod((Z_irreg_store[,i]-Y),(Z_irreg_store[,i]-Y))
    tsq_vec[i]=(1/rgamma(1,new_h,new_k))
    
    
    #alpha
    propA=alpha_vec[i-1] +rnorm(1,0,alpha_step)
    anewaold=propA-alpha_vec[i-1]
    
    utilrat=multlogU(design_vec,propA,Z_reg_store[,i],beta_vec[i-1],distmat_RR)-multlogU(design_vec,alpha_vec[i-1],Z_reg_store[,i],beta_vec[i-1],distmat_RR)
    a_pri=-(1/(2*alphavar))*( (propA-alphamean)^2-(alpha_vec[i-1]-alphamean)^2 )
    
    
    
    
    #
    alpha_accept=utilrat+a_pri
    #
    if(log(runif(1))<alpha_accept){
      alpha_vec[i]=propA
      
    }else{alpha_vec[i]=alpha_vec[i-1]}
    
    
    
    
    beta_vec[i]=0
    #ssq
    R_C=f+0.5*N
    ZM=as.vector(Z_reg_store[,i]-mp_vec[i-1])
    R_D=g+0.5*t(ZM)%*%inv_RR%*%ZM
    ssq_vec[i]=rigamma(1,alpha=R_C,beta=R_D)
    
    
    #logphi
    phi_var=1
    prop_logphi=logphi_vec[i-1]+rnorm(1,0,(phi_var))
    
    
    prop_cor_RR=exp(-distmat_RR/exp(prop_logphi))
    prop_ch_RR=chol(prop_cor_RR)
    prop_inv_RR=chol2inv(prop_ch_RR)
    
    Zr=Z_reg_store[,i]
    Zi=Z_irreg_store[,i]
    
    bigZ=c(Zr,Zi)
    
    
    oldloglike=-sum(log(diag(ch_RR))) - (0.5/ssq_vec[i])*( t(ZM)%*%inv_RR%*%(ZM)  )- b*exp(-logphi_vec[i-1])-a*logphi_vec[i-1]
    newloglike=-sum(log(diag(prop_ch_RR)))-(0.5/ssq_vec[i])*(t(ZM)%*%prop_inv_RR%*%(ZM))-b*exp(-prop_logphi)-a*prop_logphi
    
    
    MH_rat2=newloglike-oldloglike
    
    if(log(runif(1))<MH_rat2){ logphi_vec[i]=prop_logphi
    
    
    cor_IR=exp(-distmat_IR/exp(prop_logphi))
    cor_II=exp(-distmat_II/exp(prop_logphi))
    #cor_RR=exp(-distmat_RR/exp(prop_logphi))
    cor_RR=prop_cor_RR
    ch_II=chol(cor_II)
    inv_II=chol2inv(ch_II)
    # ch_RR=chol(cor_RR)
    ch_RR=prop_ch_RR
    inv_RR=prop_inv_RR
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
    propMP=mp_vec[i-1]+rnorm(1,0,0.5)
    
    mubarZ=propMP+cor_IR%*%inv_II%*%(Z_irreg_store[,i]-propMP) 
    
    ZMB=Z_reg_store[,i]-mubarZ
    XMB=Z_irreg_store[,i]-propMP
    loglikepropMP=-0.5*(t(ZMB)%*%A_RR%*%ZMB) -0.5*(t(XMB)%*%Sig1%*%XMB) - (1/(2*mpvar))*(propMP-mpmean)^2
    
    #repeating for current value
    mubarZ=mp_vec[i-1]+cor_IR%*%inv_II%*%(Z_irreg_store[,i]-mp_vec[i-1])
    ZMB=Z_reg_store[,i]-mubarZ
    XMB=Z_irreg_store[,i]-mp_vec[i-1]
    
    loglikeMP=-0.5*(t(ZMB)%*%A_RR%*%ZMB) -0.5*(t(XMB)%*%Sig1%*%XMB) - (1/(2*mpvar))*(mp_vec[i-1]-mpmean)^2
    
    MP_MH=loglikepropMP-loglikeMP
    
    if(log(runif(1))< MP_MH ){mp_vec[i]=propMP
    mp_accept=mp_accept+1}else{
      mp_vec[i]=mp_vec[i-1]
    }
    
    
    
  # print(i)
    
    
  }
  
  #plotting the results
  PredZ=rowMeans(Z_reg_store)
  sse_vec_mult[p]=sum((PredZ-Z)^2)
  alpha_vec_mult[p]=mean(alpha_vec)
  sigma_vec_mult[p]=mean(ssq_vec)
  logphi_vec_mult[p]=mean(logphi_vec)
  tau_vec_mult[p]=mean(tsq_vec)
  mp_vec_mult[p]=mean(mp_vec)
  alpha_sd_mult[p]=sd(alpha_vec)
  sigma_sd_mult[p]=sd(ssq_vec)
  logphi_sd_mult[p]=sd(logphi_vec)
  tau_sd_mult[p]=sd(tsq_vec)
  mp_sd_mult[p]=sd(mp_vec)
  
  
  predZ_store_mult[,p]=PredZ
  
  
  alpha_store_mult[,p]=alpha_vec
  sigma_store_mult[,p]=ssq_vec
  logphi_store_mult[,p]=logphi_vec
  tau_store_mult[,p]=tsq_vec
  mp_store_mult[,p]=mp_vec
  
  Z_array_mult[,,p]=Z_reg_store
  ZI_array_mult[,,p]=Z_irreg_store
  

  print(p)
  print("mult model done")
  print(sse_vec_mult[p])
  
}





####


Z_array_full_cdat=Z_array_full
ZI_array_full_cdat=ZI_array_full
Z_array_mult_cdat=Z_array_mult
ZI_array_mult_cdat=ZI_array_mult

save(Z_array_full_cdat, file="Z_array_full_cdat.Rda")
save(ZI_array_full_cdat, file="ZI_array_full_cdat.Rda")

save(Z_array_mult_cdat, file="Z_array_mult_cdat.Rda")
save(ZI_array_mult_cdat, file="ZI_array_mult_cdat.Rda")


sse_vec_mult_cdat=sse_vec_mult
alpha_vec_mult_cdat=alpha_vec_mult
predZ_store_mult_cdat=predZ_store_mult
alpha_store_mult_cdat=alpha_store_mult
sse_vec_mult_cdat=sse_vec_mult
alpha_vec_mult_cdat=alpha_vec_mult
sigma_vec_mult_cdat=sigma_vec_mult
logphi_vec_mult_cdat=logphi_vec_mult
tau_vec_mult_cdat=tau_vec_mult
mp_vec_mult_cdat=mp_vec_mult
alpha_sd_mult_cdat=alpha_sd_mult
sigma_sd_mult_cdat=sigma_sd_mult
logphi_sd_mult_cdat=logphi_sd_mult
tau_sd_mult_cdat=tau_sd_mult
mp_sd_mult_cdat=mp_sd_mult



predZ_store_full_cdat=predZ_store_full
Z_store_cdat=Z_store
alpha_store_full_cdat=alpha_store_full
beta_store_full_cdat=beta_store_full
sse_vec_full_cdat=sse_vec_full
alpha_vec_full_cdat=alpha_vec_full
beta_vec_full_cdat=beta_vec_full
sse_vec_full_cdat=sse_vec_full
alpha_vec_full_cdat=alpha_vec_full
beta_vec_full_cdat=beta_vec_full
sigma_vec_full_cdat=sigma_vec_full
logphi_vec_full_cdat=logphi_vec_full
tau_vec_full_cdat=tau_vec_full
mp_vec_full_cdat=mp_vec_full
alpha_sd_full_cdat=alpha_sd_full
beta_sd_full_cdat=beta_sd_full
sigma_sd_full_cdat=sigma_sd_full
logphi_sd_full_cdat=logphi_sd_full
tau_sd_full_cdat=tau_sd_full
mp_sd_full_cdat=mp_sd_full


alpha_store_mult_cdat=alpha_store_mult
sigma_store_mult_cdat=sigma_store_mult
logphi_store_mult_cdat=logphi_store_mult
tau_store_mult_cdat=tau_store_mult
mp_store_mult_cdat=mp_store_mult

alpha_store_full_cdat=alpha_store_full
beta_store_full_cdat=beta_store_full
sigma_store_full_cdat=sigma_store_full
logphi_store_full_cdat=logphi_store_full
tau_store_full_cdat=tau_store_full
mp_store_full_cdat=mp_store_full



save(sse_vec_mult_cdat, file="sse_vec_mult_cdat.Rda")
save(alpha_vec_mult_cdat, file="alpha_vec_mult_cdat.Rda")
save(predZ_store_mult_cdat, file="predZ_store_mult_cdat.Rda")
save(alpha_store_mult_cdat, file="alpha_store_mult_cdat.Rda")
save(sse_vec_mult_cdat, file="sse_vec_mult_cdat.Rda")
save(alpha_vec_mult_cdat, file="alpha_vec_mult_cdat.Rda")
save(sigma_vec_mult_cdat, file="sigma_vec_mult_cdat.Rda")
save(logphi_vec_mult_cdat, file="logphi_vec_mult_cdat.Rda")
save(tau_vec_mult_cdat, file="tau_vec_mult_cdat.Rda")
save(mp_vec_mult_cdat, file="mp_vec_mult_cdat.Rda")
save(alpha_sd_mult_cdat, file="alpha_sd_mult_cdat.Rda")
save(sigma_sd_mult_cdat, file="sigma_sd_mult_cdat.Rda")
save(logphi_sd_mult_cdat, file="logphi_sd_mult_cdat.Rda")
save(tau_sd_mult_cdat, file="tau_sd_mult_cdat.Rda")
save(mp_sd_mult_cdat, file="mp_sd_mult_cdat.Rda")



save(predZ_store_full_cdat, file="predZ_store_full_cdat.Rda")
save(Z_store_cdat, file="Z_store_cdat.Rda")
save(alpha_store_full_cdat, file="alpha_store_full_cdat.Rda")
save(beta_store_full_cdat, file="beta_store_full_cdat.Rda")
save(sse_vec_full_cdat, file="sse_vec_full_cdat.Rda")
save(alpha_vec_full_cdat, file="alpha_vec_full_cdat.Rda")
save(beta_vec_full_cdat, file="beta_vec_full_cdat.Rda")
save(sse_vec_full_cdat, file="sse_vec_full_cdat.Rda")
save(alpha_vec_full_cdat, file="alpha_vec_full_cdat.Rda")
save(beta_vec_full_cdat, file="beta_vec_full_cdat.Rda")
save(sigma_vec_full_cdat, file="sigma_vec_full_cdat.Rda")
save(logphi_vec_full_cdat, file="logphi_vec_full_cdat.Rda")
save(tau_vec_full_cdat, file="tau_vec_full_cdat.Rda")
save(mp_vec_full_cdat, file="mp_vec_full_cdat.Rda")
save(alpha_sd_full_cdat, file="alpha_sd_full_cdat.Rda")
save(beta_sd_full_cdat, file="beta_sd_full_cdat.Rda")
save(sigma_sd_full_cdat, file="sigma_sd_full_cdat.Rda")
save(logphi_sd_full_cdat, file="logphi_sd_full_cdat.Rda")
save(tau_sd_full_cdat, file="tau_sd_full_cdat.Rda")
save(mp_sd_full_cdat, file="mp_sd_full_cdat.Rda")


save(alpha_store_mult_cdat, file="alpha_store_mult_cdat.Rda")
save(sigma_store_mult_cdat, file="sigma_store_mult_cdat.Rda")
save(logphi_store_mult_cdat, file="logphi_store_mult_cdat.Rda")
save(tau_store_mult_cdat, file="tau_store_mult_cdat.Rda")
save(mp_store_mult_cdat, file="mp_store_mult_cdat.Rda")

save(alpha_store_full_cdat, file="alpha_store_full_cdat.Rda")
save(beta_store_full_cdat, file="beta_store_full_cdat.Rda")
save(sigma_store_full_cdat, file="sigma_store_full_cdat.Rda")
save(logphi_store_full_cdat, file="logphi_store_full_cdat.Rda")
save(tau_store_full_cdat, file="tau_store_full_cdat.Rda")
save(mp_store_full_cdat, file="mp_store_full_cdat.Rda")


