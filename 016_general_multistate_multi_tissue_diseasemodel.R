#  multistate_occ_modl

set.seed(54321); epsln=1e-5;

ofile='multistate_occ_modlsf.out'#output file

NCPUS=as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')); NREPS=NCPUS*2

##loading the likelihood 3 function

if (is.loaded('like3')) dyn.unload('like3a.so'); sink(); sink(ofile,split=T)
#in case of windows use 'ddl' extension instead of 'so'

#Example scenario of a multi-organ tissue system, where 's' tissues can be 
#infected and are sampled and tested 'k' times
#in this scenario lets say s=4 tissues were samples k=3 times for n=20 individuals
K=K #k is the number of surveys undertaken

psi=c(.5,.4,.3,.6); p=c(.2,.8,.3,.7)#initital values

psi=psi[1:4] # selecting psi values from intiital values
p=p[1:length(psi)] #similar p
S=length(psi)        #  number of sub-levels of occupancy (4)
SS=2^S               #  number of possible states, given sub-levels, including 0000
SSm1=SS-1            #  number of possible states, given sub-levels

#creating multiple possible state combination
trustate=NULL;for (i in 0:SSm1)
  trustate=rbind(trustate,as.integer(intToBits(i)[S:1])); 
#OBNG is oral blood nasal and genital tissue in our example
#Can be any tissues
states=OBNG=apply(trustate,1,paste,collapse='') 

#Creating a sample data

###create a fake sample of detection non detection for each of the tissue
sample1<- c(0,1)
a2<- cbind(sample(sample1,20, replace=T),sample(sample1,20, replace=T),sample(sample1,20, replace=T),
           sample(sample1,20, replace=T))
#each tissue sampled three times
a3<- cbind(a2,a2,a2)
a3<-as.data.frame(a3)
#give the column names, for.e.g, t1s1 refers to tissue 1 survey 1
colnames(a3)<- c("t1s1","t1s2","t1s3","t2s1","t2s2","t2s3","t3s1","t3s2","t3s3",
                 "t4s1","t4s2","t4s3")

N=nrow(a3);  nms=names(a3)
#data manipulation stuff...

h1=cbind(a3$t1s1,a3$t2s1,a3$t3s1,a3$t4s1) # det.hist. for sample 1
h2=cbind(a3$t1s2,a3$t2s2,a3$t3s2,a3$t4s2)# det.hist. for sample 2
h3=cbind(a3$t1s3,a3$t2s3,a3$t3s3,a3$t4s3) # det.hist. for sample 3
# combine samples...
h4=h40=cbind(h1,h2,h3); h40[is.na(h4)]=0; h4[is.na(h4)]=2  #  h4=det.hist as matrix,

#in case there are missing values in the data
h1a=gsub('NA','.',apply(h1,1,paste,collapse=''))           #  h40=same thing w/ missing replaced w/ 0's
h2a=gsub('NA','.',apply(h2,1,paste,collapse=''))
h3a=gsub('NA','.',apply(h3,1,paste,collapse=''))
t1=table(c(h1a,h2a,h3a));                         #  table if histories, disregarding sampling occasion
#print(cbind(1:length(t1),t1))                    #  print table for data exploration

h4a=apply(h4,1,paste,collapse=''); t4=table(h4a); #  table of histories
#print(cbind(1:length(t4),t4))                    #  print table for data exploration

dhst=h4; dfrq=rep(1,nrow(dhst))                   #  no longer summarizing histories as we'll be using
obsstate=array(as.integer(dhst),dim=c(nrow(dhst),S,K)); #  individual covariates
obsstate0=obsstate; obsstate0[obsstate==2]=0            #  obsstate=matrix of det. histories, obsstate0=same w/ missing replaced by 0's

mlogit<-function(a) {  #' computes multinomial logit of A ( exp(A)/(1+sum(exp(A))))
  x=exp(a)/(1+sum(exp(a)))
  return(x)
}

#Like 1 computes neg loglikelihood when we assume that there is no need to believe
#that infection in one tissue is dependent on another tissue
#under the assumption ofindependence this model resembles single season occupancy model
like1<- function(theta,opt,desmat,desmat2) {  #  computes the negative log-likelihood of multi-state occupancy model 
  #with independent sub-states
  # -----------------------------------------------------------------
  #  this is original function... included in source for reference.
  # 
  # ----------------------------------------------------------------
  NN=nrow(obsstate0); loglike=0; npar1=ncol(desmat); npar2=ncol(desmat2)
  psi=p=rep(1,S)
  for (i in 1:NN) {              #   for each detection history...
    for (j in 1:S) {
      psi[j]=plogis(sum(desmat[i+(j-1)*NN,] * theta[1:npar1])) # occupancy prob
      p[j]=plogis(sum(desmat2[i+(j-1)*NN,] * theta[npar1+1:npar2])) #detection prob
    }
    ms_psi=rep(1,SS);
    for (tsi in 2:SS) for (j in 1:S) ms_psi[tsi]=ms_psi[tsi]*c(1-psi[j],psi[j])[trustate[tsi,j]+1]
    ms_psi[1]=1-sum(ms_psi[-1]);
    qp=cbind(1-p,p,1)
    prd=ms_psi                   #       compute product of psi(tru-state) * p1(obs-state|tru-state)*p2(obs-state|tru-state)...
    for (srvy in 1:K) {          #       for each survey...
      for (tsi in 1:SS) {        #          for each tru-state...   (tsi=trustate index -> 1=0000, 2=0001, 3=0010, 4=0011,... 16=1111)
        if (prod((trustate[tsi,]*obsstate0[i,,srvy])==obsstate0[i,,srvy])>0) {  #  if obsstate is possible for given trustate...
          if (tsi>1)            #              ... don't compute product of p's if tru-state='0000'
            for (tiss in 1:S)     #              for each of the 4 sample types (O,B,N,G)...
              #                  if tru-state for sample-tpe tiss == 1, then multiply product for
              #                  tru-state by either 1-p, p, or 1, depending on the
              #                  observed state for the sample-type.  So, the row
              #                  of matrix qp is for the sample-type, and the column
              #                  indicates the quantity to multiply: (1-p) if not detected,
              #                  p if detected, and 1 if missing data.
              if (trustate[tsi,tiss]>0) prd[tsi]=prd[tsi]*qp[tiss,obsstate[i,tiss,srvy]+1]
        } else                  #             otherwise (obs-state not possible for tru-state),
          prd[tsi]=0;        #                  set product for that tru-state to zero.
      }
    }
    loglike=loglike+dfrq[i]*log(sum(prd))  #  sum likelihood for det.history i
  }
  
  return(-loglike)
}
#  like 2 computes the negative log-likelihood of multi-state occupancy model 
#when no assumption of independnece can  be made
#or we know for sure there is  interaction among sub-states
like2<- function(theta,opt,desmat,desmat2) {  #  computes the negative log-likelihood of multi-state occupancy model with interaction among sub-states
  # -----------------------------------------------------------------
  #  this is original function... included in source for reference.
  #  it is replaced by C function, like3 with opt argument=1
  # ----------------------------------------------------------------
  NN=nrow(obsstate0); loglike=0; npar1=ncol(desmat); npar2=ncol(desmat2)
  psi=rep(1,SSm1); p=rep(1,S)
  for (i in 1:NN) {              #   for each detection history...
    for (j in 1:SSm1) psi[j]=exp(sum(desmat[i+(j-1)*NN,] * theta[1:npar1]))
    for (j in 1:S)  p[j]=plogis(sum(desmat2[i+(j-1)*NN,] * theta[npar1+1:npar2]))
    ms_psi=psi/(1+sum(psi));
    ms_psi=c(1-sum(ms_psi),ms_psi);
    qp=cbind(1-p,p,1)
    prd=ms_psi                   #       compute product of psi(tru-state) * p1(obs-state|tru-state)*p2(obs-state|tru-state)...
    for (srvy in 1:K) {          #       for each survey...
      for (tsi in 1:SS) {        #          for each tru-state...   (tsi=trustate index -> 1=0000, 2=0001, 3=0010, 4=0011,... 16=1111)
        if (prod((trustate[tsi,]*obsstate0[i,,srvy])==obsstate0[i,,srvy])>0) {  #  if obsstate is possible for given trustate...
          if (tsi>1)            #              ... don't compute product of p's if tru-state='0000'
            for (tiss in 1:S)     #              for each of the 4 sample types (O,B,N,G)...
              #                  if tru-state for sample-tpe tiss == 1, then multiply product for
              #                  tru-state by either 1-p, p, or 1, depending on the
              #                  observed state for the sample-type.  So, the row
              #                  of matrix qp is for the sample-type, and the column
              #                  indicates the quantity to multiply: (1-p) if not detected,
              #                  p if detected, and 1 if missing data.
              if (trustate[tsi,tiss]>0) prd[tsi]=prd[tsi]*qp[tiss,obsstate[i,tiss,srvy]+1]
        } else                  #             otherwise (obs-state not possible for tru-state),
          prd[tsi]=0;        #                  set product for that tru-state to zero.
      }
    }
    loglike=loglike+dfrq[i]*log(sum(prd))  #  sum likelihood for det.history i
  }
  
  return(-loglike)
}





m=list(); i=0
#1
i=i+1; m[[i]]=indep_mod(~tissue,~tissue); print(m[[i]]$results,quote=F)
#2

i=i+1; m[[i]]=interaction_mod(~state,~tissue); print(m[[i]]$results,quote=F);

nmods=i; aictbl=mnames=NULL
for (i in 1:nmods) { aictbl=c(aictbl,m[[i]]$aic); mnames=c(mnames,gsub(' AIC.+','',m[[i]]$results$modname))}
o=order(aictbl); minaic=min(aictbl); aictbl=cbind(aictbl,aictbl-minaic); 
rownames(aictbl)=mnames; colnames(aictbl)=c('AIC','deltaAIC'); print(aictbl[o,])

x=print_psi_indep(m[[o[1]]]); cat('\nEstimates from top model:\n'); print(x)

sfStop(); sink()
save.image(file='multistate_occ_modl6s_results.rdata')


############################





