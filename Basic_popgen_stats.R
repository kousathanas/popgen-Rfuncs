
####################################
#***polymorph.stats***
#calculates polymophism statistics.
#Input: v<-an sfs vector
#Output: a dataframe with the stats
####################################
polymorph_stats<-function(v)
{
  sites=sum(v)
  v=v[-1]
  
  n=length(v)
  sum=0;S=0;pi=0;thita=0;
  a1=0;a2=0;b1=0;b2=0;c1=0;c2=0;e1=0;e2=0;
  
  b1=(n+1)/(3*(n-1))
  b2=2*(n^2+n+3)/(9*n*(n-1))
  
  for (i in 1:(n-1))
  { 
    sum=sum+i*(n-i)*v[i]
    S=S+v[i]
    a1=a1+1/i
    a2=a2+1/i^2
  }
  
  #pi and thita
  pi=2*sum/(n*(n-1))
  thita=S/a1
  
  #Tajima's D
  c1=b1-1/a1
  c2=b2-(n+2)/(a1*n)+a2/(a1^2)
  e1=c1/a1
  e2=c2/(a1^2+a2)
  
  taj_D=(pi-S/a1)/sqrt(e1*S+e2*S*(S-1))
  
  pol_stats=data.frame(
    "singletons"=v[1],
    "S"=S/sites,
    "pi"=pi/sites,
    "thitaw"=thita/sites,
    "TajimasD"=taj_D
  )
  
  return(pol_stats)
}
####################################
#***watterson.factor***
#k number of alleles
####################################
watterson_factor<-function(k)
{
  a_k=0;
  for (i in 1:(k-1))
  { 
    a_k=a_k+1/i
  }
  
  return(a_k);
}
####################################
#***pred.sfs.neut***
#k number of alleles
# see page 273, Charlesworth and Charlesworth (2010)
####################################
pred_sfs_neut<-function(k)
{
  a_k=0;
  v<-vector(length=k-2)
  for (i in 1:(k-1))
  { 
    a_k=a_k+1/i
  }
  
  for (i in 1:(k-1))
  { 
    v[i]=1/(i*a_k)
  }
  
  return(v);
}

####################################
#***div.jukes***
#calculates Jukes-Cantor divergence.
#Input: x<-total sites,y<-site diffs
####################################
div_jukes<-function(x,y)
{
  d<-vector(length=length(x));
  for (i in 1:length(x))
  {
    if (y[i]<=0){d[i]=NA;next;}
    p=y[i]/x[i]
    if ((1-(4/3)*p)<0){d[i]=NA;next;}
    d[i]=(-3/4)*log(1-(4/3)*p)
  }
  
  return(d)
}
####################################
#***calc_ne***
#calculate weighted ne 
#
####################################
calc_nw_2ep<-function(n1,n2,t)
{
  w1=n1*((1-1/(2*n2))^t);
  w2=n2*(1-exp(-t/(2*n2)));
  n_w=(n1*w1+n2*w2)/(w1+w2);
  return(n_w);
}
####################################
calc_nw_3ep<-function(n1,n2,n3,t2,t3)
{
  w1=n1*((1-1/(2*n2))^t2)*((1-1/(2*n3))^t3);
  w2=n2*(1 - exp(-t2/(2*n2)))*((1 - 1/(2*n3))^t3);
  w3=n3*(1 - exp(-t3/(2*n3)));
  n_w=(n1*w1+n2*w2+n3*w3)/(w1+w2+w3);
  return(n_w);
}
####################################
kimura_fix_prob<-function(N,s)
{
  if ((s <= 0.0)&&(N*s > -1e-07))  
  {
    num = 1.0;
    denom = 2.0*N;
  }
  else if ((s >= 0.0)&&(N*s < 1e-07))
  {
    num = 1.0;
    denom = 2.0*N;
  }
  else
  {
    num = 1 - exp(-s);
    denom = 1 - exp(-2*N*s);
  }
  return (num/denom);
}
####################################
fold_sfs<-function(v)
{
  
  l<-length(v)
  q<-floor((l-1)/2)
  
  for (i in 1:q)
  {
    v[i]<-v[i]+v[l+1-i];
    v[l+1-i]<-0;
  }
  return(v);
}
#################################
subsample_SFS<-function(sfs,nsubsample){
  nsample=length(sfs)-1
  Vsub<-vector(length=nsubsample)
  for (i in 1:sum(sfs[-1])){
    f<-sample(1:nsample,prob=sfs[-1]/sum(sfs[-1]),size=1,replace=TRUE)
    l<-rhyper(1,f,nsample-f,nsubsample)
    Vsub[l]=Vsub[l]+1
  }
  
  Vsub2<-c(sum(sfs)-sum(Vsub),Vsub)
  
  return(Vsub2)
}
#################################
integrate_fix_prob_gamma<-function(N,beta,mean){
  x<- -rgamma(1e6,beta,beta/mean)
  x[which(x< -1)]=-1
  fix_prob<-vector(length=length(x))
  
  for (i in 1:length(x)){
    if (x[i] == 0){
      fix_prob[i] = 0.5/N;
    }
    else {
      fix_prob[i] = kimura_fix_prob(N,x[i]);
    }  
  }
  
  return(mean(fix_prob))
}
#################################
calculate_load<-function(q,h,s){
  return(s*(2*h*q+(1-2*h)*q^2))
}
