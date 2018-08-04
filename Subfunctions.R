#### AGM -- Additive Genetic Model

## Lemma 1: Consider additive Genetic Model (AGM): Y=G*b+e, 
## where G~B(2, p), b=beta, e~N(0,1), G and e are indepedent.
## Searching for La and Ua for given alpha, b and p
##################################################################

#### Compute CDF of AGM Y=G*b+e, where G~B(2, p) is indepedent of e~N(0, 1). 
#### Specificly, Pr(Y<x) for a given vector x 
####
#### INPUTS
#### x: a specific point
#### b: effect size of G on Y 
#### p: MAF at the SNP
####
#### RETURN
#### The probability
####
#### Created by Huaizhen Qin on July 26, 2014
############################################################################
CDF.AGM <- function(x, b, p) 
{
  q=1-p;
  ip = q*q*pnorm(x, lower.tail = TRUE)+2*p*q*pnorm(x-b, lower.tail = TRUE)+p*p*pnorm(x-2*b, lower.tail = TRUE);
  
  return (ip);
}

#### Compute PDF of Y=G*b+e, where G~B(2, p) is indepedent of e~N(0, 1). For a 
#### fixed y, the PDF value is E(phi(x-G*b)). 
####
#### INPUTS
#### x: a vector of specific points
#### b: a vector of effect sizes of G on Y 
#### p: MAF at the SNP
####
#### RETURN
#### A vector of PDF values 
####
#### Created by Huaizhen Qin on July 26, 2014
############################################################################
PDF.AGM <- function(x, b, p) 
{
  q=1-p;
  ip = q*q*dnorm(x)+2*p*q*dnorm(x-b)+p*p*dnorm(x-2*b);
  
  return (ip);
}


#### Compute upper tail probability (UTP) of AGM Y=G*b+e, where G~B(2, p) is 
#### indepedent of e~N(0, 1). Pr(Y>x) for a given vector x, decensing w.r.t. x
####
#### INPUTS
#### x: a specific point
#### b: effect size of G on Y 
#### p: MAF at the SNP
####
#### RETURN
#### The probability
####
#### Created by Huaizhen Qin on July 26, 2014
#### Revised by Huaizhen Qin on July 26, 2014
############################################################################

UTP.AGM <- function(x, b, p)  
{
  q=1-p;
  ip = q*q*pnorm(x, lower.tail = FALSE)+2*p*q*pnorm(x-b, lower.tail = FALSE)+p*p*pnorm(x-2*b, lower.tail = FALSE);
  
  return (ip);
}



#### compute lower quantile of Y=G*b+e, where G~B(2, p) and e~N(0, 1)
####
#### INPUTS
#### x: a specific point
#### b: effect size of G on Y 
#### p: MAF at the SNP
#### alpha: the level of the quantile
####
#### RETURN
#### lower alpha quantile of Y 
####
#### Created by Huaizhen Qin on July 26, 2014
############################################################################

LowerQ <- function(b, p, alpha)
{
  f.lower <- function(x, b, p, alpha) # difference CDF.AGM-alpha
  {
    d = CDF.AGM(x, b, p) - alpha;	
    return (d);
  }
  r = uniroot(f.lower, c(-50000, 50000), b, p, alpha, tol = .Machine$double.eps^.95, maxiter = 5000);
  
  return (r$root);
}


#### Compute upper quantile of Y=G*b+e, where G~B(2, p) and e~N(0, 1)
####
#### INPUTS
#### x: a specific point
#### b: effect size of G on Y 
#### p: MAF at the SNP
#### alpha: the level of the quantile
####
#### RETURN
#### upper alpha quantile of Y 
####
#### Created by Huaizhen Qin on July 26, 2014
############################################################################

UpperQ <- function(b, p, alpha)
{
  f <- function(x, b, p, alpha) # difference UTP.AGM-alpha
  {
    SP = UTP.AGM(x, b, p) - alpha;	
    return (SP);
  }
  r = uniroot(f, c(-50000, 50000), b=b, p=p, alpha=alpha, tol = .Machine$double.eps^.85, maxiter = 5000);
  
  return (r$root);
}



#### Compute power for testing mediate effect of the SMCVM 
####
#### INPUTS
#### LT2: [lambda2, tau2]
#### NSL: nominal significance level
#### n: sample size
#### 
#### RETURN
#### Power
####
#### NOTE: Works for both SRS and EPS
####
#### Created by Huaizhen Qin on July 22, 2014
############################################################

MCVM.Power <- function(LT2, NSL, n)
{
  d2 = n*LT2[,2];
  
  ca = qchisq(p=NSL, df=1, ncp = 0,  lower.tail = FALSE, log.p = FALSE);
  pw = pchisq(q=ca/LT2[,1], df=1, ncp = d2, lower.tail = FALSE, log.p = FALSE);
  
  return(pw);
}


#### Search sample size required to achieve 
#### a given power level (e.g., 0.8)
#### 
#### INPUTS
#### LT2: matrix [lambda2, tau2]
#### NSL: nominal significance level
#### pwr: sample size
####
#### OUTPUTS
#### Vector of sample sizes
####
#### NOTE: Works for both SRS and EPS
####
#### Created by Huaizhen Qin on July 22, 2014
############################################################

MCVM.n <- function(LT2, NSL, pwr, n.range)
{
  ca = qchisq(p=NSL, df=1, ncp = 0,  lower.tail = FALSE, log.p = FALSE);
  
  pwr.ni<-function(n, LT2i, ca, pwr)
  {
    d2 = n*LT2i[2];
    pw = pchisq(q=ca/LT2i[1], df=1, ncp = d2, lower.tail = FALSE, log.p = FALSE);
    return(pw-pwr);
  }
  
  n<-c();
  for(i in c(2:dim(LT2)[1]))
  {
    n[i]= uniroot(pwr.ni, n.range, LT2i=LT2[i, ], ca=ca, pwr=pwr, 
                  tol = .Machine$double.eps^.8, maxiter = 1000)$root;
  }
  
  return(n);
}


#### Lemma 1: The first 4 moments of X
####
#### INPUTS
#### G.Ms: The first 4 moments of G
#### u.Ms: The first 4 moments of u
#### 
#### OUTPUTS 
#### The first 4 moments of X
####
#### NOTE: Works for simple random sampling (SRS)
####
#### Created by Huaizhen Qin on July 23, 2014
############################################################

X.Moments <- function(r, G.Ms, u.Ms)
{
  x.m1 <- G.Ms[1]*r;
  x.m2 <- G.Ms[2]*r*r+u.Ms[,2];
  x.m3 <- G.Ms[3]*r*r*r+3*r*G.Ms[1]*u.Ms[,2];
  x.m4 <- G.Ms[4]*r*r*r*r+6*r*r*G.Ms[2]*u.Ms[, 2]+u.Ms[,4];
  
  x.s2 <- (x.m2 - x.m1*x.m1); 
  
  X.Ms <- cbind(x.m1, x.m2, x.m3, x.m4, x.s2);
  
  return(X.Ms);
}


#### Proposition 1: Power comparison under SRS
##############################################################

#### Compute scale and non-centrality of the 
#### standard correlation test 
####
#### INPUTS
#### X.Ms: The first 4 moments of X
#### e.Ms: The first 4 moments of e
#### 
#### OUTPUTS
#### Scale and noncentrality
####
#### NOTE: Works for SRS only
####
#### Created by Huaizhen Qin on July 23, 2014
############################################################

L2T2.SRS <- function(b, X.Ms, e.Ms)
{
  bt2 <- b*b;
  l2 <- 1+(bt2/e.Ms[,2])*(6*(X.Ms[,1])^2-X.Ms[,5]+(3*(X.Ms[,1])^4 - 4*X.Ms[,1]*X.Ms[,3]+X.Ms[,4])/X.Ms[,5]);
  t2 <- bt2*X.Ms[,5]/(l2*e.Ms[,2]);
  
  lt2 <- cbind(l2, t2);
  
  return(lt2);
}

#### Compute the expectation of G*phi(x-G*b) under the AGM Y=G*b+e, 
#### where G~B(2, p) is indepedent of e~N(0, 1)
####
#### INPUTS
#### x: a specific point
#### b: effect size of G on Y 
#### p: MAF at the SNP
####
#### RETURN
#### The expectation of G*phi(x-G*b), the numerator of the 1st order 
#### derivative of quantile 
####
#### Created by Huaizhen Qin on August 16, 2014
############################################################################

N.d1.eta.AGM <- function(x, b, p) 
{
  q=1-p;
  return(2*p*(q*dnorm(x-b)+p*dnorm(x-2*b)));
}


#### Numerator of second order derivative of quantiles
N.d2.eta.AGM<-function(eta, d.eta, b, p)
{
  q=1-p;
  a0<-(q*q)*(eta)*((d.eta)^2)*dnorm(eta);
  a1<-(2*p*q)*(eta-b)*((1-d.eta)^2)*dnorm(eta-b);
  a2<-(p*p)*(eta-2*b)*((2-d.eta)^2)*dnorm(eta-2*b);
  
  return(a0+a1+a2);
}



#### Compoute the first order derivative of eta.L as a function of b 
#### for given alpha and p
####
#### INPUTS
#### b: model coefficients;
#### alpha:truncation level;
#### p: MAF;
####
#### RETUN
#### value of the first derivative of eta_(aplha, p) at the point b
####
#### Created by Huaizhen Qin on Aug 18, 2014
##############################################################################

d1.eta.La<-function(b, alpha, p)
{
  eta.L<-LowerQ(b, p, alpha);
  E.phi.L<-PDF.AGM(eta.L, b, p);
  N.d.eta.L<-N.d1.eta.AGM(eta.L, b, p);
  d.eta.L<-N.d.eta.L/E.phi.L;
  
  return (d.eta.L)
}


#### Compute the second order derivative of eta.L as a function of b 
#### for given alpha and p
####
#### INPUTS
#### b: model coefficients;
#### alpha:truncation level;
#### p: MAF;
####
#### RETUN
#### value of the second order derivative of eta_(aplha, p) at the point b
####
#### Created by Huaizhen Qin on Aug 18, 2014
##############################################################################

d2.eta.La<-function(b, alpha, p)
{
  eta.L<-LowerQ(b, p, alpha);
  E.phi.L<-PDF.AGM(eta.L, b, p);
  d.eta.L<-d1.eta.La(b, alpha, p);
  N.d2.eta.L<-N.d2.eta.AGM(eta.L, d.eta.L, b, p);
  d2.eta.L<-N.d2.eta.L/E.phi.L;
  
  return(d2.eta.L);
}


#### Define the first order derivative of eta.U as a function of b 
#### for given alpha and p
####
#### INPUTS
#### b: model coefficients;
#### alpha:truncation level;
#### p: MAF;
####
#### RETUN
#### value of the first derivative of eta_(1-aplha, p) at the point b
####
#### Created by Huaizhen Qin on Aug 18, 2014
##############################################################################

d1.eta.Ua<-function(b, alpha, p)
{
  eta.U<-UpperQ(b, p, alpha);
  E.phi.U<-PDF.AGM(eta.U, b, p);
  N.d.eta.U<-N.d1.eta.AGM(eta.U, b, p);
  d.eta.U<-N.d.eta.U/E.phi.U;
  
  return (d.eta.U)
}


#### Compute the second order derivative of eta.U as a function of b 
#### for given alpha and p
####
#### INPUTS
#### b: model coefficients;
#### alpha:truncation level;
#### p: MAF;
####
#### RETUN
#### value of the second order derivative of eta_(1-aplha, p) at the point b
####
#### Created by Huaizhen Qin on Aug 18, 2014
##############################################################################

d2.eta.Ua<-function(b, alpha, p)
{
  eta.U<-UpperQ(b, p, alpha);
  E.phi.U<-PDF.AGM(eta.U, b, p);
  d.eta.U<-d1.eta.Ua(b, alpha, p);
  N.d2.eta.U<-N.d2.eta.AGM(eta.U, d.eta.U, b, p);
  d2.eta.U<-N.d2.eta.U/E.phi.U;
  
  return(d2.eta.U);
}


#### compute E(G^k|Y<La) under model Y=Gb+e and and EPS, 
#### where e~N(0, sd.e2), Pr(Y<La)=alpha
#### 
#### INPUTS
#### k: target order, e.g., 1, 2
#### La: Lower alpha quantile of Y
#### b: Effect size of G
#### p: MAF at the SNP
#### alpha: truncation level for EPS (0~0.5)
####
#### OUTPUTS
#### The k-th conditional moment of genotype G given Y<La 
####
#### Created by Huaizhen Qin on July 31, 2014
############################################################################

mkg.L.EPS<-function(k, La, b, p, sd.e, alpha)
{
  G=c(0, 1, 2);
  Gk=G^k; 
  
  PrG <- dbinom(G, 2, p);
  PLa <- pnorm((La-b*G)/sd.e, lower.tail = TRUE);
  mkg <- sum(Gk*PrG*PLa)/alpha;
  
  return(mkg);
}



#### compute E(G^k|Y>Ua) under model Y=Gb+e and EPS, 
#### where e~N(0, sd.e2), Pr(Y>Ua)=alpha
#### 
#### INPUTS
#### k:  Target order, e.g., 1, 2
#### Ua: Lower alpha quantile of Y
#### b: Effect size of G
#### p: MAF at the SNP
#### alpha: (0~0.5)
####
#### OUTPUTS
#### The k-th conditional moment of genotype G given Y>Ua 
####
#### Created by Huaizhen Qin on July 31, 2014
############################################################################

mkg.U.EPS<-function(k, Ua, b, p, sd.e, alpha)
{
  G=c(0, 1, 2);
  Gk=G^k; 
  
  PrG <- dbinom(G, 2, p);
  PUa <- pnorm((Ua-b*G)/sd.e, lower.tail = FALSE);
  mkg <- sum(Gk*PrG*PUa)/alpha;
  
  return(mkg);
}


#### compute noncentrality under model Y=Gb+e and EPS, 
#### where e~N(0, sd.e2), Pr(Y>Ua)=alpha
#### 
#### INPUTS
#### k:  Target order, e.g., 1, 2
#### Ua: Lower alpha quantile of Y
#### b: Effect size of G
#### p: MAF at the SNP
#### alpha: (0~0.5)
####
#### OUTPUTS
#### The k-th conditional moment of genotype G given Y>Ua 
####
#### Created by Huaizhen Qin on July 31, 2014
############################################################################

delta2.G.EPS<-function(La, Ua, b, p, sd.e, alpha)
{
  
  m1x.0 <- mkg.L.EPS(k=1, La, b, p, sd.e, alpha);
  m2x.0 <- mkg.L.EPS(k=2, La, b, p, sd.e, alpha);
  
  m1x.1 <- mkg.U.EPS(k=1, Ua, b, p, sd.e, alpha);
  m2x.1 <- mkg.U.EPS(k=2, Ua, b, p, sd.e, alpha);
  
  dlt2<-(m1x.1-m1x.0)^2/(m2x.0-m1x.0^2+m2x.1-m1x.1^2);
  
  return(dlt2);
}


#### Compute the noncentrality for SNP t test under EPS
####
#### INPUTS 
#### tau.LU: The lower and upper alpha-quantiles of Y
#### beta1,...: vectors of model coefficients in the MCM model
#### p: minor allele frequency
#### alpha: EPS level for each tail
#### 
#### OUTPUTS
#### Matrix = [1, noncetrality values]
####
#### Created by Huaizhen Qin on July 31, 2014
############################################################################

L2T2.SNP.EPS<-function(tau.LU, beta1, beta2, beta3, p, alpha)
{
  ## Under the MCM, Y=G*b+e, where G=X1~B(2, p), 
  ## b=beta1*beta2*beta3, e=beta2*beta3*e2+beta3*e3+e4
  
  dlt2<-c();
  La=tau.LU[,1];
  Ua=tau.LU[,2];
  b=beta1*beta2*beta3;
  sd.e=sqrt(1+beta3^2+(beta2*beta3)^2);
  
  lb=length(beta1);
  for(i in c(1:lb))
  {
    
    dlt2[i]<-delta2.G.EPS(La[i], Ua[i], b[i], p, sd.e[i], alpha);
    
  }
  
  LT2<-cbind(1, dlt2);
  
  return(LT2);
}



#### Under the SMCVM, compute the fold change in X for given b
#### 
#### INPUTS 
#### b: effect size of X on Y=X*b+e
#### p: minor allele frequency in model X=G*r+u
#### r: effect size of G on X = G*r + u
#### sd.e: standard deviation of e
#### sd.u: standard deviation of u
#### alpha: truncation parameter 
#### 
#### OUTPUT
#### nfold: number of fold change for given parameters
####
#### RATIONALE
#### Consider the SMCVM 
#### Y=X*b+e, e~N(0, sd.e^2)
#### X=G*r+u, G~B(2, p), u~N(0, sd.u^2)
#### 
#### First, I compute La, Ua such that (Pr(Y<La)=Pr(Y>Ua)=a
#### Y=G*(b*r)+(e+b*u)
#### (Y/se)=G(b*r/se)+z, z~N(0,1)
####
#### Then, I compute the fold change as:
#### fc=E[exp(cX)|Y>Ua)/E[exp(cX)|Y<La), where c=log(2)
####
#### Created by Huaizhen Qin on August 1, 2014
############################################################################

foldchange<-function(b, p, r, sd.e, sd.u, alpha, c)
{
  se=sqrt(sd.e^2+(b*sd.u)^2);
  be=b*r/se;
  
  La<-se*LowerQ(be, p, alpha);
  Ua<-se*UpperQ(be, p, alpha);
  
  m1x.0<-integrate(intgrand.expx.0, La, b, p, r, sd.e, sd.u, alpha, c, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.85, subdivisions = 2000L)[[1]];
  m1x.1<-integrate(intgrand.expx.1, Ua, b, p, r, sd.e, sd.u, alpha, c, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.85, subdivisions = 2000L)[[1]];
  
  fc<-(m1x.1/m1x.0);
  
  return(c(fc, m1x.1, m1x.0));
}


##foldchange(b, p, r, sd.e, sd.u, alpha, c)

#### Under the SMCVM, compute the difference of fold change in X for given b
#### and a given nfold
#### 
#### INPUTS 
#### b: effect size of X on Y=X*b+e
#### p: minor allele frequency in model X=G*r+u
#### r: effect size of G on X = G*r + u
#### sd.e: standard deviation of e
#### sd.u: standard deviation of u
#### alpha: truncation parameter 
#### nfold: a given number of fold change
#### 
#### OUTPUT
#### fold.change-nfold
####
#### Created by Huaizhen Qin on Aug 2nd, 2014
##############################################################################

foldchange.eq<-function(b, p, r, sd.e, sd.u, alpha, c, nfold)
{
  
  return (foldchange(b, p, r, sd.e, sd.u, alpha, c)[1]-nfold);
}


#### Created by Huaizhen Qin on Aug 2nd, 2014
##############################################################################

nfold2beta<-function(nfold, p, r, sd.e, sd.u, alpha, c)
{
  beta = uniroot(foldchange.eq, c(0, 50), p, r, sd.e, sd.u, alpha, c, nfold, tol = .Machine$double.eps^.85, maxiter = 1000)$root;
  
  return (beta);
}


##nfold2beta(nfold, p, r, sd.e, sd.u, alpha, c)


## Search for the bt3 for given bt1, bt2 and nfold in protein expression (e.g., 3). 
nfold2beta.PRT<-function(nfold, alpha, bt1, bt2, p)
{
  #### Search beta3 given other parameters  
  var.x1 = 2*p*(1-p);		#### var(X1) = var(G);
  var.x2 = 1+(bt1*bt1)*var.x1;	#### var(X2) = var(e2+X1*bt1);
  var.x3 = 1+(bt2*bt2)*var.x2;	#### var(X3) = var(e3+X2*bt2);
  
  r=bt1*bt2; 
  sd.e=1;
  sd.u=sqrt(1+bt2^2);
  bt3<-nfold2beta(nfold, p, r, sd.e, sd.u, alpha, c);
  
  return (c(bt3, var.x3));
}


#### Determine the effect size beta3 of protein expression X3=log2(T=intensity)
#### on trait value
####
#### INPUTS
#### nfold: number of fold change in protein intensity
#### alpha: truncation parameter
#### p: MAF of the causal SNP, e.g., p =0.25
#### beta1: effect size of SNP on RNA
#### beta2: effect size of RNA on PRT
#### var.x3: varaince of X3
####
#### RETURN 
#### matrix [beta3, r2.3]
####	
#### RATIONALE
#### For this, I rewrite the MCM as a SMCVM:
####
#### Y=X*b+e: X=X3, b = beta3, e=e4~N(0,1)
#### X=G*r+u: G=X1~B(2,p), r=beta1*beta2, u=beta2*e2+e3~N(0, 1+beta2^2)
#### 
#### Created by Huaizhen Qin on Aug 2nd, 2014
#### Revised by Huaizhen Qin on Aug 11, 2014
###############################################################################

effect.sizes.PRT<-function(nfold, alpha, p, beta1, beta2, var.x3)
{
  ## First, I search the beta3.max according to beta1.max, beta2.max 
  ## and given nfold in protein expression (e.g., 3). 
  
  lb=length(beta1);
  #npts=lb-1;
  #r=beta1[lb]*beta2[lb]; 
  #sd.e=1;
  #sd.u=sqrt(1+beta2[lb]^2);
  #beta3.max<-nfold2beta(nfold, p, r, sd.e, sd.u, alpha);
  
  out.max=nfold2beta.PRT(nfold, alpha, bt1=beta1[lb], bt2=beta2[lb], p);
  
  
  ## Next, I compute the r2.3.max and r2.3 as below:
  
  #var.x3.max = 1+(beta2[lb])^2*var.x2[lb];
  r2.3.max<-out.max[1]^2*out.max[2]/(1+out.max[1]^2*out.max[2]);
  r2.3 <- 100*r2.3.max*rpts;
  
  beta3 = sqrt(r2.3/( (1-r2.3)*var.x3 ) );
  
  return (cbind(beta3, r2.3));
}


#### Power of SNP test
#### G.Ms: The first 4 moments of G
#### beta1, beta2, beta3: The model coefficients: 
#### X2 = X1*beta1+e2, 
#### X3 = X2*beta2+e3,
#### Y  = X3*beta3+e4,
#### e2, e3, e4 iid N(0,1)
####
#### Rewrite the MCM as SIMPLE MCV MODEL
#### Y=X*b+e: X=X1, b=beta1*beta2*beta3, e=beta2*beta3*e2 + beta3*e3 + e4
#### X=Gr+u:  G=X1, r=1, u=0: Pr(X=G)=Pr(u=0)=1 
####
#### Created by Huaizhen Qin on July 23, 2014
#########################################################################

SNP.Geno.pwr.SRS<-function(G.Ms, beta1, beta2, beta3, NSL, n)
{	
  r = rep(1, length(beta1));#gamma in the simple model (Fig S1)
  X.Ms <- X.Moments(r, G.Ms = G.Ms, u.Ms = cbind(0*r, 0, 0, 0));
  
  
  #epsilon1
  e.s2 <- 1+(beta3^2)+(beta2*beta3)^2; 
  e.Ms <- cbind(0, e.s2, 0, 3*(e.s2^2)); ## Normal distribution
  
  #b1=beta1*beta2*beta3
  b= beta1*beta2*beta3;
  
  LT2 <- L2T2.SRS(b, X.Ms, e.Ms);
  SNP.pwr <- MCVM.Power (LT2, NSL, n);
  
  return(SNP.pwr);
}


#### Find sample size of SNP test to reach a certain power under the MCM
####
#### X2 = X1*beta1+e2, 
#### X3 = X2*beta2+e3,
#### Y  = X3*beta3+e4,
#### e2, e3, e4 iid N(0,1)
####
#### RATIONALE
#### Rewrite the MCM as SIMPLE MCV MODEL
#### Y=X*b+e: X=X1, b=beta1*beta2*beta3, e=beta2*beta3*e2 + beta3*e3 + e4
#### X=Gr+u:  G=X1, r=1, u=0: Pr(X=G)=Pr(u=0)=1 
####
#### INPUTS
#### G.Ms: The first 4 moments of G
#### beta1, beta2, beta3: The model coefficients: 
####
#### RETURN
#### Sample size required
####
#### Created by Huaizhen Qin on July 23, 2014
##########################################################################

SNP.Geno.n.SRS<-function(G.Ms, beta1, beta2, beta3, NSL, pwr, n.range)
{	
  r = rep(1, length(beta1));#gamma in the simple model (Fig S1)
  X.Ms <- X.Moments(r, G.Ms = G.Ms, u.Ms = cbind(0*r, 0, 0, 0));
  
  
  #epsilon1
  e.s2 <- 1+(beta3^2)+(beta2*beta3)^2; 
  e.Ms <- cbind(0, e.s2, 0, 3*(e.s2^2)); ## Normal distribution
  
  #b1=beta1*beta2*beta3
  b= beta1*beta2*beta3;
  
  LT2 <- L2T2.SRS(b, X.Ms, e.Ms);
  SNP.n <- MCVM.n(LT2, NSL, pwr, n.range);
  
  return(SNP.n);
}


#### Compute power of RNA test under the MCM
####
#### X2 = X1*beta1+e2, 
#### X3 = X2*beta2+e3,
#### Y  = X3*beta3+e4,
#### e2, e3, e4 iid N(0,1)
####
#### RATIONALE
#### Rewrite the MCM as SIMPLE MCV MODEL
#### Y=X*b+e: X=X2, b=beta2*beta3, e=beta3*e3+e4
#### X=Gr+u:  G=X1, r=beta1, u=e2
####
#### INPUTS
#### G.Ms: The first 4 moments of G
#### beta1, beta2, beta3: The model coefficients: 
#### NSL: nominal significance level
#### n: sample size 
#### 
#### RETURN
#### Power
####
#### Created by Huaizhen Qin on July 24, 2014
##########################################################################

RNA.Expr.pwr.SRS <- function(G.Ms, beta1, beta2, beta3, NSL, n)
{
  #gamma in the simple model (Fig S1)
  r = beta1; 
  
  u.s2 <- 1;
  u.Ms <- cbind(0, u.s2, 0, 3*(u.s2)^2);
  X.Ms <- X.Moments(r, G.Ms, u.Ms);
  
  #epsilon2 = beta3*e3 + e4
  e.s2 <- 1+(beta3^2);
  e.Ms <- cbind(0, e.s2, 0, 3*(e.s2)^2);
  
  #b.2=beta2*beta3
  b= beta2*beta3; 
  
  LT2 <- L2T2.SRS(b, X.Ms, e.Ms);
  RNA.pwr <- MCVM.Power (LT2, NSL, n);
  
  return(RNA.pwr);
}


#### Find sample size of RNA test to reach a certain power under the MCM
####
#### X2 = X1*beta1+e2, 
#### X3 = X2*beta2+e3,
#### Y  = X3*beta3+e4,
#### e2, e3, e4 iid N(0,1)
####
#### RATIONALE
#### Rewrite the MCM as SIMPLE MCV MODEL
#### Y=X*b+e: X=X2, b=beta2*beta3, e=beta3*e3+e4
#### X=Gr+u:  G=X1, r=beta1, u=e2
####
#### INPUTS
#### G.Ms: The first 4 moments of G
#### beta1, beta2, beta3: The model coefficients: 
#### NSL: nominal significance level
#### pwr: power level to be achieved
#### n.range: range of integers from which to search the sample size
#### 
#### RETURN
#### Sample size required
####
#### Created by Huaizhen Qin on July 24, 2014
##########################################################################

RNA.Expr.n.SRS <- function(G.Ms, beta1, beta2, beta3, NSL, pwr, n.range)
{
  #gamma in the simple model (Fig S1)
  r = beta1; 
  
  u.s2 <- 1;
  u.Ms <- cbind(0, u.s2, 0, 3*(u.s2)^2);
  X.Ms <- X.Moments(r, G.Ms, u.Ms);
  
  #epsilon2 = beta3*e3 + e4
  e.s2 <- 1+(beta3^2);
  e.Ms <- cbind(0, e.s2, 0, 3*(e.s2)^2);
  
  #b.2=beta2*beta3
  b= beta2*beta3; 
  
  LT2 <- L2T2.SRS(b, X.Ms, e.Ms);
  RNA.n <- MCVM.n(LT2, NSL, pwr, n.range);
  
  return(RNA.n);
}


#### Compute power of PRT test under the MCM
####
#### X2 = X1*beta1+e2, 
#### X3 = X2*beta2+e3,
#### Y  = X3*beta3+e4,
#### e2, e3, e4 iid N(0,1)
####
#### RATIONALE
#### Rewrite the MCM as SIMPLE MCV MODEL
#### Y=X*b+e: X=X3, b=beta3, e=e4
#### X=Gr+u:  G=X1, r=beta1*beta2, u=e2*beta2+e3
####
#### INPUTS
#### G.Ms: The first 4 moments of G
#### beta1, beta2, beta3: The model coefficients: 
#### NSL: nominal significance level
#### n: sample size 
#### 
#### RETURN
#### Power
####
#### Created by Huaizhen Qin on July 25, 2014
##########################################################################

PRT.Expr.pwr.SRS <- function(G.Ms, beta1, beta2, beta3, NSL, n)
{
  #gamma in the simple model (Fig S1)
  r = beta1*beta2; 
  
  u.s2 <- 1 + (beta2^2);
  u.Ms <- cbind(0, u.s2, 0, 3*(u.s2)^2);
  X.Ms <- X.Moments(r, G.Ms, u.Ms);
  
  #epsilon3 = e4
  e.s2 <- 1;
  e.Ms <- cbind(0, e.s2, 0, 3*(e.s2)^2);
  
  #b3=beta3
  b = beta3;
  
  LT2 <- L2T2.SRS(b, X.Ms, e.Ms);
  PRT.pwr <- MCVM.Power (LT2, NSL, n);
  
  return(PRT.pwr);
}


#### Find sample size of PRT test to reach a certain power under the MCM
####
#### X2 = X1*beta1+e2, 
#### X3 = X2*beta2+e3,
#### Y  = X3*beta3+e4,
#### e2, e3, e4 iid N(0,1)
####
#### RATIONALE
#### Rewrite the MCM as SIMPLE MCV MODEL
#### Y=X*b+e: X=X3, b=beta3, e=e4
#### X=Gr+u:  G=X1, r=beta1*beta2, u=e2*beta2+e3
####
#### INPUTS
#### G.Ms: The first 4 moments of G
#### beta1, beta2, beta3: The model coefficients: 
#### NSL: nominal significance level
#### pwr: power level to be achieved
#### n.range: range of integers from which to search the sample size
#### 
#### RETURN
#### Sample size required
####
#### Created by Huaizhen Qin on July 25, 2014
##########################################################################

PRT.Expr.n.SRS <- function(G.Ms, beta1, beta2, beta3, NSL, pwr, n.range)
{
  #gamma in the simple model (Fig S1)
  r = beta1*beta2; 
  
  u.s2 <- 1 + (beta2^2);
  u.Ms <- cbind(0, u.s2, 0, 3*(u.s2)^2);
  X.Ms <- X.Moments(r, G.Ms, u.Ms);
  
  #epsilon3 = e4
  e.s2 <- 1;
  e.Ms <- cbind(0, e.s2, 0, 3*(e.s2)^2);
  
  #b3=beta3
  b = beta3;
  
  LT2 <- L2T2.SRS(b, X.Ms, e.Ms);
  PRT.n <- MCVM.n(LT2, NSL, pwr, n.range);
  
  return(PRT.n);
}


#### Searching for lower and upper quantiles of phenotype under
#### the MCM model
#### 
#### INPUTS 
#### beta1,...: vectors of model coefficients in the MCM model
#### MAF: minor allele frequency
#### alpha: EPS level for each tail
#### 
#### OUTPUTS
#### The lower and upper alpha-quantiles of Y
####
#### RATIONALE: Under the MCM model, Y = G*b+e, where G~B(2,p), 
#### b=beta1*beta2*beta3, and e~N(0, se^2). Therefore,
#### (Y/se)=G*(b/se)+e/se, where e/se ~N(0,1).
#### Find the quantiles of (Y/se) first, 
#### and then times them by se to derive quantiles of Y
####
#### Created by Huaizhen Qin on July 26, 2014
############################################################################

LU.MCM<-function(beta1, beta2, beta3, p, alpha)
{
  tau.L<-c();
  tau.U<-c();
  
  se=sqrt(1+beta2^2+(beta2*beta3)^2);
  b=beta1*beta2*beta3/se;
  for(i in c(1:length(beta1)))
  {
    tau.L[i]<-LowerQ(b=b[i], p, alpha);
    tau.U[i]<-UpperQ(b=b[i], p, alpha);
  }
  
  return (cbind(se*tau.L, se*tau.U));
}



#### Lemma 3: Compute conditonal moments of PRT expr.
#########################################################################

#### Compute f.X.x: the PDF of X. This function depends on speficif X variable
#### under the SMCVM 
####
#### Y=X*b+e, e~N(0, sd.e^2)
#### X=G*r+u, G~B(2, p), u~N(0, sd.u^2)
####
#### INPUTS
#### x: a specific point in the support of density
#### p:  MAF at the SNP
#### r: Effect size of G on X
#### sd.u: standard deviation of u
####
#### RETURN
#### The value of f.X(x)
####
#### Created by Huaizhen Qin on July 27, 2014
############################################################################

f.X.x<-function(x, p, r, sd.u)
{
  dx<-((1-p)^2*dnorm(x/sd.u)+2*p*(1-p)*dnorm((x-r)/sd.u)+p*p*dnorm((x-2*r)/sd.u))/sd.u;
  
  return(dx);
}

##fXx=f.X.x(x=0.5, p=0.25, r=2, sd.u=1);

#### Compute the integrand of the kth moment of lower truncation 
#### of X in the SMCVM 
#### Y=X*b+e, e~N(0, sd.e^2)
#### X=G*r+u, G~B(2, p), u~N(0, sd.u^2)
####
#### INPUTS
#### x: a specific point in the support of density
#### k: order of the moment
#### La: Lower alpha quantile of Y
#### b:  Effect size of X on Y
#### p:  MAF at the SNP
#### r: Effect size of G on X
#### sd.e: standard deviation of e
#### sd.u: standard deviation of u
#### alpha: truncation level for EPS (0~0.5)
####
#### RETURN
#### Value of the integrand at x 
####
#### Created by Huaizhen Qin on July 27, 2014
############################################################################

intgrand.mkx.0<-function(x, k, La, b,  p, r, sd.e, sd.u, alpha)
{
  int0<-(x^k)*f.X.x(x, p, r, sd.u)*pnorm((La-b*x)/sd.e, lower.tail = TRUE)/alpha;
  return(int0);
}

#intgrand.mkx.0(x=0.5, k=1, La=tau.LU[1,1], b=beta2[1]*beta3[1], p=0.25, r=beta1[2], sd.e=sqrt(1+beta3[2]^2), sd.u=1, alpha);

#### Compute the integrand of the kth moment of upper truncation 
#### of X in the SMCVM 
#### Y=X*b+e, e~N(0, sd.e^2)
#### X=G*r+u, G~B(2, p), u~N(0, sd.u^2)
####
#### INPUTS
#### x: a specific point in the support of density
#### k: order of the moment
#### Ua: upper alpha quantile of Y
#### b:  Effect size of X on Y
#### p:  MAF at the SNP
#### r: Effect size of G on X
#### sd.e: standard deviation of e
#### sd.u: standard deviation of u
#### alpha: truncation level for EPS (0~0.5)
####
#### RETURN
#### Value of the integrand at x 
####
#### Created by Huaizhen Qin on July 27, 2014
############################################################################

intgrand.mkx.1<-function(x, k, Ua, b, p, r, sd.e, sd.u, alpha)
{
  int1<-(x^k)*f.X.x(x, p, r, sd.u)*pnorm((Ua-b*x)/sd.e, lower.tail = FALSE)/alpha;
  return(int1);
}

#intgrand.mkx.1(x=0.5, k=1, Ua=tau.LU[1,2], b=beta2[1]*beta3[1], p=0.25, r=beta1[1], sd.e=sqrt(1+beta3[2]^2), sd.u=1, alpha);


#### Compute the integrand K(x) = exp(c*x)*f_X(x) under the SMCVM 
#### Y=X*b+e, e~N(0, sd.e^2)
#### X=G*r+u, G~B(2, p), u~N(0, sd.u^2)
####
#### INPUTS
#### x: a specific point in the support of density
#### p:  MAF at the SNP
#### r: Effect size of G on X
#### sd.u: standard deviation of u
####
#### RETURN
#### Value of the integrand at x 
####
#### Created by Huaizhen Qin on July 28, 2014
############################################################################

ecxf.X.x<-function(x, p, r, sd.u, c)
{
  dx<-exp(.5*c*c*sd.u^2)*((1-p)^2*dnorm(x/sd.u-c*sd.u)+2*p*(1-p)*dnorm((x-r)/sd.u-c*sd.u)*exp(r*c)+p*p*dnorm((x-2*r)/sd.u-c*sd.u)*exp(2*r*c))/sd.u;
  
  return(dx);
}

#### Compute the integrand of lower truncation mean of exp(cX) under the SMCVM 
#### Y=X*b+e, e~N(0, sd.e^2)
#### X=G*r+u, G~B(2, p), u~N(0, sd.u^2)
####
#### INPUTS
#### x: a specific point in the support of density
#### k: order of the moment
#### La: lower alpha quantile of Y
#### b:  Effect size of X on Y
#### p:  MAF at the SNP
#### r: Effect size of G on X
#### sd.e: standard deviation of e
#### sd.u: standard deviation of u
#### alpha: truncation level for EPS (0~0.5)
####
#### RETURN
#### Value of the integrand at x 
####
#### Created by Huaizhen Qin on July 28, 2014
############################################################################

intgrand.expx.0<-function(x, La, b,  p, r, sd.e, sd.u, alpha, c)
{
  int0<-ecxf.X.x(x, p, r, sd.u, c)*pnorm((La-b*x)/sd.e, lower.tail = TRUE)/alpha;
  return(int0);
}

#intgrand.expx.0(x=0, La, b,  p, r, sd.e, sd.u, alpha) 

#### Compute the integrand of upper truncation mean of exp(cX) under the SMCVM 
#### Y=X*b+e, e~N(0, sd.e^2)
#### X=G*r+u, G~B(2, p), u~N(0, sd.u^2)
####
#### INPUTS
#### x: a specific point in the support of density
#### k: order of the moment
#### Ua: upper alpha quantile of Y
#### b:  Effect size of X on Y
#### p:  MAF at the SNP
#### r: Effect size of G on X
#### sd.e: standard deviation of e
#### sd.u: standard deviation of u
#### alpha: truncation level for EPS (0~0.5)
####
#### RETURN
#### Value of the integrand at x 
####
#### Created by Huaizhen Qin on July 28, 2014
############################################################################

intgrand.expx.1<-function(x, Ua, b, p, r, sd.e, sd.u, alpha, c)
{
  int1<-ecxf.X.x(x, p, r, sd.u, c)*pnorm((Ua-b*x)/sd.e, lower.tail = FALSE)/alpha;
  return(int1);
}

#intgrand.expx.1(x=0, Ua, b,  p, r, sd.e, sd.u, alpha) 



#### Proposition 2: Noncentrality and fold change under EPS
############################################################################

#### Compute the noncentrality of t test for mediate effect under the SMCVM 
#### Y=X*b+e, e~N(0, sd.e^2)
#### X=G*r+u, u~N(0, sd.u^2), G~B(2, p)
####
#### INPUTS
#### La: Lower alpha quantile of Y
#### Ua: Upper alpha quantile of Y
#### b:  Effect size of X on Y
#### p:  MAF at the SNP
#### r: Effect size of G on X
#### sd.e: standard deviation of e
#### sd.u: standard deviation of u
#### alpha: truncation level for EPS (0~0.5)
####
#### RETURN
#### Noncentrality of t test for mediate effect
####
#### Created by Huaizhen Qin on July 29, 2014
############################################################################

delta2<-function(La, Ua, b, p, r, sd.e, sd.u, alpha)
{
  m1x.0<-integrate(intgrand.mkx.0, k=1, La, b, p, r, sd.e, sd.u, alpha, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.85, subdivisions = 2000L)[[1]];
  m2x.0<-integrate(intgrand.mkx.0, k=2, La, b, p, r, sd.e, sd.u, alpha, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.85, subdivisions = 2000L)[[1]];
  
  m1x.1<-integrate(intgrand.mkx.1, k=1, Ua, b, p, r, sd.e, sd.u, alpha, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.85, subdivisions = 2000L)[[1]];
  m2x.1<-integrate(intgrand.mkx.1, k=2, Ua, b, p, r, sd.e, sd.u, alpha, lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.85, subdivisions = 2000L)[[1]];
  
  dlt2<-(m1x.1-m1x.0)^2/(m2x.0-m1x.0^2+m2x.1-m1x.1^2);
  
  return(dlt2);
}



#### Compute the noncentrality for protein t test under EPS
####
#### INPUTS 
#### tau.LU: The lower and upper alpha-quantiles of Y
#### beta1,...: vectors of model coefficients in the MCM model
#### p: minor allele frequency
#### alpha: EPS level for each tail
#### 
#### OUTPUTS
#### Matrix = [1, noncetrality values]
####
#### Created by Huaizhen Qin on July 29, 2014
############################################################################

L2T2.PRT.EPS<-function(tau.LU, beta1, beta2, beta3, p, alpha)
{
  ## Connecting the MCM with the SMCVM (Lemma 3)
  ## Y=X*b+e: X=X3, b=beta3, e=e4
  ## X=Gr+u:  G=X1~B(2, p), r=beta1*beta2, u=beta2*e2+e3
  
  dlt2<-c();
  lb=length(beta1);
  for(i in c(1:lb))
  {
    La=tau.LU[i,1];
    Ua=tau.LU[i,2];
    b=beta3[i];
    r=beta1[i]*beta2[i];
    sd.e=1;
    sd.u=sqrt(1+beta2[i]^2);
    
    dlt2[i]<-delta2(La, Ua, b, p, r, sd.e, sd.u, alpha);
    
  }
  
  LT2<-cbind(1, dlt2);
  
  return(LT2);
}



#### Compute the noncentrality for RNA t test under EPS
####
#### INPUTS 
#### tau.LU: The lower and upper alpha-quantiles of Y
#### beta1,...: vectors of model coefficients in the MCM model
#### p: minor allele frequency
#### alpha: EPS level for each tail
#### 
#### OUTPUTS
#### Matrix = [1, noncetrality values]
####
#### Created by Huaizhen Qin on July 30, 2014
############################################################################

L2T2.RNA.EPS<-function(tau.LU, beta1, beta2, beta3, p, alpha)
{
  ## Connecting the MCM with the SMCVM (Lemma 2)
  ## Y=X*b+e: X=X2, b=beta2*beta3, e=beta3*e3+e4
  ## X=Gr+u:  G=X1~B(2, p), r=beta1, u=e2
  
  dlt2<-c();
  lb=length(beta1);
  for(i in c(1:lb))
  {
    La=tau.LU[i,1];
    Ua=tau.LU[i,2];
    b=beta2[i]*beta3[i];
    r=beta1[i];
    sd.e=sqrt(1+beta3[i]^2);
    sd.u=1;
    
    dlt2[i]<-delta2(La, Ua, b, p, r, sd.e, sd.u, alpha);
    
  }
  
  LT2<-cbind(1, dlt2);
  
  return(LT2);
}



#### Under the MCM, plotting powers of different tests under SRS and EPS
####
#### INPUTS
#### h2: Heritability of the causal SNP
#### PRT.pwr.EPS: power of t test on PRT effect under EPS
#### RNA.pwr.EPS: power of t test on RNA effect under EPS
#### SNP.pwr.EPS: power of t test on SNP effect under EPS
#### PRT.pwr.SRS: power of t test on PRT effect under SRS
#### RNA.pwr.SRS: power of t test on RNA effect under SRS
#### SNP.pwr.SRS: power of t test on SNP effect under SRS
#### xlimit: range(h2)*100. Note, h2[1]=0 stands for null model 
#### fig_id: figure ID
#### ppath: path for saving the figure
####
#### OUTPUT
#### A power vs. h2 figure for power comparison
####
#### Created by Huaizhen Qin on August 1, 2014
############################################################################

power.plot2<-function(h2, PRT.pwr.EPS, RNA.pwr.EPS, SNP.pwr.EPS, PRT.pwr.SRS, RNA.pwr.SRS, SNP.pwr.SRS, xlimit=range(h2[-1])*100, fig_id, ppath)
{
  tiff(filename = paste(ppath, "Fig", fig_id, ".tiff", sep=""), width = 450, height = 400, units = "px", pointsize = 12, bg = "white")
  
  op <- par(mfrow = c(1, 1), pty = "m",bty='n', pin=c(0.5, 0.3), mai=c(0.8, 0.8, 0, 0), omi=c(0.1, 0.1, 0.1, 0.1)); 
  
  plot(  h2*100, PRT.pwr.EPS*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 1, lty = 1, lwd = 2, col = "red", font=2, font.lab=2, axes=FALSE, ann=FALSE);
  points(h2*100, RNA.pwr.EPS*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 5, lty = 5, lwd = 2, col = "red");
  points(h2*100, SNP.pwr.EPS*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 3, lty = 3, lwd = 2, col = "red");
  
  
  points(h2*100, PRT.pwr.SRS*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 1, lty = 1, lwd = 1, col = "blue");
  points(h2*100, RNA.pwr.SRS*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 5, lty = 5, lwd = 1, col = "blue");
  points(h2*100, SNP.pwr.SRS*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 3, lty = 3, lwd = 1, col = "blue");
  
  x.lab=seq.int(from=0, to=120, by=10)/100;
  y.lab=seq.int(from=0, to=100, by=20);
  
  axis(1,las=1, at=x.lab, lab = sprintf("%4.1f%%", x.lab, 2), font.axis = 2);
  axis(2,las=3, at=y.lab, lab = sprintf("%d%%", y.lab), font.axis = 2);	
  
  title(xlab="Heritability",  font.lab=1, cex.lab = 1.2, srt=0);
  title(ylab="Power",  font.lab=1, cex.lab = 1.2, srt=0);
  
  leg.txt<-c("EPS: PRT", "EPS: RNA", "EPS: SNP", "SRS: PRT", "SRS: RNA", "SRS: SNP");
  legend("bottomright", leg.txt, lty = c(1, 5, 3, 1, 5, 3), lwd=c(2, 2, 2, 1, 1, 1), col=c("red", "red", "red", "blue", "blue", "blue"), bg="white");
  
  dev.off();
}



#### Under the MCM, plotting powers of different tests vs. sample size under SRS and EPS
####
#### INPUTS
#### sample.size: series of sample sizes
#### PRT.pwr.EPS.n: power of t test on PRT effect under EPS
#### RNA.pwr.EPS.n: power of t test on RNA effect under EPS
#### SNP.pwr.EPS.n: power of t test on SNP effect under EPS
#### PRT.pwr.SRS.n: power of t test on PRT effect under SRS
#### RNA.pwr.SRS.n: power of t test on RNA effect under SRS
#### SNP.pwr.SRS.n: power of t test on SNP effect under SRS
#### xlimit: limits in sample size 
#### fig_id: figure ID
#### ppath: path for saving the figure
####
#### OUTPUT
#### A power vs. sample size figure for power comparison
####
#### Created by Huaizhen Qin on August 2, 2014
############################################################################

plot.power.sample.size<-function(sample.size, PRT.pwr.EPS.n, RNA.pwr.EPS.n, SNP.pwr.EPS.n, PRT.pwr.SRS.n, RNA.pwr.SRS.n, SNP.pwr.SRS.n, xlimit=c(0, 1000), fig_id, ppath)
{
  tiff(filename = paste(ppath, "Fig", fig_id, ".tiff", sep=""), width = 450, height = 400, units = "px", pointsize = 12, bg = "white")
  
  op <- par(mfrow = c(1, 1), pty = "m",bty='n', pin=c(0.5, 0.3), mai=c(0.8, 0.8, 0, 0), omi=c(0.1, 0.1, 0.1, 0.1)); 
  
  plot(  sample.size, PRT.pwr.EPS.n*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 1, lty = 1, lwd = 2, col = "red", font=2, font.lab=2, axes=FALSE, ann=FALSE);
  points(sample.size, RNA.pwr.EPS.n*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 5, lty = 5, lwd = 2, col = "red");
  points(sample.size, SNP.pwr.EPS.n*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 3, lty = 3, lwd = 2, col = "red");
  
  points(sample.size, PRT.pwr.SRS.n*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 1, lty = 1, lwd = 1, col = "blue");
  points(sample.size, RNA.pwr.SRS.n*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 5, lty = 5, lwd = 1, col = "blue");
  points(sample.size, SNP.pwr.SRS.n*100, type="l", ylim=c(0, 100), xlim=xlimit, pch = 3, lty = 3, lwd = 1, col = "blue");
  
  x.lab=seq.int(from=0, to=10000, by=100);
  y.lab=seq.int(from=0, to=100, by=20);
  
  axis(1,las=1, at=x.lab, lab = sprintf("%d", x.lab, 2), font.axis = 2);
  axis(2,las=3, at=y.lab, lab = sprintf("%d%%", y.lab),  font.axis = 2);	
  
  title(xlab="Sample size",  font.lab=1, cex.lab = 1.2, srt=0);
  title(ylab="Power", font.lab=1, cex.lab = 1.2, srt=0);
  
  leg.txt<-c("EPS: PRT", "EPS: RNA", "EPS: SNP", "SRS: PRT", "SRS: RNA", "SRS: SNP");
  legend("bottomright", leg.txt, lty = c(1, 5, 3, 1, 5, 3), lwd=c(2, 2, 2, 1, 1, 1), col=c("red", "red", "red", "blue", "blue", "blue"), bg="white");
  
  dev.off();
}


#### Under the MCM, plotting sample sizes of different strategies to achieve a given
#### power level 
####
#### INPUTS
#### h2: Heritability of the causal SNP
#### PRT.n: sample size of PRT test
#### RNA.n: sample size of RNA test 
#### SNP.n: sample size of SNP test 
#### cur.col: color of the power curves
#### xlimit: range(h2)*100. Note, h2[1]=0 stands for null model 
#### ylimit: limits in sample size
#### fig_id: figure ID
#### ppath: path for saving the figure
####
#### OUTPUT
#### A sample size vs. h2 figure for sample size comparison
####
#### Created by Huaizhen Qin on August 3, 2014
#####################################################################################

sample.plot<-function(h2, PRT.n, RNA.n, SNP.n, cur.col, xlimit=c(0, 1), ylimit=c(0, 1000), fig_id, ppath)
{
  #ylimit=range(SNP.n[-1]);
  tiff(filename = paste(ppath, "Fig", fig_id, ".tiff", sep=""), width = 450, height = 400, units = "px", pointsize = 12, bg = "white")
  op <- par(mfrow = c(1, 1), pty = "m",bty='n', pin=c(0.5, 0.3), mai=c(0.8, 0.8, 0, 0), omi=c(0.1, 0.1, 0.1, 0.1)); 
  
  plot(  h2*100, PRT.n, type="l", ylim=ylimit, xlim=xlimit, pch = 1, lty = 1, lwd = 2, col = cur.col, font=2, font.lab=2, axes=FALSE, ann=FALSE);
  points(h2*100, RNA.n, type="l", ylim=ylimit, xlim=xlimit, pch = 2, lty = 2, lwd = 2, col = cur.col);
  points(h2*100, SNP.n, type="l", ylim=ylimit, xlim=xlimit, pch = 3, lty = 3, lwd = 2, col = cur.col);
  
  #x.lab=seq.int(from=0, to=100, by=5)/100;
  x.lab=seq.int(from=0, to=100, by=10)/100;
  y.lab=seq.int(from=0, to=25000, by=2000);
  
  axis(1,las=1, at=x.lab, lab = sprintf("%4.1f%%", x.lab, 2), font.axis = 2);
  axis(2,las=3, at=y.lab, lab = sprintf("%d", y.lab), font.axis = 2);	
  
  title(xlab="Heritability", font.lab=1, cex.lab = 1.2, srt=0);
  title(ylab="Sample size",  font.lab=1, cex.lab = 1.2, srt=0);
  
  leg.txt<-c("PRT", "RNA", "SNP");
  legend("topright", leg.txt, lty = c(1:3), lwd=c(2, 2, 2), col=c(cur.col, cur.col, cur.col), bg="white");
  
  dev.off();
}



#### Under the MCM, plotting sample sizes of different strategies to achieve a given
#### power level 
####
#### INPUTS
#### h2: Heritability of the causal SNP
#### SNP.n.SRS: sample size of correlation test on SNP effect under SRS
#### RNA.n.SRS: sample size of correlation test on RNA effect under SRS
#### PRT.n.SRS: sample size of correlation test on PRT effect under SRS
#### SNP.n.EPS: sample size of correlation test on SNP effect under EPS
#### RNA.n.EPS: sample size of correlation test on RNA effect under EPS
#### PRT.n.EPS: sample size of correlation test on PRT effect under EPS
#### xlimit: range(h2)*100. Note, h2[1]=0 stands for null model 
#### fig_id: figure ID
#### ppath: path for saving the figure
####
#### OUTPUT
#### A sample size vs. h2 figure for sample size comparison
####
#### Created by Huaizhen Qin on August 3, 2014
#####################################################################################

sample.plot2<-function(h2, SNP.n.SRS, RNA.n.SRS, PRT.n.SRS, SNP.n.EPS, RNA.n.EPS, PRT.n.EPS, xlimit, ylimit, fig_id, ppath)
{
  tiff(filename = paste(ppath, "Fig", fig_id, ".tiff", sep=""), width = 450, height = 400, units = "px", pointsize = 12, bg = "white")
  
  op <- par(mfrow = c(1, 1), pty = "m",bty='n', pin=c(0.5, 0.3), mai=c(0.8, 0.8, 0, 0), omi=c(0.1, 0.1, 0.1, 0.1)); 
  
  plot(  h2*100, SNP.n.SRS, type="l", ylim=ylimit, xlim=xlimit, pch = 1, lty = 3, lwd = 1, col = "blue", font=2, font.lab=2, axes=FALSE, ann=FALSE);
  points(h2*100, RNA.n.SRS, type="l", ylim=ylimit, xlim=xlimit, pch = 5, lty = 5, lwd = 1, col = "blue");
  points(h2*100, PRT.n.SRS, type="l", ylim=ylimit, xlim=xlimit, pch = 3, lty = 1, lwd = 1, col = "blue");
  
  points(h2*100, SNP.n.EPS, type="l", ylim=ylimit, xlim=xlimit, pch = 1, lty = 3, lwd = 2, col = "red");
  points(h2*100, RNA.n.EPS, type="l", ylim=ylimit, xlim=xlimit, pch = 5, lty = 5, lwd = 2, col = "red");
  points(h2*100, PRT.n.EPS, type="l", ylim=ylimit, xlim=xlimit, pch = 3, lty = 1, lwd = 2, col = "red");
  
  
  x.lab=seq.int(from=0, to=120, by=10)/100;
  y.lab=seq.int(from=0, to=25000, by=2000);
  
  axis(1,las=1, at=x.lab, lab = sprintf("%4.1f%%", x.lab, 2), font.axis = 2);
  axis(2,las=3, at=y.lab, lab = sprintf("%d", y.lab), font.axis = 2);	
  
  title(xlab="Heritability", font.lab=1, cex.lab = 1.2, srt=0);
  title(ylab="Sample size",  font.lab=1, cex.lab = 1.2, srt=0);
  
  leg.txt<-c("SRS: SNP", "SRS: RNA", "SRS: PRT", "EPS: SNP", "EPS: RNA", "EPS: PRT");
  legend("topright", leg.txt, lty = c(3, 5, 1, 3, 5, 1), lwd=c(1,1,1, 2, 2, 2), col=c("blue", "blue", "blue", "red", "red", "red"), bg="white");
  
  
  dev.off();
}

