## August 2nd, 2018, 01:32 PM

rm(list=ls(all="TRUE"));

path<-"C:/MCM/";
setwd(path);

source("Subfunctions.R"); ##MCM_



#### Fig 3: Power comparison at various heritability levels  
################################################################################
n=200;
NSL = 2.5e-6;	#Nominal Significance Level
alpha = 0.2;

npts = 500;
rpts =(c(0:npts)/(npts*100));
#hi=368;
r2.1 = 15*rpts; ##0.5%~5% 	#refs. SNP to RNA
r2.2 = 50*rpts; ##30~50%	#refs. (RNA to PRT) de Sousa Abreu et al. 2009 


#### Let X1 = G~B(2, p) with p = 0.25 say, then the first 4 moments of G is:
p = 0.25; ##MAF at the causal SNP
G.Ms <- c(2*p, 2*p+2*p*p, 2*p+6*p*p, 2*p+14*p*p);

#### var(G)
var.x1 = 2*p*(1-p);
beta1 = sqrt(r2.1/( (1-r2.1)*var.x1 ) );

#### var(X2) = var(e2+X1*beta1)
var.x2 = 1+(beta1*beta1)*var.x1;
beta2 = sqrt(r2.2/( (1-r2.2)*var.x2 ) );

#### var(X3)=var(e3+X2*beta2)
var.x3 = 1+(beta2*beta2)*var.x2;

c=log(2)
nfold = 5;
beta.r2.3<-effect.sizes.PRT(nfold, alpha, p, beta1, beta2, var.x3);
beta3<-beta.r2.3[,1];
r2.3<-beta.r2.3[,2];
h2<-r2.1*r2.2*r2.3; #### SNP heritability

#### Computing powers
################################################################################

#### Under SRS design

SNP.pwr.SRS <- SNP.Geno.pwr.SRS(G.Ms, beta1, beta2, beta3, NSL, n);
RNA.pwr.SRS <- RNA.Expr.pwr.SRS(G.Ms, beta1, beta2, beta3, NSL, n);
PRT.pwr.SRS <- PRT.Expr.pwr.SRS(G.Ms, beta1, beta2, beta3, NSL, n);

#### Under EPS design

## Lemma 2: For given vectors of effect sizes, searching for lower and upper 
## quantiles of Y. These quantiles work for all tests under the EPS 

tau.LU<-LU.MCM(beta1, beta2, beta3, p, alpha);

## Lemma 3-4: Computing noncentrality parameters of the three t tests

LT2.SNP.EPS <- L2T2.SNP.EPS(tau.LU, beta1, beta2, beta3, p, alpha);
LT2.RNA.EPS <- L2T2.RNA.EPS(tau.LU, beta1, beta2, beta3, p, alpha);
LT2.PRT.EPS <- L2T2.PRT.EPS(tau.LU, beta1, beta2, beta3, p, alpha);

## Proposition 2: Compute rejection rates 
SNP.pwr.EPS<-MCVM.Power(LT2.SNP.EPS, NSL, n/2);
RNA.pwr.EPS<-MCVM.Power(LT2.RNA.EPS, NSL, n/2);
PRT.pwr.EPS<-MCVM.Power(LT2.PRT.EPS, NSL, n/2);

power.plot2(h2, PRT.pwr.EPS, RNA.pwr.EPS, SNP.pwr.EPS, PRT.pwr.SRS, RNA.pwr.SRS, SNP.pwr.SRS, xlimit=c(0.1, 1.0), fig_id=3, path);

PowerComp<-cbind(h2, beta1, beta2, beta3, PRT.pwr.EPS, RNA.pwr.EPS, SNP.pwr.EPS, PRT.pwr.SRS, RNA.pwr.SRS, SNP.pwr.SRS); 
colnames(PowerComp)<-c("Heritibility",  "beta1", "beta2", "beta3", "PRT.pwr.EPS", "RNA.pwr.EPS", "SNP.pwr.EPS", "PRT.pwr.SRS", "RNA.pwr.SRS", "SNP.pwr.SRS");
write.table(PowerComp, paste(path, "Fig3.txt"), append = FALSE, sep=" ", quote = FALSE, col.names = TRUE, row.names = FALSE)

