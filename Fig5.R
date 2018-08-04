## August 2nd, 2018, 01:32 PM

rm(list=ls(all="TRUE"));

path<-"C:/MCM/";
setwd(path);

source("Subfunctions.R");


#### Fig 5: Sample sizes required to achieve 80% power 
################################################################################
pwr = 0.8;
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


#### Under SRS design

SNP.n.SRS<-SNP.Geno.n.SRS(G.Ms, beta1, beta2, beta3, NSL, pwr, n.range=c(3, 5e15));
RNA.n.SRS<-RNA.Expr.n.SRS(G.Ms, beta1, beta2, beta3, NSL, pwr, n.range=c(3, 5e12));
PRT.n.SRS<-PRT.Expr.n.SRS(G.Ms, beta1, beta2, beta3, NSL, pwr, n.range=c(3, 5e12));

#### Under the EPS design

## Lemma 2: For given vectors of effect sizes, searching for lower and upper 
## quantiles of Y. These quantiles work for all tests under the EPS 

tau.LU<-LU.MCM(beta1, beta2, beta3, p, alpha);

## Lemma 3-4: Computing noncentrality parameters of the three t tests

LT2.SNP.EPS <- L2T2.SNP.EPS(tau.LU, beta1, beta2, beta3, p, alpha);
LT2.RNA.EPS <- L2T2.RNA.EPS(tau.LU, beta1, beta2, beta3, p, alpha);
LT2.PRT.EPS <- L2T2.PRT.EPS(tau.LU, beta1, beta2, beta3, p, alpha);

SNP.n.EPS<-2*MCVM.n(LT2.SNP.EPS, NSL, pwr, n.range=c(3, 5e12));
RNA.n.EPS<-2*MCVM.n(LT2.RNA.EPS, NSL, pwr, n.range=c(3, 5e12));
PRT.n.EPS<-2*MCVM.n(LT2.PRT.EPS, NSL, pwr, n.range=c(3, 5e12));

sample.plot2(h2, SNP.n.SRS, RNA.n.SRS, PRT.n.SRS, SNP.n.EPS, RNA.n.EPS, PRT.n.EPS, xlimit=c(0.1, 1.0), ylimit=c(0, 12000), fig_id=5, path);

SampleComp<-cbind(h2, beta1, beta2, beta3, SNP.n.SRS, RNA.n.SRS, PRT.n.SRS, SNP.n.EPS, RNA.n.SRS, PRT.n.EPS); 
colnames(SampleComp)<-c("Heritibility",  "beta1", "beta2", "beta3", "SNP.n.SRS", "RNA.n.SRS", "PRT.n.SRS", "SNP.n.EPS", "RNA.n.SRS", "PRT.n.EPS");
write.table(SampleComp, paste(path, "Fig5.txt"), append = FALSE, sep=" ", quote = FALSE, col.names = TRUE, row.names = FALSE )

