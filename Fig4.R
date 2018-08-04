## August 2rd, 2018, 09:37 PM

rm(list=ls(all="TRUE"));

path<-"C:/MCM/";
setwd(path);

source("Subfunctions.R");


#### Fig 4: Power comparison at various sample sizes.  
################################################################################
NSL = 2.5e-6;	#Nominal Significance Level
alpha = 0.2;

npts = 500;
rpts =(c(0:npts)/(npts*100));
hi=368; #h2 = 1%
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

#### Computing powers for various sample sizes
################################################################################

sample.size=c(0:18296); 

#### Under SRS design

SNP.pwr.SRS.n <- SNP.Geno.pwr.SRS(G.Ms, beta1[hi], beta2[hi], beta3[hi], NSL, sample.size);
RNA.pwr.SRS.n <- RNA.Expr.pwr.SRS(G.Ms, beta1[hi], beta2[hi], beta3[hi], NSL, sample.size);
PRT.pwr.SRS.n <- PRT.Expr.pwr.SRS(G.Ms, beta1[hi], beta2[hi], beta3[hi], NSL, sample.size);

#### Under EPS design

## Lemma 2: For given vectors of effect sizes, searching for lower and upper 
## quantiles of Y. These quantiles work for all tests under the EPS 

tau.LU<-LU.MCM(beta1, beta2, beta3, p, alpha);

## Lemma 3-4: Computing noncentrality parameters of the three t tests

LT2.SNP.EPS <- L2T2.SNP.EPS(tau.LU, beta1, beta2, beta3, p, alpha);
LT2.RNA.EPS <- L2T2.RNA.EPS(tau.LU, beta1, beta2, beta3, p, alpha);
LT2.PRT.EPS <- L2T2.PRT.EPS(tau.LU, beta1, beta2, beta3, p, alpha);

halfss=sample.size/2;
SNP.pwr.EPS.n<-MCVM.Power(matrix(LT2.SNP.EPS[hi,], 1, 2), NSL, n=halfss);
RNA.pwr.EPS.n<-MCVM.Power(matrix(LT2.RNA.EPS[hi,], 1, 2), NSL, n=halfss);
PRT.pwr.EPS.n<-MCVM.Power(matrix(LT2.PRT.EPS[hi,], 1, 2), NSL, n=halfss);

plot.power.sample.size(sample.size, PRT.pwr.EPS.n, RNA.pwr.EPS.n, SNP.pwr.EPS.n, PRT.pwr.SRS.n, RNA.pwr.SRS.n, SNP.pwr.SRS.n, xlimit=c(50, 600), fig_id=4, path)

Power.Sample.Comp<-cbind(sample.size, PRT.pwr.EPS.n, RNA.pwr.EPS.n, SNP.pwr.EPS.n, PRT.pwr.SRS.n, RNA.pwr.SRS.n, SNP.pwr.SRS.n); 
colnames(Power.Sample.Comp)<-c("Sample.size", "PRT.pwr.EPS.n", "RNA.pwr.EPS.n", "SNP.pwr.EPS.n", "PRT.pwr.SRS.n", "RNA.pwr.SRS.n", "SNP.pwr.SRS.n");
write.table(Power.Sample.Comp, paste(path, "Fig4.txt"), append = FALSE, sep=" ", quote = FALSE, col.names = TRUE, row.names = FALSE)

