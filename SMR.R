# Calculation of F values
dt$se <- paste(rep(" ", 20), collapse = " ")
dt$se=sqrt(((dt$beta)^2)/qchisq(dt$p,1,lower.tail=F))
write.table(dt,"gene_name.csv",row.names=FALSE,col.names=TRUE,sep=",")

# SMR analyses between all drug target gene expression and and venous thromboembolism outcomes
# Exposure: all antipsychotic drug target genes blood or brain tissue eqtl;
# Outcome: VTE, DVT and PE
## Packages required
library(TwoSampleMR)
## Reading exposure factors
dat1 <- as.data.frame(read.xlsx("all genes.xlsx",1))
dat2 <- format_data(dat1, type="exposure")
dat2_clumped <- clump_data(dat2, clump_r2 = 0.3, pop="EUR")
dat2_clumped$se.exposure <- sqrt(((dat2_clumped$beta.exposure)^2)/qchisq(dat2_clumped$pval.exposure,1,lower.tail=F))
## Getting the outcome factors
### VTE
outcome_dat <-extract_outcome_data(snps=dat2_clumped$SNP, outcomes="finn-b-I9_VTE")
### DVT
outcome_dat <-extract_outcome_data(snps=dat2_clumped$SNP, outcomes="finn-b-I9_PHLETHROMBDVTLOW")
### PE
outcome_dat <-extract_outcome_data(snps=dat2_clumped$SNP, outcomes="finn-b-I9_PULMEMB")
## Connecting the exposure factors and the outcome factors
dat <- harmonise_data(dat2_clumped,outcome_dat)
## Viewing the list of models for the MR
mr_method_list()
## Getting the results
method_list=c("mr_egger_regression", "mr_ivw")
res <- mr(dat)
ora<-generate_odds_ratios(res)
## Heterogeneity test
het <- mr_heterogeneity(dat)
Het
## Multipleiotropy test
pleio <- mr_pleiotropy_test(dat)
pleio
## If there is strong heterogeneity between IV (Q _ pval is much less than 0.05), the random effect model is directly used to estimate the MR effect size
mr(dat,method_list=c('mr_ivw_mre'))
## Exporting result table
write.table(res,"MR_result.csv",row.names=FALSE,col.names=TRUE,sep=",")
## Differences between individual models and individual SNP estimates
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
## Forest plots of the overall assessment
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)

#SMR analyses between each drug target gene expression and and venous thromboembolism outcomes
# Exposure: each antipsychotic drug target gene blood or brain tissue eqtl;
# Outcome: VTE, DVT and PE
## Packages required
library(TwoSampleMR)
## Reading exposure factors
dat1 <- as.data.frame(read.xlsx("gene_name.xlsx",1))
dat2 <- format_data(dat1, type="exposure")
dat2_clumped <- clump_data(dat2, clump_r2 = 0.3, pop="EUR")
dat2_clumped$se.exposure <- sqrt(((dat2_clumped$beta.exposure)^2)/qchisq(dat2_clumped$pval.exposure,1,lower.tail=F))
## Getting the outcome factors
### VTE
outcome_dat <-extract_outcome_data(snps=dat2_clumped$SNP, outcomes="finn-b-I9_VTE")
### DVT
outcome_dat <-extract_outcome_data(snps=dat2_clumped$SNP, outcomes="finn-b-I9_PHLETHROMBDVTLOW")
### PE
outcome_dat <-extract_outcome_data(snps=dat2_clumped$SNP, outcomes="finn-b-I9_PULMEMB")
## Connecting the exposure factors and the outcome factors
dat <- harmonise_data(dat2_clumped,outcome_dat)
## Viewing the list of models for the MR
mr_method_list()
## Getting the results
method_list=c("mr_egger_regression", "mr_ivw")
res <- mr(dat)
ora<-generate_odds_ratios(res)
## Heterogeneity test
het <- mr_heterogeneity(dat)
Het
## Multipleiotropy test
pleio <- mr_pleiotropy_test(dat)
pleio
## If there is strong heterogeneity between IV (Q _ pval is much less than 0.05), the random effect model is directly used to estimate the MR effect size
mr(dat,method_list=c('mr_ivw_mre'))
## Exporting result table
write.table(res,"MR_result.csv",row.names=FALSE,col.names=TRUE,sep=",")
## Differences between individual models and individual SNP estimates
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
## Forest plots of the overall assessment
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)