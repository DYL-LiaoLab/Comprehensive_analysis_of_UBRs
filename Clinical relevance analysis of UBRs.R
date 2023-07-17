# Clinical relevance analysis of UBRs in pan-cancer.
# Long deyu
# 2022.09.16

library(survival)
library(survminer)
library(dplyr)

# https://ucsc-xena.gitbook.io/project/tutorials/tutorial-viewing-your-own-data#part-a
# alive=0; dead=1
tcga_survival<-read.csv("TCGA_survival_data.txt",sep = "\t")
rownames(tcga_survival)<-tcga_survival$sample
# tcga_survival<-tcga_survival[,-1]

cancer_name<-read.csv("cancer_name.csv",check.names = FALSE)

phenotype<-read.table(gzfile("TcgaTargetGTEX_phenotype.txt"),
                      sep = "\t")
colnames(phenotype)<-phenotype[1,]
phenotype<-phenotype[-1,]
colnames(phenotype)<-c("sample","detailed_category",        
                       "primary disease or tissue","primary_site",            
                       "sample_type","gender","study")

table(phenotype$primary_site)

all_tpm<-read.table(gzfile("TcgaTargetGtex_rsem_gene_tpm.gz"),
                    sep = "\t")

all_tpm[1:3,1:4]
colnames(all_tpm)<-all_tpm[1,]
all_tpm<-all_tpm[-1,]
all_tpm[1:3,1:4]

head(all_tpm$sample)
head(gsub("\\..*","",all_tpm$sample))
rownames(all_tpm)<-gsub("\\..*","",all_tpm$sample)
all_tpm<-all_tpm[,-1]
all_tpm[1:3,1:4]

tcga_all_tpm<-all_tpm[,which(grepl("TCGA",colnames(all_tpm))==TRUE)]
tcga_all_tpm[1:3,1:4]

tcga_phen<-subset(phenotype,phenotype$study=="TCGA")
cancer_name1<-as.data.frame(t(cancer_name))
cancer_name1$detailed_category<-rownames(cancer_name1)
cancer_name1[12,2]<-"Pheochromocytoma & Paraganglioma"
cancer_name1[25,2]<-"Diffuse Large B-Cell Lymphoma"
tcga_phen1<-inner_join(tcga_phen,cancer_name1,by = "detailed_category")
colnames(tcga_phen1)[8]<-"name"

pwer_anno<-read.csv("human_ubiquitination_regulators.csv")
colnames(pwer_anno)
tcga_all_tpm1<-tcga_all_tpm
tcga_all_tpm1$ensembl_id<-rownames(tcga_all_tpm1)
pwer_tcga_tpm_all<-inner_join(pwer_anno,tcga_all_tpm1,by="ensembl_id")

pwer_sur_summary<-matrix(data = NA,nrow=852,ncol=33) 
pwer_sur_summary<-as.data.frame(pwer_sur_summary)
colnames(pwer_sur_summary)<-cancer_name1$V1
rownames(pwer_sur_summary)<-pwer_tcga_tpm_all$gene_symbol


for(i in 1:length(cancer_name1$detailed_category)){
  # i=17
  print(cancer_name1$detailed_category[i])
  temp_phen<-subset(tcga_phen1,tcga_phen1$detailed_category==cancer_name1$detailed_category[i])
  temp_phen<-subset(temp_phen,
                    temp_phen$sample_type!="Solid Tissue Normal")
  temp_tpm<-pwer_tcga_tpm_all[,temp_phen$sample]
  rownames(temp_tpm)<-pwer_tcga_tpm_all$gene_symbol
  temp_tpm<-temp_tpm[sort(rownames(temp_tpm)),] 
  print(i)
  temp_survival<-tcga_survival[temp_phen$sample,]
  temp_survival_phen<-inner_join(temp_survival,temp_phen,by="sample")
  write.csv(temp_survival_phen,paste(cancer_name1$V1[i],"survival_pheno.csv"),quote = FALSE,row.names = FALSE)
  # print(i)
  for(j in 1:852){ 
    # j=1
    print(rownames(temp_tpm)[j])
    temp_gene_expr<-as.data.frame(t(temp_tpm[j,]))
    temp_gene_expr$sample<-rownames(temp_gene_expr)
    temp_cut<-median(as.numeric(temp_gene_expr[,1])) 
    temp_gene_survival_phen<-inner_join(temp_gene_expr,temp_survival_phen,by="sample")
    temp_gene_survival_phen$state<-ifelse(as.numeric(temp_gene_survival_phen[,1])>=temp_cut,"High","ALow")
    if(length(unique(temp_gene_survival_phen$state))==1){next}
    # fit<-survfit(Surv(OS.time, OS) ~ state,
    #              data = temp_gene_survival_phen)
    fit<-coxph(Surv(OS.time, OS) ~ state,
               data = temp_gene_survival_phen)
    print(fit)
    summary(fit)
    covariates<-c("state")
    univ_formulas <- sapply(covariates,
                            function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
    univ_models <- lapply( univ_formulas, function(x){coxph(x, data = temp_gene_survival_phen)})
    univ_results <- lapply(univ_models,
                           function(x){ 
                             x <- summary(x)
                             p.value<-signif(x$wald["pvalue"], digits=2)
                             wald.test<-signif(x$wald["test"], digits=2)
                             beta<-signif(x$coef[1], digits=2);#coeficient beta
                             HR <-signif(x$coef[2], digits=2);#exp(beta)
                             HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                             HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                             HR <- paste0(HR, " (", 
                                          HR.confint.lower, "-", HR.confint.upper, ")")
                             res<-c(beta, HR, wald.test, p.value)
                             names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                           "p.value")
                             return(res)
                             #return(exp(cbind(coef(x),confint(x))))
                           })
    res <- t(as.data.frame(univ_results, check.names = FALSE))
    as.data.frame(res)
    temp_res<-as.data.frame(res)
    temp_res$cancer_type<-cancer_name1$V1[i]
    temp_res$gene<-rownames(temp_tpm)[j]
    temp_res$sample_count<-length(temp_phen$sample)
    print(j)
    if(i==1 & j==1){
      total_res<-temp_res
    }else{
      total_res<-rbind(total_res,temp_res)
    }
    
    
  }
  
}


write.csv(total_res,"Survival risk matrix_COX.csv")
