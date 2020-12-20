source("../plot_pct_change.R")
library(dplyr)
library(tidyverse)
names <- ""

for (i in 0:11) {
  file <- paste0("Kidney_sc_perturbMet_neuron_subcluster_age_", i, ".RData")
  load(file)
  
  perturbed_mets<-result_ad_sc$perturbed_mets
  met_gene_pairs<-result_ad_sc$met_gene_pairs
  head(perturbed_mets[,c("mets","combined_pval","names")],n=15)
  names <- c(names, as.vector(perturbed_mets$names)[1:15])

  met2look<-as.vector(perturbed_mets[1:15,"mets"])
  for (j in seq_along(met2look)){
    sub_met_genes<-met_gene_pairs[met_gene_pairs$mets==met2look[j],]
    sub_met_genes$exprs_ref<-round(sub_met_genes$exprs_ref,digits = 1)
    sub_met_genes$exprs_treat<-round(sub_met_genes$exprs_treat,digits = 1)
    sub_met_genes$pct_ref<-(sub_met_genes$exprs_ref+1)/sum(sub_met_genes$exprs_ref+1)
    sub_met_genes$pct_treat<-(sub_met_genes$exprs_treat+1)/sum(sub_met_genes$exprs_treat+1)
    sub_met_genes$contribution<- sub_met_genes$pct_treat*log2(sub_met_genes$pct_treat/sub_met_genes$pct_ref)
    sub_met_genes<-arrange(sub_met_genes,desc(abs(contribution)))
    # print(sub_met_genes)
    p <- plot_pct_change(met_gene_pairs,met2look[j],c(0,1),6,100)
    png(file=paste0("Plot/cluster", i, "/age_", met2look[j],".png"))
    print(p)
    dev.off()
  }
}
names <- as.data.frame(names)
names <- na.omit(names)
result <- names %>% group_by(names) %>% count()
result <- result %>% arrange(desc(n))
result <- result[1:31,]

view(result)
