all_genes$uniqueness <- as.factor(all_genes$uniqueness)
omim_genes <- all_genes[!is.na(all_genes$omim),]
lof_gene_accessions <- data.frame(unique(all_lofpergene$ensembl_acc))
names(lof_gene_accessions) <- c('ensembl_acc')
lof_genes <- merge(all_genes, lof_gene_accessions, by='ensembl_acc')
lof_omim_gene_accessions <- data.frame(intersect(omim_genes$ensembl_acc, lof_genes$ensembl_acc))
names(lof_omim_gene_accessions) <- c('ensembl_acc')

omim_genes$group <- as.factor(c('OMIM'))
lof_genes$group <- as.factor(c('LOF'))
input_genes <- rbind(omim_genes, lof_genes)

model <- glm(group ~ haploinsufficiency+uniqueness+chimp_kaks+macaque_kaks+mouse_kaks+coding_gerp+promoter_gerp+paralog_num+paralog_dist+n_exon+size+spliced_size+cds_size+utr_length+domain_num+embryo_expr+tissue_expr_spec+ppi_dgr+ppi_clcf+ppi_ctty+ppi_dist2HI+ppi_lls2HI+ppi_dist2cancer+ppi_lls2cancer+cancer+yeast_hetgr, data=input_genes, family=binomial)

model_b <- glm(group ~ haploinsufficiency*macaque_kaks*coding_gerp*promoter_gerp*paralog_num*n_exon*size*tissue_expr_spec*ppi_dgr*ppi_dist2HI*cancer*yeast_hetgr, data=input_genes, family=binomial)

summary(model)

input_genes_not_null <- na.omit(input_genes[,c(6,11:19,23,25:34)])
pca <- prcomp(input_genes_not_null)