library(RMySQL)
library(ggplot2)

vplayout <- function(x, y) {
  viewport(layout.pos.row = x, layout.pos.col = y)
}

connect <- function(db) {
  dbConnect(dbDriver("MySQL"),username="team29",password="29@test",host="ia64b",dbname=db)
}

get_consequences <- function(con) {
  dbGetQuery(con, "SELECT consequence, count(*) FROM data GROUP BY consequence")
}

get_relative_genotypes <- function(con) {
  dbGetQuery(con, "SELECT relative_genotype, count(*) FROM data GROUP BY relative_genotype")
}

get_homminorgenes <- function(con) {
  dbGetQuery(con, "SELECT DISTINCT d.ensembl_acc FROM data d WHERE d.consequence IN ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST') AND d.relative_genotype IN ('minorminor','altalt') AND (d.1kg_major IS NULL OR d.1kg_maf < 0.01)")
}

get_all_homlofgenes <- function(con) {
  dbGetQuery(con, "select ensembl_acc, min(1kg_maf) as maf from data where consequence in ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST') and relative_genotype IN ('minorminor','altalt') group by ensembl_acc")
}

get_all_lofsnps <- function(con) {
  dbGetQuery(con, "select uid, 1kg_maf as maf from data where consequence in ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST')")
}

get_all_lofsnps_pilot1 <- function(con) {
  dbGetQuery(con, "select uid, 1kg_maf as maf from data where 1kg_major is not null and consequence in ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST')")
}

get_all_lofsnps_nonpilot1 <- function(con) {
  dbGetQuery(con, "select uid, 1kg_maf as maf from data where 1kg_major is null and consequence in ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST')")
}


get_all_homminorlofsnps <- function(con) {
  dbGetQuery(con, "select uid, 1kg_maf as maf from data where consequence in ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST') and relative_genotype IN ('minorminor','altalt')")
}

get_compoundhetgenes <- function(con) {
  dbGetQuery(con, "SELECT ensembl_acc, maf FROM (SELECT d.ensembl_acc, min(1kg_maf) AS maf, count(*) AS c FROM data d WHERE d.relative_genotype IN ('majorminor','refalt') AND d.consequence IN ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST') GROUP BY d.ensembl_acc) foo WHERE c > 1")
}

get_lofpergene <- function(con) {
  dbGetQuery(con, "SELECT ensembl_acc, uid FROM data WHERE consequence IN ('ESSENTIAL_SPLICE_SITE','STOP_GAINED','STOP_LOST')")
}


x <- seq(0,0.5, 0.01)
colours <- rainbow(8)

con <- connect("exome_common")
all_genes <- dbGetQuery(con, "SELECT * FROM genes;")
dbDisconnect(con)

individuals <- c('pg1','pg3','pg4','pg5','pg6','pg7','pg8','pg9')
merged_consequences <- data.frame()
merged_relative_genotypes <- data.frame()
lof_snps <- data.frame()
lof_snps_pilot1 <- data.frame()
lof_snps_nonpilot1 <- data.frame()
all_lofpergene <- data.frame()

for ( individual in individuals ) {
  if ( individual == 'pg1' ) {
    con <- connect(paste('exome', individual, 'b',sep='_'))
  } else {
    con <- connect(paste('exome', individual, sep='_'))
  }
  assign(paste(individual, 'consequences', sep='_'), get_consequences(con))
  assign(paste(individual, 'relative_genotypes', sep='_'), get_relative_genotypes(con))
  assign(paste(individual, 'homminorgenes', sep='_'), get_homminorgenes(con))
  assign(paste(individual, 'compoundhetgenes', sep='_'), get_compoundhetgenes(con))
  assign(paste(individual, 'lofgenes', sep='_'), get_all_homlofgenes(con))
  assign(paste(individual, 'lofsnps', sep='_'), get_all_lofsnps(con))
  assign(paste(individual, 'lofsnps_pilot1', sep='_'), get_all_lofsnps_pilot1(con))
  assign(paste(individual, 'lofsnps_nonpilot1', sep='_'), get_all_lofsnps_nonpilot1(con))
  assign(paste(individual, 'lofsnps_minor', sep='_'), get_all_homminorlofsnps(con))
  assign(paste(individual, 'lofsnps_y', sep='_'), as.matrix(sapply(x, function(x) { nrow(subset(get(paste(individual, 'lofsnps', sep='_')), get(paste(individual, 'lofsnps', sep='_'))$maf <= x)) } )))
  assign(paste(individual, 'lofsnps_pilot1_y', sep='_'), as.matrix(sapply(x, function(x) { nrow(subset(get(paste(individual, 'lofsnps_pilot1', sep='_')), get(paste(individual, 'lofsnps_pilot1', sep='_'))$maf <= x)) } )))
  assign(paste(individual, 'lofsnps_nonpilot1_y', sep='_'), as.matrix(sapply(x, function(x) { nrow(subset(get(paste(individual, 'lofsnps_nonpilot1', sep='_')), get(paste(individual, 'lofsnps_nonpilot1', sep='_'))$maf <= x)) } )))
  assign(paste(individual, 'lofsnps_minor_y', sep='_'), as.matrix(sapply(x, function(x) { nrow(subset(get(paste(individual, 'lofsnps_minor', sep='_')), get(paste(individual, 'lofsnps_minor', sep='_'))$maf <= x)) } )))
  assign(paste(individual, 'lofgenes_y', sep='_'), as.matrix(sapply(x, function(x) { nrow(subset(get(paste(individual, 'lofgenes', sep='_')), get(paste(individual, 'lofgenes', sep='_'))$maf <= x)) } )))
  if ( nrow(get(paste(individual, 'compoundhetgenes', sep='_'))) == 0 ) {
    assign(paste(individual, 'comphet_y', sep='_'), as.matrix(rep(0,51)))
  } else {
    assign(paste(individual, 'comphet_y', sep='_'), as.matrix(sapply(x, function(x) { nrow(subset(get(paste(individual, 'compoundhetgenes', sep='_')), get(paste(individuals, 'compoundhetgenes', sep='_'))$maf <= x)) } )))
  }
  assign(paste(individual, 'lofpergene', sep='_'), get_lofpergene(con))

  assign(paste('merged_consequences', individual, sep='_'), cbind(get(paste(individual, 'consequences', sep='_')), individual))
  merged_consequences <- rbind(merged_consequences, get(paste('merged_consequences', individual, sep='_')))

  assign(paste('merged_relative_genotypes', individual, sep='_'), cbind(get(paste(individual, 'relative_genotypes', sep='_')), individual))
  merged_relative_genotypes <- rbind(merged_relative_genotypes, get(paste('merged_relative_genotypes', individual, sep='_')))

  assign(paste('lof_snps', individual, sep='_'), data.frame(x, get(paste(individual, 'lofsnps_y', sep='_')), individual))
  lof_snps <- rbind(lof_snps, get(paste('lof_snps', individual, sep='_')))

  assign(paste('lof_snps_pilot1', individual, sep='_'), data.frame(x, get(paste(individual, 'lofsnps_pilot1_y', sep='_')), individual))
  lof_snps_pilot1 <- rbind(lof_snps_pilot1, get(paste('lof_snps_pilot1', individual, sep='_')))

  assign(paste('lof_snps_nonpilot1', individual, sep='_'), data.frame(x, get(paste(individual, 'lofsnps_nonpilot1_y', sep='_')), individual))
  lof_snps_nonpilot1 <- rbind(lof_snps_nonpilot1, get(paste('lof_snps_nonpilot1', individual, sep='_')))

  assign(paste(individual, 'lofpergene', sep='_'), data.frame(get(paste(individual, 'lofpergene', sep='_')), individual))
  all_lofpergene <- rbind(all_lofpergene, get(paste(individual, 'lofpergene', sep='_')))

  dbDisconnect(con)
}

names(merged_consequences) <- c('consequence','value','individual')
names(merged_relative_genotypes) <- c('relative_genotype','value','individual')
names(lof_snps) <- c('maf','count','individual')
names(lof_snps_pilot1) <- c('maf','count','individual')
names(lof_snps_nonpilot1) <- c('maf','count','individual')

a <- as.data.frame(table(all_lofpergene$ensembl_acc, all_lofpergene$individual))
names(a) <- c('ensembl_acc','individual','count')
b <- subset(a, a$count > 0)
c <- data.frame(table(b$ensembl_acc))
names(c) <- c('ensembl_acc','count')
nr_individuals_per_lof_gene <- data.frame(table(c$count))
names(nr_individuals_per_lof_gene) <- c('nr_individuals','count')
