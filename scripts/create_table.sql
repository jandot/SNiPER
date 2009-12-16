CREATE TABLE data (
  uid VARCHAR(12) UNIQUE NOT NULL,
  chromosome TINYINT,
  position INT,
  ref_base CHAR,
  bases CHAR(2),
  relative_genotype VARCHAR(20) DEFAULT '',
  ensembl_acc CHAR(15) DEFAULT NULL,
  consequence CHAR(200) DEFAULT '',
  dbsnp CHAR(12) DEFAULT NULL,
  1kg_major CHAR(1) DEFAULT NULL,
  1kg_minor CHAR(1) DEFAULT NULL,
  1kg_maf FLOAT DEFAULT 0,
  hm3_major CHAR(1) DEFAULT NULL,
  hm3_minor CHAR(1) DEFAULT NULL,
  hm3_maf FLOAT DEFAULT 0,
  merged BOOLEAN DEFAULT false
);

SELECT 'Loading data' FROM data LIMIT 1;
LOAD DATA LOCAL INFILE 'data.tsv' INTO TABLE data FIELDS TERMINATED BY "\t" LINES TERMINATED BY "\n";

SELECT 'Creating indices' FROM data LIMIT 1;
CREATE INDEX idx_data_uid on data(uid);
CREATE INDEX idx_data_chromosome on data(chromosome);
CREATE INDEX idx_data_position on data(position);
CREATE INDEX idx_data_chromosome_position ON data(chromosome, position);

-- QC
SELECT 'Number of SNPs where bases are different from stored: ' FROM data limit 1;
SELECT count(*)
FROM data d, exome_common.already_found_ceu c
WHERE d.uid = c.uid
AND d.bases <> c.bases;

SELECT 'Updating from SNPs that are already found' FROM data LIMIT 1;
UPDATE data d, exome_common.already_found_ceu c
SET d.relative_genotype = c.relative_genotype,
    d.ensembl_acc = c.ensembl_acc,
    d.consequence = c.consequence,
    d.dbsnp = c.dbsnp,
    d.1kg_major = c.1kg_major,
    d.1kg_minor = c.1kg_minor,
    d.1kg_maf = c.1kg_maf,
    d.hm3_major = c.hm3_major,
    d.hm3_minor = c.hm3_minor,
    d.hm3_maf = c.hm3_maf,
    d.merged = TRUE
WHERE d.uid = c.uid
AND d.bases = c.bases;

SELECT 'Merging dbSNP' FROM data LIMIT 1;
UPDATE data d, exome_common.dbsnp r
  SET d.dbsnp = r.name
  WHERE d.uid = r.uid
  AND d.merged = FALSE;

SELECT 'Merging 1000genomes' FROM data LIMIT 1;
UPDATE data d, exome_common.1kG_ceu r
  SET d.1kg_major = r.ref_base, d.1kg_minor = r.alt_base, d.1kg_maf = r.alt_allele_freq
  WHERE d.uid = r.uid
  AND d.merged = FALSE;

SELECT 'Merging hapmap3' FROM data LIMIT 1;
UPDATE data d, exome_common.hm3_ceu r
  SET d.hm3_major = r.major_allele, d.hm3_minor = r.minor_allele, d.hm3_maf = r.maf
  WHERE d.uid = r.uid
  AND d.merged = false;

SELECT 'Merging genes for chromosome 24' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 24
  AND r.chromosome = 24
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 23' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 23
  AND r.chromosome = 23
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 22' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 22
  AND r.chromosome = 22
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 21' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 21
  AND r.chromosome = 21
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 20' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 20
  AND r.chromosome = 20
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 19' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 19
  AND r.chromosome = 19
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 18' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 18
  AND r.chromosome = 18
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 17' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 17
  AND r.chromosome = 17
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 16' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 16
  AND r.chromosome = 16
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 15' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 15
  AND r.chromosome = 15
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 14' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 14
  AND r.chromosome = 14
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 13' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 13
  AND r.chromosome = 13
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 12' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 12
  AND r.chromosome = 12
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 11' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 11
  AND r.chromosome = 11
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 10' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 10
  AND r.chromosome = 10
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 9' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 9
  AND r.chromosome = 9
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 8' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 8
  AND r.chromosome = 8
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 7' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 7
  AND r.chromosome = 7
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 6' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 6
  AND r.chromosome = 6
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 5' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 5
  AND r.chromosome = 5
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 4' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 4
  AND r.chromosome = 4
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 3' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 3
  AND r.chromosome = 3
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 2' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 2
  AND r.chromosome = 2
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Merging genes for chromosome 1' FROM data LIMIT 1;
UPDATE data d, exome_common.genes r
  SET d.ensembl_acc = r.ensembl_acc
  WHERE d.chromosome = 1
  AND r.chromosome = 1
  AND d.position BETWEEN r.start AND r.stop
  AND d.merged = FALSE;

SELECT 'Swapping major<->minor' FROM data LIMIT 1;
ALTER TABLE data ADD COLUMN swap CHAR(1);
UPDATE data
  SET swap = 1kg_major, 1kg_major = 1kg_minor, 1kg_minor = swap, 1kg_maf = (1-1kg_maf)
  WHERE 1kg_maf > 0.5
  AND merged = FALSE;
ALTER TABLE data DROP COLUMN swap;

SELECT 'Updating majormajor relative genotype' FROM data LIMIT 1;
UPDATE data
  SET relative_genotype = 'majormajor'
  WHERE 1kg_major IS NOT NULL
  AND bases = CONCAT(1kg_major,1kg_major)
  AND merged = FALSE;

SELECT 'Updating majorminor relative genotype' FROM data LIMIT 1;
UPDATE data
  SET relative_genotype = 'majorminor'
  WHERE 1kg_major IS NOT NULL
  AND 1kg_minor IS NOT NULL
  AND (bases = CONCAT(1kg_major,1kg_minor) OR bases = CONCAT(1kg_minor,1kg_major))
  AND merged = FALSE;

SELECT 'Updating minorminor relative genotype' FROM data LIMIT 1;
UPDATE data
  SET relative_genotype = 'minorminor'
  WHERE 1kg_minor IS NOT NULL
  AND bases = CONCAT(1kg_minor,1kg_minor)
  AND merged = FALSE;

SELECT 'Updating refref relative genotype' FROM data LIMIT 1;
UPDATE data
  SET relative_genotype = 'refref'
  WHERE 1kg_major IS NULL
  AND bases = CONCAT(ref_base,ref_base)
  AND merged = FALSE;

SELECT 'Updating refalt relative genotype' FROM data LIMIT 1;
UPDATE data
  SET relative_genotype = 'refalt'
  WHERE 1kg_major IS NULL
  AND relative_genotype IS NULL
  AND bases LIKE CONCAT('%',ref_base,'%')
  AND merged = FALSE; -- need previous update first

SELECT 'Updating altalt relative genotype' FROM data LIMIT 1;
UPDATE data
  SET relative_genotype = 'altalt'
  WHERE 1kg_major IS NULL
  AND relative_genotype IS NULL
  AND merged = FALSE; -- need previous update first

SELECT 'Optimizing table data' FROM data LIMIT 1;
OPTIMIZE TABLE data;
