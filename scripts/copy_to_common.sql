INSERT INTO exome_common.already_found_ceu (uid, chromosome, position, ref_base, bases, relative_genotype, ensembl_acc, consequence, dbsnp, 1kg_major, 1kg_minor, 1kg_maf, hm3_major, hm3_minor, hm3_maf)
  SELECT uid, chromosome, position, ref_base, bases, relative_genotype, ensembl_acc, consequence, dbsnp, 1kg_major, 1kg_minor, 1kg_maf, hm3_major, hm3_minor, hm3_maf
  FROM data
  WHERE merged = false;
UPDATE data SET merged = true;
