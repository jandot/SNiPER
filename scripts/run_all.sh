mysql -h xxx -u yyy -pzzz -e "CREATE DATABASE exome_www;"
cat HuRef.InternalHuRef-NCBI.gff | ruby ./prepare_upload.rb > data.tsv
mysql -h xxx -u yyy -pzzz -D exome_www < create_table.sql
mysql -h xxx -u yyy -pzzz -D exome_www -e "SELECT chromosome, position, bases FROM data WHERE ensembl_acc IS NOT NULL AND merged = FALSE;" > tmp.tsv
grep -v 'chromosome' tmp.tsv | awk '{split($3, a, ""); print $1 "\t" $2 "\t" a[1] "/" a[2]}' > to_annotate.tsv
cat to_annotate.tsv | perl ./snp_funcannot.pl > annotated.tsv
cat annotated.tsv | ruby ./reformat_consequences.rb > update_consequences.sql
mysql -h xxx -u yyy -pzzz -D exome_www < update_consequences.sql
mysql -h xxx -u yyy -pzzz -D exome_www < copy_to_common.sql
