consequence_types = {'ESSENTIAL_SPLICE_SITE' => 1,
                     'STOP_GAINED' => 2,
		     'STOP_LOST' => 3,
		     'FRAMESHIFT_CODING' => 4,
		     'NON_SYNONYMOUS_CODING' => 5,
		     'SPLICE_SITE' => 6,
		     'SYNONYMOUS_CODING' => 7,
		     'REGULATORY_REGION' => 8,
		     '5PRIME_UTR' => 9,
		     '3PRIME_UTR' => 10,
		     'WITHIN_MATURE_miRNA' => 11,
		     'INTRONIC' => 12,
		     'WITHIN_NON_CODING_GENE' => 13,
		     'UPSTREAM' => 14,
		     'DOWNSTREAM' => 15}

STDIN.each do |line|
  chromosome, position, alleles, consequences = line.chomp.split(/\t/)
  consequence = (consequences.nil?) ? '' : consequences.split(/;/).sort_by{|c| consequence_types[c]}[0]
  STDOUT.puts "UPDATE data SET consequence = '#{consequence}' WHERE uid = '#{chromosome + '_' + position}';"
end
