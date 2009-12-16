iupac = {'A' => 'AA',
         'C' => 'CC',
	 'G' => 'GG',
	 'T' => 'TT',
	 'R' => 'AG',
         'Y' => 'CT',
         'M' => 'AC',
         'K' => 'GT',
         'W' => 'AT',
         'S' => 'CG',
         'B' => 'CGT',
         'D' => 'AGT',
         'H' => 'ACT',
         'V' => 'ACG',
         'N' => 'ACGT'}

STDIN.each do |line|
  next if line =~ /^#/
  fields = line.chomp.split(/\t/)
  variant_type = fields[2]
  next unless variant_type =~ /_SNP$/
  next if fields[6] == '-'

  chromosome = fields[0]
  chromosome = '23' if chromosome == 'X'
  chromosome = '24' if chromosome == 'Y'
  position = fields[4]
  key = [chromosome, position].join("_")
  attributes_array = fields[7].split(/;/)
  bases = attributes_array.shift
  next if bases.scan('/').length > 1
  bases.sub!(/\//,'')
  ref_base = bases.slice(0,1)
  STDOUT.puts [key, chromosome, position, ref_base, bases].join("\t")
end
