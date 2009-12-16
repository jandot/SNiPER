#!/usr/local/bin/perl

use strict;
use lib '/nfs/team29/projects/exome/SNiPER/ensembl_api_56/ensembl/modules';
use lib '/nfs/team29/projects/exome/SNiPER/ensembl_api_56/ensembl-variation/modules';

use Bio::EnsEMBL::Registry;

# get registry
my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');

my $vfa = $reg->get_adaptor('Human', 'variation', 'VariationFeature');
my $sa  = $reg->get_adaptor('Human', 'core', 'Slice');

#######################################################
## map build 36 to build 37

while ( my $line = <STDIN> ) {
  chomp $line;
  my @fields = split(/\t/, $line);
  my ($build36chr, $build36start, $build36end, $build36strand, $build36alleles) =
    ($fields[0], $fields[1], $fields[1], '1', $fields[2]);

  if ( $build36chr == '23' ) {
    $build36chr = 'X';
  } elsif ( $build36chr == '24' ) {
    $build36chr = 'Y';
  }

#  ## snp
#  my($build36chr,$build36start,$build36end,$build36strand,$build36alleles) = 
#      ('13','31798639','31798639','1','C/T');
  
  ## deletion
  #my($build36chr,$build36start,$build36end,$build36strand,$build36alleles) = 
  #    ('13','31798639','31798640','1','CG/-');
  
  ## insertion
  ##my($build36chr,$build36start,$build36end,$build36strand,$build36alleles) = 
  ##    ('13','31798640','31798639','1','-/A');
  
  
  ## my build 36 coordinate slice
  my $sliceA = $sa->fetch_by_region('chromosome',$build36chr,$build36start,$build36end,$build36strand,'NCBI36');
  
  ## convert to build 37
  my $newsliceA = $sliceA->project('chromosome', 'GRCh37');
  
  unless ( defined $newsliceA->[0] ) {
    print STDERR "Could not map $build36chr position $build36start onto GRCh37\n";
    next;
  }
  
  my($build37chr,$build37start,$build37end) =
      ($newsliceA->[0]->to_Slice()->seq_region_name(),$newsliceA->[0]->to_Slice()->start,$newsliceA->[0]->to_Slice()->end);
  
#  print "NEW COORDINATES ON BUILD 37 FOR $build36chr,$build36start,$build36end,$build36strand ARE ".
#      "$build37chr,$build37start,$build37end" .",". $newsliceA->[0]->to_Slice()->strand ."\n";
  
  my @all = split(/\//, $build36alleles);
#  if ( $newsliceA->[0]->to_Slice()->seq ne $all[0] && $all[0] ne '-' ) {
#    warn "$build36chr, $build37start, $build36alleles: BASE mismatch\n";
#    next;
#  } 
#  die "BASE mismatch\n" if $newsliceA->[0]->to_Slice()->seq ne $all[0] && $all[0] ne '-';
  
#  print $newsliceA->[0]->to_Slice()->seq  ." eq $all[0]\n";
  
  my $build37alleles;
  if($newsliceA->[0]->to_Slice()->strand != 1){
      $build37alleles = $newsliceA->[0]->to_Slice()->seq ."/". revcomp($all[1]);
  }else{
      $build37alleles = $build36alleles;
  }
  
  #######################################################
  # get a slice for the new feature to be attached to
  my $slice = $sa->fetch_by_region('chromosome', $build37chr);
  
  # create a new VariationFeature object
  my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
    -start => $build37start,
    -end   => $build37end,
    -slice => $slice,                      # the variation must be attached to a slice
    -allele_string  => $build37alleles,    # the first allele should be the reference allele
    -strand         => 1,
    -map_weight     => 1,
    -adaptor        => $vfa,               # we must attach a variation feature adaptor
    -variation_name => 'newSNP',
  );
  
  # get the consequence types
  my $cons = $new_vf->get_all_TranscriptVariations();
  
  my @consequences;
  foreach my $con(@{$new_vf->get_all_TranscriptVariations}) {
    foreach my $string(@{$con->consequence_type}) {
      push @consequences, $string;
#      print
#        "Variation ", $new_vf->variation_name,
#        " at position ", $new_vf->seq_region_start,
#        " on chromosome ", $new_vf->seq_region_name,
#        " has consequence ", $string,
#        " in transcript ", $con->transcript->stable_id, "\n";
    }
  }
  my $uniq_string = join(";", unique(@consequences));
  if ( $build36chr == 'X' ) {
    $build36chr = '23';
  } elsif ( $build36chr == 'Y' ) {
    $build36chr = '24';
  }
  print(join("\t", ($build36chr, $build36start, $build36alleles, $uniq_string)), "\n");
}


sub revcomp {
#  my $val = $_[0];
#  $val =~ tr/ACGTacgt/TGCAtgca/;
#  return $val;
}

sub unique {
  my @values = @_;
  my %u_values;
  foreach my $val (@values) {
    $u_values{$val} = 1;
  }
  return keys %u_values;
}

