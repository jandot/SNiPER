require '~/LocalDocuments/Projects/biorake/lib/biorake.rb'

CONFIG = YAML::load(File.read('config.yml'))

task :default => :setup_new_individual

task :setup_new_individual, :individual_accession do
    STDERR.puts "Running setup for individual " + :individual_accession.to_s
    Rake.application.invoke_task("create_directory[:individual_accession.to_s]")
end

task :create_directory, :individual_accession do
    ind_directory = './individuals/' + :individual_accession.to_s
    STDERR.puts ind_directory
    Dir.mkdir ind_directory
end

task :copy_files do
    
end

task :create_database, :individual_accession do
    
end

task :parse_raw_data do
    
end

task :load_original_data do
    
end

task :create_list_snps_in_genes do
    
end

task :perform_functional_annotation do
    
end

task :load_functional_annotation do
    
end

task :reformat_functional_annotation do
    
end

task :copy_to_common do
    
end

