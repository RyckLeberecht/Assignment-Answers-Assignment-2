#!/usr/bin/env ruby
require 'json'



# @author Ryck Leberecht
# A class which acts complementary to other classes by providing features like import / export of data, Annotations or formatting tools.
class Complement  
  
  # Loads a list of genes
  # @param [String] file a txt file with one AGI gene code each delimiter
  # @return [Array] input_genes an Array with all the Genes found in the input-file
  # @note Aborts the process if an entry is detected who does not comply with the AGI format  
  def self.load_genes(file)
      input = File.read(file) #Read in the text file
      input_genes = []
      input_regex = /^[A-Z]{2}[0-9][a-z][0-9]{5}$/ #Define the regular expression for the input structure
      rows = input.split("\n") #Split the text into rows using the observed input separator
      rows.each do |row|
        if row.match?(input_regex)# Check if the input matches the given structure
            input_genes << row # Add the input to the array if it matches the structure
        else
            abort("ERROR: Wrong format for gene: #{row}") #Abort process if an entry doesn't comply with expectation
        end
      end
    return input_genes
  end
  
  # This will recursively flatten a hash until all of the elements have been extracted into a single, one-dimensional array
  # @param [Hash] h a hash we wish to flatten
  # @return [Array] an array containing all entries in a hash
  def self.recursive_flatten(h)
  h.map { |k, v| v.is_a?(Hash) ? recursive_flatten(v) : [k, v] }.flatten
  end
    
  
  # Finds Gene Ontologies for a given gene
  # @param [String] input_gene a gene for which we have to find GOs
  # @return [Hash] a Hash with the GO ID as the key the GO Terms as entries
  def self.go_info(input_gene)
    terms = {}
    data = RestClient.get("http://togows.dbcls.jp/entry/uniprot/#{input_gene}/dr.json")
    data = JSON.parse(data.body)[0]
    terms = data["GO"].map do |go|
      [go[0].sub("GO:", ""), go[1].delete("P:")] if go[1].start_with?("P:")
    end.compact.to_h
    return terms
  end 
  
  # Finds Kegg Pathways for a given gene
  # @param [String] input_gene a gene for which we have to find KPs
  # @return [Hash] a Hash with the KP ID as the key the KP Terms as entries 
  def self.kegg_info(input_gene)
    data = RestClient.get("http://togows.org/entry/kegg-genes/ath:#{input_gene}/pathways.json").body
    data = JSON.parse(data)[0] 
    return data
  end

  # Annotates the given gene
  # @param [String] input_gene a gene which we have to annotate
  # @return [Object] a object with the genes GO, KP, interactions and potentially additional information
  def self.gene_info(input_gene)
    go_terms = self.go_info(input_gene) # Find Gene Ontologies for the input gene
    kegg_terms = self.kegg_info(input_gene) # Find Kegg Pathways for the input gene
    gene_details = fetch_gene_details(input_gene) # Add additional information about the gene
    # Build the object with the gene's GO, KP, and interactions
    object = {
      "gene_id" => input_gene.upcase,
      "kegg_pathways" => kegg_terms,
      "gene_ontologies" => go_terms,
      "name" => gene_details["name"],
      "description" => gene_details["description"]
    }
    return object
  end

  # Fetch additional details about the gene
  # @param [String] input_gene a gene for which we want additional information
  # @return [hash] with the genes name and description
  # @note I found this option in the internet and tried to implement it but either the genes are not known enough or I did something wrong.
  def self.fetch_gene_details(input_gene)
    begin
      data = RestClient.get("http://mygene.info/v3/query?q=symbol:#{input_gene}")
      data = JSON.parse(data.body)["hits"]
      if data.nil? || data.empty?
        return { "name" => "", "description" => "" }
      else
        data = data[0]
        return { "name" => data["name"], "description" => data["description"] }
      end
    rescue RestClient::ExceptionWithResponse => e
      return { "name" => "", "description" => "" }
    end
  end

  # Generate a report file
  # @param [Hash] new_groups a hash with each group of connected genes
  # @return [File] creates a txt file filled with the report
  def self.generate_report(new_groups)
  File.open("report.txt", "w") do |file|
      new_groups.each do |group, genes|
        file.puts "Group: #{group}"
        genes.each do |gene|
          gene_annotation = self.gene_info(gene)
          file.puts "Gene: #{gene_annotation['gene_id']}"
          file.puts "Name: #{gene_annotation['name']}"
          file.puts "Description: #{gene_annotation['description']}"
          file.puts "Gene Ontologies: #{gene_annotation['gene_ontologies']}"
          file.puts "Kegg Pathways: #{gene_annotation['kegg_pathways']}"
          file.puts
        end
      end
    end
  end
  
end

