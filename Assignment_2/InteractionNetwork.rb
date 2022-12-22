#!/usr/bin/env ruby

require 'rest-client'

# @author Ryck Leberecht
# A class which generates, based on the input, an Interaction Network of given Genes
class InteractionNetwork
  
  # A hash to save all first level interactions
  @@all_gene_interactions = {}
  
  # Finds interactions for the given gene
  # @param [String] input_gene a given gene, for which we want to find its interactions
  # @return [Array] a Array with all found interactions as entries
  # @note The recorded interactions are filtered by the following criteria:
  # @note 1. The interaction has to take place in the same organism (Arabidopsis Thaliana)
  # @note 2. The encoded gene has to follow the AGI format
  # @note 3. The miscore has to be higher than 0.4
  def self.find_interaction(input_gene)
      gene_regex = /([a-z,A-Z]{2}[0-9][a-z,A-Z][0-9]{5})/ # Regular expression to search for AGI coded genes
      miscore_regex = /i\w+-\w+:(0\.\d+?)/ # Regular expression to search for miscore values
      input_gene = input_gene.join('') if input_gene.is_a?(Array)
      input_gene = input_gene.upcase
      interaction = [] 
      data = RestClient.get("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{input_gene}") # data are strings!
      rows = data.split("\n") # Dependend on the input data needs to be adjusted to for example "\t"
      rows.each do |row|
        if row.split("\t")[9..10].all? { |element| element.include?("taxid:3702") } # Only execute if both interactors are Arabidopsis
            genes = row.scan(gene_regex) # Search for AGI coded genes
            genes.each { |element| element[0].upcase! } # Convert genes to uppercase
            genes = genes.uniq # Only the first and second interactor remain
            intact_miscore = row.scan(miscore_regex) # Search for miscore
            # Only records the interaction if one interactor is the given gene and the other one is an unknown one
            if genes.length == 2  
              if (genes[0][0] != input_gene) && (intact_miscore[0][0].to_f > 0.4)
                  interaction << genes[0]
              elsif (genes[0][1] != input_gene) && (intact_miscore[0][0].to_f > 0.4)
                  interaction << genes[1]
              end
            end
        end
      end
      @@all_gene_interactions[input_gene] = interaction.uniq
      return @@all_gene_interactions[input_gene]
  end
  
  # Finds interactions of a given gene
  # @param [String] gene a given gene, for which we want to find its intereactions
  # @param [Integer] max_depth a number to controll the network depth
  # @param [Array] network a array where the results are saved
  # @param [Set] visited a empty set which functions as a blacklist so no interaction loops can develop
  # @param [Integer] depth a interger which determines the starting depth of the network
  # @return [Array] a array with gene as the first entry and all found interactions following
  def self.build_interaction_network(gene, max_depth, network = [], visited = Set.new, depth = 0)
    network << gene # Add the gene to the network
    visited.add(gene) # Add the gene to the blacklist
    return network if depth >= max_depth # Return the network if the maximum depth has been reached
    interactors = self.find_interaction(gene) # Find the direct interactors of the gene
    return network if interactors.empty? # Return the network even if there are no interactions
    interactors.each do |interactor| # Iterate over the interactors and build the interaction network for each of them
      next if visited.include?(interactor) # Skip the interactor if it has already been visited
      build_interaction_network(interactor, max_depth, network, visited, depth + 1) # Rerun the function one depth level up
    end
    return network
  end

  # Compile a list of interaction networks based on a list of input genes
  # @param [Array] gene_list a array with a filled with gene codes in the form of strings
  # @param [Integer] max_depth a number to controll the network depth
  # @return [Hash] a list of gene interaction networks
  def self.build_interaction_networks(gene_list, max_depth)
    networks = []
    gene_list.each do |gene|
      network = build_interaction_network(gene, max_depth)
      networks << { gene: gene, network: network }
    end
    return networks
  end

  # Find the interaction groups of genes in the given gene list
  # @param [Hash] networks a Hash containing all recorded interaction-networks of each gene in the list
  # @return [Hash] a list of grouped genes based on interaction
  def self.group_genes_by_overlap(networks)
    groups = {}
    networks.each_with_index do |gene_data, i| # Loop over each provided network
      gene_name = gene_data[:gene]
      network = gene_data[:network]
      group_name = nil
      groups.each do |name, genes| # Check if the gene is already in a group
        if genes.include?(gene_name)
          group_name = name
          break
        end
      end
      if group_name.nil? # If the gene is not in a group, create a new group
        group_name = "Group #{groups.size + 1}" # Assign a group number
        groups[group_name] = []
      end
      groups[group_name] << gene_name
      (i + 1...networks.size).each do |j| # Compare the current network to the other networks
        other_gene_data = networks[j]
        other_gene_name = other_gene_data[:gene]
        other_network = other_gene_data[:network]
        if (network & other_network).any? # If there is an overlap, add the other gene to the same group
          groups[group_name] << other_gene_name
        end
      end
    end
    groups.delete_if { |_, genes| genes.size == 1 } # Remove groups with only one member
    new_groups = {}
    groups.each do |group_name, genes| # Renumber the remaining groups
      new_group_name = "Group #{new_groups.size + 1}"
      new_groups[new_group_name] = genes
    end

    return new_groups
  end
end
