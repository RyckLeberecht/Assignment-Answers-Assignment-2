#!/usr/bin/env ruby
require "./Complement.rb"
require "./InteractionNetwork.rb"

input_genes = Complement.load_genes("ArabidopsisSubNetwork_GeneList.txt")

data = InteractionNetwork.build_interaction_networks(input_genes, 3)

groups = InteractionNetwork.group_genes_by_overlap(data)

Complement.generate_report(groups)
