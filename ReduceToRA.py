from collections import defaultdict
import random

MolNumber = raw_input("Please insert the number of molecules produced.")
ReactionNumber = raw_input("Please insert the number of reactions produced.")
molecule_list = []
reaction_list = []

for i in range(int(MolNumber) + 1):
	molecule = "M" + str(i)
	molecule_list.append(molecule)

for i in range(int(ReactionNumber) + 1):
	reaction = "R" + str(i) + " "
	reaction_list.append(reaction)

def assign_rand_catalysis(molecules, reactions):

	catalysations = []
	threshold = 1.0 - 1.0/int(ReactionNumber) #Each molecule has the same probability of 		catalysing any one reaction; i.e. p = 1/(ReactionNumber)

	for mol in molecules:
		for react in reactions:
			if random.random() > threshold:
				catalysations.append([mol, react])

	return catalysations #returns a list of random molecule-reaction catalysation pairs

def compute_support_dict():

	edge_list = []

	with open("network_output1.txt", 'r') as network_file:
		for line in network_file:
			if "STOP HERE" in line:
				break
			elif "->" in line:
				line_no_empty_space = line.strip()
				line_no_arrow = line_no_empty_space.replace("->", "")
				clean_line = line_no_arrow.replace(";", " ")
				edge_list.append(clean_line)
	#print edge_list
	support = defaultdict(list)
	#print support

	for reaction in reaction_list:
		for reaction_edge in edge_list:
			if reaction in reaction_edge:
				edge_no_empty_space = reaction_edge.replace(reaction, "")
				clean_edge = edge_no_empty_space.strip()
				#print reaction_edge
				#print clean_edge
				#print reaction
				if clean_edge not in support[reaction]:
					support[reaction].append(clean_edge)

	#print "Support values:", support.values()
	#print "Support dictionary:", support
	return support

def compute_support_list(dict_support):
	list_supp = []
	for mol_group in dict_support.values():
		for mol in mol_group:
			if mol not in list_supp:
				list_supp.append(mol)
	return list_supp

#print "Listed support:", list_supp
	
###############################

def ReduceToRA(reactions, catalyses, LISTED_SUPP): # Eliminates all reactions that are not catalysed from the set "reactions"
	#print reactions
	for react in reactions:
		remove_flag = 1
		for mol in LISTED_SUPP:
			if [mol, react] in catalyses:
				remove_flag = 0
				#print [mol, react]
		if remove_flag == 1:
			reactions.remove(react)
	#print reactions
	return reactions

supp = compute_support_dict()
supplist = compute_support_list(supp)
catalysations = assign_rand_catalysis(molecule_list, reaction_list)
#print catalysations
ReduceToRA(reaction_list, catalysations, supplist)
