import random

def assign_rand_catalysis():

	molecule_list = []
	reaction_list = []
	catalysations = []
	MolNumber = raw_input("Please insert the number of molecules produced.")
	ReactionNumber = raw_input("Please insert the number of reactions produced.")

	for i in range(int(MolNumber) + 1):
		molecule = "M" + str(i)
		molecule_list.append(molecule)

	for i in range(int(ReactionNumber) + 1):
		reaction = "R" + str(i)
		reaction_list.append(reaction)

	threshold = 1.0 - 1.0/int(ReactionNumber) #Each molecule has the same probability of 		catalysing any one reaction; i.e. p = 1/(ReactionNumber)

	for mol in molecule_list:
		for react in reaction_list:
			if random.random() > threshold:
				catalysations.append([mol, react])

	return catalysations #returns a list of molecule-reaction catalysation pairs

dotfile = raw_input("Please insert the name of the DOT output file.")

with open(dotfile, 'a') as network_file:
	network_file.write("\n")
	for pair in catalysations:
		cat_edge = "  " + str(pair[0]) + " -> " + str(pair[1]) + " [color=red]; // Random Catalysation!;\n"
		network_file.write(cat_edge)
	network_file.write("\n}")
