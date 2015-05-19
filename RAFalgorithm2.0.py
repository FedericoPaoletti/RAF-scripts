from collections import defaultdict
import random
import copy
import linecache

MolNumber = raw_input("Please insert the number of molecules produced.")
ReactionNumber = raw_input("Please insert the number of reactions produced.")
#dotfile = raw_input("Please insert the name of the DOT output file.")
molecule_list = []
reaction_list = []

for i in range(int(MolNumber)): #changed from MolNumber + 1 to MolNumber
	molecule = "M" + str(i)
	molecule_list.append(molecule)

for i in range(int(ReactionNumber)): #changed from ReactionNumber + 1 to ReactionNumber
	reaction = "R" + str(i) + " "
	reaction_list.append(reaction)

##################DEFINE FOOD, REACTANTS AND PRODUCTS#######################################################################
	
def FindFood():
	print "\n"
	print "###############  DEFINE FOOD  ###################"
	print "\n"
	smile_list = []
	food_list = []

	
	foodsmile = raw_input('Please insert a food SMILE within \" \".')
	while foodsmile != "":
		smile_list.append(foodsmile)
		foodsmile = raw_input('Please insert a food SMILE within \" \".')
	#print smile_list
		
	with open("network_output1.txt", 'r') as dot_file:
		for smile in smile_list:
			#print smile
			dot_file.seek(0)
			for line in dot_file:
				#print line
				if smile in line:
					#print line[2:5]
					food_list.append(line[2:5]) #Be careful in changing indentation of input file (.txt file)!!!
		#print food_list
	clean_list = []
	for food in food_list:
		cleanfood = food.strip()
		clean_list.append(cleanfood)
	#print clean_list
	return set(clean_list)

def FindReactantsProducts(foodset):
	print "\n"
	print "###############  EXTRACT REACTANTS AND PRODUCTS FROM DOT FILE  ###################"
	print "\n"
	ReactsAndProds = []
	Reactants = []
	Products = []
	with open("network_output1.txt", 'r') as dot_file:
		for line in dot_file:
			if "STOP HERE" in line:
				ReactsAndProds.append((set(Reactants), set(Products))) #Added the line of code because the function was NOT adding the last reaction in from the DOT output file!
				break
			elif "shape=box" in line: #Have reached a new reaction description
				if not not Reactants and not not Products:
					#print line
					#print Reactants
					#print Products
					ReactsAndProds.append((set(Reactants), set(Products)))
					Reactants = []
					Products = []
			elif "->" in line: #Have reached a reaction edge
				if line[2] == "M":
					react = line[2:5].strip()#Position of the reactant within the text line!
					Reactants.append(react) 
				elif line[2] == "R":
					prod = (line[8:12].strip()).replace(";", "")#Position of the product within the text line!
					Products.append(prod)
	print "Reactions displayed as reactant and product tuples:"
	print ReactsAndProds
	#print "Length of ReactsAndProds:", len(ReactsAndProds)
	print "------------------------------------------------"
	return ReactsAndProds

##################ASSIGN RANDOM CATALYSATION PAIRS BASED ON A UNIFORM DISTRIBUTION#######################################################################

def assign_rand_catalysis(molecules, RandP_sets): #Changed the argument name from "reactions" to "RandP_sets", as it will now accept the list of tuples of sets produced by FindReactionProducts!
	print "\n"
	print "###############  ASSIGN RANDOM CATALYSATION PAIRS  ###################"
	print "\n"
	catalysations = []
	threshold = 1.0 - 1.0/int(ReactionNumber) #Each molecule has the same probability of catalysing any one reaction; i.e. p = 1/(ReactionNumber)

	for react in RandP_sets:
		catalyst_list = []
		for mol in molecules:
			if random.random() > threshold:
				catalyst_list.append(mol)
		catalysations.append([catalyst_list, react])

	print "Writing catalysation edges to dotfile..."
	with open("network_output1.txt", 'a') as network_file:
		network_file.write("\n")
		for pair in catalysations:
			for catalyst in pair[0]:
				cat_edge = "  " + str(catalyst) + " -> " + "R%d" % RandP_sets.index(pair[1]) + " [color=blue, style=dotted]; // Random Catalysation!;\n"
				print cat_edge
				network_file.write(cat_edge)
		network_file.write("\n}")	

	print "Catalysation pairs:"
	print catalysations
	print "------------------------------------------------"
	return catalysations #returns a list of random molecule-reaction catalysation pairs

##################COMPUTE SUPPORT OF REACTIONS#######################################################################

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

##################FIRST SUBROUTINE#######################################################################

def ReduceToRA(RandPsets, catalyses, LISTED_SUPP): # Eliminates all reactions that are not catalysed from the set "RandPsets"
	print "\n"
	print "###############  SUBROUTINE 1: REDUCE TO RA  ###################"
	print "\n"
	copy_RandPsets = copy.deepcopy(RandPsets)
	for react in copy_RandPsets:
		print "Reducing to reflexively autocatalytic... Is reaction %d catalysed?" % copy_RandPsets.index(react)
		remove_flag = 1
		catalysis_counter = 0
		for mol in LISTED_SUPP:
			if mol in catalyses[copy_RandPsets.index(react)][0]:
				remove_flag = 0
				catalysis_counter += 1
				print "Yes! Has %d catalysts." % catalysis_counter
		if remove_flag == 1:
			print "No! Eliminated the following reaction:"
			print react
			RandPsets.remove(react)
	print "------------------------------------------------"
	print "Reduced Reaction list:"
	print RandPsets
	#print "Number of Reactions:", len(RandPsets)
	print "------------------------------------------------"
	return RandPsets # --> The reduced reaction set!

##################SECOND SUBROUTINE#######################################################################

def ComputeClosure(Food, RAset):
	print "\n"
	print "###############  SUBROUTINE 2: COMPUTE CLOSURE  ###################"
	print "\n"
	closure_set = Food
	
	for reactiontuple in RAset:
		#print reactiontuple[0], reactiontuple[1]
		if reactiontuple[0].issubset(closure_set) and not reactiontuple[1].issubset(closure_set):
			closure_set.update(reactiontuple[1])
			#print "The closure set is now:", closure_set
			#print "-----------------------------------"
	return closure_set

##################THIRD SUBROUTINE#######################################################################

def ReduceToFgenerated(final_closure, RAset, nonReducedSet):
	print "\n"
	print "###############  SUBROUTINE 3: REDUCE TO F-GENERATED  ###################"
	print "\n"
	for reactiontuple in RAset:
		print "reaction:", reactiontuple
		print "Reducing to F-generated... Are the reactants of reaction %d F-generated?" % nonReducedSet.index(reactiontuple)
		if not reactiontuple[0].issubset(final_closure):
			print "No! Eliminated the following reaction:" 
			print reactiontuple
			RAset.remove(reactiontuple)
		else:
			print "Yes!"
	print "-----------------------------------"
	print "The RA, F-generated set is:"
	print RAset
	return RAset

####################HIGHLIGHT RAF SET IN DOT FILE##############################################	

def Retrieve_RAF_reactionNum_and_cat(RAFset, cat_pairs, nonReducedSet):
	catalysed_RAFset = []
	for RAF_reaction in RAFset:
		for pair in cat_pairs:
			if RAF_reaction == pair[1]:
				pair.append(nonReducedSet.index(RAF_reaction)) #Add the original reaction number to the list containing catalyst(s) + reactants&products
				catalysed_RAFset.append(pair) #Add the new list containing the catalyst ID [0], the reactants and products sets [1], and the reaction number [2]
				
	print "-----------------------------------"
	print "The RAF set with catalysations is:"
	print catalysed_RAFset
	return catalysed_RAFset

def write_RAF(cat_RAFset):
	with open("myRAF.txt", 'w') as RAF_dotfile:
		print "\n"
		print "###############  WRITE RAF TO TEXT FILE  ###################"
		print "\n"
		for reaction_list in cat_RAFset:
			reaction_list[1] = list(reaction_list[1]) #convert sets back to lists!
			reaction_list[1][0] = list(reaction_list[1][0])
			reaction_list[1][1] = list(reaction_list[1][1])
			print "Reactants of reaction %d:" % reaction_list[2], reaction_list[1][0]
			print "Products of reaction %d:" % reaction_list[2], reaction_list[1][1]

			for reactant in reaction_list[1][0]:
				reactant_edge = "  " + reactant + " -> " + "R%d;\n" % reaction_list[2]
				print 'RAF reactant edge:', reactant_edge
				RAF_dotfile.write(reactant_edge)

			for product in reaction_list[1][1]:
				product_edge = "  R%d" % reaction_list[2] + " -> " + product + ";\n"
				print 'RAF product edge:', product_edge
				RAF_dotfile.write(product_edge)
			
			for catalyst in reaction_list[0]:
				cat_edge = "  " + catalyst + " -> " + "R%d" % reaction_list[2] + " [color=blue, style=dotted];\n" #Will need to strip empty spaces from edges to compute match! in modify_toyChem
				print 'RAF catalysation edge:', cat_edge
				RAF_dotfile.write(cat_edge)

################################### USE RAF TEXT FILE TO MODIFY TOYCHEM NETWORK DOT FILE ###################################

def modify_toyChem():
	with open("myRAF.txt", 'r') as RAF_dotfile, open("RAFhighlight.txt", 'w') as finalfile, open("network_output1.txt", 'r') as toyChem_dotfile:
		for line in toyChem_dotfile: #Stoichiometry is maintained due to order of for loops!!!
			print 'networkline loop 1:', line
			match_flag = 0
			RAF_dotfile.seek(0)
			for RAFline in RAF_dotfile:
				stripped_RAFline = RAFline.strip()
				if stripped_RAFline in line:
					match_flag = 1
					print 'RAFline:', RAFline
					print 'networkline loop 2:', line
					print 'Found match:'
					print line
					if "color=blue" in RAFline:
						catline = RAFline.replace("blue", "red") #Change colour of RAF catalysation edges to RED!
						finalfile.write(catline)
					else:
						non_catline = RAFline.replace(";", " [color=green];") # RAF reaction edges are in GREEN!
						finalfile.write(non_catline)
					
			if match_flag == 0:
				finalfile.write(line)


food = FindFood()
RPsets = FindReactantsProducts(food)
catalysations = assign_rand_catalysis(molecule_list, RPsets)

support_dict = compute_support_dict()
support_list = compute_support_list(support_dict)
RAreactions = ReduceToRA(RPsets, catalysations, support_list)

closure = ComputeClosure(food, RAreactions)

RPsetsAGAIN = FindReactantsProducts(food)
myRAFset = ReduceToFgenerated(closure, RAreactions, RPsetsAGAIN)
finalRAF = Retrieve_RAF_reactionNum_and_cat(myRAFset, catalysations, RPsetsAGAIN)

write_RAF(finalRAF)
modify_toyChem()
