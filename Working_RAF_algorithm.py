from collections import defaultdict
import random
import copy
import linecache

random.seed("Federico seed8") #Uncomment this in order to obtain a fixed series of random numbers! Useful for debugging.

MolNumber = raw_input("Please insert the number of molecules produced.")
ReactionNumber = raw_input("Please insert the number of reactions produced.")
cat_probability = raw_input("Please insert the probability of a molecule to catalyse any given reaction.")
molecule_list = []
reaction_list = []

for i in range(int(MolNumber)): #changed from MolNumber + 1 to MolNumber
	molecule = "M" + str(i)
	molecule_list.append(molecule)

for i in range(int(ReactionNumber)): #changed from ReactionNumber + 1 to ReactionNumber
	reaction = "R" + str(i) + " "
	reaction_list.append(reaction)

print molecule_list
print reaction_list

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
		
	with open("network_output1.txt", 'r') as dot_file: #Extract the unique identifiers "M + integer" for each food molecule
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
	print "food list", clean_list
	print "food set", set(clean_list)
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
	print "Length of ReactsAndProds:", len(ReactsAndProds)
	print "------------------------------------------------"
	return ReactsAndProds

##################ASSIGN RANDOM CATALYSATION PAIRS BASED ON A UNIFORM DISTRIBUTION#######################################################################

def assign_rand_catalysis(molecules, RandP_sets): #Changed the argument name from "reactions" to "RandP_sets", as it will now accept the list of tuples of sets produced by FindReactionProducts!
	print "\n"
	print "###############  ASSIGN RANDOM CATALYSATION PAIRS  ###################"
	print "\n"
	catalysations = []
	threshold = 1.0 - float(cat_probability) #Here the threshold is given by a user-defined probability.

	#threshold = 1.0 - 1.0/int(ReactionNumber) #Each molecule has the same probability of catalysing any one reaction; i.e. p = 1/(ReactionNumber)

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
				cat_edge = "  " + str(catalyst) + " -> " + "R%d" % RandP_sets.index(pair[1]) + " [style=dotted]; // Random Catalysation!;\n"
				print cat_edge
				network_file.write(cat_edge)
		network_file.write("\n}")	

	print "Catalysation pairs:"
	print catalysations
	print "length of catalysation pair list", len(catalysations)
	print "------------------------------------------------"
	return catalysations #returns a list of random molecule-reaction catalysation pairs

##################COMPUTE SUPPORT OF REACTIONS#######################################################################

def compute_support_list(RandP_sets):
	list_supp = []
	for reaction in RandP_sets:
		for reactant in reaction[0]:
			#print "reactant:", reactant
			if reactant not in list_supp:
				#print "adding to list_supp..."
				list_supp.append(reactant)
		for product in reaction[1]:
			#print "product:", product
			if product not in list_supp:
				#print "adding to list_supp..."
				list_supp.append(product)
	print "support of RPlist:", list_supp
	print "Length of list_supp:", len(list_supp)
	return list_supp

##################FIRST SUBROUTINE#######################################################################

def ReduceToRA(RandPsets, catalyses, foodset): # Eliminates all reactions that are not catalysed from the set "RandPsets"
	print "\n"
	print "###############  SUBROUTINE 1: REDUCE TO RA  ###################"
	print "\n"
	counter = 0
	while True:
		counter += 1
		closure = ComputeClosure(foodset, RandPsets) #Computes the closure again each time the set of reactions is reduced!
		copy_RandPsets = copy.deepcopy(RandPsets)
		copy_catalyses = copy.deepcopy(catalyses)
		for react in RandPsets:
			print "Reducing to reflexively autocatalytic... Is reaction %d catalysed by molecules within the closure?" % RandPsets.index(react) 

#According to Hordjik et al. (2011) and Steel et al. (2013), a set of reactions "R'" is reflexively-autocatalytic if for all "r" in "R'", there is at least one molecule "x", element of closure[of FOOD with respect to "R'"] such that ('x', 'r') is part of the set of catalysations...

			print "catalyst list:", catalyses[RandPsets.index(react)][0]

			catalysis_counter = 0
			for mol in catalyses[RandPsets.index(react)][0]:
				if mol in closure:

					catalysis_counter += 1
					print "Yes! Has %d catalysts." % catalysis_counter
				else:
					print "removing catalyst %s because it is not within the closure!" % mol
					copy_catalyses[RandPsets.index(react)][0].remove(mol) #remove the catalyst that is not within the closure set!
			if not catalyses[RandPsets.index(react)][0]:
				print "No! Eliminated the following reaction from the RP list and from the catalysation list:"
				print react
				copy_RandPsets.remove(react)
				copy_catalyses.remove(catalyses[RandPsets.index(react)])
		print "Went through ReducetoRA while loop %d times." % counter

		if len(RandPsets) == len(copy_RandPsets) and all([element[0] for element in copy_catalyses]): #the all statement returns TRUE if all catalyst lists are non-empty!
			break
		else:
			RandPsets = copy_RandPsets
			catalyses = copy_catalyses
	print "------------------------------------------------"
	print "Reduced reaction list:"
	print copy_RandPsets
	print "Reduced RPlist length:", len(copy_RandPsets)
	print "Reduced catalysation list:"
	print copy_catalyses
	print "Reduced catalysation list length:", len(copy_catalyses)
	print "------------------------------------------------"
	return copy_RandPsets, copy_catalyses # --> The reduced reaction set and catalysation set!

##################SECOND SUBROUTINE#######################################################################

def ComputeClosure(Food, RAset):
	print "\n"
	print "###############  SUBROUTINE 2: COMPUTE CLOSURE  ###################"
	print "\n"
	closure_set = Food
	print "Starting closure set (food):", closure_set
	counter = 0
	while True:
		counter += 1
		copy_closure_set = copy.deepcopy(closure_set)
		for reactiontuple in RAset:
			#print reactiontuple[0], reactiontuple[1]
			if reactiontuple[0].issubset(closure_set) and not reactiontuple[1].issubset(closure_set):
				copy_closure_set.update(reactiontuple[1])
				#print "The closure set is now:", copy_closure_set
				#print "-----------------------------------"
		print "Went through closure while loop %d times." % counter
		if set(copy_closure_set) == set(closure_set):
			break
		else:
			closure_set = copy_closure_set
	print copy_closure_set
	return copy_closure_set

##################THIRD SUBROUTINE#######################################################################

def ReduceToFgenerated(final_closure, RAset, RA_catset, nonReducedSet):
	print "\n"
	print "###############  SUBROUTINE 3: REDUCE TO F-GENERATED  ###################"
	print "\n"
	print "RAset:"
	print "\n"
	print RAset

	print "\n"
	print "RA catalysation set:"
	print "\n"
	print RA_catset
	
	copy_RAset = copy.deepcopy(RAset) #ADDED A COPY LINE
	copy_RA_catset = copy.deepcopy(RA_catset)
	for reactiontuple in RAset:
		print "reaction:", reactiontuple
		print "Reducing to F-generated... Are the reactants of reaction %d F-generated?" % nonReducedSet.index(reactiontuple)
		if not reactiontuple[0].issubset(final_closure): #Are the reactants of this reaction a subset of the closure?
			print "No! Eliminated the following reaction:" 
			print reactiontuple
			copy_RAset.remove(reactiontuple)
			copy_RA_catset.remove(RA_catset[RAset.index(reactiontuple)])
		else:
			print "Yes!"
	print "-----------------------------------"
	print "The RA, F-generated set is:"
	print copy_RAset

	print "The RA, F-generated set with catalysations is:"
	print copy_RA_catset

	return copy_RA_catset

###########################RETRIEVE REACTION NUMBERS##########################

def Retrieve_RAF_reactionNum(RAFset, nonReducedSet):
	print "\n"
	print "###############  RETRIEVING REACTION NUMBER...  ###################"
	print "\n"
	print "Non-reduced RPset:", nonReducedSet
	print "\n"
	NumberedRAFset = []
	for reaction in RAFset:
		print nonReducedSet.index(reaction[1])
		reaction.append(nonReducedSet.index(reaction[1]))
		print reaction			 
		NumberedRAFset.append(reaction)
	print NumberedRAFset			
	return NumberedRAFset

####################HIGHLIGHT RAF SET IN DOT FILE##############################################	

def write_RAF(numbered_RAFset):
	with open("myRAF.txt", 'w') as RAF_dotfile:
		print "\n"
		print "###############  WRITE RAF TO TEXT FILE  ###################"
		print "\n"
		for reaction_list in numbered_RAFset:
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
				cat_edge = "  " + catalyst + " -> " + "R%d" % reaction_list[2] + " [style=dotted];\n" #Will need to strip empty spaces from edges to compute match! in modify_toyChem
				print 'RAF catalysation edge:', cat_edge
				RAF_dotfile.write(cat_edge)


def modify_toyChem():
	with open("myRAF.txt", 'r') as RAF_dotfile, open("RAFhighlight.txt", 'w') as finalfile, open("network_output1.txt", 'r') as toyChem_dotfile:
		for line in toyChem_dotfile: #Stoichiometry is maintained due to order of for loops!!!
			print "writing RAFlighlight.txt file..."
			match_flag = 0
			RAF_dotfile.seek(0)
			for RAFline in RAF_dotfile:
				stripped_RAFline = RAFline.strip()
				if stripped_RAFline in line:
					match_flag = 1
					if "style=dotted" in RAFline:
						catline = RAFline.replace("]", ", color=red]") #Change colour of RAF catalysation edges to RED!
						finalfile.write(catline)
					else:
						red_RAFline = RAFline.replace(";", " [color=red];") # All RAF edges are in RED!
						finalfile.write(red_RAFline)
					
			if match_flag == 0:
				finalfile.write(line)



food = FindFood()
RPlist = FindReactantsProducts(food)
catalysations = assign_rand_catalysis(molecule_list, RPlist)
RA_set = ReduceToRA(RPlist, catalysations, food)

RA_reactions = RA_set[0]
RA_catalysations = RA_set[1]

RPlist_2 = FindReactantsProducts(food)
finalclosure = ComputeClosure(food, RA_reactions)
RAF = ReduceToFgenerated(finalclosure, RA_reactions, RA_catalysations, RPlist_2)

if RAF:
	numberedRAF = Retrieve_RAF_reactionNum(RAF, RPlist_2) 
	write_RAF(numberedRAF)
	modify_toyChem()


