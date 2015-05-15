from collections import defaultdict

def FindFood():
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
	ReactsAndProds = []
	Reactants = []
	Products = []
	with open("network_output1.txt", 'r') as dot_file:
		for line in dot_file:
			if "STOP HERE" in line:
				break
			elif "shape=box" in line: #Have reached a new reaction description
				if not not Reactants and not not Products:
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
	print ReactsAndProds
	return ReactsAndProds

def ComputeClosure(Food, RandP_sets):
	closure_set = Food
	
	for reactiontuple in RandP_sets:
		#print reactiontuple[0], reactiontuple[1]
		if reactiontuple[0].issubset(closure_set) and not reactiontuple[1].issubset(closure_set):
			closure_set.update(reactiontuple[1])
			#print "The closure set is now:", closure_set
			#print "-----------------------------------"
	return closure_set

food = FindFood()
RPsets = FindReactantsProducts(food)
ComputeClosure(food, RPsets)
