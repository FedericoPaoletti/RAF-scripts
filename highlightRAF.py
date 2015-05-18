def highlight_RAF(cat_RAFset):
	
	with open("network_output1.txt", 'r') as network_file, open("RAF_output1.txt", 'w') as RAF_dotfile:
		for reaction_list in cat_RAFset:
			reaction_list[1] = list(reaction_list[1]) #convert sets back to lists!
			reaction_list[1][0] = list(reaction_list[1][0])
			#print "Reactants of this reaction:"
			#print reaction_list[1][0]
			reaction_list[1][1] = list(reaction_list[1][1])
			cat_edge = reaction_list[0] + " -> " + "R%d" % reaction_list[2] + " [color=red];\n"
			print cat_edge
			for line in network_file:
				#print cat_edge in line
				if cat_edge in line:
					new_line = line.replace("red", "blue") #change edge from red to blue because it is now part of the RAF set
					RAF_dotfile.write(new_line)
				#else:
					#RAF_dotfile.write(line)					

			for line in network_file:
				for reactant in reaction_list[1][0]:
					reactant_edge = reactant + " -> " + "R%d;" % reaction_list[2]
					#print reactant_edge
					if reactant_edge in line:
						clean_line = line.replace(";", "")
						new_line = clean_line + " [color=blue];\n"
						RAF_dotfile.write(new_line)
				  
			for line in network_file:
				for product in reaction_list[1][1]:
					product_edge = "R%d" % reaction_list[2] + " -> " + product + ";"
					print product_edge
					print line
					print product_edge in line
					if product_edge in line:
						clean_line = line.replace(";", "")
						new_line = clean_line + " [color=blue];\n"
						RAF_dotfile.write(new_line)
					else:
						print "Line written:"
						print line
						RAF_dotfile.write(line)
