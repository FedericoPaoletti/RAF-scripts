with open('RAFhighlight.txt', 'r') as RAF_highlighted_Network, open('Catalysed_Reactions.txt', 'w') as cat_reactions, open('Food.txt', 'w') as FOODlist:

    cat_set = []
    for line in RAF_highlighted_Network:
        if "style=dotted" in line:
            if "R" == line[8] and int(line[9:11]) not in cat_set:
                cat_reactions.write(line[9:11] + '\n')
                cat_set.append(int(line[9:11]))

            if "R" == line[9] and int(line[10:12]) not in cat_set:
                cat_reactions.write(line[10:12] + '\n')
                cat_set.append(int(line[10:12]))

    smile_list = []
    food_list = []

    foodsmile = raw_input('Please insert a food SMILE within \" \".')
    while foodsmile != "":
        smile_list.append(foodsmile)
        foodsmile = raw_input('Please insert a food SMILE within \" \".')
    with open("network_output1.txt", 'r') as dot_file: #Extract the unique identifiers "M + integer" for each food molecule
        for smile in smile_list:
            dot_file.seek(0)
            for line in dot_file:
                if smile in line:
                    food_list.append(line[2:5]) #Be careful in changing indentation of input file (.txt file)!!!
    clean_list = []
    for food in food_list:
        cleanfood = food.strip()
        clean_list.append(cleanfood)
    
    #print clean_list

    for food in clean_list:
        foodnumber = food[1:]
        FOODlist.write(foodnumber + '\n')

with open('RAFhighlight.txt', 'r') as RAF_highlighted_Network, open('RAF_molecules.txt', 'w') as RAF_molecules:
    
    for line in RAF_highlighted_Network:
        if "[color=red]" in line:
            if "M" == line[2] and line[2:5].strip() not in clean_list:
                RAFmolecule_number = line[3:5]
                RAF_molecules.write(RAFmolecule_number + '\n')
                clean_list.append(line[2:5].strip())

            if "M" == line[8] and line[8:11].strip() not in clean_list:
                RAFmolecule_number = line[9:11]
                RAF_molecules.write(RAFmolecule_number + '\n')
                clean_list.append(line[8:11].strip())

            if "M" == line[9] and line[9:12].strip() not in clean_list:
                RAFmolecule_number = line[10:12]
                RAF_molecules.write(RAFmolecule_number + '\n')
                clean_list.append(line[9:12].strip())


    #print RAF_molecules
