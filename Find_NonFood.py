with open('network_output1.txt', 'r') as Network, open('non_food_molecules.txt', 'w') as non_food:
    
    food = ['"N"', '"O"', '"C=O"', '"C#N"', '"C#CC#N"']
    non_food_list = []
    for line in Network:

        if "R" == line[2]:
            break

        if "M" == line[2]:
            okflag = 1
            for mol in food:
                if mol in line:
                    okflag = 0

            if okflag == 1:
                molecule_number = line[3:5].strip()
                if molecule_number not in non_food_list:

                    non_food.write(molecule_number + '\n')
                    non_food_list.append(molecule_number)
    #print RAF_molecules
