import numpy as N

with open('network_output1.txt', 'r') as non_highlighted_Network, open('RAFhighlight.txt', 'r') as RAF_highlighted_Network, open('Catalysts.txt', 'w') as catalysts_file:

    catalyst_list = []
    cat_react_list = []
    for line in RAF_highlighted_Network:
        if "style=dotted" in line:
            catalystNumber = line[3:5].strip()

            if "R" == line[8]:
                reactionNumber = line[9:11].strip()
            elif "R" == line[9]:
                reactionNumber = line[10:12].strip()

            if catalystNumber not in catalyst_list:
                catalysts_file.write(catalystNumber + '\n')
                catalyst_list.append(catalystNumber)

            cat_react_list.append([int(catalystNumber), int(reactionNumber)])

    #print catalyst_list

    reactions = 0
    molecules = 0

    reaction_string = "; // "
    non_highlighted_Network.seek(0)
    for line in non_highlighted_Network:
        if reaction_string in line:
            reactions += 1 #count reactions
    #print reactions

    non_highlighted_Network.seek(0)
    for line in non_highlighted_Network:

        if "M" == line[2]: #count molecules
            molecules += 1
        if "R" == line[2]:
            break
            
    catalysis_matrix = N.zeros([molecules, reactions])
    
    for i in range(reactions):
        for cat_pair in cat_react_list:
            if i == cat_pair[1]:
                catalysis_matrix[cat_pair[0], i] = 1

    integer_catalysis_matrix = catalysis_matrix.astype(int)
    #print integer_catalysis_matrix.shape
    N.savetxt("catalysis_matrix.txt", integer_catalysis_matrix, fmt='%d')
        
        
