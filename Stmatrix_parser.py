import numpy as N

with open('network_output1.txt', 'r') as reactionNetwork:
    molecules = 0
    reactions = 0
    for line in reactionNetwork:

        if "M" == line[2]: #count molecules
            molecules += 1
        if "R" == line[2]:
            break
    
    reaction_string = "; // "
    reactionNetwork.seek(0)
    for line in reactionNetwork:
        if reaction_string in line:
            reactions += 1 #count reactions

    #print molecules
    #print reactions

    stoich_matrix = N.zeros([molecules, reactions]) #generate an empty m*n matrix, where m = molecules and n = reactions
    #print stoich_matrix
    #print stoich_matrix.shape
    reactionNetwork.seek(0)

    for line in reactionNetwork:
        if line == "\n": #If you arrive at EOF, break the loop!
            break
        #print line
        #print reaction_string
        #print reaction_string in line

        if reaction_string in line:
            #print reaction_string in line
            reaction_index = line[3:5] #identify the reaction number
            #print reaction_index

        if "M" == line[2] and "shape=oval" not in line: #need not in statement to skip the molecule declaration lines
            stoich_matrix[int(line[3:5]), reaction_index] -= 1 #If you find a reactant, decrease its stoichiometric matrix value by 1 at the corresponding reaction column

        if "M" == line[8]:
            clean_line = line.replace(";", " ")
            stoich_matrix[int(clean_line[9:11]), reaction_index] += 1 #If you find a product, increase its stoichiometric matrix value by 1 at the corresponding reaction column

        if "M" == line[9]:
            clean_line = line.replace(";", " ")
            stoich_matrix[int(clean_line[10:12]), reaction_index] += 1 

    integer_matrix = stoich_matrix.astype(int)
    #print integer_matrix
    N.savetxt("stoichiometric_matrix.txt", integer_matrix, fmt='%d')
    #print stoich_matrix[:,:]
        
