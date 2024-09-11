import numpy as np
import scipy as sp


################################# Selected atom pair-wise distances #######################################

###########################################################################################################

def SelectedAtomsPairwise_distance(struc, idxs_1, idxs_2, eigenvectors):

    # Atom positions 
    atom_1_pos = struc[idxs_1]
    atom_2_pos = struc[idxs_2]

    # Feature vector 
    feature = np.diag(sp.spatial.distance.cdist(atom_1_pos, atom_2_pos))

    # Calculate coordinates on projection planes
    coords = np.zeros(len(eigenvectors))
    for i in range(len(eigenvectors)):
        coords[i] = np.dot(eigenvectors[i], feature)

    return coords

#############################################################################################################




######################################## All CA-CA distance #################################################

#############################################################################################################

def All_CACAdistance(struc, idxs, eigenvectors): 

    # Atom positions
    atom_pos = struc[idxs]

    # Feature vector
    feature = []
    for i in range(len(idxs) - 3):
        for j in range(i+3, len(idxs)):
            feature.append(np.linalg.norm(atom_pos[i] - atom_pos[j])) 
    feature = np.array(feature)

    # Calculate coordinates on projection planes
    coords = np.zeros(len(eigenvectors))
    for i in range(len(eigenvectors)):
        coords[i] = np.dot(eigenvectors[i], feature)

    return coords

#############################################################################################################
