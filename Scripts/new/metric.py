import numpy as np
import scipy as sp
from wepy.resampling.distances.distance import Distance
from new.features import SelectedAtomsPairwise_distance, All_CACAdistance



class Projection_Coord(Distance):
    
    def __init__(self, idxs, ref_pos, eigenvecs, feat):
        
        self._idxs = idxs
        self._pos = ref_pos
        self._eigenvectors = eigenvecs
        self._feat = feat
    
    def image(self, state):
        
        """
        Compute the 'image' of a walker state which should be some
        transformation of the walker state that is more
        convenient.
        
        Here 'image' is a projection of the walker state
        """
       
        # Get projections
        if self._feat == "sel_pair":
            proj_coord = SelectedAtomsPairwise_distance(state['positions'], self._idxs[0], self._idxs[1], self._eigenvectors) 
        elif self._feat == "allCA":
            proj_coord = All_CACAdistance(state['positions'], self._idxs,  self._eigenvectors)
        else:
            raise ValueError("Unrecognized feature selection")

        return proj_coord


class InvDist_FromRef(Projection_Coord):
    """
    Compute the inverse of distance between walker states from reference point.
    Also compute the distance between two walker states.
    'Images' are produced using the ProteinDistance.image method.
    """

    def image_distance(self, image1, image2):

        p1 = image1
        p2 = image2

        inv_dist1 = 1.0 / np.linalg.norm(p1 - self._pos)
        inv_dist2 = 1.0 / np.linalg.norm(p2 - self._pos)

        dist12 = np.linalg.norm(p1 - p2)


        return inv_dist1, inv_dist2, dist12

