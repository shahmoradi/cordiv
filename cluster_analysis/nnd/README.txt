Last updated by:   Amir Shahmoradi, 3:31 PM, April 31 2015, iCMB, UT Austin
 
    This folder contains codes and analysis for finding the nearest neighbour distances (nnd) of each of amino acids in all proteins.
    
    The algorithm searches for the nearest neighbor of each site, by simply calculating the distances to all sites in individual proteins.
    
    Note that for larger datasets, if computational time becomes relevant, this information (nnd) could be instead extracted from the list of nnd for individual sites using Voronoi tessellation.
    
    The algorithm also calculates the sequential neighbor distance (snd) of each site in protein from the next site in protein sequence. Since the last site in sequence does not have a next neighbor, it is assigned an snd value of -1. These negative values can be then discarded in R postprocessing.