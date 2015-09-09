Last updated by:   Amir Shahmoradi, 9:36 AM, September 9 2015, iCMB, UT Austin
 
    This folder contains data and structural analysis for 213 monomer enzyme proteins. The ultimate goals of the project are:

        -- 
        
        finding the potential variables that create the observed diversity in the correlations of the structural variables with sequence entropy, for example, the wide range of correlations observed in RSA-seqent or WCN-seqent among different proteins.

Contents/

    analysis/
        contains the extensive and comprehensive analysis of relationships of different structural properties with each other and with sequence variability (sequence evolutionary rates and sequence entropy)
        
    cluster_analysis/
        contains data and analysis for finding the nearest neighbour distances (nnd) of amino acids in individual proteins, using Voronoi partitioning, also the distribution of sequential side-chain and backbone C_alpha distances in individual proteins.
        
    contact_order_definition/
        will eventually contain data and analysis for finding the best and optimal parameter-free definition of Contact Order.
    
	ddG_calculations/
		contains all of the files, scripts, and resulting data from the FoldX ddG calculations

    structures/ 
    	contains all repaired pdb structures in the data set of Echave et al 2014

    manuscript/
        contains the preprints of the manuscripts for publication. The complete versions are currently available on OverLeaf website.
        
    rate_calculations/
        contains data and analysis for the calculation of the site-wise evolutionary rates and sequence entropy for all enzyme proteins. This includes all of the files, scripts, and resulting data from the rate4site evolutionary rate calculations
        
    wcn_best_definition/
        contains analysis and results for obtaining the best definition of WCN kernel in relation to other structural and sequence quantities, in particular the evolutionary rates and B factors.
        
    wcn_definition/
       same as in the directory "wcn_best_definition", but contains older, incomplete, unstructured analysis and results.
    
    properties/
    	contains the pdb-specific & residue-specific properties of all pdb structures, each in a single file.
