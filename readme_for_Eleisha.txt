Last updated by:   Amir Shahmoradi, 11:23 AM, September 10 2015, iCMB, UT Austin
 
Relative Solvent Accesibility:

    -- 213 enzymes:
       available at "cordiv/properties/res_prop_dssp.out". Columns have the following meaning: 1. pbd name   2. residue name   3. residue number as they appear in the original PDB file.  4. residue Accessible Surface Area (asa)  5. Residue Relative Solvent Accessibility (RSA)   6. mean Hydrogen Bond Energy of the residue (hbe)   7. The type of the local (secondary) structure in which the residue resides (refer to DSSP manual for meaning of different letters).
    
    -- viral proteins (contains also hyperthermophilic structures that are not relevant to ASAP project, ignore them):
       available at cordiv/properties/res_prop_dssp_asap.out". The header has exact same format as in the enzymes file.
        
Weighted Contact Number and B factors:
    
    -- 213 enzymes:
        available at "cordiv/properties/res_prop_wcn_bf.out". Columns have the following meaning: 1. pbd name   2. residue name   3. residue number as they appear in the original PDB file.  4. the number of the residue side chain atoms in PDB file (sizeSC)  5. the number of the amino acid side chain atoms in PDB file  (sizeAA). This should be typically always greater than sizeSC by 4 counts, unless a backbone atom had not been resolved in the crystall structure.   For columns 6-19: "wcn" stands for Weighted Contact Number, "bf" stands for the B factor, and the capital letters at the end of each column name, indicate the type of representative atoms for which WCN and B factors were calculated:
            -- CA: backbone C_alpha carbons
            -- CB: sidechain C_beta carbons
            -- N: backbone N (Nitrogen)
            -- O: backbone O (Oxygen)
            -- C: backbone C (Carbon)
            -- wcnSC: WCN calculated using the average of the side-chain coordinates (THIS IS THE ONE YOU SHOULD USE FOR ALL YOUR WCN CALCULATIONS)
            -- bfSC: bfSC calculated using the average of the side-chain B factors (THIS IS THE ONE YOU SHOULD USE FOR ALL YOUR B factor CALCULATIONS)
            -- wcnAA: WCN calculated using the average of all Amino Acid (AA) atoms coordinates.
            -- bfAA: bfSC calculated using the average of the side-chain B factors.
            
    -- viral proteins (contains also hyperthermophilic structures that are not relevant to ASAP project, ignore them):
       available at cordiv/properties/res_prop_wcn_bf_asap.out". The header has exact same format as in the enzymes file.
       
Voronoi Cell Properties:

    -- 213 enzymes:
       available at "cordiv/properties/res_prop_voronoiSC.out". Note that all data in this file were calculated using the average Side Chain (SC) coordinates, which is also indicated in the naming of the file by "SC". Data for other atomic coordinates is also available in the same directory in separate files.
       Columns have the following naming convention: 
       1.  pdb name
       2.  residue name
       3.  resname: residue name as it appears in PDB file.
       4.  resnum: residue number as it appears in PDB file.
       5.  sizeSC: the total number of Side Chain heavy atoms in the residue
       6.  sizeAA: the total number of Amino Acid (AA) heavy atoms including backbone
       7.  resvol: Residue Volume (excluding backbone)
       8.  VSCnvertices: the total number of vertices (corners) of the Voronoi cell for the corresponding residue.
       9.  VSCnedges: the total number of edges of the Voronoi cell for the corresponding residue.
       10. VSCedge_length_total: the total length of all edges of the Voronoi cell for the corresponding residue.
       11. VSCedge_length_total: the total length of all edges of the Voronoi cell for the corresponding residue.
       12. VSCnfaces: the total number of faces of the Voronoi cell for the corresponding residue.
       13. VSCarea: the total surface area of all faces of the Voronoi cell for the corresponding residue.
       13. VSCvolume: the total volume of the Voronoi cell for the corresponding residue.
       14. VSCeccentricity: The eccentricity of the Voronoi cell. (a rough measure of how assymetrical the cell looks like)
       15-20. Ignore the rest of columns: VSCfree_volume	VSCvolume_change_diff	VSCvolume_change_ratio	VSCarea_change_diff	VSCarea_change_ratio

    -- viral proteins (contains also hyperthermophilic structures that are not relevant to ASAP project, ignore them):
       available at cordiv/properties/res_prop_voronoiSC_asap.out". The header has exact same format as in the enzymes file.

            

        
    
    