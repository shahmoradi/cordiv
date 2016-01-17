! Amir Shahmoradi, Sunday 5:58 PM, March 29 2015, iCMB, UT Austin

wcn naming convention:
  wcnXY --> Weighted Coordination Number (Contact Number) using model X (e: exponential, h: Heaviside step function, p: Power-law, g: Gaussian) usinf representative coordinates Y (N, C, CA, CB, O, SC, AA)


OUTPUT FILES:
    
    There are two output files from each of the analyis in the subdirectories:
    
    -- exp_wcnXY_Z:  The rule for the naming of the file is the folowing: 
        -- exp: stands for what type of data the file contains: the free parameter of WCN (exponent), for example exponent of power-law kernel for WCN, although in general the free parameter is not necessarily an exponent. Bad naming!).
        -- wcnXY:  stands for the type of wcn kernel analyzed here in this file. X stands for the type of kernel: e --> exponential kernel,  g --> Gaussian kernel, h --> Heaviside step function kernel (simple cutoff model, equivalent to Contact Number (CN)), p --> power-law kernel (the most widely used definition of WCN).  Y stands for the set of residue-representative coordinates used in the analysis. SC --> average Side-Chain coordinates, AA --> average Amino Acid coordinates, CB --> CB (beta carbon) atomic coordinates, CA --> CA (alpha carbon) atomic coordinates, and so on for other backbone atoms.
        -- Z stands for the parameter that was corrolated with WCN definition of interest.  r4s_JC stands for the evolutionary rates calculated according to Juke-Cantor model, vorvol stands for The Voronoi cell volumes.  ddgent stands for ddG entropies calculated according Echave, Jackson, Wilke (2015) work. and so on for other definitions.
        
        The structure of the file:
        first line contains different values of the free-parameter of WCN kernel used. For example, for power-law kernel, the exponent ranges from -30 to 30.
        second column and the followings contain the Spearman correlation coefficient of the specific deifinition of WCN used with the quantity of interest (evolutionary rates, b factor, ...).
        
        first column of the file gives the protein ID for the data given on that row.
        
    -- sum_wcneSC_r4sJC:  this is the summary file, containing the most important information from exp_wcnXY_Z files, summarized for each protein.  The headings of this file are sufficiently self-explanatory.
    
    
    For any further information please contact amir@physics.utexas.edu,  amir@ices.utexas.edu,  or shahmoradi@utexas.edu, or a.shahmoradi@gmail.com  if any of the former addresses are useless in reaching me.
        
                    
        
  
