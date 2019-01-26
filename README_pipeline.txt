#PIPELINE of project "Predicting protein-protein interactions in the Solanum lycopersicum - Phytophthora infestans pathosystem"
#by Daniel Del Hoyo Gomez. Supervisors: Dick de Ridder, Sander Rodenburg
#Explanation of scripts and processes

1) Preparation of interactome
#Extraction of S. lycopersicum PPI network from STRING database
Input file (1): 4081.protein.actions.v10.5.txt.gz (predicted interactions)
Input file (2): 4081.protein.sequences.v10.5.fa.gz (protein sequences)

Scripts:    
    deredundant_actions.py : Extract the non redundant interactions from (1), specifying a score threshold and a interaction type
        output file (3): nr_4081_<interactio_type>_<score_threshold>.txt
    filter_fasta.py : Extract the protein sequences of proteins in (3) from (2)
        output file (4): filt_<score_threshold>_4081.protein.sequences.v10.5.fa
        
2) Domain annotation
#Domain annotation of protein sets with InterProScan software
#Requires InterProScan script and associated files
Input file (4): Fasta of pitg secretome and soly interactome

Scripts:
    domain_annotator.py:    Run the InterProScan software over a Fasta file with the specified databases
                            Creates a histogram for the number of domains per protein
        Output file (5):    (4)_domains.tsv (Domain annotation)
        Output file (6):    (4)_histo.csv (Histogram) 

3) Motif database creation
#Creation of a motif database using iELM and 3did databases
Input file (7): 3did_dmi_flat.gz (3did motifs involved in DMIs information)
Input file (8): elms_index.tsv (ielm motifs involved in DMIs information)

Scripts:
    3did2motif.py: Extract motifs from 3did database file (7)
        output_file (9): 3did_motifs.txt    (Motif id followed by regular expression)
    #ielm were extracted in bash and '"' replaced: cut -f 2,5 (8) > elm_motifs.txt 
        output_file (10): elm_motifs.txt    (Motif id followed by regular expression)
    #Final motif database was built concatenating (9) and (10): cat (9) (10) > motifs_database.txt
        output_file (11): motifs_database.txt

4) Surface exposure calculation
#Surface exposure prediction of protein sets 
#Requires NetSurfP software and associate files (slow)
Input file (4): Fasta of pitg secretome and soly interactome

Scripts:
        netsurfp: Run the netsurfp software (./netsurfp -i (4) -o (12))
            output_file (12): netsurfp output file
        redo_netsurf.py: Run the netsurfP software only for those proteins not found in a current netsurfp output
            #Slow program, this script avoid repeating the analysis in proteins already found
            output_file (12*): netsurfp output file

5) Motif annotation
#Motifs annotation of protein sets searching their regular expressions
Input file (4):  Fasta of pitg secretome and soly interactome
Input file (11): Motif database
Input file (12): NetSurfP output

Scripts:
    motif_annotator.py: Annotate the motifs in protein sets (4) and creates a histogram. 
        #Able to use exposure filter ((12) required) and frequency filter
        #Accept a (13) file to perform over it further filters
        output_file (13): (4)_motifs.fa[_filt][_<frequency threshold>].tsv (Motif annotation file)
        output_file (14): (13)_histo.csv (Histogram of motifs per protein)
    motif_frequency.py: Create a file with the frequency of each motif in a (13) file
        output_file (15): (13)_freq.csv (Frequncy of motifs file)
        
6) Domain Domain databases creation
#Creation of a DDI database from 3did files
Input file (16): DDI_3did.txt (Uncompressed: gzip -cd 3did_flat.gz > DDI_3did.txt )

Scripts:
    build_DDI_database.py: Adapt the format of the 3did database
        output_file (17): DDI_database.txt
        
7) Domain Motif database creation
#Creation of a DMI database from ELM and 3did files. Also DMDI database
Input file (18): elm_interaction_domains.tsv (DMI_ELM.txt)
Input file (19): DMI_3did.txt (Uncompressed: gzip -cd 3did_dmi_flat.gz > DMI_3did.txt )
Input file (20): Pfam_dictionary.txt (Mapping of 3did-Pfam identifiers)

Scripts:
    build_DMI_database.py: Build the DMI database with DMI found in (18), (19)
        output_file (21): DMI_database.txt 
    mix_DMDI.py: Build the DMDI database mixing (17) and (21), specifies score for DMI
        output_file (22): DMDI_<DMI_score>_database.txt

8) Training Dyer method
#Uses the Dyer method to score DDI(DMDI) from (3) files
Input file (3): Host (tomato) PPI network
Input file (5): Domain annotation of host proteins
Input file (13): Motif annotation of host proteins

Scripts:
    train_dyer.py: Uses the Dyer method to calculate DDI scores
    #Addition of (13) file or not for training without or with motifs
        output_file (23): trained_dyer_host_...tsv (Dyer DDIs)

9) Running Dyer method
#Uses (23) to calculate the PPI scores of (4) sets
Input file (23): Dyer DDI scores
Input file (5): Host AND pathogen domain annotations
Input file (13): Host AND pathogen motif annotations

Scripts:
    run_dyer.py: Uses Dyer DDI scores to calculate PPI scores
        output_file (24): dyer_interact_...tsv
        
10) Running Zhang (Xiaopan) method
#Uses (17) or (22) to calculate PPI scores of (4) sets
Input file (17): DDI scores database
Input file (23): DMDI scores database
Input file (5): Host AND pathogen domain annotations
Input file (13): Host AND pathogen motif annotation

Script:
    run_xiaopan.py : Uses database DDI scores to calculate PPI scores
        output_file (25): xiaopan_interact_...tsv

11) Location filter
#LOCALIZER software is used to filter for colocalized pairs of proteins
#LOCALIZER scripts and associated files are needed
Input file (4): Protein sets
Input file (24)/(25): PPI predicted file

Scripts:
    LOCALIZER.py: Run the localizer software (python LOCALIZER.py -e/-p -i (4)) #-e: effector; -p: plant
        output_folder (26*): Folder with several files saying the location
        output_file (26): Results.txt in (26*)
    ids_location.py: Use (26) to create a file with sum of location of each protein
        output_file (27): ids_location(_cyto).txt
    location_filter.py: Actual filter that uses (27) to filter location of (24)/(25)
        output_file (28): (24)/(25)_locfilt(cyt).tsv (PPI scored file filtered by location)

12) Dyer-Zhang comparations
#Compare DDI databases and PPI results
Input file (17) or (22): DDI/DMDI databases scores
Input file (23): DDI/DMDI dyer scores
Input file (24): dyer PPI scores    #or (28)
Input file (25): zhang PPI scores   #or (28)

Scripts:
    compare_traineds.py: compare (17)/(22) with (23)
        output_file (29): trcomp_....csv    (comparation of DDI scores)
    compare_interactomes.py: compare (24) with (25)
        output_file (30): comp_.....csv     (comparation of PPI scores)

13) Creating best subsets
#Take the best interactions with the sum of Dyer and Zhang scores
Input file (30): comparation of PPI scores
Input file (27): location of proteins by localizer

Scripts:
    get_x_better.py: take the x greater sum of scores froma a (30) file
        output_file (31): better_x_(30).tsv (Best subset of predicted interactions)
    add_location.py: add the location predicted to the pairs in the (31) file
        output_file (32): (31)_loc.tsv  (best subset with localization)
    make_pathointeractome.py: concatenates (3) with (31) 
        output_file (33): inter_(31)    (Best interactions over host interactome) 

14) Translating STRING ids to UNIPROT
#Translate the identifiers of STRING proteins to UNIPROT for posterior GO annotation
Input file (3): Host interactome
Input file (31): Best subset of predicted interactions
Input file (34): translations file from UNIPROT
Input file (34*): more translations from other sources (different format)

Scripts:
    translate_ids.py: translate the ids of the (3) or (31) files using (34)
        output_file (35)/(36): uni_(3)/(31)

15) Proximity of targets
#Perform the evaluation of proximity of common targtes in the host PPI
Input file (3): Host interactome
Input file (31): Best subset of predicted interactions

Scripts:
    proximity_targets.py: perform the analysis
        output_file (37): proximity_(31)    (P-value for each effector with several targets and fisher combined p-value)

16) Gene expression correlation
#Perform the analysis of gene expression correlation 
Input file (31): Best subset of predicted interactions
Input files (38): Gene expression data of host AND pathogen

Scripts:
    correlation_common_effectors.py: perform the analysis
        output_file (39): expression_(31)   
        #Pearson coefficients, Pearson pvalues and adjusted p-values for each target-effector and effector-effector

17) Known interactions
#Check the score of known interactions in predicted PPIs (28)
Input file (40): Known interactions from Whisson et al 2016 translated to UNIPROT

Script:
    compare_known.py: check scores of known interactions in (28) files
        output_file (41): Results_sp.tsv 
        #Score of each known interaction in different (28) files

18) GO terms related with infection
Input file (40): Known interactions
Input file (42): GOA files from UNIPROT (for GO annotation of UNIPROT proteins)
Input file (43): Obo file from GO (parent-child,branches relations between GO terms)
Input file (45): GO functions from GO

Scripts:
    get_go_terms.py: get go terms from a (40) file (or other interaction files)
        output_file (46): infection_gos_sp_known.tsv (Infection related GO terms)
        #GO terms found in known targets
    parse_obo.py: parse the obo file and creates a parent-child file and a branches file
        output_file (47): children_gos.tsv
        #Parent-child relations between GO terms
        output_file (48): GO_branhces.tsv 
        #Branch of GO for each GO term
    filter_go_type.py: filter the (46) file for maintain only GOs in a branch (biological process)
        output_file (49): bp_(46) (Infection related GO terms filtered by branch)
    add_all_children.py: add children of all go terms in a (49)
        output_file (50): (49)_+children.tsv

19) GO terms enrichment in the predicted targets
Input file (35): Uniprot host interactome 
Input file (36): Uniprot best subset of predicted interactions
Input file (42): GOA files from UNIPROT (for GO annotation of UNIPROT proteins) 
Input file (49): Infection related GO terms filtered by branch
Input file (50): Infection related GO terms filtered by branch + children

Scripts:
    go_enrich.R: Performs the Fisher exact test in Python script
    go_enrichment_infection.py: perform the GO enrichment analysis
        output_file (51): bp_mf_goenrich_(36)(_+children.tsv)
        #Enriched GO terms with function, pvalue and yes or no related
        output_file (52): bp_mf_goenrich_(36)(_+children)_stats.tsv
        #Statistics for the infection related classification


