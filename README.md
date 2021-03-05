# Phytochemical-Diversity-Experiment

These data and code are to accompany the manuscript: “Interaction diversity explains the maintenance of phytochemical diversity”
by Susan R. Whitehead, Ethan Bass, Alexsandra Corrigan, André Kessler, and Katja Poveda

Publication information: Accepted for publication in Ecology Letters, 02/26/21
DOI: TBD

The main purpose of this study was to test different hypotheses for how phytochemical diversity benefits plants in interactions with herbivores and pathogens. We conducted a series of bioassays with insects and fungi in which we manipulated different aspects of the chemical diversity of the growth media (including the richness, evenness, and structural complexity of the mixture) and recorded the performance of the organisms.  All data were collected at Cornell University in the Poveda Lab or at Virginia Tech in the Whitehead Lab between 2016-2019.

## Data files

There are five total data files:

### 1) Whitehead_et_al_CIDs.csv 
This is a full list of 21 compounds we considered for use in the study based on their co-occurrence in apple fruits. 

Compound: the name of the compound
CID: the PubChem Compound ID (CID) which is a a unique identifier for each molecule and can be used to retrieve structure data on PubChem

### 2) Whitehead_et_al_structures.sdf.gz
This is an output file from PubChem that contains the structural information for the compounds in "Whitehead_et_al_CIDs.csv". Additional notes on retrieving these data are provided in the script "Script_MixtureSelection.R"

### 3) Whitehead_et_al_Insects.csv
This is the main experimental data for the insect bioassays reported in the paper. Each row (replicate) represents a single individual insect reared as part of the experiment. Variables include:

Species: The species of insect; Hz=Helicoverpa zea, Sf=Spodoptera frugiperda, Cp=Cydia pomonella, Px=Plutella xylostella	

Treatment: A code indicating the specific chemical composition of the insect's diet. Diet treatments are coded either as an abbreviation of the compound name (see below) or with a code indicating the structural diversity (L=low, M=med, H=high), the richness (2-10 compounds), and the replicate randomly selected mixture at each level of structural diversity and richness (A, B, C). See the main text of the paper, Appendix S1a, and Table S1 for full details on these mixtures and how they were selected.

SD:	The structural diversity of the mixture; L=low, M=medium, H=high

Richness: The richness (number of compounds) of the mixture

TrayID: An identifier for the tray that held the cup where the insect was reared. Trays were 11" x 21" greenhouse seedling flats. Identifiers are given as the experimental round (6, 7, or 8) followed by a tray number for that round.

CupID: An identifier for the 1oz plastic deli cup in which the insect was reared. Hz and Sf were reared in individual cups, but Cp and Px were reared with multiple individuals per cup

Pupal.weight: For insects that survived to pupation, the final weight of the pupae (mg)

Sex: For insects that survived to pupation, the sex of the insect (sex could only be determined once insects pupated)

Surv: A binomial variable indicating whether the insect survived to pupation (1) or died (0). Note insects are only included in this spreadsheet if they survived the initial transfer to the diet and successfully began feeding, growing at least large enough to be easily seen as dead in the cup. So insects marked as dead here were insects that had begun feeding and development on the diet but then later died (i.e. we were able to locate a dead larvae in the cup)

Days.to.pupation: The time from initial hatching to pupation (days)

exp: The experimental round (6, 7, or 8). Earlier "rounds" (1-5) were conducted for method development and preliminary trials and were not included in the final dataset.

ChA: The % fresh mass of chlorogenic acid in the diet
CA:	The % fresh mass of caffeic acid in the diet
pCA: The % fresh mass of p-coumaric acid in the diet	
FA: The % fresh mass of ferulic acid in the diet	
GA: The % fresh mass of gallic acid in the diet	
SA: The % fresh mass of syringic acid in the diet	
GeA: The % fresh mass of gentistic acid in the diet	
Ct: The % fresh mass of catechin in the diet	
eCt: The % fresh mass of epicatechin in the diet	
R: The % fresh mass of rutin in the diet
H: The % fresh mass of hyperin in the diet	
Q: The % fresh mass of quercitin in the diet	
Phz: The % fresh mass of phloridzin in the diet	
Pht: The % fresh mass of phloretin in the diet

### Whitehead_et_al_Fungi.csv
This is the main experimental data for the fungal bioassays reported in the paper. Each row (replicate) represents a single individual fungal assay (one cell in a 96-well plate) included in the experiment. Variables include:

Fungi: The fungal species, which were TENTATIVELY identified based on morphology and disease symptoms. Botrys=Botryosphaeria dothidea (Botryosphaeriaceae; white rot), Collet=Colletotrichum sp. (Glomerellaceae; bitter rot), Penicillium=Penicillium expansum (Trichocomaceae; blue mold), and Sclerotinia=Sclerotinia sclerotiorum (Sclerotiniaceae; calyx end rot). 

Plate:	An identifier for the 96-well plate in which the assay was conducted, numbered sequentially for each species. 

Cell:	The cell number on the 96-well plate

CellID:	A unique identifier for the cell ID that included the species, the plate number, and the cell number

Treatment: A code indicating the specific chemical composition of the growth media. Treatments are coded either as an abbreviation of the compound name (see below) or with a code indicating the structural diversity (L=low, M=med, H=high), the richness (2-10 compounds), and the replicate randomly selected mixture at each level of structural diversity and richness (A, B, C). See the main text of the paper, Appendix S1a, and Table S1 for full details on how these mixtures were selected. These are the same treatment compositions used in the insect experiments.

Time:	The time at which the absorbance was measured on the plate. Each plate was measured four times 1=0 hrs, 2= ~24 hours, 3= ~48 hours, and 4= ~72 hours from the start of the assay.

Reading1:	The first of two duplicate absorbance readings taken on a plate reader at 600nm

Reading2: The second of two duplicate absorbance readings taken on a plate reader at 600nm


### Whitehead_et_al_Insects_Evenness.csv
These are data from a supplemental experiment that was conducted with insects only to test how compound evenness (i.e. the extent to which compounds were present in the mixture in equal abundance) affects insect performance. Each row (replicate) represents a single individual insect reared as part of the experiment. Variable descriptions are identical to those in "Whitehead_et_al_Insects.csv" with the exception of: 

Treatment: A code indicating the specific chemical composition of the insect's diet. Diet treatments are coded with a code indicating the structural diversity (L=low, M=med, H=high), the evenness (0.2-1), and the replicate randomly selected mixture at each level of structural diversity and evenness (A, B, C). See the Appendix S1c for full details on how these mixtures were selected.

Evenness: The evenness of the mixture, calculated using the Simpson’s Equitability index, where 0.2 indicates very low evenness and 1.0 indicates perfect evenness


## Analysis Scripts

There are three scripts:

### Script_Mixture Selection.R
This script was used to calculate the structural similarly of mixtures using ChemmineR and then semi-randomly select mixtures for inclusion in the study.

### Script_FinalAnalysesforMS.R
This script includes all statistical analyses. It is divided into sections (based on different predictions tested) that can be run independently.

### Script_FinalFigsforMS.R
This script includes code for all the figures in the manuscript and supplement, based on ggplot.
