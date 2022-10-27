This repository accompanies the manuscript "Assisted tree migration can reduce but not avert the decline of forest ecosystem services in Europe" by Mauri et al. and provides all necessary python scripts and data to reproduce included results.

DATA
1) tree_species_list.csv
This is a list of the tree species used to run the analysis.
2) traits.csv
List of traits used in the multi optimization strategy
3) services_new.csv
List of tree species with associated forest ecosystem services presented in a binary format

RASTERS input

4) EU-Trees4F_dataset.zip
Tree species datasets as taken from the EU-Trees4F
5) mask_eu.tif
Mask used to trim the data to the chosen spatial domain

SCRIPTS
6) preproccessing_tiles.py
This script rearranges the initial raster files in a proper format
7) script_single.py
This script is used to perform the following adaptation strategies: Rand, Maxi optimization and No AM
8) script_multy.py
This script is used to perform the multi optimization strategy:
-To run the scripts from terminal in multiple cores in tiles:   for i in $(seq 25); do    echo $i;    python2 script_multi.py rcp85 & done
-to run the single scripts:	python2 script_single.py rcp85
To use the scripts in Spyder, the RCP variables need to be defined within the script.
9) merge_tiles.py
This script merges the results of the script_multi.py, which was run on a high-computing cluster
