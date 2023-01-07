#terminal commands for motif analysis
#RBP_database (Ray2013_rbp_All_Species.meme) from http://cisbp-rna.ccbr.utoronto.ca/index.php

streme --verbosity 1 --rna --totallength 4000000 --time 14400 --minw 7 --maxw 15 --thresh 0.05 --align center --p input.fa -oc outputdir

tomtom -no-ssc -oc C9ALSnoFTLD_Microglia_tomto_res -verbosity 1 -min-overlap 5 -dist pearson  -eps -evalue -thresh 10..0 input/streme.xml 