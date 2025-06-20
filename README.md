# DILS-Input
This pipeline aims to help people generate the input file required for running DILS (Demographic Inferences with Linked Selection by using ABC). <br>

Refer to the Guide_DILS_input_generation.pdf file for a detailed description of each required step. Basically, the entire pipeline:<br>

1- Generates a list of locus coordinates, given a genome, a set number of loci, and a set locus length <br>
2- Generates, fo each individual, two pseudo-haplotypes for each locus, with bases with low base quality, mapping quality and read quality coded as Ns.<br>

The pipeline provides some scripts to check the correct formatting of the output fasta file, which is in the exact format required by DILS. 
