# Katzlab_phylopipeline_v3_scripts
third version of the Katzlab phylogenomic pipeline

This is the scripts folder of the pipeline. It needs to be in same path as the DataFiles folder. Here is short descrition of scripts:

* pipeline_parameter_file.txt: Here the user sets parameters for filters, guidance and overall pipeline before starting analysis
* phylopipe_3.py: The wrapper script for the pipeline
* Pipeline/__init__.py: The pipeline itself 
* Taxon/__init__.py: The TAXON class. Here is where sequences for each taxon are explored to extract OGs. Here is where filters (overlap and similarity) also occur. 
* iterUblast.py: Runs filters (similarity and overlap) during taxon step
* Gene/__init__.py : The GENE class. Once all sequences per taxon are already collected, this class continues to produce alignments and trees
* 03-guidance_v1.2c.pl: Runs guidance iterations during GENE step
* exeGuidance_2.2.sh Edits results per guidance iteration during GENE step 
* Utilities.py: This is a module that contains functions needed for many parts of the pipeline


