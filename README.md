# PhyloToL v4.1

PhyloToL: Phylogenomic Tree of Life

Manual: https://github.com/Katzlab/PhyloTOL/blob/master/manual_phyloToLv4.1.pdf

Questions, suggestions and bug reports can be sent to Mario Cerón-Romero, mceronromero@umass.edu

**Versions control**

v4.1:
- Fixed a bug in Class Taxon. This bug cased an indexing error when parsing "added taxa" (the ones in the folders ncbiFiles and BlastFiles). Because of this a few sequences were ignored each run despite having good e-values. 
- Some minor changes in structure and logs.
- New option for cleaning intermediary files (see manual) 

v4:
- Translated to python 3 
- Major changes in databases
- Two new methods for contamination cleaning were added. Previously we only performed contamination cleaning in added taxa (ncbiFiles folder). The two new methods allow cleaning contamination in only or also in the seed dataset (allOG5Files folder)
- Minor bugs were fixed

v3.1: 
- Neddle step and ingroup paralogs removal were replaced. Now we have an overlap filter (OF) and a similarity filter (SF). Both are performed using Usearch-Ublast. We removed this step from the taxon class and added as a separated script called 'iterUblast.py'. This script is called from the module 'Utilities'.
- The logic now is... all sequences per taxa should be 1.5 times smaller tha the average OG legth and pass the overlap filter. Then, the users specify if they want to run a similarity filter for the sequences.

v3.0:
- A new Guidance method: This is a bash script made by Miguel Fonseca
- A "helper" for the new Guidance method. This helper is a perl script that edits the intermediary files of each guidance run in order to allow the looping.
- A new module called "Utilities" with a bunch of extra function including functions to integrate the two components described in the previous points.  

Other version, here: https://github.com/Katzlab/PhyloToL_archive
