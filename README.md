# RADseq tutorial

This tutorial demonstrates a few way to get and filter genotypes using restriction-site associated DNA. This readme is currently skeletal, but the code is more or less complete and demonstrates how to:

1. Get genotypes using `stacks` in _de novo_ mode, step by step. 
2. Get genotypes using `stacks` in reference mapping mode. 
3. Get genotypes using the variant caller `Freebayes` with reference-mapped sequence. 
4. Filter and compare variant call sets. 
5. Do some basic QC-oriented analysis in `R`. 

Each script should be run from the `/scripts` directory and will create its own output directories as necessary. 

The dataset is from a fish, the Arctic grayling, and includes several hundred individuals distributed across three watersheds on Alaska's North slope. It has some fun issues. A bunch of samples are actually from another species. Graylings are salmonids and so have an old-ish genome duplication with residual tetrasomy which leads to some issues with paralogs. 

Unfortunately, the data are unpublished at the moment, and so are only available to users of UConn's Xanadu cluster. We will update the code to retrieve the data from NCBI once they are deposited there. 

The scripts are written to use SLURM, as configured for UConn's Xanadu cluster, but if you are at all comfortable with bash coding or SLURM you should be able to adapt them to your home server or your laptop, for your own dataset. If you're uncomfortable with those things, picking apart and figuring out code like this is a great way to learn!

If you are using UConn's Xanadu cluster, you should be able to clone this repository and simply run the scripts sequentially to complete the analysis. If you are a workshop participant from outside the University, you can also clone and run the scripts sequentially, but you will need to update the partition to "mcbstudent" instead of "general". 