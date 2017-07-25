*simpoolTE*
=======

SimpoolTE simulates the insertion and deletion of transposable elements in a population of chromosomes according to the neutral allele frequency distribution. SimpoolTE then uses pirs to simulate paired-end sequencing on these chromosomes. These simulated data can be used to test the accuracy of TE calling packages.

## Required packages

* [python](www.python.org) v2.7
* [pirs](https://github.com/galaxy001/pirs)
* [seqtk](https://github.com/lh3/seqtk)

## Testing 
Test files are provided to ensure simpoolTE and its dependencies are running correctly.

### Using simpoolTE
SimpoolTE uses the provided TE annotation to identify TE to insert and delete. If you do not have a TE annotation for your chromosome you can use RepeatMasker to generate one. TE insertions and deletions will be simulated at allele frequencies corresponding to the neutral frequency distribution. The identity of the TEs that will be inserted and deleted is chosen at random, as is the position of new insertions, and the length of the new target site duplication. One caveat to the random selection of TEs is that simpoolTE will not simulate the insertion or deletion of nested TEs. In this case, nested TEs are considered any TEs at positions within one read length (-rlen) of an annotated reference TE. One additional caveat is that since pirs does not simulate reads for any sequences containing N's, simpoolTE replaces all N's in the original chromosome with bases randomly drawn from the probability density of A's, C's, G's, and T's. This new chromosome is generated as /usr/local/teSim_simData/refChrom_ATCG_replaced.fa.  

```
usage: python /usr/local/simpoolTE_V3.2.py <required> [optional] 
    -wd [full path to working directory]
    -pirs <full path to pirs executable>
    -seqtk <full path to seqtk executable>
    -c <chromosome to simulate in fasta format>
    -b <annotation of TEs on chromosome in BED6 format>
    -ex [newline separated list of TEs to exclude from simulation (must match name in bed file)]
    -mnlen [minimum length of TEs to insert and delete]
    -mxlen [maximum length of TEs to insert and delete]
    -nchr [number of chromosomes in the simulated population]
    -nte [number of insertions and deletion to simulate in lowest frequency class] #NOTE: the total number of unique insertions and deletions will be greater than this
    -r [seed for random number generator]
    -x [coverage to simulate]
    -rlen [read length to simulate]
    -insz [insert size to simulate]
    -t [number of threads]
```

### Output
The important output data are a set of paired-end reads in fastq format and information about the positions and allele frequencies of the simulated insertions and deletions:

* /usr/local/teSim_pooledReads/simPool_1.fq
* /usr/local/teSim_pooledReads/simPool_2.fq
* /usr/local/teSim_simData/simulated_insertion_list.txt
* /usr/local/teSim_simData/simulated_deletion_list.txt


### Citation
If you use simpoolTE in your work please cite:
```
Adrion, J.R., M.J. Song, D.R. Schrider, M.W. Hahn, and S. Schaack. 2017. Genome-wide estimates of transposable element insertion and deletion rates in *Drosophila melanogaster*. Genome Biology and Evolution. https://doi.org/10.1093/gbe/evx050.
```
