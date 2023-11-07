//Parameters for the coalescence simulation program : fastsimcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
920000
720000
820000
//Haploid samples sizes   ### gli stessi di easySFS
14
14
14
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
// migration matrix
0.000 0.000 0.000
0.000 0.000 $MIG1$
0.000 $MIG2$ 0.000
// migration matrix
0.000 0.000 0.000
0.000 0.000 0.000
0.000 0.000 0.000
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
9 historical event
2333 0 0 0 1 0 1
19167 0 0 0 1 0 0
21667 0 0 0 1 0 1
50000 0 0 0 1 0 0
56167 0 0 0 1 0 1
62333 0 0 0 1 0 0
70667 0 0 0 1 0 1
87167 0 1 1 0.086 0 1
151500 1 2 1 0.139 0 1
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 1.93e-8 OUTEXP
