# ngsPSMC
[![Build Status](https://travis-ci.org/ANGSD/ngsPSMC.svg?branch=master)](https://travis-ci.org/ANGSD/ngsPSMC)
This program implements the PSMC for genotype likelihoods. The program is under development.
## Installation
Dependencies: `htslib`.
```
git clone https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsPSMC
cd ngsPSMC;make HTSSRC=../htslib/;cd ..
```
## Generate Input Files from BAM/CRAM
Before running ngsPSMC you will need to generate the input files. This is done using ANGSD with the command below
```./angsd/angsd -i input.bam -dopsmc 1  -out psmcinput -gl 1 -minq 20 -minmapq 30```
Here -gl 1 indicates the SAMtools genotype likelihood model.
## Generate Input Files from simulated data
For validating the program we simulated data like below:
```
msHOT 2 1000 -t 8196 -r 1355 3000000 -l -eN 0.01 0.05 -eN 0.0375 0.5 -eN 1.25 1.0 -T >msraw
./angsd/misc/msHOT2glf -in msraw -out sim.d8 -win 100 -psmcfa 1 -psmc2 1 -do_seq_glf 0 -err 0.005 -depth 8
```

# Run examples
```./ngsPSMC input.psmc.idx -p "1*4+25*2+1*4+1*6\" -dospline 0 -nthreads 8 -nIter 20 -init 1  -theta 0.000233095 -rho 0.005357 ```

Where input are the output from either sim.d8.psmc.idx or psmcinput.psmc.idx
You can also specify the output from the original PSMC as input:

```./ngsPSMC input.psmc.idx -infile output.from.psmc -dospline 0 -nthreads 8 -nIter 20 ```


There are two options: to initialise all the parameters similar to the original PSMC, or to read them from another PSMC file with `-infile` option (e.g. if you already have another high-coverage sample). Notice that theta (mutation rate) and maxT (last non-inifinite time point in the time discritisation) are not optimised in our version. If you do not use `-infile` option, it is highly recommended to specify theta if you run ngsPSMC with an organism other than humans. One can calculate it by the simple formula theta=2*N_0*mu*binsize (N_0 is a reference effective population size, e.g. 10000, mu is the mutation rate per basepair per generation, binsize is the window length TODO CHECK if binsize is part of theta, also is it 2 or 4?).

If `-infile` is used, there is a possibility to choose which round of optimisation is used for initialisation. To choose the last round use `-rd -1`.

You can specify the intitial population size with -init popsize


# Details and potential issues
 - Optimization of rho is currently disabled.
 - Be carefull with the binsize (default 100) when comparing results between different winsizes.
