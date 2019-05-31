# ngsPSMC
[![Build Status](https://travis-ci.org/ANGSD/ngsPSMC.svg?branch=master)](https://travis-ci.org/ANGSD/ngsPSMC)

Dependencies: `htslib`.

First one need to run angsd on a sorted bam file. `-out` is the output filename. `-gl 1` means that the standard GLs (as in SamTools) will be used.
```./angsd/angsd -i /willerslev/datasets/public/1000genomes_2015_nature/bam/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam -dopsmc 1  -out data.real1/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211 -gl 1 -minq 20 -minmap 30```

Now it's time to run ngsPSMC. There are two options: to initialise all the parameters similar to the original PSMC, or to read them from another PSMC file with `-infile` option (e.g. if you already have another high-coverage sample). Notice that theta (mutation rate) and maxT (last non-inifinite time point in the time discritisation) are not optimised in our version. If you do not use `-infile` option, it is highly recommended to specify theta if you run ngsPSMC with an organism other than humans. One can calculate it by the simple formula theta=2*N_0*mu*binsize (N_0 is a reference effective population size, e.g. 10000, mu is the mutation rate per basepair per generation, binsize is the window length TODO CHECK if binsize is part of theta, also is it 2 or 4?).

If `-infile` is used, there is a possibility to choose which round of optimisation is used for initialisation. To choose the last round use `-rd -1`.

The following command reads data from `data.psmc.idx`, and initialises all the parameters from the last round from the file `another_data.psmc`. Then it resets initial values of effective population size to 1 (`-init`). It runs for 25 rounds of optimisation (`-nIter`).

```./ngsPSMC/ngsPSMC data.psmc.idx -infile another_data.psmc -dospline 0 -dolinear 1 -rd -1 -nthreads 25 -nIter 25 -init 1  1>data.real1/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.psmc.idx.1.stdout 2>data.real1/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.psmc.idx.1.stderr```