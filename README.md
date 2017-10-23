# ngsPSMC
[![Build Status](https://travis-ci.org/ANGSD/ngsPSMC.svg?branch=master)](https://travis-ci.org/ANGSD/ngsPSMC)


//    fprintf(stderr,"%s\n",*argv);
    if(!strcasecmp(*argv,"-tole"))
      p->tole = atof(*(++argv));
    else  if(!strcasecmp(*argv,"-maxIter"))
      p->maxIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-winSize"))
      p->block = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-RD"))
      p->RD = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-nThreads"))
      p->nThreads = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-nIter"))
      p->nIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-p"))
      p->par->pattern =  strdup(*(++argv));
    else  if(!strcasecmp(*argv,"-tkfile"))
      p->tkfile =  strdup(*(++argv));
    else  if(!strcasecmp(*argv,"-nSites"))
      p->nSites = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-seed"))
      p->seed = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-infile"))
      inffilename = strdup(*++argv);
    else  if(!strcasecmp(*argv,"-doQuad"))
      p->doQuad = atoi(*++argv);
    else  if(!strcasecmp(*argv,"-r")){
      p->chooseChr = get_region(*(++argv),p->start,p->stop);
      if(!p->chooseChr)
        return NULL;
    }
