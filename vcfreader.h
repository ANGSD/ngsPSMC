#pragma once
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <map>
#include <htslib/bgzf.h>
#include "fastas.h"
#include "header.h"
#include "psmcreader.h"

rawdata readvcf(char* filename, char* chr);