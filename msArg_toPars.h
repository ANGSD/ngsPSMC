#pragma once
#include "main_psmc.h"

msarg parse_msStr(char *str);
void msarg_toPars(msarg &ms,psmc_par *par,int winsize);
