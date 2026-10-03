#include "PunctureExecution.h"
#include <string.h>
#ifndef PUNCTURES_KOKKOS
const char *Puncture_execution_name(void){return "reference";}
int Puncture_execution_initialize(int threads){(void)threads;return -1;}
int Puncture_execution_concurrency(void){return 1;}
unsigned long long Puncture_execution_device_free_bytes(void){return 0;}
int Puncture_execution_device_description(char*buffer,int capacity){(void)buffer;(void)capacity;return -1;}
int Puncture_execution_statistics(PunctureExecutionStats*s){if(!s)return -1;memset(s,0,sizeof(*s));return 0;}
void Puncture_execution_reset_statistics(void){}
#endif
