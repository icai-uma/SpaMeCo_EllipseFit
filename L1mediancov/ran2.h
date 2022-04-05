#ifndef _RAN2_H

#define _RAN2_H

/* Long period (> 2 * 10^18) random number generator of L'Ecuyer with Bays-Durham shue
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.*/

float ran2(long *idum);


#endif

