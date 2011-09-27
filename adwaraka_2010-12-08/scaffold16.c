#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include "scaffold16.h"

/* You are suppose to change the routine Product16 here to your own routine 
 * The mpz calls in the scaffolded Product16 below are the normal GMP function
 * calls and should be neglected. By casting the void pointers as normal unsigned
 * integers, you should be able to access the data values as normal 4 bytes words.
 */

void Product16(void *a, void *b, void *c, unsigned int sa,
    unsigned int bitsa, unsigned int sb, unsigned int bitsb, unsigned
        int *sc, unsigned int *bitsc){

    mpz_t x,y,z;
    mpz_init(x);
    mpz_init(y);
    mpz_init(z); 

    /* Cast a and b into short integers of size 16 bits */
    unsigned short int *sint_a = (unsigned short int *) a;
    unsigned short int *sint_b = (unsigned short int *) b;

    /* Now sint_a can be view as an array of words of size 16 bits */
    /* Similarly for sint_b */
    //printf("%u %u \n", sint_a[0], sint_a[*sa - 1]);
    //printf("%u %u \n", sint_b[0], sint_b[*sb - 1]);

    mpz_import(x, sa, ORDER, WORDBYTES, ENDIAN, NAILS, a);
    mpz_import(y, sb, ORDER, WORDBYTES, ENDIAN, NAILS, b);
    mpz_mul(z,x,y);

    c = mpz_export(c, sc, ORDER, WORDBYTES, ENDIAN, NAILS, z);
}

