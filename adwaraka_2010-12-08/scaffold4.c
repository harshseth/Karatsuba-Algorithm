#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include "scaffold4.h"

/* You are suppose to change the routine Product4 here to your own routine 
 * The mpz calls in the scaffolded Product4 below are the normal GMP function
 * calls and should be neglected. By casting the void pointers as normal unsigned
 * integers, you should be able to access the data values as normal 4 bytes words.
 */

void Product4(void *a, void *b, void *c, unsigned int sa,
    unsigned int bitsa, unsigned int sb, unsigned int bitsb, unsigned
        int *sc, unsigned int *bitsc){

    mpz_t x,y,z;
    mpz_init(x);
    mpz_init(y);
    mpz_init(z); 

    /* Cast a and b into short integers of size 8 bits */
    unsigned char *cint_a = (unsigned char *) a;
    unsigned char *cint_b = (unsigned char *) b;

    /* Now cint_a can be view as an array of words of size 8
     * bits
     * with 4 padded bits in each words.
     * Similarly for int_b
     */

    //printf("%d %d \n", cint_a[0], cint_a[*sa - 1]);
    //printf("%d %d \n", cint_b[0], cint_b[*sb - 1]);

    mpz_import(x, sa, ORDER, WORDBYTES, ENDIAN, NAILS, a);
    mpz_import(y, sb, ORDER, WORDBYTES, ENDIAN, NAILS, b);
    mpz_mul(z,x,y);

    c = mpz_export(c, sc, ORDER, WORDBYTES, ENDIAN, NAILS, z);
}

