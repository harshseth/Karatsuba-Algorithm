#ifndef SCAFFOLD4_H
#define SCAFFOLD4_H

#include <time.h>

#define WORDBYTES 1 // bytes per word
#define WORDBITS 4 // bits per word

//mpz import/export parameters
#define NAILS 0 //bit padding, usually 0
#define ORDER -1 //most significant word first, -1 for least first
#define ENDIAN 0 // native endianness

extern void Product4(void *a, void *b, void *c, unsigned int sa,
    unsigned int bitsa, unsigned int sb, unsigned int bitsb, unsigned
        int *sc, unsigned int *bitsc);

extern double mytime (void);

#endif
