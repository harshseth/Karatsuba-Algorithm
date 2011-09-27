#ifndef SCAFFOLD16_H
#define SCAFFOLD16_H

#define WORDBYTES 2 // bytes per word
#define WORDBITS 16 // bits per word

//mpz import/export parameters
#define NAILS 0 //bit padding, usually 0
#define ORDER -1 //most significant word first, -1 for least first
#define ENDIAN 0 // native endianness

extern void Product16(void *a, void *b, void *c, unsigned int sa,
    unsigned int bitsa, unsigned int sb, unsigned int bitsb, unsigned
    int *sc, unsigned int *bitsc);

extern double mytime (void);

#endif
