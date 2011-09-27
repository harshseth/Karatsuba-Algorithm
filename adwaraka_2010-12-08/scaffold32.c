#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <string.h> 
#include "scaffold32.h"
 
/* You are suppose to change the routine Product32 here to your own routine
 * The mpz calls in the scaffolded Product32 below are the normal GMP function
 * calls and should be neglected. By casting the void pointers as normal unsigned
 * integers, you should be able to access the data values as normal 4 bytes words.
 */
 
/* wa is word length of a, ba is bit length of a */

/* add32 USES ALGO 1.4 */
void add32(unsigned int *int_a, unsigned int *int_b, unsigned int *int_c,
                        unsigned int wa, unsigned int wb, unsigned int *wc)
{
     int i;
     unsigned int carry = 0, *temp;
     unsigned long long p, wordsize = 4294967296llu;

     if (wb > wa)
     {
	   temp = int_a;
	   int_a = int_b;
	   int_b = temp;
	   *temp = wa;
	   wa = wb;
	   wb = *temp;
     }

     for (i = 0; i < wa; i++)
     {
         if (i < wb)
            p = (unsigned long long)int_a[i] + int_b[i] + carry;
         else
            p = (unsigned long long)int_a[i] + carry;
         int_c[i] = p%wordsize;
         carry = p/wordsize;
     }

     while (carry > 0)
     {
	p = (unsigned long long)carry; 
	int_c[i] = p%wordsize;
	carry=p/wordsize;
	i++;
     }

     *wc = i;
}


void subtract32(unsigned int *int_a, unsigned int *int_b, unsigned int *int_c,
                        unsigned int wa, unsigned int wb, unsigned int *wc)
{
    int i = 0;
    unsigned int carry = 0, *temp;
    unsigned long long s = 0, wordsize = 4294967296llu;

    /* FOR FOLLOWING CONDITION; SWAP THE VARIABLES */
    if (wb > wa)
    {
	   temp = int_a;
	   int_a = int_b;
	   int_b = temp;
	   *temp = wa;
	   wa = wb;
	   wb = *temp;
     }
 
 
     for (i = 0; i < wa; i++)
     {
	 if (int_a[i]<int_b[i])
	 {
	    s = int_a[i] - int_b[i] + wordsize;
	    int_a[i+1]--;
	 }
	 else
	    s = int_a[i] - int_b[i];
		
	int_c[i] = s;	
     }

    *wc = wa;
}


/* USES ALGO 1.8 */
void multiply32(unsigned int *int_a, unsigned int *int_b,
    unsigned int *int_c, unsigned int wa, unsigned int wb, unsigned int *wc)
{
    int i,j;
    unsigned int carry=0;
    unsigned long long p, wordsize=4294967296llu;

    /*CATERING FOR (wa <= wb/2) OR (wa >= wb/2)*/
    if ((wa <= wb/2) || (wa >= wb/2))
       for(i = 0; i < wa; i++)   
          for(j = 0 ;j < wb; j++)
             int_c[i+j] = 0;


    for(i=0;i<wa;i++)
    {
        carry=0;
        for(j=0;j<wb;j++)
        {
	    p=(unsigned long long)int_a[i]*int_b[j]+int_c[i+j]+carry;
            int_c[i+j]=p%wordsize;
            carry=p>>32;
        }
	int_c[i+wb]=carry;	
    }
    *wc = wa + wb;
    
}

/* USES ALGO 5.2 */
void karatsuba32(unsigned int *int_a, unsigned int *int_b,
    unsigned int *int_c, unsigned int wa, unsigned int wb, unsigned int *wc)
{
    int i, j;
    unsigned int N;

    /*CHECK FOR WORDSIZES*/
    if (wa <= 35 || wb <= 35 || wa <= wb/2 || wb <= wa/2)
    {
       multiply32(int_a, int_b, int_c, wa, wb, wc);
       return;
    }


    if (wb > wa)
       N = wa>(wb/2)?wa:(wb/2);
    else
       N = wb>(wa/2)?wb:(wa/2);


    unsigned int *u2 = (unsigned int *)int_a;
    unsigned int *u1 = (unsigned int *)int_a + N/2;

    unsigned int *v2 = (unsigned int *)int_b;
    unsigned int *v1 = (unsigned int *)int_b + N/2;


    /*PRODUCT U1V1*/
    unsigned int *produ1v1 = (unsigned int *)calloc(((wa - N/2) + (wb - N/2)), sizeof(unsigned int));
    unsigned int s1;
    karatsuba32(u1, v1, produ1v1, wa - N/2, wb - N/2, &s1);

    /*SET produ1v1 INTO int_c*/
    for (i = 2*(N/2); i < 2*(N/2) + (s1); i++)
        int_c[i]+= produ1v1[i-2*(N/2)];


    /*PRODUCT U2V2*/
    unsigned int *produ2v2 = (unsigned int *)calloc(N, sizeof(unsigned int));
    unsigned int s2;
    karatsuba32(u2, v2, produ2v2, N/2 , N/2, &s2);

    /*SET produ2v2 INTO int_c*/
    for (i = 0; i < s2; i++)
	    int_c[i]+= produ2v2[i];


    /*SUM U1U2*/
    unsigned int *sumu1u2 = (unsigned int *)calloc(N/2 + 2, sizeof(unsigned int));
    unsigned int s3;
    add32(u1, u2, sumu1u2, wa - N/2, N/2, &s3);
	
	   
    /*SUM V1V2*/
    unsigned int *sumv1v2 = (unsigned int *)calloc(N/2 + 2, sizeof(unsigned int));
    unsigned int s4;
    add32(v1, v2, sumv1v2, wb - N/2, N/2, &s4);


    /*PRODUCT T1T2*/
    unsigned int *t1t2 = (unsigned int *)calloc(s3 + s4, sizeof(unsigned int));
    unsigned int s5;
    karatsuba32(sumu1u2, sumv1v2, t1t2, s3, s4, &s5);


    /*SUM U1U2 + V1V2*/
    unsigned int *tempsum = (unsigned int *)calloc(s1 + s2, sizeof(unsigned int));
    unsigned int s6;
    add32(produ1v1, produ2v2, tempsum, s1, s2, &s6);


    /*SUBTRACT T1T2 - (SUM U1U2 + V1V2)*/
    /*STORING RESULT IN sub*/
    unsigned int *sub = (unsigned int *)calloc(s3 + s4 + N/2, sizeof(unsigned int));
    memset(sub,0,(s3+s4+N/2)*sizeof(unsigned int));
    unsigned int s7;
    subtract32(t1t2, tempsum, sub+N/2, s5, s6, &s7);
	

    /*ADDING TO FINAL*/
    add32(int_c, sub, int_c, wa + wb, s7 + N/2, wc);


    /*FREEING ALL MEMORY*/
    free(produ1v1);
    free(produ2v2);
    free(sumu1u2);
    free(sumv1v2);
    free(t1t2);
    free(tempsum);
    free(sub);
}
 

void Product32(void *a, void *b, void *c, unsigned int wa,
    unsigned int ba, unsigned int wb, unsigned int bb, unsigned int
    *wc, unsigned int *bc)
{
    int i;
    unsigned int w;
 
    /* Cast a and b into short integers of size 32 bits */
    unsigned int *int_a = (unsigned int *) a;
    unsigned int *int_b = (unsigned int *) b;
    unsigned int *int_c = (unsigned int *) c;

    memset(int_c,0,(wa+wb)*sizeof(unsigned int));	

    /*Invoke Karasuba*/
    karatsuba32(int_a, int_b, int_c, wa, wb, wc);
} 
