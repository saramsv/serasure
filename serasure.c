/* *
 * Copyright (c) 2014, James S. Plank and Kevin Greenan
 * All rights reserved.
 *
 * Jerasure - A C/C++ Library for a Variety of Reed-Solomon and RAID-6 Erasure
 * Coding Techniques
 *
 * Revision 2.0: Galois Field backend now links to GF-Complete
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *  - Neither the name of the University of Tennessee nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/* Jerasure's authors:

   Revision 2.x - 2014: James S. Plank and Kevin M. Greenan.
   Revision 1.2 - 2008: James S. Plank, Scott Simmerman and Catherine D. Schuman.
   Revision 1.0 - 2007: James S. Plank.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <gf_rand.h>
#include "jerasure.h"
#include "reed_sol.h"
#include "galois.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

static void usage(char *s)
{
  fprintf(stderr, "usage: reed_sol_01 k m w seed - Does a simple Reed-Solomon coding example in GF(2^w).\n");
  fprintf(stderr, "       \n");
  fprintf(stderr, "w must be 8, 16 or 32.  k+m must be <= 2^w.  It sets up a classic\n");
  fprintf(stderr, "Vandermonde-based generator matrix and encodes k devices of\n");
  fprintf(stderr, "%ld bytes each with it.  Then it decodes.\n", sizeof(long));
  fprintf(stderr, "       \n");
  fprintf(stderr, "This demonstrates: jerasure_matrix_encode()\n");
  fprintf(stderr, "                   jerasure_matrix_decode()\n");
  fprintf(stderr, "                   jerasure_print_matrix()\n");
  fprintf(stderr, "                   reed_sol_vandermonde_coding_matrix()\n");
  if (s != NULL) fprintf(stderr, "%s\n", s);
  exit(1);
}

static void print_data_and_coding(int k, int m, int w, int size, 
		char **data, char **coding) 
{
  int i, j, x;
  int n, sp;

  if(k > m) n = k;
  else n = m;
  sp = size * 2 + size/(w/8) + 8;

  printf("%-*sCoding\n", sp, "Data");
  for(i = 0; i < n; i++) {
	  if(i < k) {
		  printf("D%-2d:", i);
		  for(j=0;j< size; j+=(w/8)) { 
			  printf(" ");
			  for(x=0;x < w/8;x++){
				printf("%02x", (unsigned char)data[i][j+x]);
			  }
		  }
		  printf("    ");
	  }
	  else printf("%*s", sp, "");
	  if(i < m) {
		  printf("C%-2d:", i);
		  for(j=0;j< size; j+=(w/8)) { 
			  printf(" ");
			  for(x=0;x < w/8;x++){
				printf("%02x", (unsigned char)coding[i][j+x]);
			  }
		  }
	  }
	  printf("\n");
  }
	printf("\n");
}

void negative_maker(char *parity, int nbyte, int w)
{
    if (w == 8)
        reed_sol_galois_w08_region_multby_2(parity, nbyte);
    else if (w == 16)
        reed_sol_galois_w16_region_multby_2(parity, nbyte);
    else if (w == 32)
        reed_sol_galois_w32_region_multby_2(parity, nbyte);
    else
        printf("Invalid field");
}

char *** step3_permutation(char ***coding_total, char ***data_total, int m, int k, int w)
{
    char *** step3_h;
    int i, j, r, x;
    step3_h = talloc(char * , m);
    for (r = 0; r<m ;r++)
    {
        step3_h[r] = talloc(char *, m);
    }
    for(i = 0; i< m; i++)
    {
        for (j = 0; j< m; j++)
        {
            step3_h[j][i] = coding_total[i][(i+j)%m];
        }
    }
    return step3_h;
}

void step4(char *** code, char *** code_copy)
{

}

void create_a_copy(char *** dest, char *** src, int m)
{
    int i,j;
    for (i = 0; i < m ; i++)
    {
        for(j = 0; j < m; j++)
            memcpy(dest[i][j], src[i][j], sizeof(long));
    }
}

int main(int argc, char **argv)
{
  long l;
  int k, w, i, j, m, instance;
  int *matrix;
  char **data, **coding, **dcopy, **ccopy;
  char ***data_total, ***coding_total,***datacopy, ***codingcopy, ***perm;
  unsigned char uc;
  int *erasures, *erased;
  uint32_t seed;
 
  if (argc != 5) usage(NULL);
  if (sscanf(argv[1], "%d", &k) == 0 || k <= 0) usage("Bad k");
  if (sscanf(argv[2], "%d", &m) == 0 || m <= 0) usage("Bad m");
  if (sscanf(argv[3], "%d", &w) == 0 || (w != 8 && w != 16 && w != 32)) usage("Bad w");
  if (sscanf(argv[4], "%u", &seed) == 0) usage("Bad seed");
  if (w <= 16 && k + m > (1 << w)) usage("k + m is too big");

  matrix = reed_sol_vandermonde_coding_matrix(k, m, w);

  printf("<HTML><TITLE>reed_sol_01 %d %d %d %d</title>\n", k, m, w, seed);
  printf("<h3>reed_sol_01 %d %d %d %d</h3>\n", k, m, w, seed);
  printf("<pre>\n");
  printf("Last m rows of the generator Matrix (G^T):\n\n");
  jerasure_print_matrix(matrix, m, k, w);

  data_total = talloc(char**, m);
  datacopy = talloc(char **, m);
  coding_total = talloc(char **, m);
  codingcopy = talloc(char **, m);
  perm = talloc(char **, m);
  
  for (instance = 0; instance < m; instance++)
  {
      seed = seed + instance;
      MOA_Seed(seed);
      data = talloc(char *, k);
      dcopy = talloc(char *, k);
      for (i = 0; i < k; i++) {
        data[i] = talloc(char, sizeof(long));
        dcopy[i] = talloc(char, sizeof(long));
        for (j = 0; j < sizeof(long); j++) {
          uc = MOA_Random_W(8, 1);
          data[i][j] = (char) uc;
        }
        memcpy(dcopy[i], data[i], sizeof(long));
      }
      data_total[instance] = data;
      datacopy[instance]= dcopy;
      
      coding = talloc(char *, m);
      ccopy = talloc(char *, m);

      for (i = 0; i < m; i++) {
        coding[i] = talloc(char, sizeof(long));
        ccopy[i] = talloc(char, sizeof(long));
      }

      jerasure_matrix_encode(k, m, w, matrix, data, coding, sizeof(long));

      for (i = 0; i < m; i++) {
        memcpy(ccopy[i], coding[i], sizeof(long));
      }
      coding_total[instance] = coding;
      codingcopy[instance] = ccopy;

      printf("Encoding Complete %d:\n\n", instance);
      print_data_and_coding(k, m, w, sizeof(long), data_total[instance], coding_total[instance]);
  }
  perm = step3_permutation(coding_total, data_total, m, k, w);
  printf("Data and parities after permutations: \n\n");
  for(i = 0; i < m; i++)
  {
      print_data_and_coding(k, m, w, sizeof(long), data_total[i], perm[i]);
  }
  int nbyte = 4;
  negative_maker(perm[0][0], nbyte, w);

  printf("After multiplying by 2: \n\n");
  for(i = 0; i < m; i++)
  {
      print_data_and_coding(k, m, w, sizeof(long), data_total[i], perm[i]);
  }
  create_a_copy(coding_total, perm, m);
  printf("Test: \n\n");
  for(i = 0; i < m; i++)
  {
      print_data_and_coding(k, m, w, sizeof(long), data_total[i], coding_total[i]);
  }
  //step4(coding_total, perm)

  /*
  erasures = talloc(int, (m+1));
  erased = talloc(int, (k+m));
  for (i = 0; i < m+k; i++) erased[i] = 0;
  l = 0;
  for (i = 0; i < m; ) {
    erasures[i] = MOA_Random_W(31, 0)%(k+m);
    if (erased[erasures[i]] == 0) {
      erased[erasures[i]] = 1;
      memcpy((erasures[i] < k) ? data[erasures[i]] : coding[erasures[i]-k], &l, sizeof(long));
      i++;
    }
  }
  erasures[i] = -1;

  printf("Erased %d random devices:\n\n", m);
  print_data_and_coding(k, m, w, sizeof(long), data, coding);
  
  i = jerasure_matrix_decode(k, m, w, matrix, 1, erasures, data, coding, sizeof(long));

  printf("State of the system after decoding:\n\n");
  print_data_and_coding(k, m, w, sizeof(long), data, coding);
  
  for (i = 0; i < k; i++) if (memcmp(data[i], dcopy[i], sizeof(long)) != 0) {
    printf("ERROR: D%x after decoding does not match its state before decoding!<br>\n", i);
  }
  for (i = 0; i < m; i++) if (memcmp(coding[i], ccopy[i], sizeof(long)) != 0) {
    printf("ERROR: C%x after decoding does not match its state before decoding!<br>\n", i);
  }*/
  return 0;
}
