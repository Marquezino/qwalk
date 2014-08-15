/* QWalk (qwcoin_io.c) 
 * Copyright (C) 2008  Franklin Marquezino
 *                                                                                                                           
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *                                                                                                                           
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA
 */

#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<string.h>
#include<math.h>
#include<limits.h>
#include<time.h>
#include "qwmem_complex.h"
#include "qwcoin_io.h"
#include "qwconsts.h"


double complex **readCoinFile1D(const char *filename){
  FILE *in;
  double complex **coin;
  int j, k;
  char keyword[100];

  in = fopen(filename,"rt");
  if(!in)
    return NULL;

  /* First we search the BEGINCOIN keyword, ... */
  do{
    fscanf(in,"%s",keyword);
  }while(STRNEQ(keyword,"BEGINCOIN") && !feof(in));

  if(feof(in))
    return NULL;
  
  coin = allocComplex2D(2, 2);
  if(!coin)
    return NULL;

  /* ...when we find it, we start reading complex numbers and we 
   * assign each of them to an entry of the matrix. The complex 
   * numbers must be given in the input file with the real and 
   * imaginary parts separated by a space.
   */
  for(j=0; j<2; j++){
    for(k=0; k<2; k++){
      double real, imag;

      fscanf(in,"%lf", &real);
      fscanf(in,"%lf", &imag);
      coin[j][k] = real+ I*imag;
    }
  }

  fclose(in);
  return coin;
}


double complex ****readCoinFile2D(const char *filename){
  FILE *in;
  double complex ****coin;
  int j, k;
  char keyword[100];

  in = fopen(filename,"rt");
  if(!in)
    return NULL;

  /* First we search the BEGINCOIN keyword... */
  do{
    fscanf(in,"%s",keyword);
  }while(STRNEQ(keyword,"BEGINCOIN") && !feof(in));

  if(feof(in))
    return NULL;
  
  coin = allocComplex4D(2, 2, 2, 2);
  if(!coin)
    return NULL;
  
  /* ...when we find it, we start reading complex numbers and we 
   * assign each of them to an entry of the matrix. The complex 
   * numbers must be given in the input file with the real and 
   * imaginary parts separated by a space.
   */
  for(j=0; j<2; j++){
    for(k=0; k<2; k++){
      int jprime, kprime;

      for(jprime=0; jprime<2; jprime++){
	for(kprime=0; kprime<2; kprime++){
	  double real, imag;

	  fscanf(in,"%lf", &real);
	  fscanf(in,"%lf", &imag);
	  coin[j][k][jprime][kprime] = real+ I*imag;
	}
      }
    }
  }
  
  return coin;
}
