/* QWalk (qwlinks.c) 
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
#include<math.h>
#include "qwmem_int.h"
#include "qwlinks.h"
#include "qwconsts.h"


int initBrokenLink1D(int ***L, int max, unsigned char type){
  int j;
  static char previousAlloc=0;
  const int rbound = (type == LINE_LATT) ? 2*max+1 : max;

  if(max<1)
    return 1;

  if(!previousAlloc){
    *L = (type == LINE_LATT) ? 
      allocInt2D(2,2*max+1) : allocInt2D(2,max);
  }

  if(*L){
    previousAlloc = 1;
  }
  else{
    previousAlloc = 0;
    return 2;
  }

  for(j=0; j<2; j++){
    int m;

    for(m=0; m<rbound; m++){
      (*L)[j][m] = (int)pow(-1.0, j);
    }
  }

  return 0;
}


int initBrokenLink2D(int *****L1, int *****L2, int max, unsigned char type){
  int j, k;
  static char previousAlloc1 = 0, previousAlloc2 = 0;
  const int rbound = (type == CYCLE_LATT) ? max : 2*max+1;

  if(max<1)
    return 1;

  if(!previousAlloc1)
    *L1 = allocInt4D(2, 2, rbound, rbound);

  if(*L1)
    previousAlloc1 = 1;
  else{
    previousAlloc1 = 0;
    return 2;
  }

  if(!previousAlloc2 && type==DIAG_LATT)
    *L2 = allocInt4D(2, 2, 2*max+1, 2*max+1);
  else if(!previousAlloc2 && type!=DIAG_LATT)
    *L2 = allocInt4D(1, 1, 1, 1);

  if(*L2)
    previousAlloc2 = 1;
  else{
    previousAlloc2 = 0;
    return 2;
  }

  /* Initializing all the links closed */
  for(j=0; j<2; j++){
    int auxj, auxk, m, n;
    auxj = (int)pow(-1.0, j);

    for(k=0; k<2; k++){
      auxk = (int)pow(-1.0,k);

      for(m=0; m<rbound; m++){
	for(n=0; n<rbound; n++){
	  (*L1)[j][k][m][n] = auxj;
	  if( type==DIAG_LATT )
	    (*L2)[j][k][m][n] = auxk;
	}
      }
    }
  }

  return 0;

}
