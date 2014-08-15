/* QWalk (qwstate.c) 
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
#include "qwmem_complex.h"
#include "qwstate.h"
#include "qwconsts.h"


double complex ****createGroverState2D(int max, unsigned int lattType){
  double complex ****matrix;
  const int size = (lattType == CYCLE_LATT) ? max : 2*max+1;

  if(max<1)
    return NULL;

  matrix = allocComplex4D(2, 2, size, size);
  if(!matrix)
    return NULL;
  cleanComplex4D(matrix, 2, 2, size, size);

  /* Note that the mathematical lattice ranges from -max to max, while the
   * computational representation of this lattice ranges from 0 to 2*max in C.
   * Therefore, if we want to describe a particle in the site (0,0) we must
   * add max to each of its component in order to perform the appropriate
   * conversion.
   */
  if(lattType == CYCLE_LATT){
    matrix[0][0][max/2][max/2] = 0.5;
    matrix[0][1][max/2][max/2] =-0.5; 
    matrix[1][0][max/2][max/2] =-0.5;
    matrix[1][1][max/2][max/2] = 0.5;
  }
  else{
    matrix[0][0][max][max] = 0.5;
    matrix[0][1][max][max] =-0.5; 
    matrix[1][0][max][max] =-0.5;
    matrix[1][1][max][max] = 0.5;
  }

  return matrix;
}



double complex ****createFourierState2D(int max, unsigned int lattType){
  double complex ****matrix;
  const int size = (lattType == CYCLE_LATT) ? max : 2*max+1;

  if(max<1)
    return NULL;

  matrix = allocComplex4D(2, 2, size, size);
  if(!matrix)
    return NULL;
  cleanComplex4D(matrix, 2, 2, size, size);

  /* Note that the mathematical lattice ranges from -max to max, while the
   * computational representation of this lattice ranges from 0 to 2*max in C.
   * Therefore, if we want to describe a particle in the site (0,0) we must
   * add max to each of its component in order to perform the appropriate
   * conversion.
   */
  if(lattType == CYCLE_LATT){
    matrix[0][0][max/2][max/2] = 0.5;
    matrix[0][1][max/2][max/2] = (1.0-I)/(2.0*sqrt(2.0)); 
    matrix[1][0][max/2][max/2] = 0.5;
    matrix[1][1][max/2][max/2] =-(1.0-I)/(2.0*sqrt(2.0));
  }
  else{
    matrix[0][0][max][max] = 0.5;
    matrix[0][1][max][max] = (1.0-I)/(2.0*sqrt(2.0)); 
    matrix[1][0][max][max] = 0.5;
    matrix[1][1][max][max] =-(1.0-I)/(2.0*sqrt(2.0));
  }


  return matrix;
}



double complex ****createHadamardState2D(int max, unsigned int lattType){
  double complex ****matrix;
  const int size = (lattType == CYCLE_LATT) ? max : 2*max+1;

  if(max<1)
    return NULL;

  matrix = allocComplex4D(2, 2, size, size);
  if(!matrix)
    return NULL;
  cleanComplex4D(matrix, 2, 2, size, size);

  /* Note that the mathematical lattice ranges from -max to max, while the
   * computational representation of this lattice ranges from 0 to 2*max in C.
   * Therefore, if we want to describe a particle in the site (0,0) we must
   * add max to each of its component in order to perform the appropriate
   * conversion.
   */
  if(lattType == CYCLE_LATT){
    matrix[0][0][max/2][max/2] = 0.5;
    matrix[0][1][max/2][max/2] = 0.5*I; 
    matrix[1][0][max/2][max/2] = 0.5*I;
    matrix[1][1][max/2][max/2] =-0.5 ;
  }
  else{
    matrix[0][0][max][max] = 0.5;
    matrix[0][1][max][max] = 0.5*I; 
    matrix[1][0][max][max] = 0.5*I;
    matrix[1][1][max][max] =-0.5 ;
  }

  return matrix;
}


double complex **createHadamardState1D(int max, unsigned int lattType){
  double complex **matrix;

  if(max<1)
    return NULL;

  matrix = (lattType == LINE_LATT) ? 
    allocComplex2D(2, 2*max+1) : allocComplex2D(2, max);
  if(!matrix)
    return NULL;
  
  if(lattType == LINE_LATT)
    cleanComplex2D(matrix, 2, 2*max+1);
  else
    cleanComplex2D(matrix, 2, max);

  /* Note that the mathematical lattice ranges from -max to max in the LINE lattice, 
   * while the computational representation of this lattice ranges from 0 to 2*max 
   * in C. Therefore, if we want to describe a particle in the site (0,0) we must
   * add max to each of its component in order to perform the appropriate
   * conversion. In lattices that range from 0 to max-1 we place the particle 
   * initially in the "middle" of the lattice.
   */
  if(lattType == LINE_LATT){
    matrix[0][max] = 1.0/sqrt(2.0);
    matrix[1][max] = 1.0/sqrt(2.0)*I;
  }
  else{
    matrix[0][max/2] = 1.0/sqrt(2.0);
    matrix[1][max/2] = 1.0/sqrt(2.0)*I;
  }
    
  return matrix;
}



int checkState2D(double complex ****matrix, options2D_t options){
  int m,n,j,k;
  double totalprob;

  totalprob = 0.0;

  if(options.lattType == CYCLE_LATT){
    for(j=0; j<2; j++)
      for(k=0; k<2; k++)
	for(m=0; m<options.max; m++)
	  for(n=0; n<options.max; n++)
	    totalprob += matrix[j][k][m][n]*conj(matrix[j][k][m][n]);
  }
  else{
    for(j=0; j<2; j++)
      for(k=0; k<2; k++)
	for(m=0; m<=2*options.max; m++)
	  for(n=0; n<=2*options.max; n++)
	    totalprob += matrix[j][k][m][n]*conj(matrix[j][k][m][n]);
  }

  if( fabs(totalprob-1.0) > WALK_TOL )
    return 0;

  return 1;
}



int checkState1D(double complex **matrix, options1D_t options){
  int m,j;
  double totalprob;

  totalprob = 0.0;

  if(options.lattType == LINE_LATT){
    for(j=0; j<2; j++)
      for(m=0; m<=2*options.max; m++)
	totalprob += matrix[j][m]*conj(matrix[j][m]);
  }
  else{
    for(j=0; j<2; j++)
      for(m=0; m<options.max; m++)
	totalprob += matrix[j][m]*conj(matrix[j][m]);
  }

  if( fabs(totalprob-1.0) > WALK_TOL )
    return 0;

  return 1;
}



int checkSymmetry1D(double complex **matrix, int max){
  int m;

  for(m=-max; m<=max; m++){
    int j;
    double probA, probB;

    probA = 0.0;
    probB = 0.0;

    for(j=0; j<2; j++){
      probA += matrix[j][max-m]*conj(matrix[j][max-m]);
      probB += matrix[j][max+m]*conj(matrix[j][max+m]);
    }
    if(fabs(probA - probB) > WALK_TOL)
      return 0;
  }

  return 1;
}



int checkXSymmetry2D(double complex ****matrix, int max){
  int m, n;

  for(m=-max; m<=max; m++){
    for(n=0; n<=max; n++){
      int j, k;
      double probA, probB;

      probA = 0.0;
      probB = 0.0;

      for(j=0; j<2; j++){
	for(k=0; k<2; k++){
	  probA += matrix[j][k][max-m][max+n]*conj(matrix[j][k][max-m][max+n]);
	  probB += matrix[j][k][max+m][max+n]*conj(matrix[j][k][max+m][max+n]);
	}
      }
      if(fabs(probA - probB) > WALK_TOL)
	return 0;
    }
  }

  return 1;
}



int checkYSymmetry2D(double complex ****matrix, int max){
  int m, n;

  for(m=-max; m<=max; m++){
    for(n=0; n<=max; n++){
      int j, k;
      double probA, probB;

      probA = 0.0;
      probB = 0.0;

      for(j=0; j<2; j++){
	for(k=0; k<2; k++){
	  probA += matrix[j][k][max+m][max-n]*conj(matrix[j][k][max+m][max-n]);
	  probB += matrix[j][k][max+m][max+n]*conj(matrix[j][k][max+m][max+n]);
	}
      }
      if(fabs(probA - probB) > WALK_TOL)
	return 0;
    }
  }

  return 1;
}

