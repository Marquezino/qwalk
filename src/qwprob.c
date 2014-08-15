/* QWalk (qwprob.c) 
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
#include "qwoptions_io.h"
#include "qwmem_real.h"
#include "qwprob.h"
#include "qwconsts.h"

double **getProbArray2D(double complex ****matrix, options2D_t opts){
  int m,n,j,k;
  double **probMatrix, totalprob;
  const int rbound = (opts.lattType == CYCLE_LATT) ?
    opts.max : 2*opts.max+1;


  probMatrix = allocReal2D(rbound, rbound);
  if(!probMatrix)
    return NULL;

  totalprob = 0.0;
  for(m=0; m<rbound; m++){
    for(n=0; n<rbound; n++){
      double aux;

      aux = 0.0;
      for(j=0; j<2; j++){
	for(k=0; k<2; k++){

	  aux += matrix[j][k][m][n]*conj(matrix[j][k][m][n]);

	}/* end-for k */
      }/* end-for j */
      probMatrix[m][n] = aux;
      totalprob += probMatrix[m][n];

    }/* end-for n */
  }/* end-for m */

  if( fabs(totalprob-1.0) > WALK_TOL ){
    printf("Warning: probability of finding the particle in the lattice is %e\n",
	   totalprob);
  }

  return probMatrix;
}

double *getProbArray1D(double complex **matrix, options1D_t opts){
  int m,j;
  double *probMatrix, totalprob;
  const int rbound = (opts.lattType == LINE_LATT) ?
    2*opts.max+1 : opts.max;

  probMatrix = (opts.lattType == LINE_LATT) ?
    allocReal1D(2*opts.max+1) : allocReal1D(opts.max);
  if(!probMatrix)
    return NULL;

  totalprob = 0.0;
  for(m=0; m<rbound; m++){
    double aux;

    aux = 0.0;
    for(j=0; j<2; j++)
      aux += matrix[j][m]*conj(matrix[j][m]);
    probMatrix[m] = aux;
    totalprob += probMatrix[m];
    
  }

  if( fabs(totalprob-1.0) > WALK_TOL ){
    printf("Warning: probability of finding the particle in the lattice is %e\n",
	   totalprob);
  }

  return probMatrix;
}

int checkProb2D(double **matrix, int max, unsigned int lattType){
  int m,n;
  double totalprob;
  const int rbound = (lattType == CYCLE_LATT) ? max : 2*max+1;


  totalprob = 0.0;

  for(m=0; m<rbound; m++)
    for(n=0; n<rbound; n++)
      totalprob += matrix[m][n];

  if( fabs(totalprob-1.0) > WALK_TOL )
    return 0;

  return 1;
}

int checkProb1D(double *array, int max, unsigned int lattType){
  int m;
  double totalprob;
  const int rbound = (lattType == LINE_LATT) ? 2*max+1 : max;

  totalprob = 0.0;

  for(m=0; m<rbound; m++)
    totalprob += array[m];

  if( fabs(totalprob-1.0) > WALK_TOL )
    return 0;

  return 1;
}

int averageProbFromState1D(double **aveMatrix, double complex **state, 
			   options1D_t options){

  int m;
  static int flag=0;
  const int rbound = (options.lattType == LINE_LATT) ?
    2*options.max+1 : options.max;

  if(options.numOfExperiments<1 || flag>options.numOfExperiments )
    return 1;

  /* If the function is called for the first time */
  if(!flag){

    if(!aveMatrix)
      return 2;
    *aveMatrix = getProbArray1D(state, options);
    if(!*aveMatrix)
      return 3;
    for(m=0; m<rbound; m++)
      (*aveMatrix)[m] /= options.numOfExperiments;

    flag++;
    return 0;
  }

  /* Now, subsequent calls are treated */
  if(!*aveMatrix)
    return 4;

  for(m=0; m<rbound; m++){
    int j;
    double prob;

    prob = 0.0;
    for(j=0; j<2; j++)
      prob += state[j][m]*conj(state[j][m]);
      
    (*aveMatrix)[m] += (prob/options.numOfExperiments);
  }

  flag++;
  return 0;
}

int averageProbFromState2D(double ***aveMatrix, double complex ****state, 
			   options2D_t opts){

  int m,n;
  static int flag=0;
  const int rbound = (opts.lattType == CYCLE_LATT) ?
    opts.max : 2*opts.max+1;


  if(opts.numOfExperiments<1 || flag>opts.numOfExperiments)
    return 1;

  /* If the function is called for the first time */
  if(!flag){

    if(!aveMatrix)
      return 2;
    *aveMatrix = getProbArray2D(state, opts);
    if(!*aveMatrix)
      return 3;
    for(m=0; m<rbound; m++)
      for(n=0; n<rbound; n++)
	(*aveMatrix)[m][n] /= opts.numOfExperiments;
    
    flag++;
    return 0;
  }

  /* Now, subsequent calls are treated */
  if(!*aveMatrix)
    return 4;

  for(m=0; m<rbound; m++){
    for(n=0; n<rbound; n++){
      int j,k;
      double prob;
      
      prob = 0.0;
      for(j=0; j<2; j++)
	for(k=0; k<2; k++)
	  prob += state[j][k][m][n]*conj(state[j][k][m][n]);
      
      (*aveMatrix)[m][n] += (prob/opts.numOfExperiments);
    }
  }

  flag++;
  return 0;
}
