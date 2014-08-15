/* QWalk (qwmem_real.c) 
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
#include "qwmem_real.h"


double **allocReal2D(int dim1, int dim2){
  double **array;
  int i;

  if(dim1<1 || dim2<1)
    return NULL;

  array = (double **)malloc(dim1*sizeof(double *));
  if(!array)
    return NULL;

  for(i=0; i<dim1; i++){
    array[i] = (double *)malloc(dim2*sizeof(double));
    if(!array[i])
      return NULL;
  }

  return array;
}



double *allocReal1D(int dim){
  double *vector;
  
  if(dim<1)
    return NULL;

  vector = (double *)malloc((2*dim+1)*sizeof(double));
  return vector;
}



int freeReal2D(double **matrix, int dim1){
  int i;

  if(dim1<1)
    return 1;
  if(!matrix)
    return 2;

  for(i=0; i<dim1; i++){
    if(!matrix[i])
      return 2;

    free(matrix[i]);
  }
  free(matrix);

  return 0;
}



int copyReal2D(double **dest, double **src, int dim1, int dim2){
  register int i,j;

  if(dim1<1 || dim2<1)
    return 1;
  if(!dest || !src)
    return 2;

  for(i=0; i<dim1; i++){
    if(!dest[i] || !src[i])
      return 2;

    for(j=0; j<dim2; j++)
      dest[i][j] = src[i][j];

  }

  return 0;
}



int cleanReal2D(double **matrix, int dim1, int dim2){
  register int i,j;

  if(dim1<1 || dim2<1)
    return 1;
  if(!matrix)
    return 2;

  for(i=0; i<dim1; i++){
    if(!matrix[i])
      return 2;

    for(j=0; j<dim2; j++)
      matrix[i][j] = 0.0;
  }

  return 0;
}



int cleanReal1D(double *vector, int dim){
  register int i;

  if(dim<1)
    return 1;
  if(!vector)
    return 2;

  for(i=0; i<dim; i++)
    vector[i] = 0.0;

  return 0;
}


