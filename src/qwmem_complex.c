/* QWalk (qwmem_complex.c) 
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


double complex ****allocComplex4D(int dim1, int dim2, int dim3, int dim4){
  double complex ****array;
  int i,j,k;

  if(dim1<1 || dim2<1 || dim3<1 || dim4<1)
    return NULL;

  array = (double complex ****)malloc(dim1*sizeof(double complex ***));
  if(!array)
    return NULL;

  for(i=0; i<dim1; i++){
    array[i] = (double complex ***)malloc(dim2*sizeof(double complex **));
    if(!array[i])
      return NULL;

    for(j=0; j<dim2; j++){
      array[i][j] = (double complex **)malloc(dim3*sizeof(double complex*));
      if(!array[i][j])
	return NULL;

      for(k=0; k<dim3; k++){
	array[i][j][k] = (double complex *)malloc(dim4*sizeof(double complex));
	if(!array[i][j][k])
	  return NULL;
      }/* end-for k */
    }/* end-for j */
  }/* end-for i */
  
  return array;
}



double complex ***allocComplex3D(int dim1, int dim2, int dim3){
  double complex ***array;
  int i,j;

  if(dim1<1 || dim2<1 || dim3<1)
    return NULL;

  array = (double complex ***)malloc(dim1*sizeof(double complex **));
  if(!array)
    return NULL;

  for(i=0; i<dim1; i++){
    array[i] = (double complex **)malloc(dim2*sizeof(double complex *));
    if(!array[i])
      return NULL;

    for(j=0; j<dim2; j++){
      array[i][j] = (double complex *)malloc(dim3*sizeof(double complex));
      if(!array[i][j])
	return NULL;
    }/* end-for j */
  }/* end-for i */
  
  return array;
}


double complex **allocComplex2D(int dim1, int dim2){
  double complex **array;
  int i;

  if(dim1<1 || dim2<1)
    return NULL;

  array = (double complex **)malloc(dim1*sizeof(double complex *));
  if(!array)
    return NULL;

  for(i=0; i<dim1; i++){
    array[i] = (double complex *)malloc(dim2*sizeof(double complex));
    if(!array[i])
      return NULL;
  }
  
  return array;
}



int freeComplex4D(double complex ****array, int dim1, int dim2, int dim3){
  int i,j,k;

  if(dim1<1 || dim2<1 || dim3<1)
    return 1;
  if(!array)
    return 2;

  for(i=0; i<dim1; i++){
    if(!array[i])
      return 2;

    for(j=0; j<dim2; j++){
      if(!array[i][j])
	return 2;

      for(k=0; k<dim3; k++){
	if(!array[i][j][k])
	  return 2;

	free(array[i][j][k]);
      }/* end-for k */

      free(array[i][j]);
    }/* end-for j */

    free(array[i]);
  }/* end-for i */

  free(array);

  return 0;
}


int freeComplex3D(double complex ***array, int dim1, int dim2){
  int i,j;

  if(dim1<1 || dim2<1)
    return 1;
  if(!array)
    return 2;

  for(i=0; i<dim1; i++){
    if(!array[i])
      return 2;

    for(j=0; j<dim2; j++){
      if(!array[i][j])
	return 2;

      free(array[i][j]);
    }/* end-for j */

    free(array[i]);
  }/* end-for i */

  free(array);

  return 0;
}


int freeComplex2D(double complex **matrix, int dim1){
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


int copyComplex4D(double complex ****dest, double complex ****src, 
		  int dim1, int dim2, int dim3, int dim4){

  register int i,j,k,m;

  if(dim1<1 || dim2<1 || dim3<1 || dim4<1)
    return 1;
  if(!dest || !src)
    return 2;

  for(i=0; i<dim1; i++){
    if(!dest[i] || !src[i])
      return 2;

    for(j=0; j<dim2; j++){
      if(!dest[i][j] || !src[i][j])
	return 2;      

      for(k=0; k<dim3; k++){
	if(!dest[i][j][k] || !src[i][j][k])
	  return 2;

	for(m=0; m<dim4; m++)
	  dest[i][j][k][m] = src[i][j][k][m];

      }/* end-for k */
    }/* end-for j */
  }/* end-for i */

  return 0;
}



int copyComplex3D(double complex ***dest, double complex ***src, int dim1, int dim2, int dim3){
  register int i,j,k;

  if(dim1<1 || dim2<1 || dim3<1)
    return 1;
  if(!dest || !src)
    return 2;

  for(k=0; k<dim1; k++){
    if(!dest[k] || !src[k])
      return 2;

    for(i=0; i<dim2; i++){
      if(!dest[k][i] || !src[k][i])
	return 2;

      for(j=0; j<dim3; j++)
	dest[k][i][j] = src[k][i][j];

    }/* end-for i */
  }/* end-for k */

  return 0;
}


int copyComplex2D(double complex **dest, double complex **src, int dim1, int dim2){
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


int cleanComplex4D(double complex ****matrix, int dim1, int dim2, int dim3, int dim4){
  register int i,j,k,m;

  if(dim1<1 || dim2<1 || dim3<1 || dim4<1)
    return 1;
  if(!matrix)
    return 2;

  for(i=0; i<dim1; i++){
    if(!matrix[i])
      return 2;

    for(j=0; j<dim2; j++){
      if(!matrix[i][j])
	return 2;

      for(k=0; k<dim3; k++){
	if(!matrix[i][j][k])
	  return 2;

	for(m=0; m<dim4; m++)
	  matrix[i][j][k][m] = 0.0;

      }/* end-for k */
    }/* end-for j */
  }/* end-for i */

  return 0;
}



int cleanComplex3D(double complex ***matrix, int dim1, int dim2, int dim3){
  register int i,j,k;

  if(dim1<1 || dim2<1 || dim3<1)
    return 1;
  if(!matrix)
    return 2;

  for(k=0; k<dim1; k++){
    if(!matrix[k])
      return 2;

    for(i=0; i<dim2; i++){
      if(!matrix[k][i])
	return 2;

      for(j=0; j<dim3; j++)
	matrix[k][i][j] = 0.0;

    }
  }

  return 0;
}



int cleanComplex2D(double complex **vector, int dim1, int dim2){
  register int i, j;

  if(dim1<1 || dim2<1)
    return 1;
  if(!vector)
    return 2;

  for(i=0; i<dim1; i++){
    if(!vector[i])
      return 2;

    for(j=0; j<dim2; j++)
      vector[i][j] = 0.0;
  }

  return 0;
}
