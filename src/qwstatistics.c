/* QWalk (qwstatistics.c) 
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
#include "qwstatistics.h"
#include "qwconsts.h"
#include "qwoptions_io.h"
#include "qwmem_real.h"


statistics_t getStatisticsFromState1D(double complex **matrix,  double *StatProb,
				      options1D_t opts, int iteration){
  statistics_t stat;
  double fstMoment, secMoment;
  int m;
  static double *SumProb=NULL;

  if(opts.max<1){
    stat.iteration=-1;
    return stat;
  }
  if(!matrix){
    stat.iteration=-2;
    return stat;
  }
  if(iteration<0){
    stat.iteration=-3;
    return stat;
  }

  if((iteration==1) && (opts.calcMix)){
    SumProb = (opts.lattType == LINE_LATT) ?
      allocReal1D(2*opts.max+1) : allocReal1D(opts.max);
    cleanReal1D(SumProb, (opts.lattType == LINE_LATT) ?
		2*opts.max+1: opts.max);
  }


  fstMoment = 0.0;
  secMoment = 0.0;

  if(opts.lattType == LINE_LATT){
    for(m=-opts.max; m<=opts.max; m++){
      double prob;
      int j;

      prob = 0.0;
      for(j=0; j<2; j++)
	prob += matrix[j][opts.max+m]*conj(matrix[j][opts.max+m]);

      fstMoment += m*prob;
      secMoment += m*m*prob;
    }
  }
  else{
    for(m=0; m<opts.max; m++){
      double prob;
      int j;

      prob = 0.0;
      for(j=0; j<2; j++)
	prob += matrix[j][m]*conj(matrix[j][m]);

      fstMoment += m*prob;
      secMoment += m*m*prob;
    }
  }

  stat.iteration = iteration;
  stat.meanX = fstMoment;
  stat.meanY = 0.0;
  stat.variance = secMoment - (fstMoment)*(fstMoment);

  /*
   *
   */
  stat.tvd = stat.tvdu = 0.0;

  if(opts.calcMix){
      const int rbound = (opts.lattType == LINE_LATT) ? 
	2*opts.max+1 : opts.max;
      const int auxSize = (opts.lattType == LINE_LATT) ? 
	2*(opts.max-opts.lattextra)+1 : opts.max;
      const double UnifProb = 1.0/auxSize;

      for(m=0; m<rbound; m++){
	int j;
	double prob = 0.0;
	
	for(j=0; j<2; j++)
	  prob += matrix[j][m]*conj(matrix[j][m]);
      
	SumProb[m] += prob;
	stat.tvd += fabs( StatProb[m] - SumProb[m]/(double)iteration );
	stat.tvdu += fabs( UnifProb - SumProb[m]/(double)iteration );
      }
  }

  if(iteration == opts.steps && opts.calcMix)
    free(SumProb);


  return stat;
}


statistics_t getStatisticsFromState2D(double complex ****matrix, double **StatProb,
				      options2D_t opts, int iteration){
  statistics_t stat;
  double fstMomentX, secMomentX, varianceX;
  double fstMomentY, secMomentY, varianceY;
  int m,n;
  static double **SumProb = NULL;
  const int auxsize = (opts.lattType == CYCLE_LATT) ? opts.max : 2*opts.max+1;
  const int lbound = (opts.lattType == CYCLE_LATT) ? 0 : -opts.max;
  const int rbound = (opts.lattType == CYCLE_LATT) ? opts.max-1 : opts.max;
  
  if(opts.max<1){
    stat.iteration=-1;
    return stat;
  }
  if(!matrix){
    stat.iteration=-2;
    return stat;
  }
  if(iteration<0){
    stat.iteration=-3;
    return stat;
  }

  if((iteration==1) && (opts.calcMix)){
    SumProb = allocReal2D(auxsize, auxsize);
    cleanReal2D(SumProb, auxsize, auxsize);
  }

  fstMomentX = 0.0;
  secMomentX = 0.0;
  fstMomentY = 0.0;
  secMomentY = 0.0;

  for(m=lbound; m<=rbound; m++){
    for(n=lbound; n<=rbound; n++){

      double prob;
      int j,k;
      int auxm = (opts.lattType == CYCLE_LATT) ? m : opts.max+m;
      int auxn = (opts.lattType == CYCLE_LATT) ? n : opts.max+n;

      prob = 0.0;
      for(j=0; j<2; j++){
	for(k=0; k<2; k++){
	  prob += matrix[j][k][auxm][auxn]*
	    conj(matrix[j][k][auxm][auxn]);
	}
      }

      fstMomentX += m*prob;
      secMomentX += m*m*prob;
      fstMomentY += n*prob;
      secMomentY += n*n*prob;     
    }
  }
  
  varianceX = secMomentX - fstMomentX*fstMomentX;
  varianceY = secMomentY - fstMomentY*fstMomentY;  

  stat.iteration = iteration;
  stat.meanX = fstMomentX;
  stat.meanY = fstMomentY;

  /* To calculate the variance in the 2D simulation we calculate the 
   * variances of both X and Y position and add them. It simplifies
   * the calculation and gives a good approximation of the expected
   * behaviour.
   */
  stat.variance = varianceX + varianceY;

  /*
   *
   */
  stat.tvd = stat.tvdu = 0.0;
  if(opts.calcMix){
    double UnifProb;
    const double diagArea = 2.0*pow((double)opts.max,2.0) - 2*opts.max + 1;
    const double natArea  = (2.0*opts.max-1.0)*(2.0*opts.max-1.0);
    const double cycArea  = opts.max*opts.max;
    if(opts.lattType == DIAG_LATT)
      UnifProb = pow(diagArea, -1.0);
    else if(opts.lattType == NATURAL_LATT)
      UnifProb = pow(natArea, -1.0);
    else
      UnifProb = pow(cycArea, -1.0);
    
    for(m=0; m<auxsize; m++){
      for(n=0; n<auxsize; n++){
	int j,k;
	double prob = 0.0;
	
	for(j=0; j<2; j++)
	  for(k=0; k<2; k++)
	    prob += matrix[j][k][m][n]*conj(matrix[j][k][m][n]);
      
	SumProb[m][n] += prob;
	stat.tvd += fabs( StatProb[m][n] - SumProb[m][n]/(double)iteration );
	if(opts.lattType != DIAG_LATT) 
	  stat.tvdu += fabs( UnifProb - SumProb[m][n]/(double)iteration );
	else if((m+n+1)%2)
	  stat.tvdu += fabs( UnifProb - SumProb[m][n]/(double)iteration );
	/* in the diagonal lattice, we only have sites (m,n) such
	 * that m+n is even.
	 */
      }
    }
  }

  if(iteration == opts.steps && opts.calcMix)
    freeReal2D(SumProb, auxsize);

  return stat;
}


/** CUIDADO! TEM QUE CONFERIR ESSA FUNCAO **/
statistics_t getStatisticsFromProb1D(double *matrix, int max, 
				     int iteration){
  statistics_t stat;
  double fstMoment, secMoment;
  int m;
  
  printf("Warning: Function getStatisticsFromProb1D should be tested. Contact franklin@lncc.br.\n");

  if(max<1){
    stat.iteration=-1;
    return stat;
  }
  if(!matrix){
    stat.iteration=-2;
    return stat;
  }
  if(iteration<0){
    stat.iteration=-3;
    return stat;
  }

  fstMoment = 0.0;
  secMoment = 0.0;

  for(m=-max; m<=max; m++){
    fstMoment += m*matrix[max+m];
    secMoment += m*m*matrix[max+m];
  }

  stat.iteration = iteration;
  stat.meanX = fstMoment;
  stat.meanY = 0.0;
  stat.variance = secMoment - (fstMoment)*(fstMoment);

  return stat;
}


/** CUIDADO! TEM QUE CONFERIR ESSA FUNCAO **/
statistics_t getStatisticsFromProb2D(double **matrix, int max, int iteration){
  statistics_t stat;
  double fstMomentX, secMomentX, varianceX;
  double fstMomentY, secMomentY, varianceY;
  int m,n;
  
  printf("Warning: Function getStatisticsFromProb2D should be tested. Contact franklin@lncc.br.\n");

  if(max<1){
    stat.iteration=-1;
    return stat;
  }
  if(!matrix){
    stat.iteration=-2;
    return stat;
  }
  if(iteration<0){
    stat.iteration=iteration;
    return stat;
  }


  fstMomentX = 0.0;
  secMomentX = 0.0;
  fstMomentY = 0.0;
  secMomentY = 0.0;

  for(m=-max; m<=max; m++){
    for(n=-max; n<=max; n++){
      fstMomentX += m*matrix[max+m][max+n];
      secMomentX += m*m*matrix[max+m][max+n];
      fstMomentY += n*matrix[max+m][max+n];
      secMomentY += n*n*matrix[max+m][max+n];     
    }
  }
  varianceX = secMomentX - fstMomentX*fstMomentX;
  varianceY = secMomentY - fstMomentY*fstMomentY;  

  stat.iteration = iteration;
  stat.meanX = fstMomentX;
  stat.meanY = fstMomentY;
  stat.variance = varianceX + varianceY;

  return stat;
}
