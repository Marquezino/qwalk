/* QWalk (qw1d_sub.c) Sub-routines for qw1d
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
#include "qwmem_complex.h"
#include "qwcoin.h"
#include "qwstate.h"
#include "qwstatistics.h"
#include "qwstatistics_io.h"
#include "qwcoin_io.h"
#include "qwstate_io.h"
#include "qwconsts.h"
#include "qw1d_sub.h"
#include "qwprob.h"


void setCoin1D(double complex ***C, int coinType, const char *filename){
  static char previousAlloc = 0;

  if(previousAlloc)
    freeComplex2D(*C, 2);

  *C = NULL;

  switch(coinType){
  case HADAMARD_COIN:
    *C = createHadamardCoin1D();
    break;
  case CUSTOM_COIN:
    *C = readCoinFile1D(filename);
    break;
  }

  if(*C)
    previousAlloc = 1;
  else
    previousAlloc = 0;

  return;
}



void setState1D(double complex ***A, options1D_t options, 
		const char *filename){
  static char previousAlloc = 0;

  if(previousAlloc)
    freeComplex2D(*A, 2);

  *A = NULL;

  switch(options.stateType){
  case HADAMARD_STATE:
    *A = createHadamardState1D(options.max, options.lattType);
    break;
  case CUSTOM_STATE:
    *A = readStateFile1D(filename, options.max, options.lattType);
    break;
  }

  if(*A)
    previousAlloc = 1;
  else
    previousAlloc = 0;

  return;
}



void randomBrokenLink1D(int ***B, options1D_t options){
  int m;

  if(!*B)
    return;

  if(options.max<1)
    return;

  if(options.blProb < 0.0)
    return;

  switch(options.lattType){
  case LINE_LATT:
    for(m=0; m<2*options.max; m++){
      if(rand() < options.blProb*RAND_MAX){
	(*B)[0][m]=0;
	(*B)[1][m+1]=0;
      }
    }
    break;

  case SEGMENT_LATT:
    for(m=0; m<options.max-1; m++){
      if(rand() < options.blProb*RAND_MAX){
	(*B)[0][m]=0;
	(*B)[1][m+1]=0;
      }
    }
    break;

  case CYCLE_LATT:
    for(m=0; m<options.max; m++){
      if(rand() < options.blProb*RAND_MAX){
	(*B)[0][m]=0;
	(*B)[1][(m+1)%(options.max)]=0;
      }
    }
    break;

  default:
    printf("Error: invalid lattice type in broken link initialization.\n");
    exit(EXIT_FAILURE);
    break;
    
  }
  
  return;
}



void iterate1D(double complex ***A, double complex ***Atemp, double complex **C, 
	       int **BrokenLinks, options1D_t options, int iteration){
  int m,j,error;
  double complex **aux;
  
  /* We define constants MAX and LATTEXTRA as shorts for options.max and
   * options.lattextra, respectively.
   */
  const int MAX = options.max;
  const int LATTEXTRA = options.lattextra;

  /* In the n-th iteration the walker cannot be farther than n sites from 
   * its initial position. Therefore, we don't need to update the entire 
   * lattice, but only an interval (lbound,rbound). Constant MAX represents 
   * the size of the lattice and constant LATTEXTRA represents an extra space 
   * in the lattice, used to avoid "bound errors" when accessing the array.
   */
  const int lbound = (options.lattType == LINE_LATT) ? 
    MAXIMUM(MAX-LATTEXTRA-iteration, 1) : 0;
  const int rbound = (options.lattType == LINE_LATT) ? 
    MINIMUM(MAX+LATTEXTRA+iteration, 2*MAX-1) : MAX-1;

  if(options.lattType == LINE_LATT)
    error = cleanComplex2D(*Atemp, 2, 2*MAX+1);
  else
    error = cleanComplex2D(*Atemp, 2, MAX);
  if(error){
    printf("Error: could not clean temporary matrix in iteration %d.\n", 
	   iteration);
    exit(EXIT_FAILURE);
  }

  switch(options.lattType){
    /* This "switch" looks ugly, but it is better for performance. It 
     * is a good programming practice to avoid "if"s inside loops.
     */
    
  case LINE_LATT:
    for(m=lbound; m<=rbound; m++){
      for(j=0; j<2; j++){
	int L, k;

	/* Further information on the matrix of broken links
	 * can be found in Physical Review A, 74, 012312 (2006) 
	 */      
	L = BrokenLinks[j][m];
	for(k=0; k<2; k++){
	  (*Atemp)[1-j][m] += C[j+L][k]* (*A)[k][m+L];
	}/* End-for k */
      }/* End-for j */
    }/* End-for m*/
    break;
  
  case CYCLE_LATT:
    for(m=0; m<MAX; m++){
      for(j=0; j<2; j++){
	int L, k;

	L = BrokenLinks[j][m];
	for(k=0; k<2; k++){
	  (*Atemp)[1-j][m] += C[j+L][k]* (*A)[k][(MAX+m+L)%MAX];
	}/* End-for k */
      }/* End-for j */
    }/* End-for m*/
    break;

  case SEGMENT_LATT:
    BrokenLinks[0][MAX-1] = 0;
    BrokenLinks[1][0]     = 0;
    for(m=0; m<MAX; m++){
      for(j=0; j<2; j++){
	int L, k;
	L = BrokenLinks[j][m];
	for(k=0; k<2; k++){
	  (*Atemp)[1-j][m] += C[j+L][k]* (*A)[k][m+L];
	}/* End-for k */
      }/* End-for j */
    }/* End-for m*/
    break;

  default:
    printf("Error: invalid lattice type for one-dimensional simulation");
    exit(EXIT_FAILURE);

  }/* end-switch*/

  /* Now we quicky exchange vectors A and Atemp */
  aux = *A;
  *A = *Atemp;
  *Atemp = aux;

  return;
}



void doStatistics1D(double complex **A, double *StatProb, options1D_t options, 
		    filenames_t fnames, int iteration, int experiment){
  int error;
  statistics_t stat;
  static statistics_t *vStat = NULL;

  if(iteration<1 || iteration>options.steps){
    printf("Error: unexpected number of steps when calculating statistics.\n");
    exit(EXIT_FAILURE);
  }
  if(experiment<1 || experiment>options.numOfExperiments){
    printf("Error: unexpected number of experiments when calculating statistics.\n");
    exit(EXIT_FAILURE);
  }

  if((experiment == 1) && (iteration == 1)){
    /* The first call of this subroutine, i.e., first iteration of
     * first experiment
     */
    register int t;

    vStat = (statistics_t *)malloc((options.steps+1)*sizeof(statistics_t));
    for(t=1; t<=options.steps; t++){
      /* Initializing the array of structures */
      vStat[t].iteration = t;
      vStat[t].variance = 0.0;
      vStat[t].meanX = 0.0;
      vStat[t].meanY = 0.0;
      vStat[t].tvd = 0.0;
      vStat[t].tvdu = 0.0;
    }
  }

  if(!vStat){
    printf("Error: invalid statistics vector\n");
    exit(EXIT_FAILURE);
  }
  
  /* Here we get the statistics for this step */
  stat = getStatisticsFromState1D(A, StatProb, options, iteration);
  if(stat.iteration<0){
    printf("Error: could not generate statistics.\n");
    exit(EXIT_FAILURE);
  }

  /* Here we accumulate the results of this iteration in 
   * the structure 
   */
  vStat[iteration].variance += stat.variance;
  vStat[iteration].meanX += stat.meanX;
  vStat[iteration].meanY += stat.meanY;
  vStat[iteration].tvd += stat.tvd;
  vStat[iteration].tvdu += stat.tvdu;

  if(experiment == options.numOfExperiments){
    /* Here we divide the previously accumulated values by
     * the number of experiments in order to obtain the
     * averages.
     */
    vStat[iteration].variance /= options.numOfExperiments;
    vStat[iteration].meanX /= options.numOfExperiments;
    vStat[iteration].meanY /= options.numOfExperiments;    
    vStat[iteration].tvd /= options.numOfExperiments;
    vStat[iteration].tvdu /= options.numOfExperiments;

    /* Here we write the statistics in the appropriate file */
    error = writeStatistics(fnames.sta_file, vStat[iteration]);
    if(error){
      printf("Error: could not write statistics.\n");
      exit(EXIT_FAILURE);
    }
  }

  if((experiment == options.numOfExperiments) && (iteration == options.steps)){
    /* The last call of this subroutine, i.e., last iteration 
     * of last experiment 
     */
    free(vStat);
    vStat = NULL;
  }
  
  return;
}



void check1D(double complex **A, options1D_t options, int iteration){

  if(options.checkState && !checkState1D(A, options)){
    printf("Error: state norm is not unitary in iteration %d.\n", iteration);
    exit(EXIT_FAILURE);
  }
  if(options.checkSymmetry && !checkSymmetry1D(A, options.max)){
    printf("Error: wave-equation is not symmetrical in iteration %d.\n", iteration);
    exit(EXIT_FAILURE);
  }
  
  return;
}


double *getStationary1D(double complex **A, double complex **C, 
			int **BLinks, options1D_t options){
  int m,t;
  double complex **Atemp;
  double *stationary;
  const int MAX = options.max;
  const int rbound = (options.lattType == LINE_LATT) ? 2*MAX+1 : MAX;

  stationary = (options.lattType == LINE_LATT) ?
    allocReal1D(2*MAX+1) : allocReal1D(MAX);

  if(!stationary)
    return NULL;

  if(options.lattType == LINE_LATT)
    cleanReal1D(stationary, 2*MAX+1);
  else
    cleanReal1D(stationary, MAX);

  Atemp = (options.lattType == LINE_LATT) ?
    allocComplex2D(2, 2*MAX+1) : allocComplex2D(2, MAX);

  if(!Atemp)
    return NULL;

  for(t=0; t<options.stepsMix; t++){
    iterate1D(&A, &Atemp, C, BLinks, options, t);
    for(m=0; m<rbound; m++){
      int j;
      double prob = 0.0;

      for(j=0; j<2; j++)
	prob += A[j][m]*conj(A[j][m]);

      stationary[m] += prob;
    }
  }

  freeComplex2D(Atemp, 2);

  for(m=0; m<rbound; m++)
    stationary[m] /= options.stepsMix;

  if(!checkProb1D(stationary,MAX, options.lattType))
    return NULL;

  return stationary;
}
