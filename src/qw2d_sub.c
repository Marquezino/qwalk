/* QWalk (qw2d_sub.c) Subroutines for qw2d
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
#include "qwmem_real.h"
#include "qwcoin.h"
#include "qwstate.h"
#include "qwstatistics.h"
#include "qwprob.h"
#include "qwcoin_io.h"
#include "qwstate_io.h"
#include "qwstatistics_io.h"
#include "qwconsts.h"
#include "qwmeasure.h"
#include "qw2d_sub.h"

void setCoin2D(double complex *****C,int coinType,const char *filename){
  static char previousAlloc = 0;

  if(previousAlloc)
    freeComplex4D(*C, 2, 2, 2);

  *C = NULL;

  switch(coinType){
  case CUSTOM_COIN:
    *C = readCoinFile2D(filename);
    break;
  case FOURIER_COIN:
    *C = createFourierCoin2D();
    break;
  case GROVER_COIN:
    *C = createGroverCoin2D();
    break;
  case HADAMARD_COIN:
    *C = createHadamardCoin2D();
    break;
  }

  if(*C)
    previousAlloc = 1;
  else
    previousAlloc = 0;

  return;
}


void setState2D(double complex *****A,options2D_t opts,const char *filename){
  static char previousAlloc = 0;

  if(previousAlloc){
    if(opts.lattType == CYCLE_LATT)
      freeComplex4D(*A, 2, 2, opts.max);
    else
      freeComplex4D(*A, 2, 2, 2*opts.max+1);
  }

  *A = NULL;

  switch(opts.stateType){
  case CUSTOM_STATE:
    *A = readStateFile2D(filename,opts.max,opts.lattType);
    break;
  case FOURIER_STATE:
    *A = createFourierState2D(opts.max,opts.lattType);
    break;
  case GROVER_STATE:
    *A = createGroverState2D(opts.max,opts.lattType);
    break;
  case HADAMARD_STATE:
    *A = createHadamardState2D(opts.max,opts.lattType);
    break;
  }

  if(*A)
    previousAlloc = 1;
  else
    previousAlloc = 0;

  return;
}




void check2D(double complex ****A, options2D_t options, int iteration){

  if(options.checkState && !checkState2D(A, options)){
    printf("Error: state norm is not unitary in iteration %d.\n", iteration);
    exit(EXIT_FAILURE);
  }
  if(options.checkXSymmetry && !checkXSymmetry2D(A, options.max)){
    printf("Error: wave-equation is not symmetrical in iteration %d.\n", iteration);
    exit(EXIT_FAILURE);
  }
  if(options.checkYSymmetry && !checkYSymmetry2D(A, options.max)){
    printf("Error: wave-equation is not symmetrical in iteration %d.\n", iteration);
    exit(EXIT_FAILURE);
  }
  
  return;
}




void iterate2D(double complex *****A,double complex *****Atemp,double complex ****C, 
	       int ****BLinks1, int ****BLinks2, options2D_t opts, int iteration){
  int m, n, error; 
  double complex ****aux;


  /* We define constants MAX and LATTEXTRA as shorts for options.max and
   * options.lattextra, respectively.
   */
  const int MAX = opts.max;
  const int LATTEXTRA = opts.lattextra;
  const int auxsize = (opts.lattType == CYCLE_LATT) ? opts.max : 2*opts.max+1;

  /* In the n-th iteration the walker cannot be farther than n sites from 
   * its initial position. Therefore, we don't need to update the entire 
   * lattice, but only a square region (lbound,rbound)X(lbound,rbound).
   * Constant MAX represents the size of the lattice and constant LATTEXTRA 
   * represents an extra space in the lattice, used to avoid "bound errors" 
   * when accessing the array.
   */
  const int lbound = (opts.lattType == CYCLE_LATT) ? 
    0 : MAXIMUM(MAX-LATTEXTRA-iteration, 1);
  const int rbound = (opts.lattType == CYCLE_LATT) ? 
    MAX-1 : MINIMUM(MAX+LATTEXTRA+iteration, 2*MAX-1);

  error = cleanComplex4D(*Atemp, 2, 2, auxsize, auxsize);
  if(error){
    printf("Error: could not clean temporary matrix in iteration %d.\n", iteration);
    exit(EXIT_FAILURE);
  }

  switch(opts.lattType){
    /* This "switch" looks ugly, but it is better for performance. It 
     * is a good programming practice to avoid "if"s inside loops.
     */

  case DIAG_LATT:
    for(m = lbound; m <= rbound; m++){
      for(n = lbound; n <= rbound; n++){
	int j,k;

	for(j=0; j<2; j++){
	  for(k=0; k<2; k++){
	    int L1,L2, jprime,kprime;
	    double complex newValue; 

	    /* Further information on the matrices of broken links
	     * can be found in Physical Review A, 74, 012312 (2006)
	     */
	    L1 = BLinks1[j][k][m][n];
	    L2 = BLinks2[j][k][m][n];

	    newValue = 0.0;
	    for(jprime=0; jprime<2; jprime++){
	      for(kprime=0; kprime<2; kprime++){

		newValue += C[j+L1][k+L2][jprime][kprime]*
		  (*A)[jprime][kprime][m+L1][n+L2];

	      }/* End-for kprime */
	    }/* End-for jprime */

	    (*Atemp)[1-j][1-k][m][n] = newValue;

	  }/* End-for k */
	}/* End-for j */

      }/* End-for n */
    }/* End-for m */
    break;

  case NATURAL_LATT:
    for(m = lbound; m <= rbound; m++){
      for(n = lbound; n <= rbound; n++){
	int j,d;

	for(j=0; j<2; j++){
	  for(d=0; d<2; d++){
	    int L,jprime,dprime;
	    double complex newValue; 

	    L = BLinks1[j][d][m][n];
	    newValue = 0.0;
	    for(jprime=0; jprime<2; jprime++){
	      for(dprime=0; dprime<2; dprime++){

		newValue += C[j+L][abs(d+L)%2][jprime][dprime]*
		  (*A)[jprime][dprime][m + L*(1-DELTA(j,d))][n + L*DELTA(j,d)];

	      }/* End-for kprime */
	    }/* End-for jprime */
	    (*Atemp)[1-j][1-d][m][n] = newValue;
	    

	  }/* End-for k */
	}/* End-for j */

      }/* End-for n */
    }/* End-for m */
    break;

  case CYCLE_LATT:
    for(m = lbound; m <= rbound; m++){
      for(n = lbound; n <= rbound; n++){
	int j,d;

	for(j=0; j<2; j++){
	  for(d=0; d<2; d++){
	    int L,jprime,dprime;
	    double complex newValue; 

	    L = BLinks1[j][d][m][n];
	    newValue = 0.0;
	    for(jprime=0; jprime<2; jprime++){
	      for(dprime=0; dprime<2; dprime++){

		newValue += C[j][d][jprime][dprime]*(*A)[jprime][dprime][m][n];

	      }/* End-for kprime */
	    }/* End-for jprime */
	    (*Atemp)[1-(j+L)][1-abs(d+L)%2][(MAX + m + L*(1-DELTA(j,d)))%MAX][(MAX + n + L*DELTA(j,d))%MAX] = newValue;

	  }/* End-for k */
	}/* End-for j */

      }/* End-for n */
    }/* End-for m */
    break;

  default:
    printf("Error: invalid lattice type for two-dimensional simulation");
    exit(EXIT_FAILURE);
        
  }/* end-switch */
  
  /* Now we quicky exchange matrices A and Atemp */
  aux = *A;
  *A = *Atemp;
  *Atemp = aux;

  return;

}



double **getStationary2D(double complex ****A, double complex ****C, 
			 int ****BLinks1, int ****BLinks2, 
			 options2D_t options){
  int m,n,t;
  double complex ****Atemp;
  double **stationary;
  const int MAX = options.max;
  const int rbound = (options.lattType == CYCLE_LATT) ? MAX : 2*MAX+1;

  stationary =  allocReal2D(rbound, rbound);
  if(!stationary)
    return NULL;

  cleanReal2D(stationary, rbound, rbound);

  Atemp = allocComplex4D(2, 2, rbound, rbound);
  if(!Atemp)
    return NULL;

  for(t=0; t<options.stepsMix; t++){

    iterate2D(&A, &Atemp, C, BLinks1, BLinks2, options, t);
    for(m=0; m<rbound; m++){
      for(n=0; n<rbound; n++){
	int j,k;
	double prob = 0.0;

	for(j=0; j<2; j++)
	  for(k=0; k<2; k++)
	    prob += A[j][k][m][n]*conj(A[j][k][m][n]);

	stationary[m][n] += prob;
      }
    }

  }

  freeComplex4D(Atemp, 2, 2, rbound);

  for(m=0; m<rbound; m++)
    for(n=0; n<rbound; n++)
      stationary[m][n] /= options.stepsMix;

  if(!checkProb2D(stationary,MAX,options.lattType))
    return NULL;

  return stationary;
}




void doStatistics2D(double complex ****A, double **StatProb, options2D_t options, 
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
  stat = getStatisticsFromState2D(A, StatProb, options, iteration);
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
    error = writeStatistics(fnames.sta_file,vStat[iteration]);
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


void randomBrokenLink2D(int *****BLinks1, int *****BLinks2, options2D_t options){
  int m,n;

  /* We define constants MAX and LATTEXTRA as shorts for options.max and
   * options.lattextra, respectively.
   */
  const int MAX = options.max;

  if(!*BLinks1)
    return;
  if(!*BLinks2)
    return;

  if(MAX<1)
    return;

  if((options.blProbA < 0.0) || options.blProbB < 0.0)
    return;

  if(options.lattType == DIAG_LATT){

    for(m=0; m<2*MAX; m++){
      for(n=0; n<2*MAX; n++){
	if(rand() < options.blProbA*RAND_MAX){
	  (*BLinks1)[0][0][ m ][ n ] = 0;
	  (*BLinks1)[1][1][m+1][n+1] = 0;
	  (*BLinks2)[0][0][ m ][ n ] = 0;
	  (*BLinks2)[1][1][m+1][n+1] = 0;
	}
      }
    }
    for(m=1; m<=2*MAX; m++){
      for(n=0; n<2*MAX; n++){
	if(rand() < options.blProbB*RAND_MAX){
	  (*BLinks1)[1][0][ m ][ n ] = 0;
	  (*BLinks1)[0][1][m-1][n+1] = 0;
	  (*BLinks2)[1][0][ m ][ n ] = 0;
	  (*BLinks2)[0][1][m-1][n+1] = 0;
	}
      }
    }
    
  }
  else if(options.lattType == NATURAL_LATT){

    for(m=0; m<2*MAX; m++){
      for(n=0; n<2*MAX; n++){
	if(rand() < options.blProbA*RAND_MAX){
	  (*BLinks1)[0][1][ m ][ n ] = 0;
	  (*BLinks1)[1][0][m+1][ n ] = 0;
	}
	if(rand() < options.blProbB*RAND_MAX){
	  (*BLinks1)[0][0][ m ][ n ] = 0;
	  (*BLinks1)[1][1][ m ][n+1] = 0;
	}
      }
    }
    for(m=0; m<2*MAX; m++){
      if(rand() < options.blProbA*RAND_MAX){
	(*BLinks1)[0][1][ m ][2*MAX] = 0;
	(*BLinks1)[1][0][m+1][2*MAX] = 0;
      }
    }
    for(n=0; n<2*MAX; n++){
      if(rand() < options.blProbB*RAND_MAX){
	(*BLinks1)[0][0][2*MAX][ n ] = 0;
	(*BLinks1)[1][1][2*MAX][n+1] = 0;
      }
    }

  }
  else{  /* if(options.lattType == CYCLE_LATT) */

    for(m=0; m<MAX; m++){
      for(n=0; n<MAX; n++){
	if(rand() < options.blProbA*RAND_MAX){
	  (*BLinks1)[0][1][ m ][ n ] = 0;
	  (*BLinks1)[1][0][(m+1)%MAX][ n ] = 0;
	}
	if(rand() < options.blProbB*RAND_MAX){
	  (*BLinks1)[0][0][ m ][ n ] = 0;
	  (*BLinks1)[1][1][ m ][(n+1)%MAX] = 0;
	}
      }
    }
    
  }

  
  return;
}
