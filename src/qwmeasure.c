/* QWalk (qwmeasure.c) 
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
#include "qwmeasure.h"
#include "qwoptions_io.h"
#include "qwconsts.h" 


int measureState2D(double complex *****A, double complex *****Atemp, 
		   options2D_t opts){
  int j, k, m, n;
  int det, dice, result, error;
  double *p, *sp;
  double complex ****aux;

  int max = opts.max;
  int detectors = opts.detectors;
  int **detector_pts = opts.detector_pts;

  const int lbound   = (opts.lattType == CYCLE_LATT) ? 0 : -max;
  const int rbound   = (opts.lattType == CYCLE_LATT) ? max-1 : max;
  const int aux_size = (opts.lattType == CYCLE_LATT) ? max : 2*max+1;

  /* Detectors are described by a collection M_m of measurement 
   * operators which satisfy the completeness equation. The index m
   * indicates the result of the measurement. Before state |Psi> is 
   * measured the probability that result m occurs is given by 
   * p(m)= <Psi| Mdagger_m M_m |Psi>, and after the measurement the state 
   * collapses to |Psi'>= 1/sqrt(p(m)) M_m|Psi>. In qw2d we may define an 
   * arbitrary number of detectors specifying a list of coordinates 
   * (m_1, n_1), ..., (m_N,n_N) [that is, int **detector_pts]. The measurement 
   * operators are M_0 = I_4 tensor I_\infty - sum(M_i, i from 1 to N) and 
   * M_i = I_4 tensor |m_i, n_i><m_i, n_i|, for 1 <= i <= N. 
   *
   * To perform the measurement we (1) generate a random number r; and 
   * (2) determine an integer m, 0 <= m <= N, such that
   * sum(p(j), j from 0 to m-1) <= m < sum(p(j), j from 0 to m)
   */


  if( (!*A) || (!*Atemp))
    return -1;
  if(max<1)
    return -2;
  if(detectors<0)
    return -3;
  if(!detector_pts)
    return -4;

  p = allocReal1D(detectors+1);
  if(!p)
    return -5;
  error = cleanReal1D(p, detectors+1);
  if(error)
    return -6;

  sp = allocReal1D(detectors+1);
  if(!sp)
    return -5;
  error = cleanReal1D(sp, detectors+1);
  if(error)
    return -6;

  /* Here we evaluate the probability of each detection... */
  for(det=1; det<=detectors; det++ ){  
    if(opts.lattType == CYCLE_LATT){
      m = detector_pts[det][0];
      n = detector_pts[det][1];
    }else{
      m = max + detector_pts[det][0];
      n = max + detector_pts[det][1];
    }

    for(j=0; j<2; j++)
      for(k=0; k<2; k++)
	p[det] += (*A)[j][k][m][n]*conj((*A)[j][k][m][n]);
  }

  /* ...and the probability of measuring the complement */
  p[0] = 1.0;
  for(j=1; j<=detectors; j++)
    p[0] -= p[j];

  /* ...  */
  sp[0] = p[0];
  for(det=1; det<=detectors; det++)
    sp[det] = sp[det-1] + p[det];

  /* We get a new random number */
  dice = rand();

  /* Based on this number we identify the corresponding result */
  for(result=0; result<=detectors; result++)
    if(dice < sp[result]*RAND_MAX) break;

  for(m=lbound; m<rbound; m++){
    for(n=lbound; n<rbound; n++){
      int delta;
      const int auxm = (opts.lattType == CYCLE_LATT) ? m : max+m;
      const int auxn = (opts.lattType == CYCLE_LATT) ? n : max+n;

      if(result==0){
	delta = 1;
	for(det=1; det<=detectors; det++)
	  delta *= (1-DELTA(m, detector_pts[det][0])*
		    DELTA(n, detector_pts[det][1]));
      }
      else{
	delta = DELTA(m, detector_pts[result][0])*
	  DELTA(n, detector_pts[result][1]);
      }
      
      for(j=0; j<2; j++)
	for(k=0; k<2; k++)
	  (*Atemp)[j][k][auxm][auxn] = 
	    delta*((*A)[j][k][auxm][auxn])/sqrt(p[result]);
    }
  }

  free(p);
  free(sp);

  /* Now we quicky exchange matrices A and Atemp */
  aux = *A;
  *A = *Atemp;
  *Atemp = aux;

  error = cleanComplex4D(*Atemp, 2, 2, aux_size, aux_size);
  if(error)
    return -6;

  return result;
}



int randMeasure2D(double complex *****A, options2D_t opts){

  int markX, markY, marked;
  int m, n, diceA, diceB;
  float sp;

  const int rbound = (opts.lattType == CYCLE_LATT) ? 
    opts.max-1 : 2*opts.max;

  /* We get a new random number */
  diceA = rand();

  /* Now we scan the lattice checking if any point is measured */
  markX = markY = -1;
  marked=0;
  sp=0.0;
  for(m=0; m<=rbound; m++){
    for(n=0; n<=rbound; n++){
      int j,k;
      float prob;

      prob=0.0;
      for(j=0; j<2; j++)
	for(k=0; k<2; k++)
	  prob += (*A)[j][k][m][n]*conj((*A)[j][k][m][n]);

      diceB = rand();
      if(diceB < opts.dtProb*RAND_MAX){ /* if site should be measured */
	sp += prob;
	if(diceA < sp*RAND_MAX && !marked){ /* if state collapsed to this site */
	  markX = m;
	  markY = n;
	  marked = 1;
	}
	else{
	  for(j=0; j<2; j++)
	    for(k=0; k<2; k++)
	      (*A)[j][k][m][n] = 0.0;
	}
      }

    }/* end-for n */
  }/* end-for m */

  for(m=0; m<=rbound; m++){
    for(n=0; n<=rbound; n++){
      int j,k;

      if(markX==-1 && markY==-1){ /* measured the complement */
	for(j=0; j<2; j++)
	  for(k=0; k<2; k++)
	    (*A)[j][k][m][n] *= pow(1.0-sp, -0.5);

      }else if(markX==m && markY==n){ /* measured this site */
	float prob;

	prob = 0.0;
	for(j=0; j<2; j++)
	  for(k=0; k<2; k++)
	    prob += (*A)[j][k][m][n]*conj((*A)[j][k][m][n]);

	for(j=0; j<2; j++)
	  for(k=0; k<2; k++)
	    (*A)[j][k][m][n] *= pow(prob, -0.5);
	

      }else{ /* measured some other site */

	for(j=0; j<2; j++)
	  for(k=0; k<2; k++)
	    (*A)[j][k][m][n] = 0.0;

      }

    }/* end-for n */
  }/* end-for m */

  return 0;
}



int measureState1D(double complex ***A, double complex ***Atemp, 
		   options1D_t opts){
  int j, m;
  int det, dice, result, error;
  double *p, *sp;
  double complex **aux;

  int MAX = opts.max; 
  int detectors = opts.detectors;
  int *detector_pts = opts.detector_pts;

  const int lbound   = (opts.lattType == LINE_LATT) ? -MAX : 0;
  const int rbound   = (opts.lattType == LINE_LATT) ? MAX : MAX-1;
  const int aux_size = (opts.lattType == LINE_LATT) ? 2*MAX+1 : MAX;




  /* Detectors are described by a collection M_m of measurement 
   * operators which satisfy the completeness equation. The index m
   * indicates the result of the measurement. Before state |Psi> is 
   * measured the probability that result m occurs is given by 
   * p(m)= <Psi| Mdagger_m M_m |Psi>, and after the measurement the state 
   * collapses to |Psi'>= 1/sqrt(p(m)) M_m|Psi>. In qw1d we may define an 
   * arbitrary number of detectors specifying a list of points 
   * (m_1, ..., m_N) [that is, int *detector_pts]. The measurement 
   * operators are M_0 = I_2 tensor I_\infty - sum(M_i, i from 1 to N) and 
   * M_i = I_2 tensor |m_i><m_i|, for 1 <= i <= N. 
   *
   * To perform the measurement we (1) generate a random number r; and 
   * (2) determine an integer m, 0 <= m <= N, such that
   * sum(p(j), j from 0 to m-1) <= m < sum(p(j), j from 0 to m)
   */


  if( (!*A) || (!*Atemp))
    return -1;
  if(MAX<1)
    return -2;
  if(detectors<0)
    return -3;
  if(!detector_pts)
    return -4;

  p = allocReal1D(detectors+1);
  if(!p)
    return -5;
  error = cleanReal1D(p, detectors+1);
  if(error)
    return -6;

  sp = allocReal1D(detectors+1);
  if(!sp)
    return -5;
  error = cleanReal1D(sp, detectors+1);
  if(error)
    return -6;

  /* Here we evaluate the probability of each detection... */
  for(det=1; det<=detectors; det++ ){  
    m = (opts.lattType == LINE_LATT) ? 
      MAX + detector_pts[det] : detector_pts[det];

    for(j=0; j<2; j++)
      p[det] += (*A)[j][m]*conj((*A)[j][m]);
  }

  /* ...and the probability of measuring the complement */
  p[0] = 1.0;
  for(j=1; j<=detectors; j++)
    p[0] -= p[j];

  /* ...  */
  sp[0] = p[0];
  for(det=1; det<=detectors; det++)
    sp[det] = sp[det-1] + p[det];

  /* We get a new random number */
  dice = rand();

  /* Based on this number we identify the corresponding result */
  for(result=0; result<=detectors; result++)
    if(dice < sp[result]*RAND_MAX) break;

  for(m=lbound; m<=rbound; m++){
    int delta;
    const int auxm = (opts.lattType == LINE_LATT) ? MAX+m : m;
    
    if(result==0){
      delta = 1;
      for(det=1; det<=detectors; det++)
	delta *= (1-DELTA(m, detector_pts[det]));
    }
    else{
      delta = DELTA(m, detector_pts[result]);
    }
      
    for(j=0; j<2; j++)
      (*Atemp)[j][auxm] = delta*((*A)[j][auxm])/sqrt(p[result]);
  }

  free(p);
  free(sp);

  /* Now we quicky exchange matrices A and Atemp */
  aux = *A;
  *A = *Atemp;
  *Atemp = aux;

  error = cleanComplex2D(*Atemp, 2, aux_size);
  if(error)
    return -6;

  return result;
}




int randMeasure1D(double complex ***A, options1D_t opts){

  int mark, marked;
  int m, diceA, diceB;
  float sp;

  const int rbound = (opts.lattType == LINE_LATT) ? 
    2*opts.max : opts.max-1;


  /* We get a new random number */
  diceA = rand();

  /* Now we scan the lattice checking if any point is measured */
  mark = -1;
  marked=0;
  sp=0.0;
  for(m=0; m<=rbound; m++){
    int j;
    float prob;

    prob=0.0;
    for(j=0; j<2; j++)
      prob += (*A)[j][m]*conj((*A)[j][m]);

    diceB = rand();
    if(diceB < opts.dtProb*RAND_MAX){ /* if site should be measured */
      sp += prob;
      if(diceA < sp*RAND_MAX && !marked){ /* if state collapsed to this site */
	mark = m;
	marked = 1;
      }
      else{
	for(j=0; j<2; j++)
	  (*A)[j][m] = 0.0;
      }
    }
    
  }/* end-for m */

  for(m=0; m<=rbound; m++){
      int j;

      if(mark==-1){ /* measured the complement */
	for(j=0; j<2; j++)
	  (*A)[j][m] *= pow(1.0-sp, -0.5);

      }else if(mark==m){ /* measured this site */
	float prob;

	prob = 0.0;
	for(j=0; j<2; j++)
	  prob += (*A)[j][m]*conj((*A)[j][m]);

	for(j=0; j<2; j++)
	  (*A)[j][m] *= pow(prob, -0.5);
	

      }else{ /* measured some other site */
	
	for(j=0; j<2; j++)
	  (*A)[j][m] = 0.0;

      }

  }/* end-for m */

  return 0;
}


