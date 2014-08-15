/* QWalk (qwprob_io.c) 
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
#include "qwoptions_io.h"
#include "qwprob_io.h"
#include "qwconsts.h"


int writeData1D(const char *filename, double *probs, options1D_t options){
  FILE *out;
  int m;
  time_t lt;
  const int MAX = options.max;

  out=fopen(filename,"wt");
  if(!out)
    return 1;
  if(!probs)
    return 2;

  fprintf(out, "#Output generated by Quantum Walk Simulator (1D)\n");

  lt = time(NULL);
  fprintf(out, "# %s\n\n",ctime(&lt));

  fprintf(out,"# Simulation options (1D):\n");

  switch(options.coinType){
  case HADAMARD_COIN: 
    fprintf(out,"#  Hadamard coin.\n");
    break;
  case CUSTOM_COIN:
    fprintf(out,"#  Customized coin (given by input file).\n");
  }

  switch(options.stateType){
  case HADAMARD_STATE: 
    fprintf(out,"#  Using a state that gives maximum spread for Hadamard coin.\n");
    break;
  case CUSTOM_STATE:
    fprintf(out,"#  Customized state (given by input file).\n");
  }

  fprintf(out,"#  Lattice type: ");
  switch(options.lattType){
  case LINE_LATT: 
    fprintf(out,"LINE.\n");
    break;
  case CYCLE_LATT:
    fprintf(out,"CYCLE.\n");
    break;
  case SEGMENT_LATT:
    fprintf(out,"SEGMENT.\n");
    break;
  }

  fprintf(out,"#  Probability of broken links: %f\n", options.blProb);
  fprintf(out,"#  Probability of measurement: %f\n", options.dtProb);
  fprintf(out,"#  Number of experiments: %d\n", options.numOfExperiments);
  fprintf(out,"#  Number of steps: %d\n", options.steps);
  if(options.calcMix)
    fprintf(out,"#  Steps to approximate stationary distribution: %d\n", 
	    options.stepsMix);

  if(options.lattType == LINE_LATT)
    fprintf(out,"#  Lattice size: -%d..%d in X axis.\n",MAX,MAX);
  else
    fprintf(out,"#  Lattice size: 0..%d in X axis (%d sites).\n",MAX-1, MAX);

  if(options.checkState)
    fprintf(out,"#  Sum of probabilities checked in each iteration.\n");
  if(options.checkSymmetry)
    fprintf(out,"#  Symmetry of probability array checked in each iteration.\n");
  fprintf(out,"#  Random seed: %d\n\n", options.seed);

  if(options.lattType == LINE_LATT){
    for(m=0; m<=2*MAX; m++){
      if(probs[m] > 0.0)
	fprintf(out,"%d\t%e\n",m-MAX,probs[m]);
    }
  }
  else{
    for(m=0; m<MAX; m++){
      if(probs[m] > 0.0)
	fprintf(out,"%d\t%e\n",m,probs[m]);
    }
  }
  fclose(out);

  return 0;
}



int writeData2D(const char *filename, double **probs, options2D_t options){
  FILE *out;
  int m, n;
  time_t lt;
  double probSum = 0.0;
  const int MAX = options.max;

  out=fopen(filename,"wt");
  if(!out)
    return 1;
  if(!*probs)
    return 2;

  fprintf(out, "#Output generated by Quantum Walk Simulator (2D)\n");

  lt = time(NULL);
  fprintf(out, "# %s\n\n",ctime(&lt));

  fprintf(out,"# Simulation options (2D):\n");

  switch(options.coinType){
  case HADAMARD_COIN: 
    fprintf(out,"#  Hadamard coin.\n");
    break;
  case GROVER_COIN:
    fprintf(out,"#  Grover coin.\n");
    break;
  case FOURIER_COIN:
    fprintf(out,"#  Fourier coin.\n");
    break;
  case CUSTOM_COIN:
    fprintf(out,"#  Customized coin (given by input file).\n");
  }

  switch(options.stateType){
  case HADAMARD_STATE: 
    fprintf(out,"#  Using the state that gives maximum spread for Hadamard coin.\n");
    break;
  case GROVER_STATE:
    fprintf(out,"#  Using the state that gives maximum spread for Grover coin.\n");
    break;
  case FOURIER_STATE:
    fprintf(out,"#  Using the state that gives maximum spread for Fourier coin.\n");
    break;
  case CUSTOM_STATE:
    fprintf(out,"#  Customized state (given by input file).\n");
  }

  fprintf(out,"#  Lattice type: ");
  switch(options.lattType){
  case DIAG_LATT: 
    fprintf(out,"DIAGONAL.\n");
    break;
  case NATURAL_LATT:
    fprintf(out,"NATURAL (default).\n");
    break;
  case CYCLE_LATT:
    fprintf(out,"CYCLE.\n");
    break;
  }


  if(options.blType==PERMANENT_BROKENLINKS)
    fprintf(out,"#  Using permanent broken links (topology given by input file).\n");
  if(options.checkState)
    fprintf(out,"#  Sum of probabilities checked in each iteration.\n");
  if(options.checkXSymmetry)
    fprintf(out,"#  x-Symmetry of probability matrix checked in each iteration.\n");
  if(options.checkYSymmetry)
    fprintf(out,"#  y-Symmetry of probability matrix checked in each iteration.\n");

  fprintf(out,"#  Probability of broken links: %f, %f\n", options.blProbA, options.blProbB);
  fprintf(out,"#  Probability of measurement: %f\n", options.dtProb);
  fprintf(out,"#  Number of experiments: %d\n", options.numOfExperiments);
  fprintf(out,"#  Number of steps: %d\n", options.steps);
  if(options.calcMix)
    fprintf(out,"#  Steps to approximate stationary distribution: %d\n", 
	    options.stepsMix);

  if(options.lattType == CYCLE_LATT)
    fprintf(out,"#  Lattice size: 0..%d in X axis and 0..%d in Y axis.\n",MAX,MAX);
  else
    fprintf(out,"#  Lattice size: -%d..%d in X axis and -%d..%d in Y axis.\n",
	    MAX,MAX,MAX,MAX);
  fprintf(out,"#  Random seed: %d\n\n", options.seed);

  if(options.lattType == CYCLE_LATT){
    for(m=0; m<MAX; m++){
      char needsBlankLine;

      if(!probs[m])
	return 2;

      needsBlankLine = 0;
      for(n=0; n<MAX; n++){
	probSum += probs[m][n];
	if(probs[m][n]>0.0){
	  fprintf(out,"%d\t%d\t%e\n",m,n,probs[m][n]);
	  needsBlankLine = 1;
	}
      }

      if(needsBlankLine)
	fprintf(out,"\n");
    }
  }else{
    for(m=0; m<=2*MAX; m++){
      char needsBlankLine;

      if(!probs[m])
	return 2;

      needsBlankLine = 0;
      for(n=0; n<=2*MAX; n++){
	probSum += probs[m][n];
	if(probs[m][n]>0.0){
	  fprintf(out,"%d\t%d\t%e\n",m-MAX,n-MAX,probs[m][n]);
	  needsBlankLine = 1;
	}
      }

      if(needsBlankLine)
	fprintf(out,"\n");
    }
  }
  fclose(out);

  if( fabs(probSum-1.0) > WALK_TOL )
    return 1;

  return 0;
}








