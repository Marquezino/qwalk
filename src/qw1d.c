/* QWalk (qw1d) This software simulates a 1D quantum walk
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
#include "qwmem_complex.h"
#include "qwmem_real.h"
#include "qwcoin.h"
#include "qwstate.h"
#include "qwprob.h"
#include "qwstatistics.h"
#include "qwlinks.h"
#include "qwscreen.h"
#include "qwmeasure.h"
#include "qwcoin_io.h"
#include "qwstate_io.h"
#include "qwprob_io.h"
#include "qwstatistics_io.h"
#include "qwoptions_io.h"
#include "qwextra_io.h"
#include "qwconsts.h"
#include "qw1d_sub.h"

int main(int argc, char **argv){
  int experiment, MAX;
  int **BrokenLinks, error;
  double complex **A, **Atemp;
  double complex **C;
  double *AverageProb, *StatProb=NULL;
  options1D_t options;
  filenames_t fnames;

  printf("QWalk 1D, version 1.4 (qw1d).\n");
  printf("Copyright (C) 2008 Franklin Marquezino.\n");
  printf("This is free software; see the source code for copying conditions.\n");
  printf("There is ABSOLUTELY NO WARRANTY; not even for MERCHANTIBILITY or\n");
  printf("FITNESS FOR A PARTICULAR PURPOSE. See the source file for details.\n\n");

  printf("Email bug reports to <franklin@lncc.br>\n\n");

  if(argc<2){
    printf("Missing arguments.\nCorrect usage:\n  qw1d inputfile\n\n");
    exit(EXIT_FAILURE);
  }

  /**************************************************
   * Reading options and performing initializations *
   **************************************************/

  options = readOptionsFile1D(argv[1]);
  switch(options.error){
  case 0:
    break;
  case 1:
    printf("Error: could not open input file.\n");
    exit(EXIT_FAILURE);
  case 2:
    printf("Error: invalid coin.\n");
    exit(EXIT_FAILURE);
  case 3:
    printf("Error: invalid state.\n");
    exit(EXIT_FAILURE);
  case 4:
    printf("Error: invalid number of steps.\n");
    exit(EXIT_FAILURE);
  case 5:
    printf("Error: could not understand CHECK option.\n");
    exit(EXIT_FAILURE);
  case 6:
    printf("Error: invalid probability of broken links.\n");
    exit(EXIT_FAILURE);
  case 7:
    printf("Error: invalid number of experiments.\n");
    exit(EXIT_FAILURE);   
  case 8: 
    printf("Error: invalid BEGIN-END structure.\n");
    exit(EXIT_FAILURE);
  case 9: 
    printf("Error: invalid lattice size.\n");
    exit(EXIT_FAILURE);
  case 10:
    printf("Error: invalid LATTEXTRA.\n");
    exit(EXIT_FAILURE);
  case 11:
    printf("Error: invalid lattice type.\n");
    exit(EXIT_FAILURE);
  case 12:
    printf("Error: invalid number of steps in mixing time calculation.\n");
    exit(EXIT_FAILURE);
  case 13:
    printf("Error: invalid probability (measurements of broken links).\n");
    exit(EXIT_FAILURE);
  case 14:
    printf("Error: invalid number of detectors.\n");
    exit(EXIT_FAILURE);
  }

  /* We define MAX as a short for options.max. It means that 
   * the lattice ranges from -MAX to MAX.
   */
  MAX = options.max;

  /* Here we set the seed for the pseudorandom number generator */
  srand(options.seed);

  /* Here we define the names of output files based on the 
   * name of input file 
   */
  fnames = createFilenames1D(argv[1], options);
  printFilenames(stdout, argv[1], fnames);

  /* Here we initialize all the links closed */
  error = initBrokenLink1D(&BrokenLinks, MAX, options.lattType);
  if(error){
    printf("Error: could not initialize all links closed.\n");
    exit(EXIT_FAILURE);
  }

  setCoin1D(&C, options.coinType, argv[1]);
  if(!C){
    printf("Error: could not allocate matrix for coin");
    exit(EXIT_FAILURE);
  }


  /* If the user requested the calculation of mixing time, we need to
   * calculate the approximate stationary distribution here.
   */
  if(options.calcMix){
    printf("Calculating (approximate) stationary distribution with %d steps.\n",
	   options.stepsMix);
    if(options.blProb > 0.0 || options.dtProb > 0.0){
      printf("Warning: this version of qwalk should not be used to calculate or to plot\n");
      printf("  the approximate stationary distribution for decoherent one-dimensional\n");
      printf("  quantum walks (it is the uniform distribution, always).\n");
    }
    printf("This calculation may take a really long time...\n");
    setState1D(&A, options, argv[1]);
    if(!A){
      printf("Error: could not allocate initial state.\n");
      exit(EXIT_FAILURE);
    } 
    StatProb = getStationary1D(A, C, BrokenLinks, options);
    if(!StatProb){
      printf("Error: could not obtain (approximate) stationary distribution.\n");
      exit(EXIT_FAILURE);
    }
  }

  Atemp = (options.lattType == LINE_LATT) ?
    allocComplex2D(2, 2*MAX+1) : allocComplex2D(2, MAX);
  if(!Atemp){
    printf("Error: could not allocate memory for temporary matrix.\n");
    exit(EXIT_FAILURE);
  }

  /***************************
   * Running the experiments *
   ***************************/
  for(experiment=1; experiment <= options.numOfExperiments; experiment++){
    int steps,t;

    printf("Starting experiment %d of %d...\n", 
	   experiment, options.numOfExperiments);

    setState1D(&A, options, argv[1]);
    if(!A){
      printf("Error: could not allocate initial state.\n");
      exit(EXIT_FAILURE);
    }

    /********************************
     * Performing a full simulation *
     ********************************/
    steps = options.steps;
    for(t=0; t<steps; t++){
      initBrokenLink1D(&BrokenLinks, MAX, options.lattType);
      randomBrokenLink1D(&BrokenLinks, options);

      check1D(A, options, t); 
      iterate1D(&A, &Atemp, C, BrokenLinks, options, t);

      if(options.detectors){
	/* If the user requested the simulation of a detector, we enter here
	 */
	int result;

	result = measureState1D(&A, &Atemp, options);

	if(result < 0){
	  printf("Error: could not measure state.");
	  exit(EXIT_FAILURE);
	}
	else if(result>0) /* non-trivial result */
	  steps = t+options.stepsAfterMeasure; /* we run some additional steps */
      }

      if(options.dtProb>0)
	randMeasure1D(&A, options);

      doStatistics1D(A, StatProb, options, fnames, t+1, experiment);
    }/* End-for t */

    error=averageProbFromState1D(&AverageProb,A,options);
    if(error){
      printf("Error: could not update average probability matrix.\n");
      exit(EXIT_FAILURE);
    }

  } /* End-for experiment */



  /**********************
   * Freeing memory (I) *
   **********************/
  freeComplex2D(Atemp, 2);
  freeComplex2D(C, 2);
  freeInt2D(BrokenLinks, 2);

  /******************* 
   * Writing results *
   *******************/

  printf("Writing probabilities file...\n");
  error = writeData1D(fnames.dat_file, AverageProb, options);
  if(error){
    printf("Error: could not write output data file.\n");
    exit(EXIT_FAILURE);
  }

  printf("Writing wave-function file...\n");
  error = writeState1D(fnames.datwav_file, A, options);
  if(error){
    printf("Error: could not write wave-function file.\n");
    exit(EXIT_FAILURE);
  }

  if(options.calcMix){
    printf("Writing approximate stationary distribution...\n");
    error = writeData1D(fnames.datpb_file, StatProb, options);
    if(error){
      printf("Error: could not write stationary data.\n");
      exit(EXIT_FAILURE);
    }
  }

  printf("Writing gnuplot script...\n");
  error = writeScript1D(fnames, options);
  if(error)
    printf("Warning: could not write script file for gnuplot.\n");


  /***********************
   * Freeing memory (II) *
   ***********************/
  freeComplex2D(A, 2);
  if(options.calcMix)
    free(StatProb);
  free(AverageProb);

  printf("\nSimulation finished.\n");
  printf("Please, report bug reports to franklin@lncc.br.\n\n");

  exit(EXIT_SUCCESS);
}

