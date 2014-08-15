/* QWalk (qw2d) This software simulates a 2D quantum walk
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
#include "qw2d_sub.h"

int main(int argc, char **argv){
  int MAX, error, experiment;
  int ****BrokenLinks1, ****BrokenLinks2;
  double complex ****A, ****Atemp, ****C;
  double **AverageProb = NULL, **StatProb = NULL;
  options2D_t options;
  filenames_t fnames;
  screen_t screen;

  printf("QWalk 2D, version 1.4 (qw2d).\n");
  printf("Copyright (C) 2008 Franklin Marquezino.\n");
  printf("This is free software; see the source code for copying conditions.\n");
  printf("There is ABSOLUTELY NO WARRANTY; not even for MERCHANTIBILITY or\n");
  printf("FITNESS FOR A PARTICULAR PURPOSE. See the source code for details.\n\n");

  printf("Email bug reports to <franklin@lncc.br>\n\n");

  if(argc<2){
    printf("Missing arguments.\nCorrect usage:\n  qw2d inputfile\n\n");
    exit(EXIT_FAILURE);
  }

  /**************************************************
   * Reading options and performing initializations *
   **************************************************/
  /* First we read the options file (usually with extension .in) */
  options = readOptionsFile2D(argv[1]);
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
    printf("Error: invalid number os experiments.\n");
    exit(EXIT_FAILURE);
  case 7:
    printf("Error: invalid lattice size.\n");
    exit(EXIT_FAILURE);
  case 8:
    printf("Error: invalid number of detectors.\n");
    exit(EXIT_FAILURE);
  case 9:
    printf("Error: error in the BEGIN-END structure.\n");
    exit(EXIT_FAILURE);
  case 10:
    printf("Error: invalid LATTEXTRA option\n");
    exit(EXIT_FAILURE);
  case 11:
    printf("Error: invalid probability (measurements of broken links)\n");
    exit(EXIT_FAILURE);
  case 12:
    printf("Error: invalid number of steps in mixing time calculation\n");
    exit(EXIT_FAILURE);
  case 13:
    printf("Error: invalid lattice type\n");
    exit(EXIT_FAILURE);
  }

  /* Here we define MAX as a short for options.max. It means 
   * that the lattice ranges from -MAX to MAX in X axys, and
   * from -MAX to MAX in Y axys.
   */
  MAX = options.max;

  /* Here we set the seed for the pseudorandom number generator */
  srand(options.seed);

  /* Here we define the names of output files based on the 
   * name of input file 
   */
  fnames = createFilenames2D(argv[1],options);
  printFilenames(stdout, argv[1],fnames);

  /* First we initialize all the links closed */
  error = initBrokenLink2D(&BrokenLinks1, &BrokenLinks2, MAX, options.lattType);
  if(error){
    printf("Error: could not initialize all links closed.\n");
    exit(EXIT_FAILURE);
  }

  /* In case of broken links, we read the broken link file in 
   * order to open the appropriate links. Note that the broken 
   * link information can be written in the same file used for 
   * the general options (usually a file with extension .in)     
   */
  if(options.blType == PERMANENT_BROKENLINKS){
    error = readBrokenLinkFile2D(argv[1], MAX, &BrokenLinks1, &BrokenLinks2,
				 options.lattType);
    if(error){
      printf("Error: could not read broken links from file.\n");
      exit(EXIT_FAILURE);
    }
  }

  /* If an observation screen was required in the input file then this
   * screen is initialized here. The screen is represented by a straight
   * line from (a0,a1) to (b0,b1) 
   */
  if(options.screen){
    error = initScreen(&screen, options.screen_pta[0], options.screen_pta[1], 
		       options.screen_ptb[0], options.screen_ptb[1]);
    if(error){
      printf("Error: could not initialize screen detector.\n");
      exit(EXIT_FAILURE);
    }
  }

  setCoin2D(&C, options.coinType, argv[1]);
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
    if((options.blProbA > 0.0) || (options.blProbB > 0.0) || (options.dtProb > 0.0)){
      printf("Warning: this version of qwalk should not be used to calculate or to plot\n");
      printf("  the approximate stationary distribution for decoherent two-dimensional\n");
      printf("  quantum walks (it is the uniform distribution, always).\n");
    }
    printf("This calculation may take a really long time...\n");
    setState2D(&A, options, argv[1]);
    if(!A){
      printf("Error: could not allocate initial state.\n");
      exit(EXIT_FAILURE);
    } 
    StatProb = getStationary2D(A, C, BrokenLinks1, BrokenLinks2, options);
    if(!StatProb){
      printf("Error: could not obtain (approximate) stationary distribution.\n");
      exit(EXIT_FAILURE);
    }
  }

  Atemp = (options.lattType == CYCLE_LATT) ?
    allocComplex4D(2, 2, MAX, MAX) : allocComplex4D(2, 2, 2*MAX+1, 2*MAX+1);
  if(!Atemp){
    printf("Error: could not allocate memory for temporary matrix.\n");
    exit(EXIT_FAILURE);
  }

  /*********************************************
   * Running the experiments                   *
   *********************************************/
  for(experiment=1; experiment <= options.numOfExperiments; experiment++){
    int steps, t;

    printf("Starting experiment %d of %d...\n", 
	   experiment, options.numOfExperiments);

    setState2D(&A, options, argv[1]);
    if(!A){
      printf("Error: could not allocate initial state.\n");
      exit(EXIT_FAILURE);
    } 

    /********************************
     * Performing a full simulation *
     ********************************/
    steps = options.steps;
    for(t=0; t<steps; t++){
      
      if((options.blProbA > 0.0) || (options.blProbB > 0.0)){
	/* If the user requested simulation of random broken links we enter here at 
	 * every step. We start by initializing all broken links closed.
	 */
	error = initBrokenLink2D(&BrokenLinks1, &BrokenLinks2, MAX, options.lattType);
	if(error){
	  printf("Error: could not initialize all links closed.\n");
	  exit(EXIT_FAILURE);
	}
	/* If besides random broken links we have also permanent broken links, then
	 * we need to read those permanent broken links from a file before breaking 
	 * the random ones.
	 */
	if(options.blType == PERMANENT_BROKENLINKS){
	  error = readBrokenLinkFile2D(argv[1], MAX, &BrokenLinks1, &BrokenLinks2,
				       options.lattType);
	  if(error){
	    printf("Error: could not read broken links from file.\n");
	    exit(EXIT_FAILURE);
	  }
	}
	/* Finally, we may break the random broken links.
	 */
	randomBrokenLink2D(&BrokenLinks1, &BrokenLinks2, options);
		
      }

      check2D(A, options, t);
      iterate2D(&A, &Atemp, C, BrokenLinks1, BrokenLinks2, options, t);

      if(options.detectors){
	/* If the user requested the simulation of a detector, we enter here
	 */
	int result;

	result = measureState2D(&A, &Atemp, options);
	if(result < 0){
	  printf("Error: could not measure state.");
	  exit(EXIT_FAILURE);
	}
	else if(result>0) /* non-trivial result */
	  steps = t+options.stepsAfterMeasure; /* we run some additional steps */
      }

      if(options.dtProb>0)
	randMeasure2D(&A, options);

      /* Now we calculate expectation, variance, standard deviation, etc, 
       * and save in a file.
       */
      doStatistics2D(A, StatProb, options, fnames, t+1, experiment);

      if(options.screen){
	/* If the user requested an observation screen, then we do it here */
	error = updateScreen(&screen, A, MAX);
	if(error){
	  printf("Error: could not update screen.\n");
	  exit(EXIT_FAILURE);
	}
      }

    }/* End-for t */

    /* Here we take the probability distribution obtained after of one experiment
     * and accumulate it to obtain, at the end, an average probability distribution
     * (which should not be confused with the average distribution used to define
     * the stationary distribution of a quantum Markov chain)
     */
    error = averageProbFromState2D(&AverageProb, A, options);
    if(error){
      printf("Error: could not update average probability matrix.\n");
      exit(EXIT_FAILURE);
    }

  }/* End-for experiments */


  /**********************
   * Freeing memory (I) *
   **********************/
  if(options.lattType == CYCLE_LATT){
    freeComplex4D(Atemp, 2, 2, MAX);
    freeInt4D(BrokenLinks1, 2, 2, MAX);
  }
  else{
    freeComplex4D(Atemp, 2, 2, 2*MAX+1);
    freeInt4D(BrokenLinks1, 2, 2, 2*MAX+1);
  }
  if(options.lattType == DIAG_LATT) 
    freeInt4D(BrokenLinks2, 2, 2, 2*MAX+1);
  else  
    freeInt4D(BrokenLinks2, 1, 1, 1);
  freeComplex4D(C, 2, 2, 2);

  /******************* 
   * Writing results *
   *******************/
  printf("Writing probabilities file...\n");
  error = writeData2D(fnames.dat_file, AverageProb, options);
  if(error){
    printf("Error: could not write output data file.\n");
    exit(EXIT_FAILURE);
  }

  printf("Writing wave-function file...\n");
  error = writeState2D(fnames.datwav_file, A, options);
  if(error){
    printf("Error: could not write wave-function file.\n");
    exit(EXIT_FAILURE);
  }


  if(options.calcMix){
    printf("Writing approximate stationary distribution...\n");
    error = writeData2D(fnames.datpb_file, StatProb, options);
    if(error){
      printf("Error: could not write stationary data.\n");
      exit(EXIT_FAILURE);
    }
  }


  if(options.screen){
    printf("Writing observation screen file...\n");
    error = writeScreen(fnames.datscr_file, screen);
    if(error){
      printf("Error: could not write screen detector file.\n");
      exit(EXIT_FAILURE);
    }
  }

  printf("Writing gnuplot script...\n");
  error = writeScript2D(fnames, options);
  if(error)
    printf("Warning: could not write script file for gnuplot.\n");

  /***********************
   * Freeing memory (II) *
   ***********************/
  if(options.lattType == CYCLE_LATT){
    freeComplex4D(A, 2, 2, MAX);
    if(options.calcMix)
      freeReal2D(StatProb, MAX);
    freeReal2D(AverageProb, MAX);
  }
  else{
    freeComplex4D(A, 2, 2, 2*MAX+1);
    if(options.calcMix)
      freeReal2D(StatProb, 2*MAX+1);
    freeReal2D(AverageProb, 2*MAX+1);
  }
  if(options.screen)
    free(screen.values);

  exit(EXIT_SUCCESS);
}
