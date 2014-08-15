/* QWalk (qwoptions_io_read.c) 
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
#include "qwoptions_io_read.h"
#include "qwmem_int.h"
#include "qwconsts.h"

int readOptions_coin2D(FILE *in, options2D_t *options){
  /* If a COIN keyword is found, then we expect one of the keywords:
   * FOURIER, HADAMARD, GROVER or CUSTOM. The last one requires the 
   * definition of the coin in a different file or in a separate part 
   * of the same input file, using the keywords BEGINCOIN and ENDCOIN
   */

  char keyword[100];

  fscanf(in,"%s",keyword); 
   
  if(STREQ(keyword,"CUSTOM"))
    options->coinType=CUSTOM_COIN;
  else if(STREQ(keyword,"FOURIER"))
    options->coinType=FOURIER_COIN;
  else if(STREQ(keyword,"GROVER"))
    options->coinType=GROVER_COIN;
  else if(STREQ(keyword,"HADAMARD"))
    options->coinType=HADAMARD_COIN;
  else{
    options->error = 2;
    return 2;
  }

  return 0;
}


int readOptions_state2D(FILE *in, options2D_t *options){
  /* If a STATE keyword is found, then we expect one of the keywords:
   * FOURIER, HADAMARD, GROVER or CUSTOM. The last one requires the 
   * definition of the initial state in a different file of in a
   * separate part of the same input file, using the keywords 
   * BEGINSTATE and ENDSTATE
   */
  
  char keyword[100];

  fscanf(in,"%s",keyword);    

  if(STREQ(keyword,"CUSTOM"))
    options->stateType=CUSTOM_STATE;
  else if(STREQ(keyword,"FOURIER"))
    options->stateType=FOURIER_STATE;
  else if(STREQ(keyword,"GROVER"))
    options->stateType=GROVER_STATE;
  else if(STREQ(keyword,"HADAMARD"))
    options->stateType=HADAMARD_STATE;
  else{
    options->error = 3;
    return 3;
  }
  
  return 0;
}


int readOptions_steps2D(FILE *in, options2D_t *options){
  /* If a STEPS keyword is found then we expect then a positive 
   * integer containing the number of steps that will be simulated.
   * NOTE: If LATTEXTRA and LATTSIZE keywords are used, STEPS must
   * come after LATTEXTRA and before LATTSIZE.
   */

  fscanf(in,"%d",&(options->steps)); 
  if(options->steps<1){
    options->error = 4;
    return 4;
  }
  
  /* The lattice ranges from -options.max to options.max in both axes.
   * The default size of the lattice is a little bigger than the number 
   * of steps simulated. Keyword LATTEXTRA sets options.lattextra,
   * which represents this extra space in the  lattice.
   */ 
  options->max = options->steps + options->lattextra;

  return 0;
}


int readOptions_afterm2D(FILE *in, options2D_t *options){
  /* If a AFTERMEASURE keyword is found then we expect a non-negative
   * integer describing the number of steps that will be simulated
   * after a measurement returns a non-trivial result.
   */
  
  fscanf(in, "%d", &(options->stepsAfterMeasure));

  if(options->stepsAfterMeasure<0){
    options->error = 4;
    return 4;
  }

  return 0;
}


int readOptions_check2D(FILE *in, options2D_t *options){
  /* If a CHECK keyword is found, then we expect one of the 
   * keywords: STATEPROB, XSYMMETRY, YSYMMETRY or DAGGER. 
   * If STATEPROB is found then the programm will check in each
   * step if the norm of wave-function is unitary. If XSYMMETRY
   * is found, the programm will check in each step the symmetry 
   * of the probability matrix around the y axis, i.e., it will 
   * check if the probabilities for the x>0 sites are equal to the 
   * probabilities for the x<0 sites. Analogously for the YSYMMETRY 
   * keyword.
   */

  char keyword[100];

  fscanf(in, "%s", keyword);

  if(STREQ(keyword,"STATEPROB"))
    options->checkState = 1;
  else if(STREQ(keyword,"XSYMMETRY"))
    options->checkXSymmetry = 1;
  else if(STREQ(keyword,"YSYMMETRY"))
    options->checkYSymmetry = 1;
  else{
    options->error = 5;
    return 5;
  }

  return 0;
}


int readOptions_cmix2D(FILE *in, options2D_t *options){
  /* If a CALCMIX keyword is found, then we expect a non-negative integer
   * describing how many steps will be used in the approximation of the
   * stationary distribution.
   */

  fscanf(in,"%d",&(options->stepsMix));
  if(options->stepsMix < 0){
    options->error = 12;
    return 12;
  }
  if(options->stepsMix == 0)
    options->calcMix = 0;
  else
    options->calcMix = 1;


  return 0;
}

int readOptions_blprob2D(FILE *in, options2D_t *options){
  /* If a BLPROB keyword is found, then we expect two non-negative 
   * real (double precision) numbers describing the probability of 
   * broken links in each direction.
   */
  
  fscanf(in,"%f",&(options->blProbA)); 
  if((options->blProbA < 0.0) || (options->blProbA > 1.0)){
    options->error = 11;
    return 11;
  }
  fscanf(in,"%f",&(options->blProbB)); 
  if((options->blProbB < 0.0) || (options->blProbB > 1.0)){
    options->error = 11;
    return 11;
  }

  return 0;
}



int readOptions_dtprob2D(FILE *in, options2D_t *options){
  /* If a DTPROB keyword is found, then we expect a non-negative 
   * real (double precision) number describing the probability of 
   * measurement in each lattice site
   */
  
  fscanf(in,"%f",&(options->dtProb)); 
  if((options->dtProb < 0.0) || (options->dtProb > 1.0)){
    options->error = 11;
    return 11;
  }

  return 0;
}



int readOptions_exp2D(FILE *in, options2D_t *options){
  /* If a EXPERIMENTS keyword is found then we expect a positive
   * integer describing the number of experiments that will be
   * carried out.
   */

  fscanf(in, "%d", &(options->numOfExperiments));
  if(options->numOfExperiments<1){
    options->error = 6;
    return 6;
  }

  return 0;
}


int readOptions_lsize2D(FILE *in, options2D_t *options){
  /* If a LATTSIZE keyword is found then we expect a positive
   * integer describing the size of the lattice. We consider
   * that the lattice ranges from -options.max to options.max
   * in both axes. An extra space, given by options.lattextra,
   * is added to the size of the lattice in order to avoid
   * "bound errors" when accessing it during the execution of
   * the software.
   * NOTE: If LATTEXTRA and STEPS keywords are used, LATTSIZE
   * must come last.
   */

  fscanf(in, "%d", &(options->max));
  if(options->max<1){
    options->error=7;
    return 7;
  }
  if(options->lattType != CYCLE_LATT)
    options->max += options->lattextra;

  return 0;
}


int readOptions_lextra2D(FILE *in, options2D_t *options){
  /* If a LATTEXTRA keyword is found then we expect a non-negative
   * integer describing the extra space reserved for the lattice.
   * This options is very important, for example, when the initial
   * condition is not entirely localized in (0,0).
   * NOTE: If STEPS and LATTSIZE keywords are used, LATTEXTRA must 
   * come first.
   */

  fscanf(in, "%d", &(options->lattextra));
  if(options->lattextra<0){
    options->error=10;
    return 10;
  }

  return 0;
}


int readOptions_ltype2D(FILE *in, options2D_t *options){
  /* If a LATTTYPE keyword is found, then we expect one of the keywords:
   * NATURAL or DIAGONAL. If DIAGONAL keyword is found, the simulador
   * will use the evolution equation which makes the physical lattice
   * diagonal in relation to the mathematicas lattice. If NATURAL keyword
   * is found, the simulator will use the evolution equation which makes
   * the physical lattice coincides with the mathematical lattice.
   */

  char keyword[100];

  fscanf(in, "%s", keyword);

  if(STREQ(keyword,"NATURAL"))
    options->lattType = NATURAL_LATT;
  else if(STREQ(keyword,"DIAGONAL"))
    options->lattType = DIAG_LATT;
  else if(STREQ(keyword,"CYCLE")){
    options->lattType = CYCLE_LATT;

    options->checkXSymmetry = 0;
    options->checkYSymmetry = 0;

    /* just in case the LATTSIZE statement has already been used */
    options->max -= options->lattextra;

    /* in the LATTSIZE statement is used after this point, there
     * will be no adjustment via options.lattextra
     */
    options->lattextra = 0;
  }
  else{
    options->error = 13;
    return 13;
  }


  return 0;
}



int readOptions_detec2D(FILE *in, options2D_t *options){
  /* If a DETECTORS keyword is found then we expect a positive integer,
   * describing the number of detectors used in the simulation. After
   * that, for each detector we expect two integers describing the
   * positions of these detectors.
   */
  
  int i;
  
  fscanf(in, "%d", &(options->detectors));
  if(options->detectors<1){
    options->error=8;
    return 8;
  }
  
  options->detector_pts = allocInt2D(options->detectors+1, 2);
  if(!options->detector_pts){
    options->error=8;
    return 8;
  }

  for(i=1; i<=options->detectors; i++){
    fscanf(in, "%d", &(options->detector_pts[i][0]));
    fscanf(in, "%d", &(options->detector_pts[i][1]));
  }

  return 0;
}


int readOptions_seed2D(FILE *in, options2D_t *options){
  /* If a SEED keyword is found then we expect an integer describing
   * a seed for the pseudorandom number generator.
   */
    
  fscanf(in, "%d", &(options->seed));
  options->seed = abs(options->seed);

  return 0;
}


int readOptions_screen2D(FILE *in, options2D_t *options){
  /* If a SCREEN keyword is found, then we expect four integers
   * describing its position. The first two integers, say xa and
   * ya, describe the first point, (xa,ya). The next two integers,
   * say xb and yb, describe the second point, (xb,yb). With these
   * two points it is possible to define the observation screen.
   * The screen must be in horizontal, vertical ou in 45 degrees.
   */

  options->screen = 1;
  fscanf(in,"%d",&(options->screen_pta[0])); 
  fscanf(in,"%d",&(options->screen_pta[1])); 
  fscanf(in,"%d",&(options->screen_ptb[0])); 
  fscanf(in,"%d",&(options->screen_ptb[1])); 

  return 0;
}


int readOptions_blperm2D(FILE *in, options2D_t *options){
  /* If a BLPERMANENT keyword is found we set the field blType 
   * in options2D_t structure as PERMANENT_BROKENLINKS, so that 
   * the program reads the broken links file later.
   */
  
  options->blType = PERMANENT_BROKENLINKS;
  return 0;
}


int readOptions_coin1D(FILE *in, options1D_t *options){
  /* If a COIN keyword is found, we expect then one of the keywords:
   * HADAMARD or CUSTOM. The last one requires the definition
   * of the coin in a different file or in a separate part of the 
   * same input file, using the keywords BEGINCOIN and ENDCOIN
   */

  char keyword[100];

  fscanf(in,"%s",keyword);    
  if(STREQ(keyword,"CUSTOM"))
    options->coinType=CUSTOM_COIN;
  else if(STREQ(keyword,"HADAMARD"))
    options->coinType=HADAMARD_COIN;
  else{
    options->error = 2;
    return 2;
  }

  return 0;
}


int readOptions_state1D(FILE *in, options1D_t *options){
  /* If a STATE keyword is found, then we expect one of the 
   * keywords: HADAMARD or CUSTOM. The last one requires the 
   * definition of the initial state in a different file or in a
   * separate part of the same input file, using the keywords 
   * BEGINSTATE and ENDSTATE
   */
  
  char keyword[100];
  
  fscanf(in,"%s",keyword);    
  if(STREQ(keyword,"CUSTOM"))
    options->stateType=CUSTOM_STATE;
  else if(STREQ(keyword,"HADAMARD"))
    options->stateType=HADAMARD_STATE;
  else{
    options->error = 3;
    return 3;
  }

  return 0;
}


int readOptions_steps1D(FILE *in, options1D_t *options){
  /* If a STEPS keyword is found, then we expect a positive 
   * integer containing the number of steps that will be 
   * simulated.
   */
  
  fscanf(in,"%d",&(options->steps)); 
  if(options->steps<1){
    options->error = 4;
    return 4;
  }

  /* The lattice ranges from -options.max to options.max. The 
   * default size of the lattice is a little bigger than the number 
   * of steps simulated. Keyword LATTEXTRA sets options.lattextra,
   * which represents this extra space in the  lattice.
   */ 
  options->max = options->steps + options->lattextra;

  return 0;
}


int readOptions_check1D(FILE *in, options1D_t *options){
  /* If a CHECK keyword is found, then we expect one of the 
   * keywords: STATEPROB, SYMMETRY or DAGGER. 
   * If STATEPROB is found then the programm will check in each
   * step if the norm of wave-function is unitary. If SYMMETRY
   * is found, the programm will check in each step the symmetry 
   * of the probability matrix around the axis, i.e., it will 
   * check if the probabilities for the x>0 sites are equal to the 
   * probabilities for the x<0 sites. 
   */

  char keyword[100];

  fscanf(in, "%s", keyword);
  if(STREQ(keyword,"STATEPROB"))
    options->checkState = 1;
  else if(STREQ(keyword,"SYMMETRY"))
    options->checkSymmetry = 1;
  else
    options->error = 5;

  return 0;
}


int readOptions_blprob1D(FILE *in, options1D_t *options){
  /* If a BLPROB keyword is found, then we expect then a positive 
   * float containing the probability of broken links in the simulation,
   * i.e., in each simulation step each link has this probability of 
   * being open.
   */
  
  fscanf(in,"%f",&(options->blProb)); 
  if((options->blProb < 0.0) || (options->blProb > 1.0)){
    options->error = 6;
    return 6;
  }

  return 0;
}



int readOptions_dtprob1D(FILE *in, options1D_t *options){
  /* If a DTPROB keyword is found, then we expect a non-negative 
   * real (double precision) number describing the probability of 
   * measurement in each lattice site
   */
  
  fscanf(in,"%f",&(options->dtProb)); 
  if((options->dtProb < 0.0) || (options->dtProb > 1.0)){
    options->error = 13;
    return 13;
  }

  return 0;
}



int readOptions_exp1D(FILE *in, options1D_t *options){
  /* If EXPERIMENTS keyword is found, then we expect a positive 
   * integer containing the experiments that will be carried out 
   * (useful for simulations involving random broken links).
   */

  fscanf(in,"%d",&(options->numOfExperiments)); 
  if(options->numOfExperiments < 1){
    options->error = 7;
    return 7;
  }

  return 0;
}


int readOptions_seed1D(FILE *in, options1D_t *options){
  /* If a SEED keyword is found then we expect an integer describing
   * a seed for the pseudorandom number generator.
   */
    
  fscanf(in, "%d", &(options->seed));
  options->seed = abs(options->seed);

  return 0;
}


int readOptions_lsize1D(FILE *in, options1D_t *options){
  /* If a LATTSIZE keyword is found then we expect a positive
   * integer describing the size of the lattice. We consider
   * that the lattice ranges from -options.max to options.max
   * An extra space, given by options.lattextra, is added to 
   * the size of the lattice in order to avoid "bound errors" 
   * when accessing it during the execution of the software.
   * NOTE: If LATTEXTRA and STEPS keywords are used, LATTSIZE
   * must come last.
   */

  fscanf(in, "%d", &(options->max));
  if(options->max < 1){
    options->error = 9;
    return 9;
  }
  if(options->lattType == LINE_LATT)
    options->max += options->lattextra;

  return 0;
}


int readOptions_lextra1D(FILE *in, options1D_t *options){
  /* If a LATTEXTRA keyword is found then we expect a non-negative
   * integer describing the extra space reserved for the lattice.
   * This options is very important, for example, when the initial
   * condition is not entirely localized in x=0.
   * NOTE: If STEPS and LATTSIZE keywords are used, LATTEXTRA must 
   * come first.
   */

  fscanf(in, "%d", &(options->lattextra));
  if(options->lattextra < 0){
    options->error=10;
    return 10;
  }

  return 0;
}

int readOptions_ltype1D(FILE *in, options1D_t *options){
  /* If a LATTTYPE keyword is found, then we expect one of the keywords:
   * LINE, CYCLE or SEGMENT. If LINE keyword is found, the simulador
   * will use the evolution equation for an infinite one-dimensional
   * lattice. If SEGMENT keyword is found, the simulator will set 
   * permanent broken links with the previous evolution equation in
   * order to define a finite one-dimensional lattice with reflecting
   * boundaries. If CYCLE keyword is found, the simulator will use the
   * evolution equation for the cycle. 
   */

  char keyword[100];

  fscanf(in, "%s", keyword);

  if(STREQ(keyword,"LINE"))
    options->lattType = LINE_LATT;
  else if(STREQ(keyword,"CYCLE")){
    options->lattType = CYCLE_LATT;
    options->checkSymmetry = 0;

    /* just in case the LATTSIZE statement has already been used */
    options->max -= options->lattextra;

    /* in the LATTSIZE statement is used after this point, there
     * will be no adjustment via options.lattextra
     */
    options->lattextra = 0;
  }
  else if(STREQ(keyword,"SEGMENT")){
    options->lattType = SEGMENT_LATT;
    options->checkSymmetry = 0;

    /* just in case the LATTSIZE statement has aldeary been used */
    options->max -= options->lattextra;

    /* in the LATTSIZE statement is used after this point, there
     * will be no adjustment via options.lattextra
     */
    options->lattextra = 0;
  }
  else{
    options->error = 11;
    return 11;
  }

  return 0;
}


int readOptions_cmix1D(FILE *in, options1D_t *options){
  /* If a CALCMIX keyword is found, then we expect a non-negative integer
   * describing how many steps will be used in the approximation of the
   * stationary distribution.
   */

  fscanf(in,"%d",&(options->stepsMix));
  if(options->stepsMix < 0){
    options->error = 12;
    return 12;
  }
  if(options->stepsMix == 0)
    options->calcMix = 0;
  else
    options->calcMix = 1;

  return 0;
}


int readOptions_detec1D(FILE *in, options1D_t *options){
  /* If a DETECTORS keyword is found then we expect a positive integer,
   * describing the number of detectors used in the simulation. After
   * that, for each detector we expect one integers describing the
   * positions of these detectors.
   */
  
  int i;
  
  fscanf(in, "%d", &(options->detectors));
  if(options->detectors<1){
    options->error=8;
    return 14;
  }
  
  options->detector_pts = (int *)malloc((options->detectors+1)*sizeof(int));
  if(!options->detector_pts){
    options->error=14;
    return 14;
  }

  for(i=1; i<=options->detectors; i++)
    fscanf(in, "%d", &(options->detector_pts[i]));

  return 0;
}

int readOptions_afterm1D(FILE *in, options1D_t *options){
  /* If a AFTERMEASURE keyword is found then we expect a non-negative
   * integer describing the number of steps that will be simulated
   * after a measurement returns a non-trivial result.
   */
  
  fscanf(in, "%d", &(options->stepsAfterMeasure));

  if(options->stepsAfterMeasure<0){
    options->error = 4;
    return 4;
  }

  return 0;
}


