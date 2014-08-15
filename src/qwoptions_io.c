/* QWalk (qwoptions_io.c) 
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
#include "qwoptions_io_read.h"
#include "qwmem_int.h"
#include "qwconsts.h"

options2D_t readOptionsFile2D(const char *filename){
  FILE *in;
  options2D_t options;
  char keyword[100];

  /* Default values */
  options.coinType = HADAMARD_COIN;
  options.stateType = HADAMARD_STATE;
  options.blType = NO_BROKENLINKS;
  options.lattType = NATURAL_LATT;

  options.blProbA = 0.0;
  options.blProbB = 0.0;
  options.dtProb  = 0.0;

  options.steps = 100;
  options.stepsMix = 0;
  options.numOfExperiments = 1;
  options.stepsAfterMeasure = 0;
  options.lattextra = 1;
  options.max = options.steps + options.lattextra;

  options.screen = 0;
  options.screen_pta[0] = 0;
  options.screen_pta[1] = 0;
  options.screen_ptb[0] = 0;
  options.screen_ptb[1] = 0;

  options.detectors = 0;
  options.detector_pts = NULL;

  options.calcMix = 0;
  options.checkState = 0;
  options.checkXSymmetry = 0;
  options.checkYSymmetry = 0;

  options.seed = time(0);


  in = fopen(filename,"rt");
  if(!in){
    options.error = 1;
    return options;
  }

  /* In the next loop we seach the BEGIN keyword. */
  do{
    fscanf(in,"%s",keyword); 
  }while(STRNEQ(keyword,"BEGIN") && !feof(in));

  if(feof(in)){
    options.error = 9;
    return options;
  }


  while(STRNEQ(keyword,"END")){
    /* In this loop we read a sequence of keywords and interpret them. 
     * We finish when we find an END keyword.
     */

    int error = 0;

    fscanf(in,"%s",keyword);

    if(STREQ(keyword,"COIN"))
      error = readOptions_coin2D(in, &options);
    else if(STREQ(keyword,"STATE"))
      error = readOptions_state2D(in, &options);
    else if(STREQ(keyword,"STEPS"))
      error = readOptions_steps2D(in, &options);
    else if(STREQ(keyword,"AFTERMEASURE"))
      error = readOptions_afterm2D(in, &options);
    else if(STREQ(keyword,"CHECK"))
      error = readOptions_check2D(in, &options);
    else if(STREQ(keyword,"CALCMIX")){
      printf("Warning: CALCMIX keyword is deprecated. Use MIXTIME instead.\n");
      error = readOptions_cmix2D(in, &options);
    }
    else if(STREQ(keyword,"MIXTIME"))
      error = readOptions_cmix2D(in, &options);
    else if(STREQ(keyword,"BLPROB"))
      error = readOptions_blprob2D(in, &options);
    else if(STREQ(keyword,"DTPROB"))
      error = readOptions_dtprob2D(in, &options);
    else if(STREQ(keyword,"EXPERIMENTS"))
      error = readOptions_exp2D(in, &options);
    else if(STREQ(keyword,"LATTSIZE"))
      error = readOptions_lsize2D(in, &options);
    else if(STREQ(keyword,"LATTEXTRA"))
      error = readOptions_lextra2D(in, &options);
    else if(STREQ(keyword,"LATTTYPE")){
      printf("Warning: LATTTYPE keyword is deprecated. Use LATTYPE instead.\n");
      error = readOptions_ltype2D(in, &options);
    }
    else if(STREQ(keyword,"LATTYPE"))
      error = readOptions_ltype2D(in, &options);
    else if(STREQ(keyword,"DETECTORS"))
      error = readOptions_detec2D(in, &options);
    else if(STREQ(keyword,"SEED"))
      error = readOptions_seed2D(in, &options);
    else if(STREQ(keyword,"SCREEN"))
      error = readOptions_screen2D(in, &options);
    else if(STREQ(keyword,"BLPERMANENT"))
      error = readOptions_blperm2D(in, &options);

    if(error)
      return options;

  }/* end-while */

  fclose(in);

  options.error = 0;
  return options;
}


options1D_t readOptionsFile1D(const char *filename){
  FILE *in;
  options1D_t options;
  char keyword[100];

  /* Default values */
  options.coinType = HADAMARD_COIN;
  options.stateType = HADAMARD_STATE;
  options.lattType = LINE_LATT;
  options.blProb = 0.0;
  options.dtProb  = 0.0;

  options.steps = 100;
  options.lattextra = 1;
  options.max = options.steps+options.lattextra;
  options.numOfExperiments = 1;
  options.stepsAfterMeasure = 0;

  options.checkState = 0;
  options.checkSymmetry = 0;

  options.seed = time(0); 

  options.detectors = 0;
  options.detector_pts = NULL;

  options.stepsMix = 0;
  options.calcMix = 0;

  in = fopen(filename,"rt");
  if(!in){
    options.error = 1;
    return options;
  }

  /* In the next loop we seach the BEGIN keyword. */
  do{
    fscanf(in, "%s", keyword); 
  }while(STRNEQ(keyword,"BEGIN") && !feof(in));

  if(feof(in)){
    options.error = 8;
    return options;
  }

  while(STRNEQ(keyword,"END")){
    /* In this loop we read a sequence of keywords and interpret them. 
     * We finish when we find an END keyword.
     */

    int error = 0;

    fscanf(in,"%s",keyword);

    if(STREQ(keyword,"COIN"))
      error = readOptions_coin1D(in, &options);
    else if(STREQ(keyword,"STATE"))
      error = readOptions_state1D(in, &options);
    else if(STREQ(keyword,"STEPS"))
      error = readOptions_steps1D(in, &options);
    else if(STREQ(keyword,"CHECK"))
      error = readOptions_check1D(in, &options);
    else if(STREQ(keyword,"BLPROB"))
      error = readOptions_blprob1D(in, &options);
    else if(STREQ(keyword,"DTPROB"))
      error = readOptions_dtprob1D(in, &options);
    else if(STREQ(keyword,"EXPERIMENTS"))
      error = readOptions_exp1D(in, &options);
    else if(STREQ(keyword,"SEED"))
      error = readOptions_seed1D(in, &options);
    else if(STREQ(keyword,"LATTSIZE"))
      error = readOptions_lsize1D(in, &options);
    else if(STREQ(keyword,"LATTEXTRA"))
      error = readOptions_lextra1D(in, &options);
    else if(STREQ(keyword,"LATTTYPE")){
      printf("Warning: LATTTYPE keyword is deprecated. Use LATTYPE instead.\n");
      error = readOptions_ltype1D(in, &options);
    }
    else if(STREQ(keyword,"LATTYPE"))
      error = readOptions_ltype1D(in, &options);
    else if(STREQ(keyword,"CALCMIX")){
      printf("Warning: CALCMIX keyword is deprecated. Use MIXTIME instead.\n");
      error = readOptions_cmix1D(in, &options);
    }
    else if(STREQ(keyword,"MIXTIME"))
      error = readOptions_cmix1D(in, &options);
    else if(STREQ(keyword,"DETECTORS"))
      error = readOptions_detec1D(in, &options);
    else if(STREQ(keyword,"AFTERMEASURE"))
      error = readOptions_afterm1D(in, &options);


    if(error) 
      return options;

  }/* end-while */

  fclose(in);

  options.error = 0;
  return options;
}
