/* QWalk (qwextra_io.c) 
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
#include "qwextra_io.h"
#include "qwoptions_io.h"
#include "qwconsts.h"

filenames_t createFilenames1D(const char *input_filename, options1D_t options){
  filenames_t fnames;

  /* probabilities */
  changeFileEnding(&(fnames.dat_file), input_filename, ".dat");

  /* wave */
  changeFileEnding(&(fnames.datwav_file), input_filename, "-wave.dat");

  /* Stationary distribution, $\bar{P}$ */
  if(options.calcMix){
    changeFileEnding(&(fnames.datpb_file), input_filename, "-pb.dat");
    changeFileEnding(&(fnames.epspb_file), input_filename, "-pb.eps");
  }
  else{
    fnames.datpb_file = "";
    fnames.epspb_file = "";
  }

  /* statistics */
  changeFileEnding(&(fnames.sta_file), input_filename, ".sta");

  /* plot */
  changeFileEnding(&(fnames.eps2d_file), input_filename, ".eps");

  /* gnuplot script */
  changeFileEnding(&(fnames.gpt_file), input_filename, ".plt");

  /* some files are only used in 2D simulations */
  fnames.datscr_file = "";
  fnames.epsscr_file = "";
  fnames.eps3d_file = "";

  return fnames;
}


filenames_t createFilenames2D(const char *input_filename, options2D_t options){
  filenames_t fnames;

  /* probabilities */
  changeFileEnding(&(fnames.dat_file), input_filename, ".dat");

  /* wave */
  changeFileEnding(&(fnames.datwav_file), input_filename, "-wave.dat");

  /* Stationary distribution, $\bar{P}$ */
  if(options.calcMix){
    changeFileEnding(&(fnames.datpb_file), input_filename, "-pb.dat");
    changeFileEnding(&(fnames.epspb_file), input_filename, "-pb.eps");
  }
  else{
    fnames.datpb_file = "";
    fnames.epspb_file = "";
  }


  /* observation screen */
  if(options.screen){
    changeFileEnding(&(fnames.datscr_file), input_filename, "-screen.dat");
    changeFileEnding(&(fnames.epsscr_file), input_filename, "-screen.eps");
  }
  else{
    fnames.datscr_file = "";
    fnames.epsscr_file = "";
  }

  /* statistics */
  changeFileEnding(&(fnames.sta_file), input_filename, ".sta");

  /* 3D plot */
  changeFileEnding(&(fnames.eps3d_file), input_filename, "-3d.eps");

  /* contour plot */
  changeFileEnding(&(fnames.eps2d_file), input_filename, "-2d.eps");

  /* gnuplot script */
  changeFileEnding(&(fnames.gpt_file), input_filename, ".plt");

  return fnames;
}


int printFilenames(FILE *out, char *input_filename, filenames_t fnames){

  if(!out)
    return 1;
  if(!input_filename)
    return 2;

  fprintf(out,"Input filename ..........: %s\n\n", input_filename);

  fprintf(out,"Files generated by this simulation:\n");
  fprintf(out,"  Probabilities .........: %s\n", fnames.dat_file);
  fprintf(out,"  Wave function .........: %s\n", fnames.datwav_file);
  fprintf(out,"  Stationary distribution: %s\n", fnames.datpb_file);
  fprintf(out,"  Screen data ...........: %s\n", fnames.datscr_file);
  fprintf(out,"  Statistics ............: %s\n", fnames.sta_file);
  fprintf(out,"  Gnuplot script ........: %s\n\n", fnames.gpt_file);

  fprintf(out,"Files generated after running gnuplot:\n");
  fprintf(out,"  3D plot postscript.....: %s\n", fnames.eps3d_file);
  fprintf(out,"  2D plot postscript.....: %s\n", fnames.eps2d_file);
  fprintf(out,"  Screen postscript......: %s\n", fnames.epsscr_file);
  fprintf(out,"  Stationary postscript..: %s\n\n", fnames.epspb_file);


  return 0;
}


int readBrokenLinkFile2D(const char *filename, int max, int *****L1, 
			 int *****L2, int type){
  FILE *in;
  char keyword[100];

  if(max<1)
    return 1;
  if(!*L1 || !*L2)
    return 2;
  in = fopen(filename,"rt");
  if(!in)
    return 3;

  do{
    fscanf(in,"%s",keyword);
  }while(STRNEQ(keyword,"BEGINBL") && !feof(in));

  if(feof(in))
    return 6;


  while(STRNEQ(keyword,"ENDBL")){
    /* In this loop we read a sequence of keywords and interpret them. 
     * We finish when we find a ENDBL keyword.
     */

    int xi,yi,xf,yf;
    int j,k, auxj,auxk;

    fscanf(in,"%s",keyword);

    if(STREQ(keyword,"POINT")){
      /* If we find a POINT keyword we read its coordinates and, if 
       * they are valid, we break all the four links conected to this
       * point, i.e., we isolate the point.
       */

      fscanf(in,"%d",&xi);
      if(xi<-max || xi>max)
	return 4;

      fscanf(in,"%d",&yi);
      if(yi<-max || yi>max)
	return 4;

      for(j=0; j<2; j++){ 
	auxj = (int)pow(-1.0, j);

	for(k=0; k<2; k++){  
	  auxk = (int)pow(-1.0,k);

	  /* Remember that xi and yi passed by the user are ranging
	   * from -max to max while the computer understands only
	   * coordinates ranging from 0 to 2*max+1. Therefore we must
	   * add max to the coordinates in order to make the correct
	   * conversion.
	   */
	  (*L1)[j][k][max+xi][max+yi] = 0; 
	  if(type==DIAG_LATT)
	    (*L2)[j][k][max+xi][max+yi] = 0; 

	  /* If the point is on a boundary site of the lattice 
	   * then we don't have to set the complement (it would
	   * make no sense) 
	   */
	  if(xi==max && auxj==1) 
	    continue;
	  if(xi==-max && auxj==-1)
	    continue;
	  if(yi==max && auxk==1)
	    continue;
	  if(yi==-max && auxk==-1)
	    continue;

	  /* In any other case we MUST set the complement to zero */
	  if(type==DIAG_LATT){
	    (*L1)[1-j][1-k][max+xi+auxj][max+yi+auxk] = 0;
	    (*L2)[1-j][1-k][max+xi+auxj][max+yi+auxk] = 0;
	  }
	  else
	    (*L1)[1-j][1-k][max+xi+(auxj*(1-DELTA(j,k)))][max+yi+(auxj*DELTA(j,k))] = 0;

	}/* end-for k */
      }/* end-for j */
 

    }
    else if(STREQ(keyword,"LINE")){
      /* If we find a LINE keyword we read its coordinates and, if 
       * they are valid, we calculate the points intercepted by this
       * line. Each of these points are isolated according to the
       * comment given for the previous "if".
       */

      int distx, disty, xvar, yvar;
      int t;

      fscanf(in,"%d",&xi);
      if(xi<-max || xi>max)
	return 4;

      fscanf(in,"%d",&yi);
      if(yi<-max || yi>max)
	return 4;

      fscanf(in,"%d",&xf);
      if(xf<-max || xf>max)
	return 4;

      fscanf(in,"%d",&yf);
      if(yf<-max || yf>max)
	return 4;

      distx = abs(xf-xi);
      disty = abs(yf-yi);

      /* The line must be vertical, horizontal, or must have 
       * slope $\pm 1$ 
       */
      if(distx && disty && distx-disty)
	return 5;
      
      /* xvar and yvar indicates if the line is horizontal, vertical or 
       * diagonal. Note that the operator (condition ? A : B) in C language 
       * simply means A if condition is true or B if condition is false.
       * See, for example: Schildt, H. "C: the complete reference".
       */
      xvar = (xf==xi ? 0 : (xf-xi)/distx);
      yvar = (yf==yi ? 0 : (yf-yi)/disty);
      
      for(t=0; t<= MAXIMUM(distx,disty); t++){
	int m,n;

	m = xi + xvar*t;
	n = yi + yvar*t;

	for(j=0; j<2; j++){  
	  auxj = (int)pow(-1.0, j);

	  for(k=0; k<2; k++){
	    auxk = (int)pow(-1.0,k);

	    /* Remember that the coordinates passed by the user are ranging
	     * from -max to max while the computer understands only
	     * coordinates ranging from 0 to 2*max+1. Therefore we must
	     * add max to the coordinates in order to make the correct
	     * conversion.
	     */
	    (*L1)[j][k][max+m][max+n] = 0; 
	    if(type==DIAG_LATT)
	      (*L2)[j][k][max+m][max+n] = 0; 

	    /* If the point is on a boundary site of the lattice 
	     * then we don't have to set the complement (it would
	     * make no sense)
	     */
	    if(m==max && auxj==1) 
	      continue;
	    if(m==-max && auxj==-1)
	      continue;
	    if(n==max && auxk==1)
	      continue;
	    if(n==-max && auxk==-1)
	      continue;

	    /* In any other case we MUST set the complement to zero */
	    if(type==DIAG_LATT){
	      (*L1)[1-j][1-k][max+m+auxj][max+n+auxk] = 0;
	      (*L2)[1-j][1-k][max+m+auxj][max+n+auxk] = 0;
	    }
	    else
	      (*L1)[1-j][1-k][max+m+(auxj*(1-DELTA(j,k)))][max+n+(auxj*DELTA(j,k))] = 0;

	  }/* end-for k */
	}/* end-for j */
      }/* end-for t */
	      
    }/* end-if */
  }/* end-while */

  fclose(in);

  return 0;
}


int writeScreen(const char *filename, screen_t screen){
  FILE *out;
  int n;

  out = fopen(filename, "wt");
  if(!out)
    return 1;

  fprintf(out, "#Index\tX\tY\tIntensity\n");
  for(n=0; n< screen.numpts; n++)
    if(screen.values[n]>0.0)
      fprintf(out, "%d\t%d\t%d\t%e\n", n, screen.xa + n*screen.xvar, 
	      screen.ya + n*screen.yvar, screen.values[n]);

  fclose(out);
  
  return 0;
}


int writeScript1D(filenames_t fname, options1D_t options){
  FILE *out;
  time_t lt;

  out = fopen(fname.gpt_file, "wt");
  if(!out)
    return 1;

  fprintf(out, "#Script generated by Quantum Walk Simulator (1D)\n");

  lt = time(NULL);
  fprintf(out, "# %s\n\n",ctime(&lt));

  fprintf(out,"reset\n\n");

  fprintf(out,"###########################\n");
  fprintf(out,"# Final distribution plot #\n");
  fprintf(out,"###########################\n");
  
  fprintf(out,"print \"Generating plot for final distribution. Please, wait...\"\n");
  fprintf(out,"set grid\n");
  fprintf(out,"set xlabel \"x\"\n");
  fprintf(out,"set ylabel \"Probability\"\n");
  if(options.lattType == LINE_LATT)
    fprintf(out,"set xrange [-%d:%d]\n", options.max-1, options.max-1);
  else
    fprintf(out,"set xrange [0:%d]\n", options.max-1);
  fprintf(out,"set xtics %d\n", (options.max-1)/5);
  fprintf(out,"set mxtics 5\n");

  fprintf(out,"set terminal postscript enhanced monochrome \'Times\' 22\n");
  fprintf(out,"set output \"%s\"\n", fname.eps2d_file);
  fprintf(out,"plot \"%s\" title \"\" with lines lw 2\n",fname.dat_file);
  fprintf(out,"print \"Plot saved as %s.\"\n\n", fname.eps2d_file);

  if(options.calcMix){
    fprintf(out,"##############################\n");
    fprintf(out,"# Limiting distribution plot #\n");
    fprintf(out,"##############################\n");

    fprintf(out,"print \"Generating plot for limiting distribution. Please, wait...\"\n");
    fprintf(out,"set output \"%s\"\n", fname.epspb_file);
    fprintf(out,"plot \"%s\" title \"limiting\" with lines lw 2, 1.0/%d title\"uniform\"\n",
	    fname.datpb_file, options.max);
    fprintf(out,"print \"Plot saved as %s.\"\n\n", fname.epspb_file);
  }

  fprintf(out,"print \"Done. If the plot present any inconsistency, try to adjust");
  fprintf(out," the script manually and run gnuplot again.\"\n");

  fclose(out);

  return 0;
}


int writeScript2D(filenames_t fname, options2D_t options){
  FILE *out;
  time_t lt;

  out = fopen(fname.gpt_file, "wt");
  if(!out)
    return 1;

  fprintf(out, "#Script generated by Quantum Walk Simulator (2D).\n");

  lt = time(NULL);
  fprintf(out, "# %s\n\n",ctime(&lt));

  fprintf(out,"reset\n\n");
  fprintf(out,"print \"The script may need fine-tuning to meet the user preferences.\"\n");

  if(options.screen){
    fprintf(out,"###########################\n");
    fprintf(out,"# Observation screen plot #\n");
    fprintf(out,"###########################\n");
    fprintf(out,"print \"Generating observation screen plot. Please, wait...\"\n");
    fprintf(out,"set xlabel \"site number\"\n");
    fprintf(out,"set ylabel \"Intensity\"\n");
    fprintf(out,"set grid\n");
    fprintf(out,"set terminal postscript eps monochrome \'Times\' 22\n");
    fprintf(out,"set output \"%s\"\n", fname.epsscr_file);
    fprintf(out,"plot \"%s\" using 1:4 title \"\" with lines lw 2.5\n", fname.datscr_file);
    fprintf(out,"# You can replace \"using 1:4\" by \"using 2:4\" or \"using 3:4\"\n");
    fprintf(out,"# if the observation screen is respectively in horizontal or vertical\n");
    fprintf(out,"# position.\n\n");
    fprintf(out,"print \"Observation screen plot saved as %s.\"\n\n", fname.epsscr_file);
  }

  fprintf(out,"###########\n");
  fprintf(out,"# 3D plot #\n");
  fprintf(out,"###########\n");
  fprintf(out,"print \"Generating 3D plot. Please, wait...\"\n");
  fprintf(out,"set border 4095\n");
  if(options.lattType == CYCLE_LATT){
    fprintf(out,"set xrange [0:%d]", options.max-1);
    fprintf(out,"\t\t# these options may require fine-tuning!\n");
    fprintf(out,"set yrange [0:%d]\n", options.max-1);
  }
  else{
    fprintf(out,"set xrange [-%d:%d]", options.max-1, options.max-1);
    fprintf(out,"\t\t# these options may require fine-tuning!\n");
    fprintf(out,"set yrange [-%d:%d]\n", options.max-1, options.max-1);
  }
  fprintf(out,"set ztics\n");
  fprintf(out,"set xtics %d\n", (options.max-1)/4);
  fprintf(out,"set ytics %d\n", (options.max-1)/4);
  fprintf(out,"set mxtics 5\n");
  fprintf(out,"set mytics 5\n");
  fprintf(out,"set zlabel \"Probability\"\n");
  fprintf(out,"set xlabel \"x\"\n");
  fprintf(out,"set ylabel \"y\"\n");
  fprintf(out,"set ticslevel 0\n");
  fprintf(out,"set size 1,1\n");
  fprintf(out,"set hidden3d\n");
  fprintf(out,"set view 53.00,44.00\n");
  fprintf(out,"set style data lines\n");
  fprintf(out,"set dgrid3d %d,%d,16\n", options.max, options.max);
  fprintf(out,"# If this option is not well adjusted the graphic may present inconsistent results.\n");
  fprintf(out,"# Hint: See the maximum absolute value of the coordinates in %s and add one.\n",
	  fname.dat_file);
  fprintf(out,"#       If problem persists, try to define a huge grid. It will be slower but safer.\n\n");
  fprintf(out,"set terminal postscript enhanced monochrome \'Times\' 22\n");
  fprintf(out,"set output \"%s\"\n", fname.eps3d_file);
  fprintf(out,"splot \"%s\" using 1:2:3 title \"\" with lines\n",fname.dat_file);
  fprintf(out,"print \"3D plot saved as %s.\"\n\n", fname.eps3d_file);

  if(options.calcMix){
    fprintf(out,"##############################\n");
    fprintf(out,"# Limiting distribution plot #\n");
    fprintf(out,"##############################\n");

    fprintf(out,"print \"Generating plot for limiting distribution. Please, wait...\"\n");
    fprintf(out,"set output \"%s\"\n", fname.epspb_file);
    fprintf(out,"splot \"%s\" using 1:2:3 title \"\" with lines\n",fname.datpb_file);
    fprintf(out,"print \"Plot saved as %s.\"\n\n", fname.epspb_file);
  }

  fprintf(out,"###########\n");
  fprintf(out,"# 2D plot #\n");
  fprintf(out,"###########\n");
  fprintf(out,"print \"Generating 2D plot. Please, wait...\"\n");
  fprintf(out,"set samples 50, 50\n");
  fprintf(out,"set isosamples 50, 50\n");
  fprintf(out,"set contour base\n");
  fprintf(out,"unset surface\t\t#in old versions of gnuplot, replace by set nosurface\n");
  fprintf(out,"set cntrparam levels auto 500\t\t#this parameter may need to be fine-tuned by the user\n");
  fprintf(out,"set zrange [0.0:1.0] noreverse nowriteback\n");
  fprintf(out,"set size square\n");
  fprintf(out,"set table \'contour.tbl\'\n");
  fprintf(out,"splot \"%s\" with lines\n", fname.dat_file);
  fprintf(out,"unset table\n");
  fprintf(out,"set terminal postscript enhanced monochrome 'Times' 22\n");
  fprintf(out,"unset label\t\t#in old versions of gnuplot, replace by set nolabel\n");
  fprintf(out,"set grid\n");
  fprintf(out,"set output \"%s\"\n", fname.eps2d_file);
  fprintf(out,"unset key\t\t#in old versions of gnuplot, replace by set nokey\n");
  fprintf(out,"plot 'contour.tbl' with lines\n");
  fprintf(out,"print \"2D plot saved as %s.\"\n", fname.eps2d_file);

  fprintf(out,"print \"Done. If the plots present any inconsistency, try to adjust");
  fprintf(out," the script manually and run gnuplot again.\"\n");

  fclose(out);

  return 0;
}

int changeFileEnding(char **newname, const char *oldname, const char *end){
  const int aux = strlen(end);

  *newname = (char *)malloc((strlen(oldname)+aux)*sizeof(char));
  strcpy(*newname, oldname);
  strtok(*newname, ".");
  strcat(*newname, end);

  return 0;
}