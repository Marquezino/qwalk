/* QWalk (qwamplify) This software amplifies regions of wave-equations
 * Copyright (C) 2007  Franklin Marquezino
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

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include "qwconsts.h"
#include "qwamplify.h"

int main(int argc, char *argv[]){
  FILE *in, *out;
  int i, type, error;
  int xi, yi, xf, yf;
  float factor;
  char *filename;

  printf("QWalk Amplify, version 1.2 (qwamplify).\n");
  printf("Copyright (C) 2007 Franklin Marquezino.\n");
  printf("This is free software; see the source code for copying conditions.\n");
  printf("There is ABSOLUTELY NO WARRANTY; not even for MERCHANTIBILITY or\n");
  printf("FITNESS FOR A PARTICULAR PURPOSE. See the source code for details.\n\n");

  printf("Email bug reports to <franklin@lncc.br>\n\n");

  if(argc<2){
    printf("Missing arguments.\n");
    helpMessage();
    exit(EXIT_FAILURE);
  }

  type = NO;

  for(i=2; i<argc; i++){
    if(STREQ(argv[i],"-f")){
      sscanf(argv[i+1], "%f", &factor);      
      i++;
      continue;
    }
    else if(STREQ(argv[i],"-ir")){
      type = IR;
      sscanf(argv[i+1], "%d", &xi);
      sscanf(argv[i+2], "%d", &yi);
      sscanf(argv[i+3], "%d", &xf);
      sscanf(argv[i+4], "%d", &yf);

      if(xf<xi || yf<yi){
	printf("You can't enter neither xi>xf nor yi>xf.\n");
	helpMessage();
	return EXIT_FAILURE;
      }

      i+=4;
      continue;
    }
    else if(STREQ(argv[i],"-or")){
      type = OR;
      sscanf(argv[i+1], "%d", &xi);
      sscanf(argv[i+2], "%d", &yi);
      sscanf(argv[i+3], "%d", &xf);
      sscanf(argv[i+4], "%d", &yf);

      if(xf<xi || yf<yi){
	printf("You can't enter neither xi>xf nor yi>xf.\n");
	helpMessage();
	return EXIT_FAILURE;
      }

      i+=4;
      continue;
    }
    else if(STREQ(argv[i],"-ud")){
      int dx, dy;

      type = UD;
      sscanf(argv[i+1], "%d", &xi);
      sscanf(argv[i+2], "%d", &yi);
      sscanf(argv[i+3], "%d", &xf);
      sscanf(argv[i+4], "%d", &yf);

      if(xf<xi){
	printf("You can't enter xi>xf.\n");
	helpMessage();
	return EXIT_FAILURE;
      }

      dx = xf-xi;
      dy = yf-yi;
      if( abs(dx)!=abs(dy) ){
	printf("Option -ud requires a diagonal (45 degrees) line\n");
	helpMessage();
	return EXIT_FAILURE;
      }

      i+=4;
      continue;
    }
    else if(STREQ(argv[i],"-od")){
      int dx, dy;

      type = OD;
      sscanf(argv[i+1], "%d", &xi);
      sscanf(argv[i+2], "%d", &yi);
      sscanf(argv[i+3], "%d", &xf);
      sscanf(argv[i+4], "%d", &yf);

      if(xf<xi){
	printf("You can't enter xi>xf.\n");
	helpMessage();
	return EXIT_FAILURE;
      }

      dx = xf-xi;
      dy = yf-yi;
      if( abs(dx)!=abs(dy) ){
	printf("Option -od requires a diagonal (45 degrees) line\n");
	helpMessage();
	return EXIT_FAILURE;
      }

      i+=4;
      continue;
    }
    else{
      printf("Wrong arguments.\n");
      helpMessage();
      return EXIT_FAILURE;
    }/* end-if */

  }/* end-for i */


  filename = (char *)malloc((strlen(argv[1])+5)*sizeof(char));
  strcpy(filename,argv[1]);
  strcat(filename,"~~");
  error = rename(argv[1], filename);
  if(error){
    printf("Error: something wrong with the input filename.\n");
    return EXIT_FAILURE;
  }
  printf("Backup: input file renamed to %s\n",filename);

  in = fopen(filename, "rt");
  if(!in){
    printf("Error: could not open input file %s.\n", filename);
    return EXIT_FAILURE;
  }
  out= fopen(argv[1], "wt");
  if(!out){
    printf("Error: could not open output file %s.\n", argv[1]);
    return EXIT_FAILURE;
  }

  fprintf(out, "# This data file was amplified from %s using qwamplify\n", filename);
  fprintf(out, "# Comments from original input file were ignored.\n");

  /* While the end of file is not reached */
  while(!feof(in)){
    int x,y;
    float prob;
    char line[200];

    fgets(line,200,in);

    /* If the line is a gnuplot comment we ignore it */
    if(line[0]=='#')
      continue;

    /* If we read a blank line then we copy this blank line 
     * in the output file 
     */
    if(line[0]=='\n'){
      fprintf(out,"\n");
      continue;
    }

    sscanf(line, "%d%d%f", &x, &y, &prob);

    switch(type){
      int dist,slope;
    case IR:
      if(xi <= x && x <= xf && yi <= y && y <= yf)
	prob *= factor;
      break;
    case OR:
      if(xi > x || x > xf || yi > y || y > yf)
	prob *= factor;
      break;
    case UD:
      dist = (x-xi);
      slope = (yf-yi) ? (yf-yi)/abs(yf-yi) : 0;  
      /* Note that the operator (condition ? A : B) in C language 
       * simply means A if condition is true or B if condition is false.
       * See, for example: Schildt, H. "C: the complete reference".
       */

      if(xi <= x && x <= xf && y < yi + slope*dist )
	prob *= factor;
      break;
    case OD:
      dist = (x-xi);
      slope = (yf-yi) ? (yf-yi)/(yf-yi) : 0;
      if(xi <= x && x <= xf && y > yi + slope*dist )
	prob *= factor;
      break;
    }

    if(prob>0.0)
      fprintf(out,"%d\t%d\t%e\n", x, y, prob);
  }
  
  fclose(in);
  fclose(out);

  printf("Output saved in %s\n\n", argv[1]);

  return EXIT_SUCCESS;
}

void helpMessage(){
  printf("Usage: qwamplify file [options]\n\n");
  printf("Options:\n");
  printf(" -f factor\t\tspecifies the amplification factor\n");
  printf(" -ir xi yi xf yf\tamplifies inside the rectangular region (x,y)\n");
  printf("                \tsuch that xi <= x <= xf and yi <= y <= yf\n");
  printf(" -or xi yi xf yf\tamplifies outside the rectangular region (x,y)\n");
  printf("                \tsuch that xi <= x <= xf and yi <= y <= yf\n");
  printf(" -ud xi yi xf yf\tamplifies the region under the diagonal line\n");
  printf("                \tdelimited by points (xi,yi) and (xf,yf)\n");
  printf(" -od xi yi xf yf\tamplifies the region over the diagonal line\n");
  printf("                \tdelimited by points (xi,yi) and (xf,yf)\n");

  return;
}
