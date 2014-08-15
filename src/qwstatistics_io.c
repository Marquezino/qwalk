/* QWalk (qwstatistics_io.c) 
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
#include "qwstatistics.h"
#include "qwstatistics_io.h"
#include "qwconsts.h"


int writeStatistics(const char *filename, statistics_t stat){
  FILE *out;
  static char firstCall = 1;

  if(STREQ(filename,"")){
    firstCall=1;
    return 0;
  }

  if(firstCall){
    out = fopen(filename,"wt");
    if(!out)
      return 1;
    fprintf(out,"#Iter\tMean X\t\tMean Y\t\tVariance\tStd deviation\tTVD\tTVD (unif)\n\n");
    firstCall = 0;
  }
  else{
    out = fopen(filename,"at");
    if(!out)
      return 1;
  }


  fprintf(out,"%d\t%e\t%e\t%e\t%e\t%e\t%e\n",
	  stat.iteration, stat.meanX, stat.meanY,
	  stat.variance, sqrt(stat.variance),
	  stat.tvd, stat.tvdu);

  fclose(out);

  return 0;
}


