/* QWalk (qwscreen.c) 
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
#include "qwscreen.h"
#include "qwconsts.h" 


int initScreen(screen_t *screen, int xa, int ya, int xb, int yb){
  int distx, disty, j;

  if(!screen)
    return 1;

  screen->xa = xa;
  screen->xb = xb;
  screen->ya = ya;
  screen->yb = yb;

  distx = abs(xb-xa);
  disty = abs(yb-ya);

  /* The screen must be in vertical orientation, 
   * horizontal, or must have slope $\pm 1$ 
   */
  if(distx && disty && distx-disty)
    return 2;

  /* xvar and yvar indicates if the line is horizontal, vertical or 
   * diagonal. Note that the operator (condition ? A : B) in C language 
   * simply means A if condition is true or B if condition is false.
   * See, for example: Schildt, H. "C: the complete reference".
   */
  screen->xvar = (screen->xb == screen->xa ? 
		   0 : (screen->xb - screen->xa)/distx);
  screen->yvar = (screen->yb == screen->ya ? 
		   0 : (screen->yb - screen->ya)/disty);

  screen->numpts = MAXIMUM(distx, disty)+1;
  screen->values = (double *)malloc(screen->numpts*sizeof(double));
  if(!screen->values)
    return 3;

  for(j=0; j< screen->numpts; j++)
    screen->values[j] = 0.0;

  return 0;
}

int updateScreen(screen_t *screen, double complex ****state, int max){
  int t;

  if(!screen)
    return 1;
  if(!state)
    return 2;
  if(max<1)
    return 3;
 
  for(t=0; t< screen->numpts; t++){
    int m, n, j, k;

    m = max + screen->xa + screen->xvar*t;
    n = max + screen->ya + screen->yvar*t;
    
    for(j=0; j<2; j++)
      for(k=0; k<2; k++)
	screen->values[t] += state[j][k][m][n]*conj(state[j][k][m][n]);
  }
  return 0;
}
