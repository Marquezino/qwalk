/* QWalk (qwcoin.c) 
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
#include "qwcoin.h"


double complex ****createHadamardCoin2D(){
  double complex ****coin;

  coin = allocComplex4D(2, 2, 2, 2);
  if(!coin)
    return NULL;

  coin[0][0][0][0] = 0.5;
  coin[0][0][0][1] = 0.5;
  coin[0][0][1][0] = 0.5;
  coin[0][0][1][1] = 0.5;

  coin[0][1][0][0] = 0.5;
  coin[0][1][0][1] =-0.5;
  coin[0][1][1][0] = 0.5;
  coin[0][1][1][1] =-0.5;

  coin[1][0][0][0] = 0.5;
  coin[1][0][0][1] = 0.5;
  coin[1][0][1][0] =-0.5;
  coin[1][0][1][1] =-0.5;

  coin[1][1][0][0] = 0.5;
  coin[1][1][0][1] =-0.5;
  coin[1][1][1][0] =-0.5;
  coin[1][1][1][1] = 0.5;

  return coin;
}


double complex ****createFourierCoin2D(){
  double complex ****coin;

  coin = allocComplex4D(2, 2, 2, 2);
  if(!coin)
    return NULL;

  coin[0][0][0][0] = 0.5;
  coin[0][0][0][1] = 0.5;
  coin[0][0][1][0] = 0.5;
  coin[0][0][1][1] = 0.5;

  coin[0][1][0][0] = 0.5;
  coin[0][1][0][1] = 0.5*I;
  coin[0][1][1][0] =-0.5;
  coin[0][1][1][1] =-0.5*I;

  coin[1][0][0][0] = 0.5;
  coin[1][0][0][1] =-0.5;
  coin[1][0][1][0] = 0.5;
  coin[1][0][1][1] =-0.5;

  coin[1][1][0][0] = 0.5;
  coin[1][1][0][1] =-0.5*I;
  coin[1][1][1][0] =-0.5;
  coin[1][1][1][1] = 0.5*I;

  return coin;
}


double complex ****createGroverCoin2D(){
  double complex ****coin;

  coin = allocComplex4D(2, 2, 2, 2);
  if(!coin)
    return NULL;

  coin[0][0][0][0] =-0.5;
  coin[0][0][0][1] = 0.5;
  coin[0][0][1][0] = 0.5;
  coin[0][0][1][1] = 0.5;

  coin[0][1][0][0] = 0.5;
  coin[0][1][0][1] =-0.5;
  coin[0][1][1][0] = 0.5;
  coin[0][1][1][1] = 0.5;

  coin[1][0][0][0] = 0.5;
  coin[1][0][0][1] = 0.5;
  coin[1][0][1][0] =-0.5;
  coin[1][0][1][1] = 0.5;

  coin[1][1][0][0] = 0.5;
  coin[1][1][0][1] = 0.5;
  coin[1][1][1][0] = 0.5;
  coin[1][1][1][1] =-0.5;

  return coin;
}


double complex **createHadamardCoin1D(){
  double complex **coin;

  coin = allocComplex2D(2,2);
  if(!coin)
    return NULL;

  coin[0][0] = 1.0/sqrt(2.0);
  coin[0][1] = 1.0/sqrt(2.0);
  coin[1][0] = 1.0/sqrt(2.0);
  coin[1][1] =-1.0/sqrt(2.0);

  return coin;
}
