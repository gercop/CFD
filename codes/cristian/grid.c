/* -|----|----|----|----|----|----|----|----|----|----|----|----|----|*/
/*      10        20        30        40        50        60        70*/
/* -|----|----|----|----|----|----|----|----|----|----|----|----|----|*/
/*
Set up the grid and stretching factors. Save them to files.
Calculate delta (= average grid size) in case of RANS/LES.

This routine is a little bit messy, since naming conventions
are not obvious (e.g., using x,y instead of i,j for point-indeces
and different use of xi,eta in program and printout, etc).
*/
 
#include <string.h>
#include <math.h>
#include "externs.h"
#include "macros.h"
#include "functype.h"
 

/* -|----|----|----|----|----|----|----|----|----|----|----|----|----|*/
int setup_grid (void)
{
int	d, x, y, x0, i;
REAL	xi, mapy0;
#if (defined(STRETCH_X) || defined(STRETCH_Y))
  int	max_size;
  REAL	mapf0, mapf1, mapf3, mapf5;
  #ifdef IOSLIB_ON
  #ifdef THREE_D
  char	*sflags[2] = {"ggg", "sss"};
  #else /* TWO_D */
  char	*sflags[2] = {"gg", "ss"};
  #endif /* THREE_D */
  #endif /* IOSLIB_ON */
#endif /* defined(STRETCH_X) || defined(STRETCH_Y) */
#ifdef STRETCH_X
  REAL	xi_mid, xi_0;
#endif /* STRETCH_X */
#ifdef STRETCH_Y
  REAL	eta, h, h_s, K_0, mapf2, true_slope;
  int	y0;
  REAL	mapy0a;
  #ifdef ANY_WAKE
  REAL	temp, eta_mid, K_1;
  #endif /* ANY_WAKE */
#endif /* STRETCH_Y */
paramst **parpt = prmpts[0]->params,
 	*parpt0 = prmpts[0]->params[0];


/*----------------------------------------------------------------*/
/* X-GRID */
/*----------------------------------------------------------------*/
/* The x-zero is placed on the base; for a flat plate it is placed
   at the left boundary.
   Depending on the value of x_map_type, different stretching functions
   are applied:
   0 = none (automatically set if STRETCH_X is not defined)
   1 = polynomial stretching (default if STRETCH_X is defined)
*/

/*----------------------------------------------------------------*/
if (param0.x_map_type == STRETCH_NONE) { 
 	/* No grid-stretching in x (only equidistant grids) */
 printf ("x-stretching, none:      %12.8f x\n", 1.);

 for (x=0; x<param0.x_shape_size[0]; x++) {
   miscdat(0).x_grid_all[x] = (REAL)(x - parpt0->x_diff) * parpt0->dx;
   miscdat(0).x_slope_all[x] = 1.;
 }
} 

#ifdef STRETCH_X
/*----------------------------------------------------------------*/
else if (param0.x_map_type == STRETCH_POLY) {
 	/* Polynomial x-stretching:

 The seven constants are:

 x_mapf[0] - stretching factor at the base wall
 x_mapf[1] - stretching factor at the outflow (see 3)
 x_mapf[2] - stretching factor at the inflow; if this is not set or is
	< zero, the same stretching is used upstream as downstream
 x_mapf[3] - from this location downstream to the outflow, the grid is
	again constant; if this is not set or is <zero, then a 3rd-order
	stretching is used all the way to the outflow (see 1). 
	Otherwise, a 5th-order stretching is used from the base to this
	location, then constant to the base.
 x_mapf[4] - stretching value to outflow.  If this is present (1) is
	just used to the point specified by (5), and then after ramping
	from (5) to (6), (4) is used.
 x_mapf[5] - begin change to outflow value
 x_mapf[6] - end change to outflow value (if not, then 3rd to outflow)
 */

 /*
 Note: to change the domain size and keep the stretching, change map[1]
 as follows (this is if a 3rd-order mapping is used all the way
 downstream, or towards the free-stream):

	new_map[1] = (old_map[1] - map[0]) * (new_size-1)^2 + map [0]
					     --------------
					     (old_size-1)^2
 */

 /* upstream mapping */

	/* use a different stretching upstream, if set */
 if (parpt0->x_diff > 0) {
   if (param0.x_mapf[2] > 0.) {
      mapf1 = param0.x_mapf[0];
      mapf3 = (param0.x_mapf[2] - param0.x_mapf[0])
	/ (3.* pow (parpt0->x_diff*parpt0->dx, 2.));
   } else {
      mapf1 = param0.x_mapf[0];
      mapf3 = (param0.x_mapf[1] - param0.x_mapf[0])
	/ (3.* pow ((param0.x_shape_size0[0]-parpt0->x_diff)
						*parpt0->dx, 2.));
   }

   printf ("x-stretching, upstream:  %12.8f x + %12.8f x^3\n", 
   		mapf1, mapf3);
   printf ("    down to xi = %12.8f (pt %d)\n",
		0., parpt0->x_diff);

   for (x=0; x<parpt0->x_diff; x++) {
      xi = (REAL)(x - parpt0->x_diff) * parpt0->dx;
      miscdat(0).x_grid_all[x] = (mapf3 * xi*xi + mapf1) * xi;
      miscdat(0).x_slope_all[x] = 1./(3.*mapf3 * xi*xi + mapf1);
   }
 }

 /* downstream mapping */

 if (param0.x_mapf[3] > 0.) {	/* double, 5th/1st-order stretching */

	/* 5th-order zone, up to param0.x_mapf[3] */

	/* xi = 15*x2 / (8*x0' + 7*x2') */
   xi_mid = 15.* param0.x_mapf[3] /
	(8.* param0.x_mapf[0] + 7.*param0.x_mapf[1]);
   mapf1 = param0.x_mapf[0];
   mapf3 = 2. * (param0.x_mapf[1] - param0.x_mapf[0])
	/ (3.* xi_mid*xi_mid);
   mapf5 = -0.3 * mapf3
		/ (xi_mid*xi_mid);

   printf ("x-stretching, by base:   %12.8f x + %12.8f x^3 + %12.8f x^5\n",
				mapf1, mapf3, mapf5);
   printf ("    down to xi = %12.8f (pt %d)\n",
		xi_mid, (int)(xi_mid/parpt0->dx)+parpt0->x_diff);

   for (x=parpt0->x_diff;
		x <= (int)(xi_mid/parpt0->dx)+parpt0->x_diff
			 && x<param0.x_shape_size[0]; x++) {
      xi = (REAL)(x - parpt0->x_diff) * parpt0->dx;
      miscdat(0).x_grid_all[x]
		= ((mapf5 * xi*xi + mapf3) * xi*xi + mapf1) * xi;
      miscdat(0).x_slope_all[x]
		= 1./((5.*mapf5 * xi*xi +3.*mapf3) * xi*xi + mapf1);
   }

	/* 1st-order zone, param0.x_mapf[2] to outflow (maybe) */

   mapf1 = param0.x_mapf[1];
   mapf0 = param0.x_mapf[3] - xi_mid * param0.x_mapf[1];

   printf ("x-stretching, downstream:%12.8f x + %12.8f\n",
				mapf1, mapf0);

   for ( ; x<param0.x_shape_size[0]; x++) {
      xi = (REAL)(x - parpt0->x_diff) * parpt0->dx;
      miscdat(0).x_grid_all[x]
		= mapf1 * xi + mapf0;
      miscdat(0).x_slope_all[x]
		= 1./mapf1;
   }

   if (param0.x_mapf[4] > 0.
	&& param0.x_mapf[5] > param0.x_mapf[3]
	&& param0.x_mapf[6] > param0.x_mapf[5]) {
	/* Add more stretching near outflow--5th to 1st. */

	/* begin stretch location */
      xi_0 = (param0.x_shape_size0[0] - parpt0->x_diff) * parpt0->dx
	- (miscdat(0).x_grid_all[param0.x_shape_size0[0]]
				- param0.x_mapf[5]) / param0.x_mapf[1];
	/* length in xi-coords */
      xi_mid = 15.* (param0.x_mapf[5] - param0.x_mapf[6]) /
	(8.* param0.x_mapf[4] + 7.*param0.x_mapf[1]);
      mapf0 = param0.x_mapf[6];
      mapf1 = param0.x_mapf[4];
      mapf3 = 2. * (param0.x_mapf[1] - param0.x_mapf[4])
	/ (3.* xi_mid*xi_mid);
      mapf5 = -0.3 * mapf3
		/ (xi_mid*xi_mid);

		/* check to make sure within domain limits */
      if ((int)((xi_0-xi_mid)/parpt0->dx)+parpt0->x_diff
					<= param0.x_shape_size[0]) {

         printf ("x-stretching, near out:  %12.8f x + %12.8f x^3 + %12.8f x^5\n",
				mapf1, mapf3, mapf5);
         printf ("    from xi = %12.8f (pt %d) to xi = %12.8f (pt %d)\n",
		xi_0, (int)(xi_0/parpt0->dx)+parpt0->x_diff, 
		xi_0+xi_mid,
		(int)((xi_0-xi_mid)/parpt0->dx)+parpt0->x_diff);

         for (x=(int)(xi_0/parpt0->dx)+parpt0->x_diff+1;
	 	x <= (int)((xi_0-xi_mid)/parpt0->dx)+parpt0->x_diff
			 && x<param0.x_shape_size[0]; x++) {
            xi = (REAL)(x - parpt0->x_diff) * parpt0->dx - xi_0 
	    		+ xi_mid;
            miscdat(0).x_grid_all[x]
		= ((mapf5 * xi*xi + mapf3) * xi*xi + mapf1) * xi +mapf0;
            miscdat(0).x_slope_all[x]
		= 1./((5.*mapf5 * xi*xi +3.*mapf3) * xi*xi + mapf1);
         }

	/* 1st-order zone, param0.x_mapf[4] to outflow */

         mapf1 = param0.x_mapf[4];
         mapf0 = param0.x_mapf[6] - (xi_0 - xi_mid) * param0.x_mapf[4];

         printf ("x-stretching, to outflow:%12.8f x + %12.8f \n",
				mapf1, mapf0);

         for ( ; x<param0.x_shape_size[0]; x++) {
            xi = (REAL)(x - parpt0->x_diff) * parpt0->dx;
            miscdat(0).x_grid_all[x]
		= mapf1 * xi + mapf0;
            miscdat(0).x_slope_all[x]
		= 1./mapf1;
         }

      } else {
         printf ("x-stretching, near out:  over limits, %d > %d\n",
			(int)((xi_0-xi_mid)/parpt0->dx)+parpt0->x_diff,
			param0.x_shape_size[0]);
      }


   } else if (param0.x_mapf[4] > 0.
	&& param0.x_mapf[5] > param0.x_mapf[3]
	&& param0.x_mapf[6] < 0.) {
	/* Add more stretching near outflow--3rd to exit */

	/* begin stretch location */
      xi_0 = (param0.x_shape_size0[0] - parpt0->x_diff) * parpt0->dx
	- (miscdat(0).x_grid_all[param0.x_shape_size0[0]]
				- param0.x_mapf[5]) / param0.x_mapf[1];
      x0=(int)(xi_0/parpt0->dx)+parpt0->x_diff;

      mapf0 = param0.x_mapf[5];
      mapf1 = param0.x_mapf[1];
      mapf3 = (param0.x_mapf[4] - param0.x_mapf[1])
		/ (3.* pow ( 
		   (param0.x_shape_size0[0]-parpt0->x_diff)
					* parpt0->dx - xi_0, 2.) );

      printf ("x-stretching, to out:    %12.8f x + %12.8f x^3 + %12.8f\n",
				mapf1, mapf3, param0.x_mapf[5]);
      printf ("    from xi = %12.8f (pt %d)\n",
		xi_0, x0);

      for (x=x0+1; x<param0.x_shape_size[0]; x++) {
         xi = (REAL)(x-parpt0->x_diff) * parpt0->dx - xi_0;
         miscdat(0).x_grid_all[x] = (mapf3 * xi*xi + mapf1) * xi 
	 				+ mapf0;
         miscdat(0).x_slope_all[x] = 1./(3.*mapf3 * xi*xi + mapf1);
      }



   }
   printf ("\n");	/* no more terms */

 } else {			/* single, 3rd-order stretching */

   mapf1 = param0.x_mapf[0];
   mapf3 = (param0.x_mapf[1] - param0.x_mapf[0])
	/ (3.* pow ((param0.x_shape_size0[0]-parpt0->x_diff)*parpt0->dx,
								2.));

   printf ("x-stretching, downstream:%12.8f x + %12.8f x^3\n",
		mapf1, mapf3);

   for (x=parpt0->x_diff; x<param0.x_shape_size[0]; x++) {
      xi = (REAL)(x - parpt0->x_diff) * parpt0->dx;
      miscdat(0).x_grid_all[x] = (mapf3 * xi*xi + mapf1) * xi;
      miscdat(0).x_slope_all[x] = 1./(3.*mapf3 * xi*xi + mapf1);
   }

 }
}

/*----------------------------------------------------------------*/
#endif /* STRETCH_X */

else { 
	/* Invalid stretching option */

 printf ("\n\nERROR: Invalid x-stretching option x_map_type = %d!\n\n",
    param0.x_map_type);
 shutdown();

}

/*----------------------------------------------------------------*/
	/* fill in domains */
for (d=0; d<params.n_domains; d++) {
   for (x=0; x<parpt[d]->x_size; x++) {
      miscdat(d).x_grid[x] = miscdat(0).x_grid_all[x+parpt[d]->x_off];
      miscdat(d).x_slope[x] = miscdat(0).x_slope_all[x+parpt[d]->x_off];
   }
}


/*----------------------------------------------------------------*/
/* Y-GRID */
/*----------------------------------------------------------------*/
/* The y-zero is placed at the center of the base; this places 
   the sidewall at mapy0, i.e., 0.5 for wakes, 0.0 for a flat plate 
   and 1.0 for a step flow. 
   Depending on the value of y_map_type, different stretching functions
   are applied:
   0 = none (automatically set if STRETCH_Y is not defined)
   1 = polynomial stretching (default if STRETCH_Y is defined)
   2 = tanh-mapping
*/

#if defined(ANY_WAKE)

#if (SHAPE == BACK_STEP)
mapy0 = 1.0*parpt0->body_thick;		/* backward-facing step,
					   channel-step */
#else
mapy0 = 0.5*parpt0->body_thick;		/* wakes, lone corner */
#endif /* BACK_STEP */

#else
mapy0 = parpt0->body_thick;		/* flat plate, channel, 
					   free-flows */
#endif /* ANY_WAKE */


/*----------------------------------------------------------------*/
if (param0.y_map_type == STRETCH_NONE) { 
	/* No grid-stretching in y (only equidistant grids) */
 printf ("y-stretching, none:      %12.8f y\n\n", 1.);

 for (y=0; y<param0.y_shape_size[0]; y++) {
   miscdat(0).y_grid[y] = (REAL)y * parpt0->dy + mapy0;
   miscdat(0).y_slope[y] = 1.;
 }

 #if defined(ANY_WAKE)
 for (y=0; y<param0.y_shape_size[1]; y++) {
   miscdat(1).y_grid[y] = (REAL)(y - param0.y_shape_size[1])
						* parpt0->dy + mapy0;
   miscdat(1).y_slope[y] = 1.;
 }
 #endif /* ANY_WAKE */

 #if (SHAPE == FULL_WAKE)
 for (y=0; y<param0.y_shape_size[2]; y++) {
   miscdat(2).y_grid[y] = 
	(REAL)(y - param0.y_shape_size[1] - param0.y_shape_size[2])
						* parpt0->dy + mapy0;
   miscdat(2).y_slope[y] = 1.;
 }
 #endif /* SHAPE == FULL_WAKE */
}

#ifdef STRETCH_Y
/*----------------------------------------------------------------*/
else if (param0.y_map_type == STRETCH_POLY) {
	/* Polynomial y-stretching:

  The four constants are:

  y_mapf[0] - stretching factor at the side walls
  y_mapf[1] - stretching factor at the free stream (or upper wall)
  y_mapf[2] - if set, give distance from y=0 to have a uniform grid
	  across the base
  y_mapf[3] - if set, gives distance towards the free stream from the
	  side walls to have uniform grid.
 */

  mapf1 = param0.y_mapf[0];

  if (param0.y_mapf[3] > 0.) {
   y0 = min( (int)(param0.y_mapf[3]/(parpt0->dy*param0.y_mapf[0])),
		param0.y_shape_size[0]);

   printf ("y-stretching, by wall:   %12.8f y\n", mapf1);
   printf ("    up to eta = %12.8f (pt %d)\n",
				(REAL)y0*parpt0->dy, y0);

   for (y=0; y<y0; y++) {
      eta = (REAL)y * parpt0->dy;
      miscdat(0).y_grid[y] = mapf1 * eta + mapy0;
      miscdat(0).y_slope[y] = 1./mapf1;
     #if (SHAPE == FULL_WAKE)
      miscdat(2).y_grid[param0.y_shape_size0[2]-y] = -mapf1*eta - mapy0;
      miscdat(2).y_slope[param0.y_shape_size0[2]-y] = 1./mapf1;
     #endif /* SHAPE == FULL_WAKE */
   }

   mapy0a = mapy0 + (REAL)y0 * parpt0->dy * mapf1;

  } else {
   y0 = 0;
   mapy0a = mapy0;
  }

  mapf3 = (param0.y_mapf[1] - param0.y_mapf[0])
	/ (3.* pow ((param0.y_shape_size0[0]-y0) * parpt0->dy, 2.));

  printf ("y-stretching, outer:     %12.8f y + %12.8f y^3 + %12.8f\n",
		mapf1, mapf3, (REAL)y0*parpt0->dy);

  for (y=y0; y<param0.y_shape_size[0]; y++) {
   eta = (REAL)(y-y0) * parpt0->dy;
   miscdat(0).y_grid[y] = (mapf3 * eta*eta + mapf1) * eta + mapy0a;
   miscdat(0).y_slope[y] = 1./(3.*mapf3 * eta*eta + mapf1);
  }

  #if (SHAPE == FULL_WAKE)
  for (y=0; y<param0.y_shape_size[2]-y0; y++) {
   eta = (REAL)(y - param0.y_shape_size0[2] + y0) * parpt0->dy;
   miscdat(2).y_grid[y] = (mapf3 * eta*eta + mapf1) * eta - mapy0a;
   miscdat(2).y_slope[y] = 1./(3.*mapf3 * eta*eta + mapf1);
  }
  #endif /* SHAPE == FULL_WAKE */

  /* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

  #if defined(ANY_WAKE)		/* across base */

  if (param0.y_mapf[2] > 0.) {	/* double, 5th/1st-order stretching */

	/* 5th-order zone, from corner to mapy0 - param0.y_mapf[2] */

	/* mapf0 = center grid slope,
	   from y_mapf[2], which is the distance from the corner that
	   the linear begins */
   temp = 0.5 + 4./7.*(param0.y_mapf[2]/mapy0 - param0.y_mapf[0]);
   mapf0 = temp + sqrt( square(temp)
			+ 8./7. * param0.y_mapf[0]
				* (1. - param0.y_mapf[2]/mapy0) );

	/* location in eta coordinates */
   eta_mid = (param0.y_mapf[2] - (1. - mapf0) * mapy0) / mapf0;

   mapf1 = param0.y_mapf[0];
   mapf3 = 2. * (mapf0 - param0.y_mapf[0]) / (3. * eta_mid * eta_mid);
   mapf5 = -0.2 * (mapf0 - param0.y_mapf[0])
		/ (eta_mid * eta_mid * eta_mid * eta_mid);

   printf ("y-stretching, corner:    %12.8f y + %12.8f y^3 + %12.8f y^5\n",
				mapf1, mapf3, mapf5);
   printf ("    down to eta = %8.4f - %12.8f (pt -%d)\n",
		mapy0, eta_mid, (int)(eta_mid/parpt0->dy));

   for (y=1; y <= (int)(eta_mid/parpt0->dy); y++) {
      eta = (REAL)(y) * parpt0->dy;
      miscdat(1).y_grid[param0.y_shape_size[1]-y]
		= mapy0 -
		((mapf5 * eta*eta + mapf3) * eta*eta + mapf1) * eta;
      miscdat(1).y_slope[param0.y_shape_size[1]-y]
		= 1./((5.*mapf5 * eta*eta +3.*mapf3) * eta*eta + mapf1);
   }

	/* 1st-order zone, center base */

   mapf1 = mapf0;

   printf ("y-stretching, center base: %12.8f y\n\n",
				mapf1);

   for (y=parpt0->y_mid;
	y<param0.y_shape_size[1]-(int)(eta_mid/parpt0->dy); y++) {
      eta = (REAL)(y - parpt0->y_mid) * parpt0->dy;
      miscdat(1).y_grid[y]
		= mapf1 * eta;
      miscdat(1).y_slope[y]
		= 1./mapf1;
   }

	/* fill in bottom */
   for (y=0; y<parpt0->y_mid; y++) {
      miscdat(1).y_grid[y]
		= -miscdat(1).y_grid[param0.y_shape_size0[1]-y];
      miscdat(1).y_slope[y]
		= miscdat(1).y_slope[param0.y_shape_size0[1]-y];
   }

  /* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */
  } else {			/* 5th order all the way across */

   mapf1 = 1. + 0.875*(1. - param0.y_mapf[0]);
   mapf3 =     -1.250*(1. - param0.y_mapf[0])/square(mapy0);
   mapf5 =      0.375*(1. - param0.y_mapf[0])/square(square(mapy0));

   printf ("y-stretching, inner:     %12.8f y + %12.8f y^3 + %12.8f y^5\n\n",
				mapf1, mapf3, mapf5);

   for (y=0; y<param0.y_shape_size[1]; y++) {
      eta = (REAL)(y - parpt0->y_mid) * parpt0->dy;
      miscdat(1).y_grid[y]
		= ((mapf5 * eta*eta + mapf3) * eta*eta + mapf1) * eta;
      miscdat(1).y_slope[y]
		= 1./((5.*mapf5 * eta*eta +3.*mapf3) * eta*eta + mapf1);
   }

  }

  #else /* SHAPE == FLAT_PLATE */

  printf ("\n");

  #endif /* ANY_WAKE */
} 

/*----------------------------------------------------------------*/
else if (param0.y_map_type == STRETCH_TANH) {
	/* tanh-mapping in y 

  step/base height on equidistant grid:
    h_s = y_size_1 * dy		for steps, sym. wakes, and corner flows
        = mapy0	+ tanh_eta_y	otherwise

  domain height:
    h   = h_s + (y_size_0 - 1) * dy

  for 0 <= y_grid < mapy0	(behind base/step)
    eta = y * dy		with y=0,y_size_1
    y_grid = K_1 * xi_1
    	     * [ 1 - tanh(gamma_1*(xi_1-eta))/tanh(gamma_1*xi_1) ]

    K_1 = mapy0 / { xi_1 * 
    		[ 1 - tanh(gamma_1*(xi_1-h_s))/tanh(gamma_1*xi_1) ] }

  for mapy0 <= y_grid < h	(above sidewall)
    eta = h_s + y * dy		with y=0,y_size_0
    y_grid = mapy0 + K_0 * (xi_0 - h_s) 
     	     * [ 1 - tanh(gamma_0*(xi_0-eta))/tanh(gamma_0*(xi_0-h_s)) ]

    K_0 = (h - mapy0) / (h - h_s)
    
    this K_0 is only correct for xi_0 = h, otherwise it should be
    
    K_0 = (h - mapy0) / [ (xi_0 - h_s) { 1 - tanh[gamma_0(xi_0 - h)]/
    		tanh[gamma_0(xi_0 - h_s)] } ]

  The four constants are (defaults in brackets):

  y_mapf[0] = mapf0: 	gamma_0 	(1.0)
  y_mapf[1] = mapf1: 	gamma_1		(1.0)
  y_mapf[2] = mapf2: 	xi_0    	((h+h_s)/2)
  y_mapf[3] = mapf3:	xi_1    	(h_s/2)


  NOTE: gamma_0 is currently always read from the input file, but for
        ANY_WAKE set, it should be computed using Newton-Raphson by
	solving

	K_0 * theta - lambda * sinh(theta) = 0

	with theta = 2 * gamma_0 * ( xi_0 - h_s )
	and  lambda= K_1 * gamma_1 * xi_1 /
		( tanh(gamma_1*xi_1)*[cosh(gamma_1*(xi_1-h_s))]^2 )

	So currently, the user has to do this calculation on his own
	in order to determine gamma_0.
  */

 /* grid point of zero y-location */
 #if (SHAPE == FULL_WAKE)
   y0 = 1 + param0.y_shape_size[1]/2;
 #else
   y0 = 0;
 #endif /* SHAPE == FULL_WAKE */

 /* set equidistant step/base height */
 #ifdef ANY_WAKE
   #if (SHAPE == FULL_WAKE)
   h_s = (REAL)(param0.y_shape_size[1]/2) * parpt0->dy;
   #else
   h_s = (REAL)param0.y_shape_size[1] * parpt0->dy;
   #endif /* SHAPE == FULL_WAKE */
 #else
   h_s = mapy0 + param0.tanh_eta_y;
 #endif /* ANY_WAKE */

 /* set domain height */
 h   = h_s + (REAL)param0.y_shape_size0[0] * parpt0->dy; 

 /* set mapping factors */
 mapf0 = param0.y_mapf[0];			/* gamma_0 */
 if (param0.y_mapf[2]<0) /* set to default */
 	param0.y_mapf[2] = 0.5 * (h + h_s);
 mapf2 = param0.y_mapf[2];			/* xi_0 */
 #ifdef ANY_WAKE
 mapf1 = param0.y_mapf[1];			/* gamma_1 */
 if (param0.y_mapf[3]<0) /* set to default */
 	param0.y_mapf[3] = 0.5 * h_s;
 mapf3 = param0.y_mapf[3];			/* xi_1 */
 #endif /* ANY_WAKE */

 /* set coefficients K_0 and K_1 */
 K_0 = (h - mapy0) / (h - h_s);
 /* It seems to work also for (mapf2 != h), it may be required to set
 K_0 = (h - mapy0) / ( (mapf2 - h_s)
 			* ( 1. - tanh( mapf0 * (mapf2 - h  ) )
			        /tanh( mapf0 * (mapf2 - h_s) ) ) );
 */
 #ifdef ANY_WAKE
 K_1 = mapy0 / ( mapf3 * ( 1. - tanh( mapf1 * (mapf3 - h_s) )
			       /tanh( mapf1 * mapf3 ) ) );
 #endif /* ANY_WAKE */

 printf ("y-stretching, outer:     %12.8f\n",	     mapy0);
 printf ("\t+ %12.8f * (%12.8f - %12.8f)\n",	     K_0, mapf2, h_s);
 printf ("\t* [1 - tanh(%12.8f *(%12.8f - y))\n",    mapf0, mapf2);
 printf ("\t/ tanh(%12.8f *(%12.8f - %12.8f))]\n\n", mapf0, mapf2, h_s);

 for (y=0; y<param0.y_shape_size[0]; y++) {
   eta			 = h_s + (REAL)y * parpt0->dy;
   miscdat(0).y_grid[y]  = mapy0 
   	+ K_0 * ( mapf2 - h_s) * ( 1. - tanh(mapf0*(mapf2-eta))
		   		       /tanh(mapf0*(mapf2-h_s)) );
   true_slope		 = K_0 * mapf0 * (mapf2 - h_s) /
   	( tanh(mapf0*(mapf2-h_s)) * square(cosh(mapf0*(mapf2-eta))) );
   miscdat(0).y_slope[y] = 1. / true_slope;
 }

 #if defined(ANY_WAKE)
 printf ("y-stretching, inner:     ");
 printf ("%12.8f * %12.8f\n", 			  K_1, mapf3);
 printf ("\t* [1 - tanh(%12.8f *(%12.8f - y))\n", mapf1, mapf3);
 printf ("\t/ tanh(%12.8f * %12.8f)]\n\n",	  mapf1, mapf3);

 for (y=y0; y<param0.y_shape_size[1]; y++) {
   eta			 = (REAL)y * parpt0->dy;
   miscdat(1).y_grid[y]  = K_1 * mapf3 * ( 1. - tanh(mapf1*(mapf3-eta))
					      /tanh(mapf1*mapf3) );
   true_slope		 = K_1 * mapf1 * mapf3 /
   	( tanh(mapf1*mapf3) * square(cosh(mapf1*(mapf3-eta))) );
   miscdat(1).y_slope[y] = 1. / true_slope;
 }
 #endif /* ANY_WAKE */

 #if (SHAPE == FULL_WAKE)
 printf("\n\nERROR:  tanh-grid for full wake not implemented!!!\n\n");
 shutdown();
 
 /* FOR FULL_WAKE, NEED TO ADJUST y_grid AND true_slope BELOW !!! */
 for (y=0; y<y0; y++) {
   miscdat(2).y_grid[y]  = 1.;
   true_slope		 = 1.;
   miscdat(2).y_slope[y] = 1. / true_slope;
 }
 
 for (y=0; y<param0.y_shape_size[2]; y++) {
   miscdat(2).y_grid[y]  = 1.;
   true_slope		 = 1.;
   miscdat(2).y_slope[y] = 1. / true_slope;
 }
 #endif /* SHAPE == FULL_WAKE */

}

/*----------------------------------------------------------------*/
#endif /* STRETCH_Y */

else { 
 /* Invalid stretching option */

 printf ("\n\nERROR: Invalid y-stretching option y_map_type = %d!\n\n",
    param0.y_map_type);
 shutdown();

}

/*----------------------------------------------------------------*/
	/* fill in other y-grids with multiple threads */

for (d=SHAPE_DOMAINS; d<params.n_domains; d++)
for (y=0; y<parpt[d]->y_size; y++) {
   miscdat(d).y_grid[y] = miscdat(d % SHAPE_DOMAINS).y_grid[y];
   miscdat(d).y_slope[y] = miscdat(d % SHAPE_DOMAINS).y_slope[y];
}

/*----------------------------------------------------------------*/
	/* put all the y-grids in one array */

#if (SHAPE == FLAT_PLATE)
miscdat(0).y_grid_all = miscdat(0).y_grid;
miscdat(0).y_slope_all = miscdat(0).y_slope;
#endif /* FLAT_PLATE */
#if (SHAPE == SYM_WAKE || SHAPE == BACK_STEP || SHAPE == LONE_CORNER)
for (y=0; y<param0.y_shape_size[1]; y++) {
   miscdat(0).y_grid_all[y] = miscdat(1).y_grid[y];
   miscdat(0).y_slope_all[y] = miscdat(1).y_slope[y];
}
for (y=0; y<param0.y_shape_size[0]; y++) {
   miscdat(0).y_grid_all[param0.y_shape_size[1]+y]
					= miscdat(0).y_grid[y];
   miscdat(0).y_slope_all[param0.y_shape_size[1]+y]
					= miscdat(0).y_slope[y];
}
#endif /* SHAPE == ??? */
#if (SHAPE == FULL_WAKE)
for (y=0; y<param0.y_shape_size[2]; y++) {
   miscdat(0).y_grid_all[y] = miscdat(2).y_grid[y];
   miscdat(0).y_slope_all[y] = miscdat(2).y_slope[y];
}
for (y=0; y<param0.y_shape_size[1]; y++) {
   miscdat(0).y_grid_all[param0.y_shape_size[2]+y]
					= miscdat(1).y_grid[y];
   miscdat(0).y_slope_all[param0.y_shape_size[2]+y]
					= miscdat(1).y_slope[y];
}
for (y=0; y<param0.y_shape_size[0]; y++) {
   miscdat(0).y_grid_all[param0.y_shape_size[2]+param0.y_shape_size[1]+y]
					= miscdat(0).y_grid[y];
   miscdat(0).y_slope_all[param0.y_shape_size[2]+param0.y_shape_size[1]+y]
					= miscdat(0).y_slope[y];

}
#endif /* SHAPE == FULL_WAKE */


/*----------------------------------------------------------------*/
/* CALCULATE INVERSE SLOPES (actually, back to normal) */
/*----------------------------------------------------------------*/
	/* used for the outlet body forcing calculation */
for (d=0; d<params.n_domains; d++) {
   for (x=0; x<parpt[d]->x_size; x++)
      miscdat(d).ix_slope[x] = 1./miscdat(d).x_slope[x];
   for (y=0; y<parpt[d]->y_size; y++)
      miscdat(d).iy_slope[y] = 1./miscdat(d).y_slope[y];
}


/*----------------------------------------------------------------*/
/* CALCULATE REVERSE (INVERTED) X-GRID */
/*----------------------------------------------------------------*/
	/* used for free stream boundary interpolations */
	/* x_num = (x-x_zero) * x_factor */

loop_d_all (d) {
	/* mapping offset (zero) */
   parpt[d]->re_x_z = miscdat(0).x_grid[0];
   if (parpt0->M > 1.) {
      parpt[d]->re_x_z -= (miscdat(0).y_grid[parpt0->y_size0]
			- miscdat(0).y_grid[parpt0->y_size2])
	* sqrt (1. - 1/square(parpt0->M) ) * 1.5;
   }

	/* mapping scaling factor */
   parpt[d]->re_x_f = (REAL)(4*parpt0->x_total-1)
	/ (miscdat(0).x_grid_all[parpt0->x_total-1] - parpt0->re_x_z);
}

	/* get interpolation values */
x = 0;
x0 = 0;
d = 0;
for (i=0; i<4*parpt0->x_total; i++) {
   xi = (REAL)i / parpt0->re_x_f + parpt0->re_x_z;
   if (xi <= miscdat(0).x_grid_all[0]) {
   } else {
      while (xi > miscdat(0).x_grid_all[x+1] && x < parpt0->x_total-1) {
         x++;
         x0++;
         if (x0 >= parpt[d]->x_size) {
            x0 = 0;
            d += SHAPE_DOMAINS;
         }
      }
      miscdat(0).re_x_grid[i] = (REAL)(x)
		+ (xi - miscdat(0).x_grid_all[x])
		  / (miscdat(0).x_grid_all[x+1]
					- miscdat(0).x_grid_all[x]);
   }
/*
printf ("%d %d %f %f %f %f %f\n", i, x, parpt0->re_x_z, parpt0->re_x_f,
	miscdat(0).x_grid[x], xi, miscdat(0).re_x_grid[i]);
*/
}

d = 0;
x0 = 0;
for (x=0; x<parpt0->x_total; x++) {
   if (x0 >= parpt[d]->x_size) {
      d += SHAPE_DOMAINS;
      x0 = 0;
   }
   miscdat(0).re_x_gridd[x] = d;
   miscdat(0).re_x_gridx[x] = x0;
   x0++;
}

	/* copy to all domains with boundaries */
for (d=1; d<params.n_domains; d++)
if (parpt[d]->bounds.top.free_st.on
		|| parpt[d]->bounds.bottom.free_st.on) {
   memcpy (miscdat(d).re_x_grid, miscdat(0).re_x_grid,
				4*parpt0->x_total * sizeof(REAL) );
   memcpy (miscdat(d).re_x_gridd, miscdat(0).re_x_gridd,
				parpt0->x_total * sizeof(int) );
   memcpy (miscdat(d).re_x_gridx, miscdat(0).re_x_gridx,
				parpt0->x_total * sizeof(int) );
}


/*----------------------------------------------------------------*/
/* fOR RANS/LES, CALCULATE DELTA ACROSS GRID */
/*----------------------------------------------------------------*/
#ifdef USE_LES

for (d=0; d<params.n_domains; d++)
for (x=0; x<parpt[d]->x_size; x++)
for (y=0; y<parpt[d]->y_size; y++) {
/* delta = arithmetic mean of grid sizes: sqrt[(dx^2+dy^2+dz^2)/3] */
#ifdef THREE_D
   miscdat(d).delta[INDEX2(d,x,y)] = sqrt ((
	  parpt0->dx * miscdat(d).ix_slope[x]
	* parpt0->dx * miscdat(d).ix_slope[x]
	+ parpt0->dy * miscdat(d).iy_slope[y]
	* parpt0->dy * miscdat(d).iy_slope[y]
	+ parpt0->dz * parpt0->dz)/3.);
#else /* TWO_D */
   miscdat(d).delta[INDEX2(d,x,y)] = sqrt ((
	  parpt0->dx * miscdat(d).ix_slope[x]
	* parpt0->dx * miscdat(d).ix_slope[x]
	+ parpt0->dy * miscdat(d).iy_slope[y]
	* parpt0->dy * miscdat(d).iy_slope[y])*0.5);
#endif /* THREE_D */
/* delta^2 = geometric mean of grid sizes: (dx*dy*dz)^2/3 
   actually need only delta -> ()^1/3 !!! */
/*    miscdat(d).delta[INDEX2(d,x,y)] = pow (
	  parpt0->dx * miscdat(d).ix_slope[x]
	* parpt0->dy * miscdat(d).iy_slope[y]
	* parpt0->dz,
			parpt0->i23); */
}

#endif /* USE_LES */


/*----------------------------------------------------------------*/
/* FIND ZERO LOCATION FOR PROFILE */
/*----------------------------------------------------------------*/

for (x=0; x<param0.x_shape_size[0] &&
	parpt0->body_length + miscdat(0).x_grid_all[x] < 0.; x++) {

   for (d=0; d<params.n_domains; d++)
      parpt[d]->x_zero = x;
}


/*----------------------------------------------------------------*/
/* SAVE GRID AND STRETCHING FILES */
/*----------------------------------------------------------------*/

#if (defined(STRETCH_X) || defined(STRETCH_Y))
					/* save only if stretched */

	/* put the grid, stretching in phys.full, consecutively */

max_size = max (parpt0->x_total, parpt0->y_total);
for (x=0; x<6*max_size; x++)
   params.phys.full[x] = 0.;

for (x=0; x<parpt0->x_total; x++) {
   params.phys.full[           x] = miscdat(0).x_grid_all[x];
   params.phys.full[3*max_size+x] = miscdat(0).x_slope_all[x];
}

for (y=0; y<parpt0->y_total; y++) {
   params.phys.full[  max_size+y] = miscdat(0).y_grid_all[y];
   params.phys.full[4*max_size+y] = miscdat(0).y_slope_all[y];
}

#ifdef THREE_D	/* save z-direction grid also */
for (x=0; x<max_size; x++) {
   params.phys.full[2*max_size+x] = (REAL)x * parpt0->dz;
   params.phys.full[5*max_size+x] = 1.;
}
#endif /* THREE_D */   

#ifdef IOSLIB_ON
/* output a grid-file only if ios-libraries are enabled,
   otherwise grid-info is included in eas3-files anyway */
open_output_file (GRID_FILENO,    param0.grid_file,   sflags[0],1,TRUE);
open_output_file (STRETCH_FILENO, param0.stretch_file,sflags[1],1,TRUE);

save_data (GRID_FILENO, FALSE);
save_data (STRETCH_FILENO, FALSE);

close_file (GRID_FILENO);
close_file (STRETCH_FILENO);
#endif /* IOSLIB_ON */

#endif /* defined(STRETCH_X) || defined(STRETCH_Y) */

return TRUE;
}



/* -|----|----|----|----|----|----|----|----|----|----|----|----|----|*/
/*
Given an x,y-location, return the corresponding domain and x,y-points.
*/

void xy_mid (REAL x_loc, REAL y_loc, int *d_ret, int *x_ret, int *y_ret)
{
int	d, x, y;
paramst	**parpt = prmpts[0]->params;

		/* first find y-location */
d = SHAPE_DOMAINS-1;
while (d >= 0 && y_loc > miscdat(d).y_grid[y=0]) {
   while (y < parpt[d]->y_size && y_loc > miscdat(d).y_grid[y])
      y++;
   if (y < parpt[d]->y_size)
      break;
   d--;
}
if (d < 0) {
   d++;
   y = parpt[d]->y_size;
}

		/* now find x-location */
while (d < params.n_domains && x_loc > miscdat(d).x_grid[x=0]) {
   while (x < parpt[d]->x_size && x_loc > miscdat(d).x_grid[x])
      x++;
   if (x < parpt[d]->x_size)
      break;
   d += SHAPE_DOMAINS;
}
if (d > params.n_domains) {
   d -= SHAPE_DOMAINS;
   x = parpt[d]->x_size;
}

*d_ret = d;
*x_ret = x;
*y_ret = y;

return;
}
