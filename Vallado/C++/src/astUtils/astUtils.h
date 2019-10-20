#ifndef _astUtils_h_
#define _astUtils_h_


#include "ast2Body.h"
#include "astMath.h"

#pragma once


namespace astUtils {

  /*------------------------------------------------------------------------------
  %
  %                           function sight
  %
  %  this function takes the position vectors of two satellites and determines
  %    if there is line-of-sight between the two satellites.  an oblate earth
  %    with radius of 1 er is assumed.  the process forms the equation of
  %    a line between the two vectors.  differentiating and setting to zero finds
  %    the minimum value, and when plugged back into the original line equation,
  %    gives the minimum distance.  the parameter tmin is allowed to range from
  %    0.0  to 1.0 .  scale the k-component to account for oblate earth because it's
  %    the only qunatity that changes.
  %
  %  author        : david vallado                  719-573-2600   31 oct 2003
  %
  %  revisions
  %                -
  %
  %  inputs          description                    range / units
  %    r1          - position vector of the 1st sat km
  %    r2          - position vector of the 2nd sat km
  %    whichkind   - spherical or ellipsoidal earth 's', 'e'*default
  %
  %  outputs       :
  %    los         - line of sight                  true,false
  %
  %  locals        :
  %    tr1         - scaled r1 vector
  %    tr2         - scaled r2 vector
  %    adotb       - dot product of a dot b
  %    tmin        - minimum value of t from a to b
  %    distsqrd    - min distance squared to earth
  %    asqrd       - magnitude of a squared
  %    bsqrd       - magnitude of b squared
  %
  %  coupling:
  %
  %  references    :
  %    vallado       2001, 291-295, alg 35, ex 5-3
  %
  % [los] = sight ( r1,r2, whichkind );
  % ------------------------------------------------------------------------------*/
  bool sight(double r1[3], double r2[3], char whichkind='e');



  /* ------------------------------------------------------------------------------
  %
  %                           function light
  %
  %  this function determines if a spacecraft is sunlit or in the dark at a
  %    particular time.  an oblate earth and cylindrical shadow is assumed.
  %
  %  author        : david vallado                  719-573-2600   27 may 2002
  %
  %  revisions
  %                -
  %
  %  inputs          description                    range / units
  %    r           - position vector of sat         er
  %    jd          - julian date at desired time    days from 4713 bc
  %    whichkind   - spherical or ellipsoidal earth 's', 'e'*default
  %
  %  outputs       :
  %    vis         - visibility flag                true, false
  %
  %  locals        :
  %    rtasc       - suns right ascension           rad
  %    decl        - suns declination               rad
  %    rsun        - sun vector                     au
  %    auer        - conversion from au to er
  %
  %  coupling      :
  %    sun         - position vector of sun
  %    lncom1      - multiple a vector by a constant
  %    sight       - does line-of-sight exist beteen vectors
  %
  %  references    :
  %    vallado       2001, 291-295, alg 35, ex 5-6
  %
  % [lit] = light ( r, jd, whichkind );
  % ------------------------------------------------------------------------------*/
  bool light(double r[3], double jdtdb, double jdtdbF, char whichkind = 'e');

  
  /*****************************************************************************/
  // Obliquity of ecliptic
  //
  // Obliquety of ecliptic eps (in radians)
  // at time t (in Julian Centuries from J2000.0)
  // Taken from the KDE AStroLib
  /*****************************************************************************/
  double eps(double t);
  
  
  /*****************************************************************************/
  //  Ecliptic to Equatorial
  //
  // Convert position vector r1 from ecliptic into equatorial coordinates
  // at t (in Julian Centuries from J2000.0)
  /*****************************************************************************/
  void eclequ(double t, double * rIn, double * rOut);

  /* ****************** Function refraction ******************* *\
   *                                                            *
   * this function compute the correction to true elevation     *
   * due to light refraction                                    *
   * it should be used only for positive elevation              *
   * input:                                                     *
   *       el rad                                               *
   * output:                                                    *
   *     appel rad                                              *
   * we use formula 16.4 of Astronomical Algorithms p.106       *
   * we do not add the correction to recover no refraction for  *
   * elevation=90deg to be consistent with previsat             *
   *                                                            *
   \* ********************************************************** */
  void refraction(double el, double * appel);
  
  
  /* ******************** Function Rebox ********************** *\
   *                                                            *
   * this function take an angle and makes it in [0:2*pi]       *
   *                                                            *
   \* ********************************************************** */
  void rebox(double & angle);

  
}  // namespace

#endif
