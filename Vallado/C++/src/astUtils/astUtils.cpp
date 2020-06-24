#include "astUtils.h"

namespace astUtils {
  
  /* ------------------------------------------------------------------------------
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
  bool sight(double r1[3], double r2[3], char whichkind) {
      
    double eesqrd = 0.006694385000;     // eccentricity of earth sqrd
    double re      = 6378.137;           // km
      
    double tr1[3], tr2[3];
      
    // -------------------------  implementation   -----------------
    for(int i=0; i<3; ++i){
      tr1[i] = r1[i];
      tr2[i] = r2[i];
    }
      
    double magr1 = astMath::mag(tr1);
    double magr2 = astMath::mag(tr2);
    
    double temp = 1.0;;
      
    // --------------------- scale z component ---------------------
    if(whichkind == 'e')
      temp = 1.0 / sqrt(1.0 - eesqrd);
      
    tr1[2] = tr1[2] * temp;
    tr2[2] = tr2[2] * temp;
      
    double bsqrd = magr2*magr2;
    double asqrd = magr1*magr1;
    
    double adotb = astMath::dot(tr1,tr2);
      
    // ---------------------- find tmin ----------------------------
    double distsqrd = 0.0;
      
    double tmin = 0.0;
      
    if(!(fabs(asqrd + bsqrd - 2.0 *adotb) < 0.0001))
      tmin = (asqrd - adotb) / (asqrd + bsqrd - 2.0 * adotb);
      
      
    // ----------------------- check los ---------------------------
    if((tmin < 0.0 ) | (tmin > 1.0 )) {
        
      return true;
        
    } else {
      
      distsqrd = ((1.0 -tmin)*asqrd + adotb*tmin)/(re*re);
      
      if(distsqrd > 1.0)
        return true;
      else
        return false;
      
    }
      
      return false;
    
    }
        
        
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
    bool light(double r[3], double jdtdb, double jdtdbF, char whichkind) {
      
      double auer = 149597870.0 / 6378.1363;
      
      double rsun[3];
      double rtasc, decl;
      
      // -------------------------  implementation   -------------------------
      ast2Body::sun(jdtdb, jdtdbF, rsun, rtasc, decl);
      
      rsun[0] *= auer;
      rsun[1] *= auer;
      rsun[2] *= auer;
      
      // ------------ is the satellite in the shadow? ----------------
      return sight(rsun, r, whichkind);
      
    }
  
  
    /*****************************************************************************/
    // Obliquity of ecliptic
    //
    // Obliquety of ecliptic eps (in radians)
    // at time t (in Julian Centuries from J2000.0)
    // Taken from the KDE AStroLib
    /*****************************************************************************/
    double eps(double t) {

      double tp;
      
      tp = 23.43929111 - (46.815+(0.00059-0.001813*t)*t)*t/3600.0;
      tp = 1.74532925199e-2 * tp;
      
      return tp;
      
    }
  

    /*****************************************************************************/
    //  Ecliptic to Equatorial
    //
    // Convert position vector r1 from ecliptic into equatorial coordinates
    // at t (in Julian Centuries from J2000.0)
    /*****************************************************************************/
    void eclequ(double t, double * rIn, double * rOut) {
      
      std::vector< std::vector<double> > m;
      
      astMath::rot1mat(-eps(t), m); // m = xrot (-eps(t));
      
      astMath::matvecmult(m, rIn, rOut); // r2 = mxvct(m, r1);
      
    }

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
    void refraction(double el, double& appel) {
      
      const double rad2deg = 180.0 / M_PI; // rad2deg
      double eltmp;
      
      if(el>0.0) {      
        eltmp = el*rad2deg;        
        eltmp = (eltmp+10.3/(eltmp+5.11));       
        appel = el + 1.02/(60.*tan(eltmp/rad2deg))/rad2deg;        
        if(appel<0.0) appel = el;        
      }
      
    }

    /* ******************** Function Rebox ********************** *\
     *                                                            *
     * this function take an angle and makes it in [0:2*pi]       *
     *                                                            *
     \* ********************************************************** */
    void rebox(double & angle) {
      
      while(angle > 2.0*M_PI)
       angle-=2.*M_PI;
      
      while(angle < 0)
        angle+=2.*M_PI;
      
    }

}   // namespace

