/*
 * MIT License
 *
 * Copyright © 2017 S5Lab
 * Created by Leonardo Parisi (leonardo.parisi[at]gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _H_ASTRO_PHOTO_H_
#define _H_ASTRO_PHOTO_H_

#include <cstdlib>
#include <cstdio>

#include <vector>

//****************************************************************************/
// namespace astro::photo
//****************************************************************************/
namespace astro {

/* ********************* Rotation Matrix ******************** *\
 *                                                            *
 * Define the rotation matrix in 3D given the rotation angles *
 * along each cartesian axis                                  *
 * Transpose as well is computed                              *
 *                                                            *
\* ********************************************************** */
void matrot(double tex, double tey, double tez, double rot[3][3], double trot[3][3])
{
  double rx[3][3],ry[3][3],rz[3][3];
  double rxt[3][3],ryt[3][3],rzt[3][3];
  double text, teyt, tezt;
  int i,j,k,l;
  
  text=-tex;
  teyt=-tey;
  tezt=-tez;
  
  // rotation matrix along x-axis
  rx[0][0]=1.0; rx[0][1]=0.0;      rx[0][2]=0.0;
  rx[1][0]=0.0; rx[1][1]=cos(tex); rx[1][2]=-sin(tex);
  rx[2][0]=0.0; rx[2][1]=sin(tex); rx[2][2]=cos(tex);  
  // rotation matrix along y-axis
  ry[0][0]=cos(tey);  ry[0][1]=0.0; ry[0][2]=sin(tey);
  ry[1][0]=0.0;       ry[1][1]=1.0; ry[1][2]=0.0;
  ry[2][0]=-sin(tey); ry[2][1]=0.0; ry[2][2]=cos(tey);
  // rotation matrix along z-axis
  rz[0][0]=cos(tez); rz[0][1]=-sin(tez); rz[0][2]=0.0;
  rz[1][0]=sin(tez); rz[1][1]=cos(tez);  rz[1][2]=0.0;
  rz[2][0]=0.0;      rz[2][1]=0.0;       rz[2][2]=1.0;
  
  // transpose matrix
  // rotation matrix along x-axis
  rxt[0][0]=1.0; rxt[0][1]=0.0;       rxt[0][2]=0.0;
  rxt[1][0]=0.0; rxt[1][1]=cos(text); rxt[1][2]=-sin(text);
  rxt[2][0]=0.0; rxt[2][1]=sin(text); rxt[2][2]=cos(text);
  // rotation matrix along y-axis
  ryt[0][0]=cos(teyt);  ryt[0][1]=0.0; ryt[0][2]=sin(teyt);
  ryt[1][0]=0.0;        ryt[1][1]=1.0; ryt[1][2]=0.0;
  ryt[2][0]=-sin(teyt); ryt[2][1]=0.0; ryt[2][2]=cos(teyt);
  // rotation matrix along z-axis
  rzt[0][0]=cos(tezt); rzt[0][1]=-sin(tezt); rzt[0][2]=0.0;
  rzt[1][0]=sin(tezt); rzt[1][1]=cos(tezt);  rzt[1][2]=0.0;
  rzt[2][0]=0.0;       rzt[2][1]=0.0;        rzt[2][2]=1.0;
  
  // initiate
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      rot[i][j]=0.0;
      trot[i][j]=0.0;
    }
  }
  
  // loop
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	for(int l=0; l<3; l++){
	  rot[i][l]=rot[i][l]+rz[i][j]*ry[j][k]*rx[k][l];
	  trot[i][l]=trot[i][l]+rxt[i][j]*ryt[j][k]*rzt[k][l];
	}
      }
    }
  }
  
}


/* ********************** AZEL_XYZ ****************** * 
 *                                                    *
 * Transformation from topocentric azimuthal frame of *
 * reference to cartesian SEZ coordinates.            *
 * North is pointing along the x-axis                 *
 *                                                    *
 * ************************************************** */
void azel_xyz(double& f, double& az, double& el, edirection direct, double x[3])
{

  if (direct == eTo)
    {
      x[0] = f*cos(el)*cos(az);
      x[1] = -f*cos(el)*sin(az);
      x[2] = f*sin(el);
    }

  if (direct == eFrom)
    {
      f = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      if (x[0]==0.0)
        {
          az=0.0;
        }
      else
        {
          az = atan2(-x[1],x[0]);;
          astUtils::rebox(az);
        }
      if (f==0.0)
        {
          el = 0.0;
        }
      else
        {
          el = asin(x[2]/f);
          astUtils::rebox(el);
        }
    }
  
}

/* ********************** RADEC_XYZ ****************** * 
 *                                                     * 
 * Transformation from topocentric equatorial frame of *
 * reference to cartesian coordinates                  *
 *                                                     *
 * *************************************************** */
void radec_xyz(double& f, double& ra, double& dec, edirection direct, double x[3])
{

  if (direct == eTo)
    {
      x[0] = f*cos(dec)*cos(ra);
      x[1] = f*cos(dec)*sin(ra);
      x[2] = f*sin(dec);
    }

  if (direct == eFrom)
    {
      f = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      if (x[0]==0.0)
        {
          ra=0.0;
        }
      else
        {
          ra = atan2(x[1],x[0]);;
          astUtils::rebox(ra);
        }
      if (f==0.0)
        {
          dec = 0.0;
        }
      else
        {
          dec = asin(x[2]/f);
          astUtils::rebox(dec);
        }
    }

}

    
  /* *********************** Projection ************************* *	\
   *                                                              *
   * Project a unit vector r on the sensor, by knowing the vector *
   * of the center and computing the distance from it,            *
   * hence the modulus ap=f/cos(d)=f/rp.rc                        *
   *                                                              *
   \* ************************************************************ */
  void project(double xc[3], double r[3], double x[3])
{
  double rc[3];
  double prod, ar, rr;
  
  rr=sqrt(xc[0]*xc[0]+xc[1]*xc[1]+xc[2]*xc[2]);
  for(int i=0; i<3; i++){
    rc[i] = xc[i]/rr;
  }  
  prod=0.0;
  for(int i=0; i<3; i++){
    prod = prod + rc[i]*r[i];
  }
  ar=rr/prod;
  for(int i=0; i<3; i++){
    x[i] = ar*r[i];
  }
  
}
  
  
/* ******************* Rotate function ******************** *\
 *                                                          *
 * 3d rotation. It rotates the vector xp into the vector xr *
 * using the matrix mat 3x3.                                *
 *                                                          *
\* ******************************************************** */
void rotate3d(double mat[3][3], double xp[3], double xr[3])
{
  
  for(int i=0; i<3; i++){
    xr[i]=0.0;
  }
  for(int l=0; l<3; l++){
    for(int k=0; k<3; k++){
      xr[l] = xr[l] + mat[l][k]*xp[k];
    }
  }
  
}


/* ******************* Rotate function ******************** *\
 *                                                          *
 * 2D rotation. It rotates the vector xp into the vector xr *
 * using the matrix mat 2x2.                                *
 *                                                          *
\* ******************************************************** */
void rotate2d(double angle, double x[2], double xr[2])
{
  
  double A11, A12, A21, A22;
  
  A11=cos(angle);
  A12=-sin(angle);
  A21=sin(angle);
  A22=cos(angle);
  xr[0] = A11*x[0] + A12*x[1];
  xr[1] = A21*x[0] + A22*x[1];
  
}  
  
/* ******************* Orientation ******************* *\
 *                                                     *
 * Compute the angle between the axis of the photo and *
 * the equatorial axis. Orientation angle is defined   *
 * as the clockwise East rotation from North.          *
 *                                                     *
 * *************************************************** */
void orientation(double xc[3], double latobs, double lst, double trot[3][3], double& angle1, double& angle2)
{
  
  const double pi=2.0*asin(1.0);
  const double rad2deg = 180.0 / pi;
  double x[3], x0[3], rx[3];
  double y[3], y0[3], ry[3];
  double f, rr, delta, az0, el0, ra0, dec0, dra, ddec;
  double rax, decx, ray, decy, azx, elx, azy, ely;
  
  // center of image
  astro::azel_xyz(f, az0, el0, eFrom, xc);
  
  // ra/dec center
  astIOD::radec_azel(ra0, dec0, lst, latobs, eFrom, az0, el0);  

  // axis in ra/dec
  delta=0.00001;
  dra=delta*pi/180.; ddec=delta*pi/180.;  
  rax=ra0+dra;
  decx=dec0;
  astIOD::radec_azel(rax, decx, lst, latobs, eTo, azx, elx);
  astUtils::rebox(azx);
  ray=ra0;
  decy=dec0+ddec;
  astIOD::radec_azel(ray, decy, lst, latobs, eTo, azy, ely);
  astUtils::rebox(azy);
  
  // projection of ra axis
  rr=1.0;
  astro::azel_xyz(rr, azx, elx, eTo, rx);
  project(xc, rx, x);
  // projection of dec axis
  rr=1.0;
  astro::azel_xyz(rr, azy, ely, eTo, ry);
  astro::project(xc, ry, y);
  // back rotation
  astro::rotate3d(trot, x, x0);
  astro::rotate3d(trot, y, y0);
  
  // compute angle
  angle1=-atan2(x0[2],x0[1]);
  angle2=-atan2(-y0[1],y0[2]);
  astUtils::rebox(angle1);
  astUtils::rebox(angle2);
  
}


/* ******************* Orientation ******************* *\
 *                                                     *
 * Compute the angle between the axis of the photo and *
 * the equatorial axis. Orientation angle is defined   *
 * as the clockwise East rotation from North.          *
 *                                                     *
 * *************************************************** */
void orientation_j2k(double f, double az0, double el0, double jd, double latobs, double lonobs,
		     double& ra0, double& dec0, double& angle1, double& angle2)
{
  
  const double pi=2.0*asin(1.0);
  const double rad2deg = 180.0 / pi;
  double x[3], x0[3], rx[3], xc[3];
  double y[3], y0[3], ry[3];
  double rr, delta, ra1, dec1, dra, ddec;
  double deltara, deltadec, t0, t1, lst, gst;
  double rax, decx, ray, decy, azx, elx, azy, ely;
  double rax1, decx1, ray1, decy1;
  double dra_abe, ddec_abe, dra_prec, ddec_prec, dra_nut, ddec_nut;
  double tex, tey, tez, rot[3][3], trot[3][3];

  astro::azel_xyz(f, az0, el0, eTo, xc);
  tex=0.0; tey=-el0; tez=-az0;
  astro::matrot(tex, tey, tez, rot, trot);
  
  // time
  t0 = (jd - 2451545.)/36525.;
  t1 = (2451545. - jd)/36525.;
  astTime::lstime(lonobs, jd, lst, gst);
  
  // ra/dec center
  astIOD::radec_azel(ra1, dec1, lst, latobs, eFrom, az0, el0);
  astUtils::rebox(ra1);
  astro::precession2(ra1, dec1, t0, t1, ra0, dec0);
  astUtils::rebox(ra0);
  dra_prec = ra1-ra0; ddec_prec = dec1-dec0;
  astro::nutation_radec(ra0, dec0, t0, deltara, deltadec);
  dra_nut = deltara; ddec_nut = deltadec;
  astro::aberration(ra0, dec0, t0, deltara, deltadec);
  dra_abe = deltara; ddec_abe = deltadec;
  ra0 = ra1 - dra_prec - dra_nut - dra_abe;
  dec0 = dec1 - ddec_prec - ddec_nut - ddec_abe;
  
  // ra axis
  delta=0.025;
  dra=delta*pi/180.; ddec=delta*pi/180.;
  rax=ra0+dra;
  decx=dec0;
  astro::precession(rax, decx, t0, rax1, decx1);
  astUtils::rebox(rax1);
  dra_prec = rax1-rax; ddec_prec = decx1-decx;
  astro::nutation_radec(rax, decx, t0, deltara, deltadec);
  dra_nut = deltara; ddec_nut = deltadec;
  astro::aberration(rax, decx, t0, deltara, deltadec);
  dra_abe = deltara; ddec_abe = deltadec;
  rax1 = rax + dra_prec + dra_nut + dra_abe;
  decx1 = decx + ddec_prec + ddec_nut + ddec_abe;
  astIOD::radec_azel(rax1, decx1, lst, latobs, eTo, azx, elx);
  astUtils::rebox(azx);
  // dec axis
  ray=ra0;
  decy=dec0+ddec;
  astro::precession(ray, decy, t0, ray1, decy1);
  astUtils::rebox(ray1);
  dra_prec = ray1-ray; ddec_prec = decy1-decy;
  astro::nutation_radec(ray, decy, t0, deltara, deltadec);
  dra_nut = deltara; ddec_nut = deltadec;
  astro::aberration(ray, decy, t0, deltara, deltadec);
  dra_abe = deltara; ddec_abe = deltadec;
  ray1 = ray + dra_prec + dra_nut + dra_abe;
  decy1 = decy + ddec_prec + ddec_nut + ddec_abe;
  astIOD::radec_azel(ray1, decy1, lst, latobs, eTo, azy, ely);
  astUtils::rebox(azy);
  
  // projection of ra axis
  rr=1.0;
  astro::azel_xyz(rr, azx, elx, eTo, rx);
  project(xc, rx, x);
  // projection of dec axis
  rr=1.0;
  astro::azel_xyz(rr, azy, ely, eTo, ry);
  astro::project(xc, ry, y);
  // back rotation
  rotate3d(trot, x, x0);
  rotate3d(trot, y, y0);
  
  // compute angle
  angle1=-atan2(x0[2],x0[1]);
  angle2=-atan2(-y0[1],y0[2]);
  astUtils::rebox(angle1);
  astUtils::rebox(angle2);
  
}

/* ******************* Orientation ******************* *\ 
 *                                                     *
 * Compute the angle between the axis of the photo and *
 * the equatorial axis. Orientation angle is defined   *
 * as the clockwise East rotation from North.          *
 *                                                     *
 * *************************************************** */
void orientation_radec(double xc[3], double latobs, double lst, double trot[3][3],
                 double& angle1, double& angle2)
{

  const double pi=2.0*asin(1.0);
  const double rad2deg = 180.0 / pi;
  double x[3], x0[3], rx[3];
  double y[3], y0[3], ry[3];
  double f, az0, el0, ra0, dec0, dra, ddec;
  double rax, decx, ray, decy;
  double azx, elx, azy, ely, rr;
  
  // center of image
  astro::radec_xyz(f, ra0, dec0, eFrom, xc);
  
  // axis in ra/dec  
  dra=0.001*pi/180.; ddec=0.001*pi/180.;  
  rax=ra0+dra;
  decx=dec0;
  ray=ra0;
  decy=dec0+ddec;
  
  // projection of ra axis
  rr=1.0;
  astro::radec_xyz(rr, rax, decx, eTo, rx);
  astro::project(xc, rx, x);
  // projection of dec axis
  rr=1.0;
  astro::radec_xyz(rr, ray, decy, eTo, ry);
  astro::project(xc, ry, y);
  // back rotation
  astro::rotate3d(trot, x, x0);
  astro::rotate3d(trot, y, y0);
  
  // compute angle
  angle1=-atan2(x0[2],x0[1]);
  angle2=-atan2(-y0[1],y0[2]);
  astUtils::rebox(angle1);
  astUtils::rebox(angle2);

}

/* ***************** Orientation_Theo **************** *\
 *                                                     *
 * Theoretical prediction of the orientation angle     *
 *                                                     *
\* *************************************************** */
void orientation_theo(double az, double el, double latobs, double tex,
		      double& gamma)
{
  double x, y;
  
  y = cos(latobs)*sin(az);
  x = cos(el)*sin(latobs)-cos(az)*cos(latobs)*sin(el);
  gamma = tex + atan2(y,x);
  astUtils::rebox(gamma);

}

/* ************* Projective matrix (RT) ************* *\
 *                                                    *
 * Projective matrix in ECEF                          *
 *                                                    *
\* ************************************************** */

void matRT(double lon, double lat, double r_earth[3], double mat[3][4])
{
  
  double t[3], rot[3][3];
  int i,j,k,l;
  
  // compute rotation matrix R
  rot[0][0]=sin(lon);          rot[0][1]=-cos(lon);         rot[0][2]=0.0;
  rot[1][0]=sin(lat)*cos(lon); rot[1][1]=sin(lat)*sin(lon); rot[1][2]=-cos(lat);
  rot[2][0]=cos(lat)*cos(lon); rot[2][1]=cos(lat)*sin(lon); rot[2][2]=sin(lat);
  
  // compute traslation vector T
  for(i=0; i<3; i++){
    t[i]=0.0;    
  }
  for(l=0; l<3; l++){
    for(k=0; k<3; k++){
      t[l] = t[l] + rot[l][k]*(-r_earth[k]);      
    }
  }
  
  // compose RT matrix 3x4
  for(l=0; l<3; l++){
    for(k=0; k<3; k++){
      mat[l][k] = rot[l][k];
    }
  }
  mat[0][3]=t[0]; mat[1][3]=t[1]; mat[2][3]=t[2];
  
}


} /* namespace astro */

#endif /* _H_ASTRO_PHOTO_H_ */
