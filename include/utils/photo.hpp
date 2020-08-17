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
  //int i,j,k,l;
  
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
  
  //const double rad2deg = 180.0 / M_PI;
  //double rr;
  
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
          az = atan2(-x[1],x[0]);
          astUtils::rebox(az);
        }
      if (f==0.0)
        {
          el = 0.0;
        }
      else
        {
          el = asin(x[2]/f);
          //astUtils::rebox(el);
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

    
  /* *********************** Projection ************************* *\
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
  
/* *********************** Orientation ******************** *\
 *                                                          *
 * Compute the image center in ra/dec JNOW and the angle of *
 * orientation between the axis of the photo and the axis   *
 * of the equatorial frame in JNOW.                         *
 * Orientation angle is defined as the clockwise rotation   *
 * from North towards East.                                 *
 *                                                          * 
 * ******************************************************** */
void orientation_azel_jnow(double f, double az0, double el0, 
			   double jd, double latobs, double lonobs,
			   double& ra0, double& dec0, double& angle1, double& angle2)
{
  
  const double rad2deg = 180.0 / M_PI;
  double x[3], x0[3], rx[3], xc[3];
  double y[3], y0[3], ry[3];
  double rr, delta, dra, ddec;
  double lst, gst;
  double rax, decx, ray, decy, azx, elx, azy, ely;
  double tex, tey, tez, rot[3][3], trot[3][3];
  
  astro::azel_xyz(f, az0, el0, eTo, xc);
  tex=0.0; tey=-el0; tez=-az0;
  astro::matrot(tex, tey, tez, rot, trot);
  
  // sidereal time
  astTime::lstime(lonobs, jd, lst, gst);
  
  // ra/dec center JNOW
  astIOD::radec_azel(ra0, dec0, lst, latobs, eFrom, az0, el0);
  astUtils::rebox(ra0);
  
  // ra axis
  delta=0.001; // 0.001
  dra=delta/rad2deg; ddec=delta/rad2deg;
  rax=ra0+dra;
  decx=dec0;
  astIOD::radec_azel(rax, decx, lst, latobs, eTo, azx, elx);
  astUtils::rebox(azx);
  
  // dec axis
  ray=ra0;
  decy=dec0+ddec;
  astIOD::radec_azel(ray, decy, lst, latobs, eTo, azy, ely);
  astUtils::rebox(azy);
  
  // projection of ra axis
  rr=1.0;
  astro::azel_xyz(rr, azx, elx, eTo, rx);
  astro::project(xc, rx, x);
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


/* ******************* Orientation_J2K ******************* *\
 *                                                         *
 * Compute the image center in ra/dec J2K and the angle of *
 * orientation between the axis of the photo and the axis  *
 * of the equatorial frame in J2K.                         *
 * Orientation angle is defined as the clockwise rotation  *
 * from North towards East.                                *
 *                                                         *
 * ******************************************************* */
void orientation_azel_j2k(double f, double az0, double el0,
		     double jd, double latobs, double lonobs, double altobs,
		     double& ra0, double& dec0, double& angle1, double& angle2)
{
  
  const double rad2deg = 180.0 / M_PI;
  double x[3], x0[3], rx[3], xc[3];
  double y[3], y0[3], ry[3];
  double rr, ff, delta, dra, ddec; // ra1, dec1,
  double t0, t1, lst, gst; //deltara, deltadec,
  double rax, decx, ray, decy, azx, elx, azy, ely;
  //double rax1, decx1, ray1, decy1;
  //double dra_abe, ddec_abe, dra_prec, ddec_prec, dra_nut, ddec_nut;
  double tex, tey, tez, rot[3][3], trot[3][3];

  astro::azel_xyz(f, az0, el0, eTo, xc);
  tex=0.0; tey=-el0; tez=-az0;
  astro::matrot(tex, tey, tez, rot, trot);
  
  // time
  t0 = (jd - 2451545.)/36525.;
  t1 = (2451545. - jd)/36525.;
  astTime::lstime(lonobs, jd, lst, gst);
  
  // ra/dec center
  astro::j2k_azel(ff, ra0, dec0, latobs, lonobs, altobs, jd, eFrom, f, az0, el0);
  
  // increments
  delta=0.001;
  dra=delta/rad2deg; ddec=delta/rad2deg;

  // ra axis
  rax=ra0+dra;
  decx=dec0;
  astro::j2k_azel(f, rax, decx, latobs, lonobs, altobs, jd, eTo, ff, azx, elx);
  
  // dec axis
  ray=ra0;
  decy=dec0+ddec;
  astro::j2k_azel(f, ray, decy, latobs, lonobs, altobs, jd, eTo, ff, azy, ely);
  
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
void orientation_radec_j2k(double f, double ra0, double dec0, double& angle1, double& angle2)
{

  const double rad2deg = 180.0 / M_PI;
  double xc[3];
  double x[3], x0[3], rx[3];
  double y[3], y0[3], ry[3];
  double rr, dra, ddec, rax, decx, ray, decy;
  double tex, tey, tez, rot[3][3], trot[3][3];

  // center of image
  astro::radec_xyz(f, ra0, dec0, eTo, xc);
  tex=0.0; tey=-dec0; tez=ra0;
  astro::matrot(tex, tey, tez, rot, trot);
    
  // axis in ra/dec  
  dra=0.001/rad2deg; ddec=0.001/rad2deg;
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

/* ******************* Orientation ******************* *\ 
 *                                                     *
 * Compute the angle between the axis of the photo and *
 * the equatorial axis. Orientation angle is defined   *
 * as the clockwise East rotation from North.          *
 *                                                     *
 * *************************************************** */
void orientation_radec_jnow(double f, double ra0, double dec0, double jd, double latobs, double lonobs, double altobs,
		       double& ra1, double& dec1, double& angle1, double& angle2)
{

  const double rad2deg = 180.0 / M_PI;
  double xc[3];
  double x[3], x0[3], rx[3];
  double y[3], y0[3], ry[3];
  double rr, ff, lst, gst, az0, el0, dra, ddec, rax, decx, ray, decy;
  double rax1, decx1, ray1, decy1, azx, elx, azy, ely;
  double tex, tey, tez, rot[3][3], trot[3][3];
  
  // center of image
  astro::radec_xyz(f, ra0, dec0, eTo, xc);
  tex=0.0; tey=-dec0; tez=ra0;
  astro::matrot(tex, tey, tez, rot, trot);

  // sidereal time
  astTime::lstime(lonobs, jd, lst, gst);

  // center in ra/dec jnow
  astro::j2k_azel(f, ra0, dec0, latobs, lonobs, altobs, jd, eTo, ff, az0, el0);
  astIOD::radec_azel(ra1, dec1, lst, latobs, eFrom, az0, el0);
  astUtils::rebox(ra1);
  
  // increments
  dra=0.001/rad2deg; ddec=0.001/rad2deg;
  
  // axis in ra  
  rax1=ra1+dra;
  decx1=dec1;
  astIOD::radec_azel(rax1, decx1, lst, latobs, eTo, azx, elx);
  astro::j2k_azel(ff, rax, decx, latobs, lonobs, altobs, jd, eFrom, f, azx, elx);
  
  // axis in dec
  ray1=ra1;
  decy1=dec1+ddec;
  astIOD::radec_azel(ray1, decy1, lst, latobs, eTo, azy, ely);
  astro::j2k_azel(ff, ray, decy, latobs, lonobs, altobs, jd, eFrom, f, azy, ely);

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
void orientation_theo(double az, double el, double latobs, double tex, double& gamma)
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
  int i,k,l;
  
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


/* ******************** AUTOCENTER  ******************** *\
 *                                                       *
 * this function extract the center of image from the    *
 * barycenter of the satellite track                     *
 *                                                       *
\* ***************************************************** */
void autocenter_azel(char satnum[128], char nameobs[128], double& az0, double& el0, double& appel0)
{
  
  // variables                                                             
  double az, appel, el, appra, appdec, rajnw, decjnw, raj2k, decj2k, alt, range;
  int    year, mon, day, hr, min;
  double sec;
  FILE   *infile;
  char   filename[128];
  int    iframe;
  
  sprintf(filename,"viewsite.%s.%s",satnum,nameobs);
  if ((infile=fopen(filename,"r"))==NULL){
    printf(" ERROR: file %s does not exist \n",filename);
    exit(9);
  }
  
  iframe=0;
  az0=0.0; el0=0.0; appel0=0.0;
  while (feof(infile)==0){
    do
      {
        iframe=iframe+1;
        fscanf(infile," %i %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
               &year, &mon, &day, &hr, &min, &sec,
               &az, &appel, &el, &appra, &appdec, &rajnw, &decjnw, &raj2k, &decj2k, &alt, &range);
        if (iframe==1){
          az0=az;
          el0=el;
          appel0=appel;
        }	
      } while ((feof(infile) == 0));
  }// end of infile
  fclose(infile);
  
  if (fabs(az0-az)>180.) az0=az0+360.;
  
  // middle point
  az0=0.5*(az0+az);
  el0=0.5*(el0+el);
  appel0=0.5*(appel0+appel);
  
#ifdef NORIFRA
  appel0=el0;
#endif
  
}

/* ******************** AUTOCENTER  ******************** *\                                         
 *                                                       *   
 * this function extract the center of image from the    *
 * barycenter of the satellite track                     *                                             
 *                                                       * 
\* ***************************************************** */
void autocenter_radec(char satnum[128], char nameobs[128], double& ra0, double& dec0)
{
  
  double az, appel, el, appra, appdec, rajnw, decjnw, raj2k, decj2k, alt, range;
  int    year, mon, day, hr, min;
  double sec;
  FILE   *infile;
  char   filename[128];
  int    iframe;

  sprintf(filename,"viewsite.%s.%s",satnum,nameobs);
  if ((infile=fopen(filename,"r"))==NULL){
    printf(" ERROR: file %s does not exist \n",filename);
    exit(9);
  }

  iframe=0;
  ra0=0.0; dec0=0.0;
  while (feof(infile)==0){
    do
      {
        iframe=iframe+1;
        fscanf(infile," %i %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
               &year, &mon, &day, &hr, &min, &sec,
               &az, &appel, &el, &appra, &appdec, &rajnw, &decjnw, &raj2k, &decj2k, &alt, &range);
        if (iframe==1){
          ra0=raj2k;
          dec0=decj2k;
        }
      } while ((feof(infile) == 0));
  }// end of file
  fclose(infile);
  
  if (fabs(ra0-raj2k)>180.) ra0=ra0+360.;
  
  // middle point
  ra0=0.5*(ra0+raj2k);
  dec0=0.5*(dec0+decj2k);

}

// positioning radec/azel
void positioning(double te1, double te2, double xc[3], double trot[3][3], double angle,
		 double xy[2], double xyr[2], double sizex, double sizey, double pxsize)
{
  
  double rr, az, el, ra, dec;
  double rp[3], xp[3], xp0[3], xq[2], xq0[2];
  
  rr=1.0;
  
#ifdef RADEC
  ra=te1; dec=te2;
  astro::radec_xyz(rr, ra, dec, eTo, rp);
#endif
  
#ifdef AZEL
  az=te1; el=te2;
  astro::azel_xyz(rr, az, el, eTo, rp);
#endif
  
  // project and back-rotate
  astro::project(xc, rp, xp);
  astro::rotate3d(trot, xp, xp0);
  
  // invert axis
  xq0[0]=-xp0[1];
  xq0[1]=-xp0[2];
  
  // positions in pixel
  xy[0]=(xq0[0]+0.5*sizex)/pxsize;
  xy[1]=(xq0[1]+0.5*sizey)/pxsize;
  
  // rotate using orientation angle
  astro::rotate2d(angle, xq0, xq);
  
  // rotated positions in pixel
  xyr[0]=(xq[0]+0.5*sizex)/pxsize;
  xyr[1]=(xq[1]+0.5*sizey)/pxsize;  
  
}


// write photosat output 
void outsat(int iframe, double xy[2], double rasat, double decsat, double rasat_eci, double decsat_eci,
	    double rsat_ecef[3], double rsat_eci[3], int year, int mon, int day, int hr, int min, double sec,
	    double sizex, double sizey, double pxsize, FILE* photofile)
{

  const double rad2deg = 180.0 / M_PI;
  
  if ((xy[0] > 0.0) && (xy[0] < sizex/pxsize)){
    if ((xy[1] > 0.0) && (xy[1] < sizey/pxsize)){
      fprintf(photofile," %03i %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %i %02i %02i %02i %02i %04.2f \n",
	      iframe, xy[0], xy[1], rasat*rad2deg, decsat*rad2deg, rasat_eci*rad2deg, decsat_eci*rad2deg,
	      rsat_ecef[0], rsat_ecef[1], rsat_ecef[2], rsat_eci[0], rsat_eci[1], rsat_eci[2],
	      year, mon, day, hr, min, sec);
    }
  }
  

}

// write stars in ouput file
void outstar(int iframe, int istar, double xy[2], double rastar, double decstar,
	     double magstar, double rac, double decc, double angle,
	     int year, int mon, int day, int hr, int min, double sec,
	     double sizex, double sizey, double pxsize, FILE* photofile)
{

  const double rad2deg = 180.0 / M_PI;
    
  // FRAME ID_STAR X Y RA DEC MAG RA0 DEC0 ANGLE DATE
  if ((xy[0] > 0.0) && (xy[0] < sizex/pxsize)){
    if ((xy[1] > 0.0) && (xy[1] < sizey/pxsize)){
      fprintf(photofile," %03i %i %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %i %02i %02i %02i %02i %04.2f \n",
	      iframe, istar, xy[0], xy[1], rastar*rad2deg, decstar*rad2deg, magstar, rac*rad2deg, decc*rad2deg, angle*rad2deg,
	      year, mon, day, hr, min, sec);
    }
  }
  
}


} /* namespace astro */

#endif /* _H_ASTRO_PHOTO_H_ */
