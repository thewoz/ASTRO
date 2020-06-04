 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "astro/astro.hpp"

using namespace std;

/* ******************** AUTOCENTER  ******************** *	\
 *                                                       *
 * this function extract the center of image from the    *
 * barycenter of the satellite track                     *
 *                                                       *
\* ***************************************************** */
void autocenter(char satnum[128], char nameobs[128], int rifra, double& az0, double& el0)
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
  az0=0.0; el0=0.0;
  while (feof(infile)==0){
    do
      {
	iframe=iframe+1;
	fscanf(infile," %i %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
	       &year, &mon, &day, &hr, &min, &sec,
	       &az, &appel, &el, &appra, &appdec, &rajnw, &decjnw, &raj2k, &decj2k, &alt, &range);
	if (rifra != 0){
	  el = appel;
	}
	if (iframe==1){
	  az0=az;
	  el0=el;
	}
      } while ((feof(infile) == 0));
  }// end of infile
  fclose(infile);
  
  if (fabs(az0-az)>180.) az0=az0+360.;
  
  // middle point
  az0=0.5*(az0+az);
  el0=0.5*(el0+el);
  
}



/* ********************************************************** *\
 *                                                            * 
 *                            M A I N                         *            
 *                                                            *
\* ********************************************************** */


int main()
{
  
  // universal constants  
  const double rad2deg = 180.0 / M_PI;
  const double rad2hrs = 12.0 / M_PI;
  const double pxsize = 0.0065;
    
  // sat infos
  char satnum[128], satname[128];

  //photo parameters
  double rr, rp[3], xp[3], xp0[3];
  double xr[2], xr0[2], xx0[3], xx[3];
  double xc1[3], xc2[3];
  double focal, sizex, sizey;
  int rifra;
  FILE *paramfile;
  
  // time variables
  int year, mon, day, hr, min, timezone;
  int iframe, istar;
  double sec, jd, jdfrac, t, t0, t1;
  
  // i/o files
  char filename[128];
  FILE *infile, *infile1, *infile2;
  FILE *matfile1, *matfile2;
  FILE *photofile1, *photofile2;
  FILE *photofile1r, *photofile2r;
  FILE *eceffile;

  // correction
  double dra_abe, ddec_abe, dra_prec, ddec_prec, dra_nut, ddec_nut;
  double deltara, deltadec;
  
  // obs1
  double azsat1, elsat1, appelsat1, altsat1, rangesat1;
  double rasat1, decsat1, rasat11, decsat11, apprasat1, appdecsat1;
  double rajnwsat1, decjnwsat1, raj2ksat1, decj2ksat1;
  double lonobs1, latobs1, altobs1, lst1, gst1;
  double azc1, elc1, rac1, decc1, rac11, decc11;
  double tex1, tey1, tez1, rot1[3][3], trot1[3][3];
  double angle1, gamma1, angle_ra1, angle_dec1;
  double robs1_ecef[3], vobs1_ecef[3];
  double robs1_eci[3], vobs1_eci[3];
  double rsat_ecef[3], vsat_ecef[3];
  double lon1, lat1, mrt1[3][4];
  char nameobs1[128];

  // obs2
  double azsat2, elsat2, appelsat2, altsat2, rangesat2;
  double rasat2, decsat2, rasat22, decsat22, apprasat2, appdecsat2;
  double rajnwsat2, decjnwsat2,	raj2ksat2, decj2ksat2;
  double lonobs2, latobs2, altobs2, lst2, gst2;
  double azc2, elc2, rac2, decc2, rac22, decc22;
  double tex2, tey2, tez2, rot2[3][3], trot2[3][3];
  double angle2, gamma2, angle_ra2, angle_dec2;
  double robs2_ecef[3],	vobs2_ecef[3];
  double robs2_eci[3], vobs2_eci[3];
  double lon2, lat2, mrt2[3][4];
  char nameobs2[128];
  
  
  /* ********************* PARAMETERS ******************** */
  // read satellite number, coordinates of sites,
  // focal, sensor sizes, rifraction index
  if ((paramfile=fopen("param.txt","r"))==NULL){
    printf(" ERROR: file param.txt does not exist \n");
    exit(9);
  }
  fscanf(paramfile," %s %s ", satname, satnum);
  fscanf(paramfile," %s %lf %lf %lf ", nameobs1, &lonobs1, &latobs1, &altobs1);  
  fscanf(paramfile," %s %lf %lf %lf ", nameobs2, &lonobs2, &latobs2, &altobs2);
  fscanf(paramfile," %lf ", &focal);
  fscanf(paramfile," %lf %lf ", &sizex, &sizey);
  fscanf(paramfile," %d ", &rifra);
  fclose(paramfile);
  
  //Compute image center as the barycenter of the trajectory
  autocenter(satnum, nameobs1, rifra, azc1, elc1);
  autocenter(satnum, nameobs2, rifra, azc2, elc2);
  
  // write logfile (before conversion)
  printf(" \n");
  printf(" ******************* PHOTOSAT ******************* \n");
  printf(" satellite name: %s \n", satname);
  printf(" satellite number: %s \n", satnum);
  printf(" observatory:  %s \n", nameobs1);
  printf(" lon, lat, alt = %lf %lf %lf \n", lonobs1, latobs1, altobs1);
  printf(" center in az,el  = %lf, %lf \n", azc1, elc1);
  printf(" ************************************************ \n");
  printf(" observatory:  %s \n", nameobs2);
  printf(" lon, lat, alt = %lf %lf %lf \n", lonobs2, latobs2, altobs2);
  printf(" center in az,el  = %lf, %lf \n", azc2, elc2);
  printf(" ************************************************ \n");
  printf(" focal length = %lf \n", focal);
  printf(" sensor sizes = %lf %lf \n", sizex, sizey);
  if (rifra == 0)
    {
      printf(" no refraction, rifra = %d \n", rifra);
    }
  else
    {
      printf(" with refraction, rifra = %d \n", rifra);
    }
  printf(" ************************************************ \n");
  
  // degree conversion in radians
  lonobs1 /= rad2deg; latobs1 /= rad2deg; altobs1 /= 1000.;
  lonobs2 /= rad2deg; latobs2 /= rad2deg; altobs2 /= 1000.;  
  azc1 /= rad2deg; elc1 /= rad2deg;
  azc2 /= rad2deg; elc2 /= rad2deg;

  // site position in ecef
  astIOD::site(latobs1, lonobs1, altobs1, robs1_ecef, vobs1_ecef);
  astIOD::site(latobs2, lonobs2, altobs2, robs2_ecef, vobs2_ecef);  
  /* *********************************************************** */
  
  
  /* ********************* SENSOR ROTATION ********************* */
  // define rotation matrix for the sensor
  // obs 1
  tex1=0.0; tey1=-elc1; tez1=-azc1;
  astro::matrot(tex1, tey1, tez1, rot1, trot1);
  // obs 2
  tex2=0.0; tey2=-elc2; tez2=-azc2;
  astro::matrot(tex2, tey2, tez2, rot2, trot2);
  // center of the image 
  astro::azel_xyz(focal, azc1, elc1, eTo, xc1);
  astro::azel_xyz(focal, azc2, elc2, eTo, xc2);  
  /* *********************************************************** */
  
  
  /* **************** OPEN INPUT/OUTPUT FILES ****************** */
  // read satellite coordinates in viewsite
  sprintf(filename,"viewsite.%s.%s",satnum,nameobs1);  
  if ((infile1=fopen(filename,"r"))==NULL){
    printf(" ERROR: infile1 does not exist \n");
    exit(9);
  }
  sprintf(filename,"viewsite.%s.%s",satnum,nameobs2);
  if ((infile2=fopen(filename,"r"))==NULL){
    printf(" ERROR:  does not exist \n");
    exit(9);
  }
  // read 3D ECEF position of the satellite
  sprintf(filename,"ecef.%s",satnum);
  if ((eceffile=fopen(filename,"r"))==NULL){
    printf(" ERROR:  does not exist \n");
    exit(9);
  }
  // open output files
  sprintf(filename,"mat.%s.%s.dat",satnum,nameobs1);
  matfile1=fopen(filename, "w");
  sprintf(filename,"photosat.%s.%s.dat",satnum,nameobs1);
  photofile1 = fopen(filename, "w");
  sprintf(filename,"photosat.rot.%s.%s.dat",satnum,nameobs1);
  photofile1r = fopen(filename, "w");
  //
  sprintf(filename,"mat.%s.%s.dat",satnum,nameobs2);
  matfile2=fopen(filename, "w");
  sprintf(filename,"photosat.%s.%s.dat",satnum,nameobs2);
  photofile2 = fopen(filename, "w");
  sprintf(filename,"photosat.rot.%s.%s.dat",satnum,nameobs2);
  photofile2r = fopen(filename, "w");
  /* *********************************************************** */

  /* ****************** PHOTOGRAPH ******************** */
  // loop on sat
  iframe=0;
  while (feof(infile1)==0){
    do
      {
	
	iframe=iframe+1;
	// reading obs1
	fscanf(infile1," %i %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
	       &year, &mon, &day, &hr, &min, &sec, &azsat1, &appelsat1, &elsat1, &apprasat1, &appdecsat1,
	       &rajnwsat1, &decjnwsat1, &raj2ksat1, &decj2ksat1, &altsat1, &rangesat1);
	if (rifra !=0 )
	  {
	    elsat1 = appelsat1;
	    rajnwsat1 = apprasat1;
	    decjnwsat1 = appdecsat1;
	  }
	
	// reading obs2
	fscanf(infile2," %i %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
               &year, &mon, &day, &hr, &min, &sec, &azsat2, &appelsat2, &elsat2, &apprasat2, &appdecsat2,
	       &rajnwsat2, &decjnwsat2, &raj2ksat2, &decj2ksat2, &altsat2, &rangesat2);
        if (rifra !=0 )
	  {
	    elsat2 = appelsat2;
	    rajnwsat2 = apprasat2;
	    decjnwsat2 = appdecsat2;
	  }
	
        // convert input data in radians
        azsat1 = azsat1/rad2deg; elsat1 = elsat1/rad2deg;
        rajnwsat1 = rajnwsat1/rad2deg; decjnwsat1 = decjnwsat1/rad2deg;
        raj2ksat1 = raj2ksat1/rad2deg; decj2ksat1 = decj2ksat1/rad2deg;
        azsat2 = azsat2/rad2deg; elsat2 = elsat2/rad2deg;
        rajnwsat2 = rajnwsat2/rad2deg; decjnwsat2 = decjnwsat2/rad2deg;
	raj2ksat2 = raj2ksat2/rad2deg; decj2ksat2 = decj2ksat2/rad2deg;

	// reading satellite ecef coordinates
	fscanf(eceffile," %i %i %i %i %i %lf %lf %lf %lf %lf %lf %lf",
               &year, &mon, &day, &hr, &min, &sec, &rsat_ecef[0], &rsat_ecef[1], &rsat_ecef[2], &vsat_ecef[0], &vsat_ecef[1], &vsat_ecef[2]);
	
	// set epoch of JNOW
	timezone=0;
	astTime::jday(year, mon, day, hr, min, sec, jd, jdfrac);
	printf(" Frame = %03i | Date = %i/%02i/%02i %02i:%02i:%04.2f \n", iframe, year, mon, day, hr, min, sec);
	t = (jd + jdfrac - 2451545.)/36525.;
        t0 = (jd + jdfrac - 2451545.)/36525.;
        t1 = (2451545. - jd - jdfrac)/36525.;
	
	
	/* ******************************************************************* */
	/* compute center in ra/dec and orientation angle */
	// obs1
	astTime::lstime(lonobs1, jd+jdfrac, lst1, gst1);
        astro::orientation_j2k(focal, azc1, elc1, jd+jdfrac, latobs1, lonobs1, rac11, decc11, angle_ra1, angle_dec1);
	astro::orientation_theo(azc1, elc1, latobs1, tex1, gamma1);
	// obs2
	astTime::lstime(lonobs2, jd+jdfrac, lst2, gst2);
        astro::orientation_j2k(focal, azc2, elc2, jd+jdfrac, latobs2, lonobs2, rac22, decc22, angle_ra2, angle_dec2);	
	astro::orientation_theo(azc2, elc2, latobs2, tex2, gamma2);
	/* ******************************************************************* */
	
	
	/* ******************************************************************* */
	// PROJECTIVE MATRIX
	// obs1
	lon1=rac1-gst1; lat1=decc1;
	astUtils::rebox(lon1);
	astro::matRT(lon1, lat1, robs1_ecef, mrt1);
        fprintf(matfile1," %16.12f %16.12f %16.12f %16.12f \n", mrt1[0][0], mrt1[0][1], mrt1[0][2], mrt1[0][3]);
        fprintf(matfile1," %16.12f %16.12f %16.12f %16.12f \n", mrt1[1][0], mrt1[1][1], mrt1[1][2], mrt1[1][3]);
        fprintf(matfile1," %16.12f %16.12f %16.12f %16.12f \n", mrt1[2][0], mrt1[2][1], mrt1[2][2], mrt1[2][3]);
	// obs2
	lon2=rac2-gst2; lat2=decc2;
	astUtils::rebox(lon2);
        astro::matRT(lon2, lat2, robs2_ecef, mrt2);
        fprintf(matfile2," %16.12f %16.12f %16.12f %16.12f \n", mrt2[0][0], mrt2[0][1], mrt2[0][2], mrt2[0][3]);
        fprintf(matfile2," %16.12f %16.12f %16.12f %16.12f \n", mrt2[1][0], mrt2[1][1], mrt2[1][2], mrt2[1][3]);
        fprintf(matfile2," %16.12f %16.12f %16.12f %16.12f \n", mrt2[2][0], mrt2[2][1], mrt2[2][2], mrt2[2][3]);
	/* ******************************************************************* */


	/* ******************************************************************* */
	// CORRECTIONS TO SATELLITE IN J2K
	// obs1
        astIOD::radec_azel(rasat1, decsat1, lst1, latobs1, eFrom, azsat1, elsat1);
	astUtils::rebox(rasat1);
        astro::precession2(rasat1, decsat1, t0, t1, rasat11, decsat11);
        astUtils::rebox(rasat11);
        dra_prec = rasat1-rasat11; ddec_prec = decsat1-decsat11;
        astro::nutation_radec(rasat11, decsat11, t0, deltara, deltadec);
        dra_nut = deltara; ddec_nut = deltadec;
        astro::aberration(rasat11, decsat11, t0, deltara, deltadec);
        dra_abe = deltara; ddec_abe = deltadec;
        rasat11 = rasat1 - dra_prec - dra_nut - dra_abe;
        decsat11 = decsat1 - ddec_prec - ddec_nut - ddec_abe;
	// obs2
	astIOD::radec_azel(rasat2, decsat2, lst2, latobs2, eFrom, azsat2, elsat2);
	astUtils::rebox(rasat2);
        astro::precession2(rasat2, decsat2, t0, t1, rasat22, decsat22);
        astUtils::rebox(rasat22);
        dra_prec = rasat2-rasat22; ddec_prec = decsat2-decsat22;
        astro::nutation_radec(rasat22, decsat22, t0, deltara, deltadec);
        dra_nut = deltara; ddec_nut = deltadec;
        astro::aberration(rasat22, decsat22, t0, deltara, deltadec);
        dra_abe = deltara; ddec_abe = deltadec;
        rasat22 = rasat2 - dra_prec - dra_nut - dra_abe;
        decsat22 = decsat2 - ddec_prec - ddec_nut - ddec_abe;
	/* ******************************************************************* */
	
	
	/* ******************************************************************* */
	// PHOTO OF SATELLITE
	// obs 1
	rr=1.0;
	astro::azel_xyz(rr, azsat1, elsat1, eTo, rp);
	astro::project(xc1, rp, xp);
	astro::rotate3d(trot1, xp, xp0);
	// invert axis for photo with positive parity
        xr0[0]=-xp0[1];
        xr0[1]=-xp0[2];
	// in pixel
	xx0[0]=(xr0[0]+0.5*sizex)/pxsize;
	xx0[1]=(0.5*sizey+xr0[1])/pxsize;
	// Formato per Leo: FRAME X Y RA DEC X Y Z TIME
	if ((xr0[0] > -0.5*sizex) && (xr0[0] < 0.5*sizex)){
          if ((xr0[1] > -0.5*sizey) && (xr0[1] < 0.5*sizey)){
	    fprintf(photofile1," %03i %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %i %02i %02i %02i %02i %04.2f \n",
		    iframe, xx0[0], xx0[1], rasat11*rad2deg, decsat11*rad2deg, rsat_ecef[0], rsat_ecef[1], rsat_ecef[2],
		    year, mon, day, hr, min, sec);
          }
        }
	// rotate by gamma for recostruction
	angle1=angle_dec1;
	astro::rotate2d(angle1, xr0, xr);
	// in pixel
	xx[0]=(xr[0]+0.5*sizex)/pxsize;
        xx[1]=(0.5*sizey+xr[1])/pxsize;	
	if ((xr[0] > -0.5*sizex) && (xr[0] < 0.5*sizex)){
	  if ((xr[1] > -0.5*sizey) && (xr[1] < 0.5*sizey)){
	    fprintf(photofile1r," %03i %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %i %02i %02i %02i %02i %04.2f \n",                    
		    iframe, xx[0], xx[1], rasat11*rad2deg, decsat11*rad2deg, rsat_ecef[0], rsat_ecef[1], rsat_ecef[2],
		    year, mon, day, hr, min, sec);
	  }
	}
	
	// obs 2
        rr=1.0;
	astro::azel_xyz(rr, azsat2, elsat2, eTo, rp);
	astro::project(xc2, rp, xp);
	astro::rotate3d(trot2, xp, xp0);
	// invert axis for photo with positive parity
	xr0[0]=-xp0[1];
        xr0[1]=-xp0[2];
	// in pixel
        xx0[0]=(xr0[0]+0.5*sizex)/pxsize;
        xx0[1]=(0.5*sizey+xr0[1])/pxsize;
	// for LEO:  FRAME X Y RA DEC X Y Z TIME 	
        if ((xr0[0] > -0.5*sizex) && (xr0[0] < 0.5*sizex)){
          if ((xr0[1] > -0.5*sizey) && (xr0[1] < 0.5*sizey)){
	    fprintf(photofile2," %03i %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %i %02i %02i %02i %02i %04.2f \n",                    
		    iframe, xx0[0], xx0[1], rasat22*rad2deg, decsat22*rad2deg, rsat_ecef[0], rsat_ecef[1], rsat_ecef[2],
		    year, mon, day, hr, min, sec);
          }
        }
	// rotate by gamma for reconstruction
	angle2=angle_dec2;
	astro::rotate2d(angle2, xr0, xr);
	xx[0]=(xr[0]+0.5*sizex)/pxsize;
	xx[1]=(0.5*sizey+xr[1])/pxsize;
	// for LEO: FRAME X Y RA DEC X Y Z TIME 
	if ((xr[0] > -0.5*sizex) && (xr[0] < 0.5*sizex)){
          if ((xr[1] > -0.5*sizey) && (xr[1] < 0.5*sizey)){
            fprintf(photofile2r," %03i %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %i %02i %02i %02i %02i %04.2f \n",
		    iframe, xx[0], xx[1], rasat22*rad2deg, decsat22*rad2deg, rsat_ecef[0], rsat_ecef[1], rsat_ecef[2],
		    year, mon, day, hr, min, sec);
          }
        }
	/* ******************************************************************* */
	
	
      } while ((feof(infile1) == 0));
  }// end of infile
  
  // close input files
  fclose(infile1);
  fclose(infile2);
  fclose(eceffile);
  
  // close output files
  fclose(matfile1);
  fclose(photofile1);
  fclose(photofile1r);
  fclose(matfile2);
  fclose(photofile2);
  fclose(photofile2r);
  
  printf(" ************************************************ \n");
  printf(" \n ");
  
  return 0;
}
