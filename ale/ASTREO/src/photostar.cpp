
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

/* ******************** AUTOCENTER  ******************** *\  
 *                                                       *
 * this function extract the center of image from the    *
 * barycenter of the satellite track                     *
 *                                                       *
\* ***************************************************** */
void autocenter(char satnum[128], char nameobs[128], int rifra, double& az0, double& el0)
{
  
  // variables
  double az, appel, el, alt, range;
  double appra, appdec, rajnw, decjnw, raj2k, decj2k;
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
  const double pi = 2.0*asin(1.0);
  const double rad2deg = 180.0 / M_PI;
  const double rad2hrs = 12.0 / M_PI;
  const double pxsize = 0.0065;
  
  // sat infos
  char satnum[128], satname[128];

  //obs1 variables
  double robs1_ecef[3], vobs1_ecef[3];
  double lonobs1, latobs1, altobs1, lst1, gst1;
  double azc1, elc1, rac1, decc1, rac11, decc11;
  double angle1, angle_ra1, angle_dec1, gamma1;
  double tex1, tey1, tez1, rot1[3][3], trot1[3][3], xc1[3];
  double lon1, lat1, mrt1[3][4];
  char nameobs1[128];

  //obs2 variables
  double robs2_ecef[3], vobs2_ecef[3];
  double lonobs2, latobs2, altobs2, lst2, gst2;
  double azc2, elc2, rac2, decc2, rac22, decc22;
  double angle2, angle_ra2, angle_dec2, gamma2;
  double tex2, tey2, tez2, rot2[3][3], trot2[3][3], xc2[3];
  double lon2, lat2, mrt2[3][4];
  char nameobs2[128];

  //photo parameters
  double rr, rp[3], xp[3], xp0[3];
  double xr[2], xr0[2], xx[3], xx0[3];
  double focal, sizex, sizey;
  int rifra;
  FILE *paramfile;

  // time variables
  int year, mon, day, hr, min, timezone;
  int iframe, istar;
  double sec;
  double jd, jdfrac, t, t0, t1;

  // star variables
  double rastar, decstar, pmrastar, pmdecstar, magstar, plxstar;
  double ra1, dec1, deltara, deltadec;
  double dra_abe, ddec_abe, dra_prec, ddec_prec, dra_nut, ddec_nut;
  double rastar1, decstar1, azstar1, elstar1, appelstar1;
  double rastar2, decstar2, azstar2, elstar2, appelstar2;
  int tyc1, tyc2, tyc3, hip, hd;
  
  // i/o files
  char filename[128];
  FILE *infile, *timesheet;
  FILE *matfile1, *photofile1, *photofile1r;
  FILE *matfile2, *photofile2, *photofile2r;
  
  
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
  printf(" ******************* PHOTOSTAR ******************* \n");
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
  printf(" ************************************************* \n");
  
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
  // read time sheet
  if ((timesheet=fopen("timesheet.txt","r"))==NULL){
    printf(" ERROR: file timesheet.txt does not exist \n");
    exit(9);
  }
  // open output files
  sprintf(filename,"mat.%s.%s.dat",satnum,nameobs1);
  matfile1=fopen(filename, "w");
  sprintf(filename,"photostar.%s.%s.dat",satnum,nameobs1);
  photofile1 = fopen(filename, "w");
  sprintf(filename,"photostar.rot.%s.%s.dat",satnum,nameobs1);
  photofile1r = fopen(filename, "w");
  //
  sprintf(filename,"mat.%s.%s.dat",satnum,nameobs2);
  matfile2=fopen(filename, "w");
  sprintf(filename,"photostar.%s.%s.dat",satnum,nameobs2);
  photofile2 = fopen(filename, "w");
  sprintf(filename,"photostar.rot.%s.%s.dat",satnum,nameobs2);
  photofile2r = fopen(filename, "w");
  /* *********************************************************** */

  
  /* loop on timesheet */
  iframe=0;
  while (feof(timesheet)==0){
    do
      {

	/* ******************************************************************* */
	// TIME
	iframe=iframe+1;
	// read date
	fscanf(timesheet," %i %i %i %i %i %lf ", &year, &mon, &day, &hr, &min, &sec);	
	// set julian date
	timezone=0;
	astTime::jday(year, mon, day, hr, min, sec, jd, jdfrac);
	printf(" Frame = %03i | Date = %i/%02i/%02i %02i:%02i:%04.2f \n", iframe, year, mon, day, hr, min, sec);
	t = (jd + jdfrac - 2451545.)/36525.;
	t0 = (jd + jdfrac - 2451545.)/36525.;
        t1 = (2451545. - jd - jdfrac)/36525.;
	/* ******************************************************************* */

	
        /* ******************************************************************* */
        // RA/DEC CENTER & ORIENTATION ANGLE
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
	// OPEN INPUT/OUTPUT FILES
	// open Tycho-2 catalog
	if ((infile=fopen("/usr/local/include/astro/data/tycho2_catalog.dat","r"))==NULL){
	  printf(" ERROR: Tycho-2 catalog not found \n");
	  exit(9);
	}
	/* ******************************************************************* */
	
	// loop on stars	
	istar=0;
	while (feof(infile)==0){
	  do
	    {
	      // read input from Tycho-2 catalog
	      istar=istar+1;
	      fscanf(infile," %lf %lf %lf %lf %lf %lf %i %i %i %i %i ",
		     &rastar, &decstar, &pmrastar, &pmdecstar, &magstar, &plxstar, &tyc1, &tyc2, &tyc3, &hip, &hd);
	      rastar = rastar/rad2deg;
	      decstar = decstar/rad2deg;
	      
	      /* correction due to annual aberration (in J2000) */
	      astro::aberration(rastar, decstar, t, deltara, deltadec);
	      dra_abe = deltara; ddec_abe = deltadec;	      
	      /* correction due to precession (transform from J2000 to JNOW) */
	      astro::precession(rastar, decstar, t, ra1, dec1);
	      astUtils::rebox(ra1);
	      dra_prec = ra1-rastar; ddec_prec = dec1-decstar;              
	      /* correction due to nutation  */
	      astro::nutation_radec(rastar, decstar, t, deltara, deltadec);
	      dra_nut = deltara; ddec_nut = deltadec;	      
	      /* add corrections */
	      rastar1 = rastar + dra_prec + dra_nut + dra_abe;
	      decstar1 = decstar + ddec_prec + ddec_nut + ddec_abe;
              rastar2 = rastar + dra_prec + dra_nut + dra_abe;
              decstar2 = decstar + ddec_prec + ddec_nut + ddec_abe;

	      // obs1
	      /* compute azimuth and elevation */
	      astTime::lstime(lonobs1, jd+jdfrac, lst1, gst1);
	      astIOD::radec_azel(rastar1, decstar1, lst1, latobs1, eTo, azstar1, elstar1);
	      astUtils::rebox(azstar1);
	      /* correction due to atmospheric refraction */
	      if(rifra != 0)
		{
		  astUtils::refraction(elstar1, &appelstar1);
		}
	      else
		{
		  appelstar1 = elstar1;
		}
              // obs2
              /* compute azimuth and elevation */
	      astTime::lstime(lonobs2, jd+jdfrac, lst2, gst2);
	      astIOD::radec_azel(rastar2, decstar2, lst2, latobs2, eTo, azstar2, elstar2);
	      astUtils::rebox(azstar2);
              /* correction due to atmospheric refraction */
              if(rifra != 0)
		{
                  astUtils::refraction(elstar2, &appelstar2);
                }
	      else
                {
                  appelstar2 = elstar2;
                }

	      
	      /* ******************************************************************* */
	      // PHOTO
	      // obs 1
	      if (appelstar1 > 0.0)
		{
		  rr=1.0;
		  astro::azel_xyz(rr, azstar1, appelstar1, eTo, rp);
		  astro::project(xc1, rp, xp);
		  astro::rotate3d(trot1, xp, xp0);
		  // invert and rotate by gamma
		  xr0[0]=-xp0[1];
		  xr0[1]=-xp0[2];
		  // in pixel
		  xx0[0]=(xr0[0]+0.5*sizex)/pxsize;
		  xx0[1]=(0.5*sizey+xr0[1])/pxsize;
		  // for LEO: FRAME ID X Y RA DEC MAG RA0 DEC0 DAT
		  if ((xr0[0] > -0.5*sizex) && (xr0[0] < 0.5*sizex)){
                    if ((xr0[1] > -0.5*sizey) && (xr0[1] < 0.5*sizey)){
                      fprintf(photofile1," %03i %i %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  %i %02i %02i %02i %02i %04.2f \n",
                              iframe, istar, xx0[0], xx0[1], rastar*rad2deg, decstar*rad2deg,
                              magstar, rac11*rad2deg, decc11*rad2deg, angle_dec1*rad2deg,
	                      year, mon, day, hr, min, sec);
                    }
		  }
		  // rotated photo
		  angle1=angle_dec1;
		  astro::rotate2d(angle1, xr0, xr);
		  xx[0]=(xr[0]+0.5*sizex)/pxsize;
                  xx[1]=(0.5*sizey+xr[1])/pxsize;
		  // for LEO: FRAME ID X Y RA DEC MAG RA0 DEC0 DAT
		  if ((xr[0] > -0.5*sizex) && (xr[0] < 0.5*sizex)){
		    if ((xr[1] > -0.5*sizey) && (xr[1] < 0.5*sizey)){
		      fprintf(photofile1r," %03i %i %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  %i %02i %02i %02i %02i %04.2f \n",
			      iframe, istar, xx[0], xx[1], rastar*rad2deg, decstar*rad2deg,
			      magstar, rac11*rad2deg, decc11*rad2deg, 0.0,
			      year, mon, day, hr, min, sec);
		    }
		  }		  
		} // close if on appelstar1

	      // obs 2  
	      if (appelstar2 > 0.0)
		{
		  rr=1.0;
		  astro::azel_xyz(rr, azstar2, appelstar2, eTo, rp);
		  astro::project(xc2, rp, xp);
		  astro::rotate3d(trot2, xp, xp0);
		  // invert and rotate by gamma
		  xr0[0]=-xp0[1];
		  xr0[1]=-xp0[2];
		  // in pixel
		  xx0[0]=(xr0[0]+0.5*sizex)/pxsize;
                  xx0[1]=(0.5*sizey+xr0[1])/pxsize;
		  // for LEO: FRAME ID X Y RA DEC MAG RA0 DEC0 DAT
		  if ((xr0[0] > -0.5*sizex) && (xr0[0] < 0.5*sizex)){
                    if ((xr0[1] > -0.5*sizey) && (xr0[1] < 0.5*sizey)){
                      fprintf(photofile2," %03i %i %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  %i %02i %02i %02i %02i %04.2f \n",
                              iframe, istar, xx0[0], xx0[1], rastar*rad2deg, decstar*rad2deg,
                              magstar, rac22*rad2deg, decc22*rad2deg, angle_dec2*rad2deg,
                              year, mon, day, hr, min, sec);
                    }
                  }
		  // rotated photo
		  angle2=angle_dec2;
		  astro::rotate2d(angle2, xr0, xr);
		  xx[0]=(xr[0]+0.5*sizex)/pxsize;
                  xx[1]=(0.5*sizey+xr[1])/pxsize;
		  // for LEO: FRAME ID X Y RA DEC MAG RA0 DEC0 DATE
		  if ((xr[0] > -0.5*sizex) && (xr[0] < 0.5*sizex)){
		    if ((xr[1] > -0.5*sizey) && (xr[1] < 0.5*sizey)){
		      fprintf(photofile2r," %03i %i %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  %i %02i %02i %02i %02i %04.2f \n",
			      iframe, istar, xx[0], xx[1], rastar*rad2deg, decstar*rad2deg,
			      magstar, rac22*rad2deg, decc22*rad2deg, 0.0,
			      year, mon, day, hr, min, sec);
		    }
		  }
		} // close if on appelstar2
	      
	      /* ******************************************************************* */
	      
	    } while ((feof(infile) == 0));
	}// end of infile
	
      } while ((feof(timesheet) == 0));
  }// end of timesheet
  printf(" ************************************************* \n");
  printf(" \n ");
  
  fclose(timesheet);
  fclose(matfile1);
  fclose(matfile2);
  fclose(photofile1);
  fclose(photofile2);
  fclose(photofile1r);
  fclose(photofile2r);
  
  return 0;
}
