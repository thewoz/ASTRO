/* -------------------------------------------------------------------- *\
 *                                                                      *
 *  VIESAT_2OBS                                                         *
 *  Assumptions: TEME==JNOW  /  ECI==J2000                              *
 *    sun gives coordinates in MOD==TEME (see vallado book)             *
 *                                                                      *
 * -------------------------------------------------------------------- */

#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>

#include "SGP4.h"
#include "ast2Body.h"
#include "astMath.h"
#include "astTime.h"
#include "astIOD.h"
#include "astUtils.h"
#include "coordFK5.h"
#include "EopSpw.h"

//DEFINES
#define VIEW //if want to check visibility conditions
//#define ALLOUTPUT //if want to write all output files

// structure for degrees,minutes,seconds
typedef struct{ 
  int d,m;
  double s;
}dms;
// structure for hours,minutes,seconds
typedef struct{
  int h,m;
  double s;
}hms;

//Wrapper for rad to dms with sign treatment
dms rad_to_dms(double x){
  int deg,min;
  double sec;
  dms tmp;
  double appo;
  appo=fabs(x);
  astTime::dms_rad(deg,min,sec,eFrom,appo);
  deg*=(int) (x/appo);
  tmp.d=deg;tmp.m=min;tmp.s=sec;
  return tmp;
}
//Wrapper for rad to hms with sign check
hms rad_to_hms(double x){
  int hr,min;
  double sec;
  hms tmp;
  if(x<0){
    fprintf(stderr, "cannot convert %g in HMS because is negative\n",x);
    exit(9);
  }
  astTime::hms_rad(hr,min,sec,eFrom,x);
  tmp.h=hr;tmp.m=min;tmp.s=sec;
  return tmp;
}


/*  ***************** GLOBAL VARIABLES ********************* */
char typerun, typeinput, opsmode;
gravconsttype  whichconst;
int whichcon;


/*  ********************************************************* *\
 *                                                            * 
 *                            M A I N                         *            
 *                                                            *
\*  ********************************************************* */


int main(int argc, char* argv[]){
  
  // SATELLITE VARIABLES IN DIFFERENT RF
  elsetrec satrec;
  double rsat_teme[3], vsat_teme[3], asat_teme[3];
  double rsat_ecef[3], vsat_ecef[3], asat_ecef[3];
  double rsat_eci[3], vsat_eci[3], asat_eci[3];  
  double latsat_gc, latsat_gd, lonsat, altsat;
  double lst1,gst1,lst2,gst2; 
  // orbital parameters
  double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;

  // SAT VARIABLES FROM 1ST OBSERVER
  double rangemagsat_sezobs1, azsat_sezobs1, elsat_sezobs1, appelsat_sezobs1;
  double drangemagsat_sezobs1, dazsat_sezobs1, delsat_sezobs1;
  double rangesat_teqobs1, rasat_teqobs1, decsat_teqobs1;
  double apprasat_teqobs1, appdecsat_teqobs1;
  double drangesat_teqobs1, dradsat_teqobs1, ddecsat_teqobs1;
  double rangesat_teqnowobs1, rasat_teqnowobs1, decsat_teqnowobs1;
  double apprasat_teqnowobs1, appdecsat_teqnowobs1;
  double drangesat_teqnowobs1, dradsat_teqnowobs1, ddecsat_teqnowobs1;  
  double rangesat_teq2kobs1, rasat_teq2kobs1, decsat_teq2kobs1;
  double apprasat_teq2kobs1, appdecsat_teq2kobs1;
  double drangesat_teq2kobs1, dradsat_teq2kobs1, ddecsat_teq2kobs1;
  bool   seesat1, seesun1;
  int 	 elsundeg1;
  // SAT VARIABLES FROM 2ND OBSERVER
  double rangemagsat_sezobs2, azsat_sezobs2, elsat_sezobs2, appelsat_sezobs2;
  double drangemagsat_sezobs2, dazsat_sezobs2, delsat_sezobs2;
  double rangesat_teqobs2, rasat_teqobs2, decsat_teqobs2;
  double apprasat_teqobs2, appdecsat_teqobs2;
  double drangesat_teqobs2, dradsat_teqobs2, ddecsat_teqobs2;
  double rangesat_teqnowobs2, rasat_teqnowobs2, decsat_teqnowobs2;
  double apprasat_teqnowobs2, appdecsat_teqnowobs2;
  double drangesat_teqnowobs2, dradsat_teqnowobs2, ddecsat_teqnowobs2;
  double rangesat_teq2kobs2, rasat_teq2kobs2, decsat_teq2kobs2;
  double apprasat_teq2kobs2, appdecsat_teq2kobs2;
  double drangesat_teq2kobs2, dradsat_teq2kobs2, ddecsat_teq2kobs2;
  bool   seesat2, seesun2;
  int 	 elsundeg2;
  
  //SUN VARIABLES
  double rsun_eci[3],vsun_eci[3],asun_eci[3]; //Not used
  double rsun_teme[3],vsun_teme[3],asun_teme[3];
  double rsun_ecef[3],vsun_ecef[3],asun_ecef[3];
  double latsun_gc, latsun_gd, lonsun, altsun;
  dms    sunel_obs1,sunel_obs2;
  
  // SUN VARIABLES FROM 1ST OBSERVER
  double rangemagsun_sezobs1, azsun_sezobs1, elsun_sezobs1;
  double drangemagsun_sezobs1, dazsun_sezobs1, delsun_sezobs1;
  double rangesun_teqobs1, rasun_teqobs1, decsun_teqobs1;
  double drangesun_teqobs1, dradsun_teqobs1, ddecsun_teqobs1;
  // SUN VARIABLES FROM 2ND OBSERVER
  double rangemagsun_sezobs2, azsun_sezobs2, elsun_sezobs2;
  double drangemagsun_sezobs2, dazsun_sezobs2, delsun_sezobs2;
  double rangesun_teqobs2, rasun_teqobs2, decsun_teqobs2;
  double drangesun_teqobs2, dradsun_teqobs2, ddecsun_teqobs2;
  
  // 1ST OBSERVER VARIABLES
  double latobs1_gd, lonobs1, altobs1;
  double robs1_ecef[3], vobs1_ecef[3], aobs1_ecef[3];
  double robs1_teme[3], vobs1_teme[3], aobs1_teme[3];
  double robs1_eci[3], vobs1_eci[3], aobs1_eci[3];
  char nameobs1[128];
  // 2ND OBSERVER VARIABLES
  double latobs2_gd, lonobs2, altobs2;
  double robs2_ecef[3], vobs2_ecef[3], aobs2_ecef[3];
  double robs2_teme[3], vobs2_teme[3], aobs2_teme[3];
  double robs2_eci[3], vobs2_eci[3], aobs2_eci[3];
  char nameobs2[128];
  FILE *obsfile;
  
  // TLE STUFF
  char tle1[69],tle2[69];
  char longstr1[130],longstr2[130];
  char outname[64];
  char str[2];
  char infilename[75];  
  FILE *infile;
 
  // TIME VARIABLES
  double t0, t1;
  double sec, jd, jdFrac, tsince, startmfe, stopmfe, deltamin;
  int year, mon, day, hr, min, timezone, dat;
  double dut1, ut1, tut1, jdut1, jdut1Frac, utc, tai;
  double tt, ttt, jdtt, jdttFrac, tcg, tdb, ttdb, jdtdb, jdtdbFrac, tcb;
  double jdtcg, jdtcgfrac, jdtcb, jdtcbfrac; //extra variabili per versione di Leo
  typedef char str3[4];

  //VARIABLES FOR OUTPUT-MANAGEMENT 
  int    iseen1,iseen2,iwrite;
  double Sstart_jd, Sstart_jdFrac;
  double Sstop_jd, Sstop_jdFrac;
    
  // EOP AND IAU80 VARIABLES
  char eopfilein[140];
  double  ddpsi, ddeps, dx, dy, x, y, s, deltapsi, deltaeps;
  double  lod, xp, yp;
  double jdeop, jdeopFrac;
  std::vector<eopdata> eopinfo;
  char   EopLoc[128];
  iau80data iau80rec;
  
  // CONSTANTS FOR CONVERSIONS
  const double pi = 2.0*asin(1.0);
  const double au2km=149597870.7;  //astronomical units 2 Km
  const double rad2deg = 180.0 / M_PI; //rad2deg 
  
  // OUTPUT FILES
  char filename[128];
  FILE *fin;
  FILE *temefile, *ecifile, *ecefile, *filevisibility;
  FILE *geofile, *anglefilegc, *orbfile;
  FILE *sunfile1, *anglefile1, *filevisibility1;
  FILE *sunfile2, *anglefile2, *filevisibility2;
  FILE *ecefsite, *temesite1, *ecisite1, *temesite2, *ecisite2;
  
  // ------------------------  implementation   --------------------------
  
  printf("%s \n", SGP4Version);  
  //define operation modes
  opsmode = 'a';
  typerun = 'c';  //m manual -- f from file
  typeinput = 'e'; // e input date hours
  whichconst = wgs84; // wgs72 to compare with previsat
  
  /* NB we work in Greenwhich time */
  timezone=0;
  
  /* ***************** READ SITE COORDINATES ******************* */  
  if ((obsfile=fopen("obs.inp","r"))==NULL){
    printf(" ERROR: file obs.inp does not exist \n");
    exit(9);
  }
  fscanf(obsfile,"%s %lf %lf %lf", nameobs1, &lonobs1, &latobs1_gd, &altobs1);
  fscanf(obsfile,"%s %lf %lf %lf", nameobs2, &lonobs2, &latobs2_gd, &altobs2); 
  fclose(obsfile);
  printf(" \n");
  printf(" ******************* VIEWSAT ******************* \n");
  printf(" observatory:  %s \n", nameobs1);
  printf(" lon, lat, alt = %lf %lf %lf \n", lonobs1, latobs1_gd, altobs1);
  printf(" observatory:  %s \n", nameobs2);
  printf(" lon, lat, alt = %lf %lf %lf \n", lonobs2, latobs2_gd, altobs2);
  printf(" ************************************************ \n");
  //convert to proper units
  lonobs1/=rad2deg; //in rads
  latobs1_gd/=rad2deg;  //in rads
  altobs1/=1000.;   // in km
  lonobs2/=rad2deg;
  latobs2_gd/=rad2deg;
  altobs2/=1000.;
  
  //site position in ecef
  for(int k=0;k<3;k++){
    vobs2_ecef[k]=vobs1_ecef[k]=0.0;
  }
  astIOD::site(latobs1_gd, lonobs1, altobs1,robs1_ecef, vobs1_ecef);
  aobs1_ecef[0]=aobs1_ecef[1]=aobs1_ecef[2]=0.0;
  astIOD::site(latobs2_gd, lonobs2, altobs2,robs2_ecef, vobs2_ecef);
  aobs2_ecef[0]=aobs2_ecef[1]=aobs2_ecef[2]=0.0;  
  // write ecef site position
  // SITE 1
  sprintf(filename,"ecef.site.%s",nameobs1);
  ecefsite  = fopen(filename, "w");
  fprintf(ecefsite," %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
	  robs1_ecef[0], robs1_ecef[1], robs1_ecef[2], vobs1_ecef[0], vobs1_ecef[1], vobs1_ecef[2]);
  fclose(ecefsite);
  // SITE 2
  sprintf(filename,"ecef.site.%s",nameobs2);
  ecefsite  = fopen(filename, "w");
  fprintf(ecefsite," %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
          robs2_ecef[0], robs2_ecef[1], robs2_ecef[2], vobs2_ecef[0], vobs2_ecef[1], vobs2_ecef[2]);
  fclose(ecefsite);
  /* **************************************************************** */

  /*  ******************* READ IAU80 and EOP data ******************* */   
  //  set  iau80 file and compute variables
  sprintf(EopLoc,"/usr/local/include/vallado/data/nut80.dat");
  if ((fin=fopen(EopLoc,"r"))==NULL){
    printf(" ERROR: file %s does not exist \n",EopLoc);
    exit(9);
  }else{
    fclose(fin);
  }
  coordFK5::iau80in(iau80rec,EopLoc);
  
  // set EOP values
  sprintf(eopfilein, "./EOP-All.txt");
  if ((fin=fopen(eopfilein,"r"))==NULL){
    printf(" ERROR: file %s does not exist \n",eopfilein);
    exit(9);
  }else{
    fclose(fin);
  }
  
  EopSpw::initeop(eopinfo, eopfilein, jdeop, jdeopFrac);
  /* **************************************************************** */
  
  
  // ---------------- setup files for operation ------------------
  // input 2-line element set file
  //printf("input elset filename: \n");       
  //std::cin >> infilename;
  //infile = fopen(infilename, "r"); 
  
  if(argc!=2){
    printf("usage %s <TLE> \n",argv[0]);
    exit(1);
  }
  infile=fopen(argv[1],"r");
  
  /* ********************* LOOP ON TLEs in file *********************** */
  
  
  while (feof(infile) == 0){      
    /* ******************   READ TLE line 1 ****************************** */
    do
      {
	fgets(longstr1, 130, infile);
	strncpy(str, &longstr1[0], 1);
	str[1] = '\0';
      } while ((strcmp(str, "#") == 0) && (feof(infile) == 0));
    
    if (feof(infile) == 0)
       {
	 /* ****************** READ TLE line 2 ****************************** */
	 fgets(longstr2, 130, infile);
	 // convert the char string to sgp4 elements
	 // includes initialization of sgp4 and jd, jdfrac and sgp4init
	 
	 /* ****************** INITIALIZE FROM TLE ****************************** */
	 SGP4Funcs::twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst,
			       startmfe, stopmfe, deltamin, satrec);
	 printf(" NORAD OBJECT %ld \n", satrec.satnum);


	 /* ****************** READ TIMES *********************** */	 
	 int startyear, startmon, startday, starthr, startmin;
	 double startsec, jdstart, jdstartF;
	 int stopyear,	stopmon, stopday, stophr, stopmin;
         double stopsec, jdstop, jdstopF;
	 FILE* initfile;
	 if ((initfile=fopen("init.txt","r"))==NULL){
	   printf(" ERROR: file init.txt does not exist \n");
	   exit(9);
	 }
	 fscanf(initfile,"%i %i %i %i %i %lf ", &startyear, &startmon, &startday, &starthr, &startmin, &startsec);
         fscanf(initfile,"%i %i %i %i %i %lf %lf ", &stopyear, &stopmon, &stopday, &stophr, &stopmin, &stopsec, &deltamin);
	 fclose(initfile);
	 astTime::jday(startyear, startmon, startday, starthr, startmin, startsec, jdstart, jdstartF);	 
	 astTime::jday(stopyear, stopmon, stopday, stophr, stopmin, stopsec, jdstop, jdstopF);
	 startmfe = (jdstart - satrec.jdsatepoch) * 1440.0 + (jdstartF - satrec.jdsatepochF) * 1440.0;
	 stopmfe = (jdstop - satrec.jdsatepoch) * 1440.0 + (jdstopF - satrec.jdsatepochF) * 1440.0;
	 /* ****************************************************** */
	 
	 
	 /* ******************   OPEN OUTPUT FILES ****************************** */
	 //open output files (a regime alcuni di questi potrebbero non essere piu necessari)
	 
	 //FILES WITH DATA DEPENDING ON SITE 1
#ifdef ALLOUTPUT
	 sprintf(filename,"azelradec.%ld.%s",satrec.satnum,nameobs1);
	 anglefile1  = fopen(filename, "w");
	 sprintf(filename,"sun.%ld.%s",satrec.satnum,nameobs1);
	 sunfile1  = fopen(filename, "w");
	 sprintf(filename,"teme.site.%s",nameobs1);
	 temesite1  = fopen(filename, "w");
	 sprintf(filename,"eci.site.%s",nameobs1);
	 ecisite1  = fopen(filename, "w");
#endif
	 
         //FILES WITH DATA DEPENDING ON SITE 2
#ifdef ALLOUTPUT
	 sprintf(filename,"azelradec.%ld.%s",satrec.satnum,nameobs2);
	 anglefile2  = fopen(filename, "w");
	 sprintf(filename,"sun.%ld.%s",satrec.satnum,nameobs2);
	 sunfile2  = fopen(filename, "w");
	 sprintf(filename,"teme.site.%s",nameobs2);
	 temesite2  = fopen(filename, "w");
	 sprintf(filename,"eci.site.%s",nameobs2);
	 ecisite2  = fopen(filename, "w");
#endif
	 
#ifdef VIEW
         sprintf(filename,"viewsite.%ld.%s",satrec.satnum,nameobs1);
         filevisibility1  = fopen(filename, "w");
	 sprintf(filename,"viewsite.%ld.%s",satrec.satnum,nameobs2);
         filevisibility2  = fopen(filename, "w");
#endif	 

	 
	 //FILES WITH DATA NOT DEPENDING ON SITES
         sprintf(filename,"ecef.%ld",satrec.satnum);
         ecefile  = fopen(filename, "w");
         sprintf(filename,"eci.%ld",satrec.satnum);
         ecifile  = fopen(filename, "w");
#ifdef ALLOUTPUT
	 sprintf(filename,"teme.%ld",satrec.satnum);
	 temefile = fopen(filename, "w");	 
	 sprintf(filename,"orb.%ld",satrec.satnum);
	 orbfile  = fopen(filename, "w");
	 sprintf(filename,"radecgc.%ld",satrec.satnum);
	 anglefilegc  = fopen(filename, "w");
	 sprintf(filename,"geo.%ld",satrec.satnum);
         geofile  = fopen(filename, "w");
#endif
	 
	 /* ******************   SGP4 INITIALIZE  ****************************** */
	 // call the propagator to get the initial state vector value
	 // no longer need gravconst since it is assigned in sgp4init
	 SGP4Funcs::sgp4(satrec, 0.0, rsat_teme, vsat_teme);

	 /* ******************   WRITE SOME LOGS ****************************** */
	 // generate .e files for stk
	 jd = satrec.jdsatepoch;
	 jdFrac = satrec.jdsatepochF;
	 strncpy(outname, &longstr1[2], 5);
	 outname[5] = '.';
	 outname[6] = 'e';
	 outname[7] = '\0';
	 
	 SGP4Funcs::invjday(jd, jdFrac, year, mon, day, hr, min, sec);
	 printf(" Date: %i/%02i/%02i %02i:%02i:%04.2f \n", year, mon, day, hr, min, sec);
	 printf(" ************************************************ \n");
	 printf(" Ephemeris Begins ... \n");
	 printf(" Number of Points = 146 \n");
	 printf(" Interpolation Method: Lagrange \n");
	 printf(" Interpolation Order = 5 \n");
	 printf(" Central Body : Earth \n");
	 printf(" CoordinateSystem : TEME \n");
	 printf(" Distance Unit: Kilometers \n");
	 
	 tsince = startmfe;
	 // check so the first value isn't written twice
	 if (fabs(tsince) > 1.0e-8)
	   tsince = tsince - deltamin;

	 printf(" dt = %g \n", deltamin);
	 printf(" start time = %g \n", startmfe);
	 printf(" final time = %g \n", stopmfe);
	 /* ******************   PROPAGATE TLE ****************************** */ 
	 // ----------------- loop to perform the propagation ----------------

	 iseen1=iseen2=iwrite=0; //set flags observation to zero
	 stopmfe -= deltamin;
	 while ((tsince < stopmfe) && (satrec.error == 0))
	   {
	     tsince = tsince + deltamin;	      
	     //per evitare doppia stampa alla fine
	     //	     if (tsince > stopmfe)
	     //  tsince = stopmfe;	     
	     // propagate satellite by up to tsince epoch
	     SGP4Funcs::sgp4(satrec, tsince, rsat_teme, vsat_teme);
	     asat_teme[0] = 0.0;asat_teme[1] = 0.0; asat_teme[2] = 0.0;
	     
	     if (satrec.error > 0)
	       printf("# *** error: t:= %f *** code = %3d\n", satrec.t, satrec.error);
	     
	     if (satrec.error == 0) //IF NO ERROR
	       {
		 //manage julian date
		 jd = satrec.jdsatepoch;
		 jdFrac = satrec.jdsatepochF + tsince / 1440.0;
		 if (jdFrac < 0.0)
		   {
		     jd = jd - 1.0;
		     jdFrac = jdFrac + 1.0;
		   }
		 // convert jday to gregorian day and time		  		  
		 SGP4Funcs::invjday(jd, jdFrac, year, mon, day, hr, min, sec);	    
		 if(fabs(sec-60.)<1.e-8){
		   sec=0.0;
		    min++;
		 }
		 if(min==60){
		   min=0;
		   hr++;
		 }

		 /* ***************** FIND SUN POSITION ******************** */
		 //sun position in  mod from Vallado book Here we assume MOD==TEME
		 // this should be the ra/dec in  geocentric based on teme/mod
		 double ras,decs;
		 ast2Body::sun(jd,jdFrac,rsun_teme,ras,decs);
		 //convert from au to Km and initialize sun vel acc
		 for(int k=0;k<3;k++){
		   rsun_teme[k]*=au2km;
		   vsun_teme[k]=0.0;
		   asun_teme[k]=0.0;
		 }
		 
		 /* ********* TIME COVERSIONS ACCOUNTING FOR EOP  ********* */		 
		 EopSpw::findeopparam(jd, jdFrac, 'l', 'f', eopinfo, jdeop, dut1, dat,
				      lod, xp, yp, ddpsi, ddeps, dx, dy, x, y, s, deltapsi, deltaeps);
		 astTime::convtime(year, mon, day, hr, min, sec, timezone, dut1, dat,ut1,tut1,jdut1,jdut1Frac,utc,tai, tt, ttt,
				   jdtt, jdttFrac,tdb,ttdb,jdtdb,jdtdbFrac,tcg,jdtcg, jdtcgfrac,tcb,jdtcb,jdtcbfrac);
		 
                 // Per il confronto con PREVISAT ho capito che loro non usano
		 // tutte le correzioni EOP etc per ottenere gli stessi valori in ecef ad esempio
		 // basta implementare le linee che seguono, rimane la differenza sul sole perche usano una routine diversa
		 //#define PREVISAT_ECEF_OK //scommentare per fare confronto con previsat
#ifdef PREVISAT_ECEF_OK
		 xp=yp=0.0;
		 jdut1=jd;
		 jdut1Frac=jdFrac;
		 lod=0;
		 ttt = (jdut1+jdut1Frac - 2451545.0) / 36525.0;
#endif
		 
		 /* ************ FRAME CONVERSIONS ******************** */

		 //FOR SAT		  
		 //teme->ecef
		 coordFK5::teme_ecef(rsat_teme, vsat_teme, asat_teme, eTo, rsat_ecef, vsat_ecef, asat_ecef, ttt, jdut1+jdut1Frac, lod, xp, yp, 2);
		 //ecef-->lat/lon/alt
		 ast2Body::ijk2ll(rsat_ecef, latsat_gc, latsat_gd, lonsat, altsat);
		 //teme->eci
		 coordFK5::teme_eci(rsat_teme,vsat_teme,asat_teme,eTo,rsat_eci,vsat_eci,asat_eci,iau80rec,ttt,ddpsi,ddeps);

		 
		 //FOR SITE 1
		 //teme->ecef
		 coordFK5::teme_ecef(robs1_teme, vobs1_teme, aobs1_teme, eFrom, robs1_ecef, vobs1_ecef, aobs1_ecef, ttt, jdut1+jdut1Frac, lod, xp, yp, 2);
		 //teme->eci
		 coordFK5::teme_eci(robs1_teme,vobs1_teme,aobs1_teme,eTo,robs1_eci,vobs1_eci,aobs1_eci,iau80rec,ttt,ddpsi,ddeps);
		 // local sidereal time
		 astTime::lstime(lonobs1, jdut1+jdut1Frac, lst1, gst1);
		 //FOR SITE 2
		 //teme->ecef
		 coordFK5::teme_ecef(robs2_teme, vobs2_teme, aobs2_teme, eFrom, robs2_ecef, vobs2_ecef, aobs2_ecef, ttt, jdut1+jdut1Frac, lod, xp, yp, 2);
		 //teme->eci
		 coordFK5::teme_eci(robs2_teme,vobs2_teme,aobs2_teme,eTo,robs2_eci,vobs2_eci,aobs2_eci,iau80rec,ttt,ddpsi,ddeps);
		 // local sidereal time
		 astTime::lstime(lonobs2, jdut1+jdut1Frac, lst2, gst2);

		 		 
		 //FOR SUN
		 //teme->ecef
		 coordFK5::teme_ecef(rsun_teme, vsun_teme, asun_teme, eTo, rsun_ecef, vsun_ecef, asun_ecef, ttt, jdut1+jdut1Frac, lod, xp, yp, 2);
		 //ecef-->lat/lon/alt
                 ast2Body::ijk2ll(rsun_ecef, latsun_gc, latsun_gd, lonsun, altsun);
		 
		 		 
		 /* ************ COMPUTATIONS OF QUANTITIES ON SATELLITE ******************** */
		 
		 // compute satellite orbital elements
		 SGP4Funcs::rv2coe(rsat_teme, vsat_teme, satrec.mu, p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper);

		 
		 //FOR SITE 1
		 // compute azimuth and elevation of the satellite with respect to site in topocentric SEZ system
		 astIOD::rv_razel(rsat_ecef, vsat_ecef, robs1_ecef,latobs1_gd, lonobs1, eTo,
				  rangemagsat_sezobs1, azsat_sezobs1, elsat_sezobs1, drangemagsat_sezobs1, dazsat_sezobs1, delsat_sezobs1);
		 astUtils::rebox(azsat_sezobs1);
		 astUtils::refraction(elsat_sezobs1, &appelsat_sezobs1);
		 
		 // compute right ascension and declination of the satellite with respect to site in topocentric equatorial system (JNOW)
		 astIOD::rv_tradec(rsat_teme, vsat_teme, robs1_teme, eTo,
				   rangesat_teqnowobs1, rasat_teqnowobs1, decsat_teqnowobs1, drangesat_teqnowobs1, dradsat_teqnowobs1, ddecsat_teqnowobs1);
		 astUtils::rebox(rasat_teqnowobs1);
		 //the call here below should be used if we want to compute the apparent ra/dec otherwise comment the next lines
		 astTime::lstime(lonobs1, jdut1+jdut1Frac, lst1, gst1);
		 astIOD::radec_azel(apprasat_teqnowobs1, appdecsat_teqnowobs1, lst1, latobs1_gd, eFrom, azsat_sezobs1, appelsat_sezobs1);
		 astUtils::rebox(apprasat_teqnowobs1);

		 // compute right ascension and declination of the satellite with respect to site in topocentric equatorial system (J2000)
		 astIOD::rv_tradec(rsat_eci, vsat_eci, robs1_eci, eTo,
				   rangesat_teq2kobs1, rasat_teq2kobs1, decsat_teq2kobs1, drangesat_teq2kobs1, dradsat_teq2kobs1, ddecsat_teq2kobs1);
                 astUtils::rebox(rasat_teq2kobs1);
		 
		 
		 //FOR SITE 2
		 // compute azimuth and elevation of the satellite with respect to site in topocentric SEZ system
		 astIOD::rv_razel(rsat_ecef, vsat_ecef, robs2_ecef, latobs2_gd, lonobs2, eTo,
				  rangemagsat_sezobs2, azsat_sezobs2, elsat_sezobs2, drangemagsat_sezobs2, dazsat_sezobs2, delsat_sezobs2);
		 astUtils::rebox(azsat_sezobs2);
		 astUtils::refraction(elsat_sezobs2, &appelsat_sezobs2);
		 
		 //compute right ascension and declination of the satellite with respect to site in topocentric equatorial system (JNOW)		 
		 astIOD::rv_tradec(rsat_teme, vsat_teme, robs2_teme, eTo,
		 rangesat_teqnowobs2, rasat_teqnowobs2, decsat_teqnowobs2, drangesat_teqnowobs2, dradsat_teqnowobs2, ddecsat_teqnowobs2);
                 astUtils::rebox(rasat_teqnowobs2);
		 //the call here below should be used if we want to compute the apparent ra/dec otherwise comment the next lines
		 astTime::lstime(lonobs2, jdut1+jdut1Frac, lst2, gst2);
		 astIOD::radec_azel(apprasat_teqnowobs2, appdecsat_teqnowobs2, lst2, latobs2_gd, eFrom, azsat_sezobs2, appelsat_sezobs2);
		 astUtils::rebox(apprasat_teqnowobs2);

                 //compute right ascension and declination of the satellite with respect to site in topocentric equatorial system (J200)
		 astIOD::rv_tradec(rsat_eci, vsat_eci, robs2_eci, eTo,
				   rangesat_teq2kobs2, rasat_teq2kobs2, decsat_teq2kobs2, drangesat_teq2kobs2, dradsat_teq2kobs2, ddecsat_teq2kobs2);
                 astUtils::rebox(rasat_teq2kobs2);

		 
		 // test compute ra dec in geocentric ref
		 // Assume ECI==J2000 and use rv_radec
		 double rangesat_gc, rasat_gc, decsat_gc,  drangesat_gc,  drassat_gc,  ddecsat_gc;
		 astIOD::rv_radec(rsat_eci, vsat_eci, eTo,rangesat_gc, rasat_gc, decsat_gc,  drangesat_gc,  drassat_gc, ddecsat_gc);
		 astUtils::rebox(rasat_gc);
		 
		 
		 /* ************ COMPUTATIONS ON SUN ******************** */
		 
		 // FROM SITE 1
		 // compute azimuth and elevation of the sun with respect to site in topocentric SEZ system
		 // nb ecef qui è ottenuto da teme di sun ma in realta' è mod
		 astIOD::rv_razel(rsun_ecef, vsun_ecef, robs1_ecef,latobs1_gd, lonobs1,eTo,
				  rangemagsun_sezobs1, azsun_sezobs1, elsun_sezobs1, drangemagsun_sezobs1, dazsun_sezobs1, delsun_sezobs1);
		 astUtils::rebox(azsun_sezobs1);
		 
		 // compute right ascension and declination of the sun with respect to site in topocentric equatorial system
		 astIOD::rv_tradec(rsun_teme, vsun_teme, robs1_teme,eTo,
				   rangesun_teqobs1, rasun_teqobs1, decsun_teqobs1,drangesun_teqobs1, dradsun_teqobs1, ddecsun_teqobs1);
		 astUtils::rebox(rasun_teqobs1);

		 // FROM SITE 2
		 astIOD::rv_razel(rsun_ecef, vsun_ecef, robs2_ecef,latobs2_gd, lonobs2,eTo,
				  rangemagsun_sezobs2, azsun_sezobs2, elsun_sezobs2, drangemagsun_sezobs2, dazsun_sezobs2, delsun_sezobs2);
		 astUtils::rebox(azsun_sezobs2);
		 
		 astIOD::rv_tradec(rsun_teme, vsun_teme, robs2_teme,eTo,
				   rangesun_teqobs2, rasun_teqobs2, decsun_teqobs2,drangesun_teqobs2, dradsun_teqobs2, ddecsun_teqobs2);
		 astUtils::rebox(rasun_teqobs2);

		 // Conversions of angles from rad to dms/hms               
                 // Sun
                 sunel_obs1  = rad_to_dms(elsun_sezobs1);
                 sunel_obs2  = rad_to_dms(elsun_sezobs2);

		 
		 /* **************************** WRITE OUTPUT FILES ABOUT ORBIT *************************************** */

		 // SATELLITE INFO DEPENDING ON SITE

#ifdef ALLOUTPUT

		 // FOR SITE 1
		 // sat range az appel el appra appdec (topo jnow) ra dec (topo j2000) 
		 fprintf(anglefile1, "%i %02i %02i %02i %02i %09.6f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n",
		 	 year, mon, day, hr, min, sec, azsat_sezobs1*rad2deg, appelsat_sezobs1*rad2deg, elsat_sezobs1*rad2deg,
			 apprasat_teqnowobs1*rad2deg, appdecsat_teqnowobs1*rad2deg, rasat_teqnowobs1*rad2deg, decsat_teqnowobs1*rad2deg,
			 rasat_teq2kobs1*rad2deg, decsat_teq2kobs1*rad2deg, rangemagsat_sezobs1, rangesat_teqobs1);
		 		 		 		 
		 //SUN
		 fprintf(sunfile1, "%i %02i %02i %02i %02i %09.6f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n",
			 year, mon, day, hr, min, sec, rangemagsun_sezobs1/au2km, azsun_sezobs1*rad2deg, elsun_sezobs1*rad2deg,
			 rangesun_teqobs1/au2km, rasun_teqobs1*rad2deg, decsun_teqobs1*rad2deg);

                 // FOR SITE 2
		 // sat range az el ra (topo) dec (topo)
		 fprintf(anglefile2, "%i %02i %02i %02i %02i %09.6f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n",
			 year, mon, day, hr, min, sec, azsat_sezobs2*rad2deg, appelsat_sezobs2*rad2deg, elsat_sezobs2*rad2deg,
			 apprasat_teqnowobs2*rad2deg, appdecsat_teqnowobs2*rad2deg, rasat_teqnowobs2*rad2deg, decsat_teqnowobs2*rad2deg,
			 rasat_teq2kobs2*rad2deg, decsat_teq2kobs2*rad2deg, rangemagsat_sezobs2, rangesat_teqobs2);                 
		 
		 //SUN
		 fprintf(sunfile2, "%i %02i %02i %02i %02i %09.6f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n",
			 year, mon, day, hr, min, sec, rangemagsun_sezobs2/au2km, azsun_sezobs2*rad2deg, elsun_sezobs2*rad2deg,
			 rangesun_teqobs2/au2km, rasun_teqobs2*rad2deg, decsun_teqobs2*rad2deg);
		 
		 //SITE 1
		 //teme 
		 fprintf(temesite1, "%i %02i %02i %02i %02i %09.6f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
			 year, mon, day, hr, min, sec, robs1_teme[0], robs1_teme[1], robs1_teme[2], vobs1_teme[0], vobs1_teme[1], vobs1_teme[2]);

		 //eci
		 fprintf(ecisite1, "%i %02i %02i %02i %02i %09.6f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
			 year, mon, day, hr, min, sec, robs1_eci[0], robs1_eci[1], robs1_eci[2], vobs1_eci[0], vobs1_eci[1], vobs1_eci[2]);

		 //SITE 2
		 //teme 
		 fprintf(temesite2, "%i %02i %02i %02i %02i %09.6f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
			 year, mon, day, hr, min, sec, robs2_teme[0], robs2_teme[1], robs2_teme[2], vobs2_teme[0], vobs2_teme[1], vobs2_teme[2]);

		 //eci
		 fprintf(ecisite2, "%i %02i %02i %02i %02i %09.6f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
			 year, mon, day, hr, min, sec, robs2_eci[0], robs2_eci[1], robs2_eci[2], vobs2_eci[0], vobs2_eci[1], vobs2_eci[2]);


		 		 
		 // SATELLITE INFO NOT DEPENDING ON SITE
		 // sat teme 
		 fprintf(temefile, "%i %02i %02i %02i %02i %09.6f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
			 year, mon, day, hr, min, sec, rsat_teme[0], rsat_teme[1], rsat_teme[2], vsat_teme[0], vsat_teme[1], vsat_teme[2]);
		 
                 //sat lat (gc) lat (gd) lon alt + sun lat (gc) lat (gd) lon alt
                 fprintf(geofile, "%i %02i %02i %02i %02i %09.6f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f \n",
                         year, mon, day, hr, min, sec, latsat_gc*rad2deg, latsat_gd*rad2deg, lonsat*rad2deg, altsat,
			 latsun_gc*rad2deg, latsun_gd*rad2deg, lonsun*rad2deg, ttt, jdut1+jdut1Frac);
		 
		 //sat  ra (geoc) dec (geoc) 
		 fprintf(anglefilegc, "%i %02i %02i %02i %02i %09.6f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n",
			 year, mon, day, hr, min, sec, rangesat_gc, rasat_gc*rad2deg, decsat_gc*rad2deg, drangesat_gc, drassat_gc, ddecsat_gc);

		 // sat orbit parameters
		 fprintf(orbfile, "%i %02i %02i %02i %02i %09.6f %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
			 year, mon, day, hr, min, sec,
			 a, ecc, incl*rad2deg, node*rad2deg, argp*rad2deg, nu*rad2deg, m*rad2deg);
		 
#endif
		 
                 //sat eci
                 fprintf(ecifile, "%i %02i %02i %02i %02i %09.6f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
                         year, mon, day, hr, min, sec, rsat_eci[0], rsat_eci[1], rsat_eci[2], vsat_eci[0], vsat_eci[1], vsat_eci[2]);
		 
                 //sat ecef
                 fprintf(ecefile, "%i %02i %02i %02i %02i %09.6f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
                         year, mon, day, hr, min, sec, rsat_ecef[0], rsat_ecef[1], rsat_ecef[2], vsat_ecef[0], vsat_ecef[1], vsat_ecef[2]);
		 
		 
#ifdef VIEW
		 //detect visibility condition (we should add min elevation for  satellite and max elevation of sun
		 seesat1=astUtils::sight(robs1_ecef,rsat_ecef,'e');
		 seesat2=astUtils::sight(robs2_ecef,rsat_ecef,'e');
		 seesun1=astUtils::sight(robs1_ecef,rsun_ecef,'e');
		 seesun2=astUtils::sight(robs2_ecef,rsun_ecef,'e');
		 if(seesat1 && seesat2 && (!seesun1) && (!seesun2) && astUtils::sight(rsat_ecef,rsun_ecef,'e') && (sunel_obs1.d<-10) && (sunel_obs2.d<-10)){
		   // SAT INFO FROM SITE 1
		   fprintf(filevisibility1," %i %2i %2i %2i %2i %9.6f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f \n",
			   year, mon, day, hr, min, sec, azsat_sezobs1*rad2deg, appelsat_sezobs1*rad2deg, elsat_sezobs1*rad2deg,
			   apprasat_teqnowobs1*rad2deg, appdecsat_teqnowobs1*rad2deg, rasat_teqnowobs1*rad2deg, decsat_teqnowobs1*rad2deg,
			   rasat_teq2kobs1*rad2deg, decsat_teq2kobs1*rad2deg, altsat, rangemagsat_sezobs1);
		   // SAT INFO FROM SITE 2
		   fprintf(filevisibility2," %i %2i %2i %2i %2i %9.6f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f \n",
			   year, mon, day, hr, min, sec, azsat_sezobs2*rad2deg, appelsat_sezobs2*rad2deg, elsat_sezobs2*rad2deg,
			   apprasat_teqnowobs2*rad2deg, appdecsat_teqnowobs2*rad2deg, rasat_teqnowobs2*rad2deg, decsat_teqnowobs2*rad2deg,
			   rasat_teq2kobs2*rad2deg, decsat_teq2kobs2*rad2deg, altsat, rangemagsat_sezobs2);
		   iseen1=1;
		 }else{
		   iseen1=0;
		 } 
		 

		 if(iseen1 && !iseen2){//if seen now and not before store start time (from 1 min before)
		   Sstart_jd=jd;   
		   Sstart_jdFrac=jdFrac-1./1440.0;
		 }else if(!iseen1 && iseen2){//if seen before and not now store stop time
		   Sstop_jd=jd;
		   Sstop_jdFrac=jdFrac;
		   iwrite=1; //and set write flag
		   fprintf(filevisibility1,"\n\n");
		   fprintf(filevisibility2,"\n\n");
		 }

		 iseen2=iseen1;

		 // you have to write start/end for refined integration
		 if(iwrite){
		   SGP4Funcs::invjday(Sstart_jd, Sstart_jdFrac, year, mon, day, hr, min, sec);	    
		   if(fabs(sec-60.)<1.e-8){
		     sec=0.0;
		     min++;
		   }
		   if(min==60){
		     min=0;
		     hr++;
		   }
		   // loading TLE strings
 		   for(int kk=0;kk<69;kk++){
		     tle1[kk]=longstr1[kk];
		     tle2[kk]=longstr2[kk];
		   }
		   // print out start/end
		   fprintf(stdout,"%s %i %i %i %i %i %lf\n",tle1,year,mon,day,hr,min,sec);
		   SGP4Funcs::invjday(Sstop_jd, Sstop_jdFrac, year, mon, day, hr, min, sec);	    
		   if(fabs(sec-60.)<1.e-8){
		     sec=0.0;
		     min++;
		   }
		   if(min==60){
		     min=0;
		     hr++;
		   }		   
		   double deltaminphoto= 0.5/60.; // rate of the photocamera. 1/100s is enough? 
		   fprintf(stdout,"%s %i %i %i %i %i %lf %lf\n",tle2,year,mon,day,hr,min,sec,deltaminphoto);
		   iwrite=0;
		 }
		 
#endif
		 
	       } // if satrec.error == 0
	     
	   } // while propagating the orbit
	 
	 printf(" ... Ephemeris ends \n");
	 printf(" ************************************************ \n");
	 printf(" \n");
	  
	 //close output files	 
#ifdef ALLOUTPUT
	 fclose(anglefile1);
	 fclose(sunfile1);
	 fclose(anglefile2);
	 fclose(sunfile2);
	 fclose(ecefile);
         fclose(geofile);
         fclose(temefile);
         fclose(ecifile);
         fclose(orbfile);
         fclose(anglefilegc);
#endif
	 
#ifdef VIEW
	 fclose(filevisibility1);
	 fclose(filevisibility2);	 
#endif

	 
       } // if not eof
    
  } // while through the input file
 
 return 0;
}  // testSGP4

