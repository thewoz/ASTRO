
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>


/* ********************************************************** *\
 *                                                            * 
 *                            M A I N                         *            
 *                                                            *
\* ********************************************************** */

int main(int argc, char* argv[]){
  
  int nframe, nstar, nlabel;
  const int nmax=1000;
  
  // output files
  FILE *infile, *outfile;
  char filename[128];
  char filename1[128];
  int i, j, k, ifr, icount;
 
  int frame[nmax], lab[nmax], nlab[nmax], labu[nmax];
  double x[nmax], y[nmax], ra[nmax], dec[nmax], mag[nmax];
  double ra0[nmax], dec0[nmax], angle[nmax];
  int yr[nmax], mon[nmax], day[nmax], hr[nmax], min[nmax];
  double sec[nmax];  
  
  /* read frames and filename from command line */
  if(argc!=3){
    printf("usage %s <#Frames> <FILE> \n",argv[0]);
    exit(1);
  }
  nframe=atoi(argv[1]);
  sprintf(filename,"%s",argv[2]);
  
  /* reading file */
  icount=0;
  infile=fopen(filename,"r");
  while (feof(infile)==0){
    do
      {
	icount=icount+1;
	fscanf(infile," %d %i %lf %lf %lf %lf %lf %lf %lf %lf %i %d %d %d %d %lf ", &frame[icount], &lab[icount],
	       &x[icount], &y[icount], &ra[icount], &dec[icount], &mag[icount], &ra0[icount], &dec0[icount], &angle[icount],
	       &yr[icount], &mon[icount], &day[icount], &hr[icount], &min[icount], &sec[icount]);
      } while ((feof(infile)==0));
  }
  fclose(infile);
  nstar=icount;
  
  /* detect unique labels */
  icount=1;
  for (i=1; i<nstar+1; i++) // loop on stars
    {
      if(i == 1) // first label
	{
	  labu[icount]=lab[i];
	  icount = icount +1;
	}
      if(icount != 1) // other labels
	{
	  for (j=1; j<icount; j++)
	    { // loop on previous labels
	      if(lab[i] == labu[j])
		{
		  // label is not new
		  goto endloop;
		}
	    }
	  // label is brand new
	  icount=icount+1;
	  labu[icount]=lab[i];
	}
    endloop:;
    }
  nlabel=icount;
  
  /* count labels up to the number of frames */
  for (i=1; i<nlabel+1; i++) nlab[i]=0;
  for (i=1; i<nlabel+1; i++) // loop on labels
    {
      for (j=1; j<nstar; j++) // loop on stars
	{
	  if (lab[j]==labu[i]) nlab[i]=nlab[i]+1;
	}
    }
  
  /* order label number */
  icount = 0;
  for (i=1; i<nlabel+1; i++) // loop on labels
    {
      if (nlab[i] == nframe)
	{
	  icount = icount + 1;
	  labu[icount]=labu[i];
	}
    }
  nlabel=icount;
  
  /* write sorted star */
  sprintf(filename1,"%s_sorted",argv[2]);
  outfile = fopen(filename1,"w");
  for(i=1; i<nstar; i++){
    for(j=1; j<nlabel; j++){
      if (lab[i] == labu[j]){
	fprintf(outfile,"%03i %03i %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %i %02i %02i %02i %02i %4.2f \n",	    
		j, frame[i], x[i], y[i], ra[i], dec[i], mag[i], ra0[i], dec0[i], angle[i], yr[i], mon[i], day[i], hr[i], min[i], sec[i]);
      }
    }
  }
  fclose(outfile);  
  
  return 0;
}
