
#define MAIN
#include "/P/hmark/headers/mysac.h"
#include "/P/hmark/headers/sac_db.h"
#include "/P/hmark/headers/koftan.h"
#define SLEN 400000

// This version does not need prior minV and maxV; the values are set to 1.7 and 4 km/s
// Period range is ~4-50 seconds

/* Modified version of spectral_snr.c code to include Misha's filter for faster narrow-band filtering for larger data sets. This version does not write out the narrow-band filtered waveforms (as SAC files) */
// Note that this version computes maximum values ***based on the envelope function*** (filter4_cv)

// script neatened 06.22.2020 (hfm)


/* FUNCTION PROTOTYPES */
SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, long nmax);
void write_sac (char *fname, float *sig, SAC_HD *SHD);
void filter4cv_(double *f1, double *f2, double *f3, double *f4, int *npow, double *dt, 
	      int *n, float seis_in[], float seis_out[]);
int get_snr(char *fname);


/*--------------------------------------------------------------------------*/
SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, long nmax)
/*--------------------------------------------------------------------------*/
/* function to read sac files given the name, fname. The function outputs the time signal to the pointer sig, fills the header SHD, if the signal has fewer than nmax points */
{
  FILE *fsac;

  if((fsac = fopen(fname, "rb")) == NULL) {
    printf("could not open sac file to read%s \n", fname);
    exit(1);
  }

  if ( !fsac ) {
    /*fprintf(stderr,"file %s not find\n", fname);*/
    return NULL;
  }

  if ( !SHD ) SHD = &SAC_HEADER;

  fread(SHD,sizeof(SAC_HD),1,fsac);
  if ( SHD->npts > nmax ) {
    fprintf(stderr,"ATTENTION !!! %s npts is limited to %d.\n",fname,nmax);
    SHD->npts = nmax;
  }

  fread(sig,sizeof(float),(int)(SHD->npts),fsac);
  fclose (fsac);

/*-------------  calcule de t0  ----------------*/
   {
     int eh, em ,i;
     float fes;
     char koo[9];

     for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
     koo[8] = '\0';

     SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
     SHD->nzsec + SHD->nzmsec*.001;

     sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

     SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

   return SHD;
}


/*--------------------------------------------------------------------------*/
void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*--------------------------------------------------------------------------*/
{
  FILE *fsac;
  int i;
  if((fsac = fopen(fname, "wb"))==NULL) {
    printf("could not open sac file to write\n");
    exit(1);
  }

  if ( !SHD ) {
    SHD = &SAC_HEADER;
  }

  SHD->iftype = (long)ITIME;
  SHD->leven = (long)TRUE;
  SHD->lovrok = (long)TRUE;
  SHD->internal4 = 6L;
  SHD->depmin = sig[0];
  SHD->depmax = sig[0];

  for ( i = 0; i < SHD->npts ; i++ ) {
    if ( SHD->depmin > sig[i] ) {
      SHD->depmin = sig[i];
    }
    if ( SHD->depmax < sig[i] ) {
      SHD->depmax = sig[i];
    }
   }

  fwrite(SHD,sizeof(SAC_HD),1,fsac);
  fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);

  fclose (fsac);
}


float sig0[SLEN];
float sig1[SLEN];
char fname[300];
/*--------------------------------------------------------------*/
int get_snr(char *fname )
/*--------------------------------------------------------------*/
{
  FILE *fp1, *fp2, *fd, *fsp, *fstep;
  int k,i, iswap, cnt, *ih, *nsamples;
  int nf = 20;
  int npow;
  float f_cent;
  static int n;
  double dt, delta, multfac, addfac;
  double maxP= 50.0;
  double minP= 4.0  ;
  double per[nf],f[nf],snr[nf],snr1[nf];
  double b, e, fb, fe, step,cve;
  double dist, minV, maxV, minT, maxT, window, signalmax, noisemax,noisemax1;
  double maxT1,minT1;
  double f1, f2, f3, f4;
  double sum,sum1;
  static float seis_in[300000],seis_out[4000000];
  float hed[158],sig[nf],noise[nf],noise1[nf];
  char name[160], fname1[300], fname2[300];
  int kk,kk1;

  fstep = fopen("step.txt","w");

  if ( read_sac (fname, sig0, &SAC_HEADER, SLEN) == NULL )
    {
      fprintf(stderr,"File %s not found\n", fname);
      return 0;
    }

  fb=1.0/maxP;
  fe=1.0/minP;
  step=(log(fb)-log(fe))/(nf-1);

/* Calculate the frequency steps */
  for(k=0;k<nf;k++)
    {
      f[k]=exp(log(fe)+k*step);
      per[k]=1/f[k];
    }

  b=SAC_HEADER.b;
  e=SAC_HEADER.e;
  cve = e;
  dist = SAC_HEADER.dist;
  if (dist!=dist) dist=-1;

  if (dist < 0)
    {
    fprintf(stderr,"dist: %g. There is problem with the dist.\n",dist);
    return 0;
    }

  minV = 1.7;
  maxV = 4.0;

/* Call to Misha's filter, filter4_cv.f */
/* ==========================================================
   Parameters for filter4_cv function:
   Input parameters:
   f1,f2   - low corner frequences, f2 > f1, Hz, (double)
   f3,f4   - high corner frequences, f4 > f3, Hz, (double)
   npow    - power of cosine tapering,  (int)
   dt      - sampling rate of seismogram in seconds, (double)
   n       - number of input samples, (int)
   seis_in - input array length of n, (float)
   Output parameters:
   seis_out - output array length of n, (float)
   output is the envelope function!!!!!!!!!!!!!!!!!!!!!!!!
   ==========================================================
*/

  for(k=1;k<nf-1;k++)
    {
      minT1 = dist/maxV - 6*per[k];
      maxT1 = dist/maxV - 1*per[k];
      if (minT1 < 0) minT1 = 0;
      if (maxT1 < minT1) maxT1 = minT1 + 3*per[k];

      minT = dist/maxV - per[k];
      maxT = dist/minV + 2*per[k];
      if(minT<b+10.)
          minT=b+10.;
      if( maxT> ( e - 400) )
          maxT= (e - 400);
      window=maxT-minT;

      f2=f[k+1];
      f3=f[k-1];

      f_cent = 2/(f2+f3); 

      addfac= .01;
      f1=f2 - addfac;
      f4=f3 + addfac; 

      fprintf(fstep,"%d %g %g %g %g %g %g %g\n",k,f1,f2,f3,f4,1/f2,1/f[k],1/f3);

      dt = SAC_HEADER.delta;
      delta = dist;
      n = SAC_HEADER.npts;
      for (i=0;i<n;i++) {
          seis_in[i] = sig0[i];  
          }
      npow=1;
      filter4cv_(&f1,&f2,&f3,&f4,&npow,&dt,&n,seis_in,seis_out);

      for(cnt=0; cnt<n; cnt++) {
        if(seis_out[cnt] > 0) {
          sig1[cnt]=seis_out[cnt];
        }
        else {
          sig1[cnt]=(-1)*seis_out[cnt];
        }
      }

/* loop over sig1 to find maxima within signal and noise time periods */
      signalmax=0;
      noisemax=0.;
      noisemax1=0.;

      for(i=(int)minT/dt;i<maxT/dt;i++)
	{
	  if((double)sig1[i]>signalmax) 
	    signalmax=(double)sig1[i];
	}

      sum = 0.;
      sum1 = 0.;
      kk = 0;
      for(i=(int)((maxT+500)/dt);i<cve;i++)
	{
	   sum = sum+(double)sig1[i]*(double)sig1[i];
//           sum = sum + sig0[i];
           kk ++ ;
	}
   
      kk1 = 0;
      for(i=(int)(minT1/dt);i<(int)(maxT1/dt);i++)
        {
           sum1 = sum1 + (double)sig1[i]*(double)sig1[i];
           kk1++;
        }

      noisemax = pow(sum/kk,0.5);
      noisemax1 = pow(sum1/kk1,0.5);
      snr[k]=signalmax/noisemax;      
      snr1[k]=signalmax/noisemax1;
      sig[k] = signalmax;
      noise[k] = noisemax;
      noise1[k] = noisemax1;
    }

  sprintf(fname1,"%s_snr.cv.p.txt",fname);
  fp2=fopen(fname1,"w");
  for(i=1;i<nf-1;i++)
    {
      fprintf(fp2,"%lf %lf %lf %lf %lf %lf\n",per[i],snr[i],snr1[i],sig[i],noise[i],noise1[i]);
    }
  fclose(fp2);
  fclose(fstep);
  return 1;
}


SAC_DB sdb;
/*--------------------------------------------------------------*/
int main (int argc, char *argv[])
/*--------------------------------------------------------------*/
{
  FILE *ff,*ff1;
  char filename[20];
  int i;

//check input command line arguments
  if(argc != 2) {
    fprintf(stderr,"USAGE: spectral_snr_f [list file]\n");
    fprintf(stderr,"Req.s a list of all cut files for spectral SNR values.\n");
    exit(-1);
  }

  if((ff = fopen(argv[1], "r"))==NULL) {
    fprintf(stderr,"Cannot open file list: %s\n", argv[1]);
    exit(1);
  }
  for(i=0;;i++)
    {
      if(fscanf(ff,"%s",&filename)==EOF)
	break;
      //if (fmod(i,50)==1) 
      fprintf(stderr,"working on %s %d\n",filename,i);
      if(get_snr(filename)==0) continue;
    }
  return 1;
}
