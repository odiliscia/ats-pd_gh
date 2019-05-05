/* ------------------------ ats_sinnoi~ 0.1 ----------------------------- */
/*
ats_sinnoi~ external by Oscar Pablo Di Liscia

odiliscia@unq.edu.ar pdiliscia@iuna.edu.ar
Programa de Investigacion "Sistemas Temporales y Sintesis Espacial en el Arte Sonoro"
http://stseas.web.unq.edu.ar/
odiliscia@unq.edu.ar
Escuela Universitaria de Artes, Universidad Nacional de Quilmes, Argentina
  
  This external was programed based on the osckbank~ external by Richie Eakin  reakinator@gmail.com 10-15-2007
  https://github.com/pd-l2ork/pd/tree/master/externals/oscbank~
  
  The modifications made by Pablo Di Liscia, to connect this external with the atsread external 
  (by Alex Norman) and to allow the deterministic plus the residual parts synthesis of ATS files.
  For a full explanation of the ATS analysis technique, see:
  http://wiki.dxarts.washington.edu/groups/general/wiki/39f07/attachments/55bd6/ATS_theory.pdf
  See also the PD examples included in this package.
 
  Basically, the modifications by Pablo Di Liscia include the addition of a bank with as many "randi" units
  (which produce linearly interpolated random values in the range of 1 to -1 at a frequeny set by the user) 
  and an inlet to retrieve the residual RMS power data for each partial, which must be sent by atsread as a 
  succession of float values.  This must be computed previously and sent by the atsread external 
  once a file having both residual and deterministic part is loaded by it. The residual part, in this case, 
  is synthesized using randi units as sources, but the output of each one is multiplied by the output 
  of each deterministic trajectory and scaled by the RMS of the noise computed for each partial. 
  The deterministic part is synthesized as varying frequency and amplitude sine waves. 
  An extra audio outlet was added to the external as well, in order to have individual outputs of both, 
  the residual and the deterministic parts. This allows the user to further mix them with individual amplitude 
  scaling, if desired.

*/
/* ----------------------------------------------------------------*/
#include "../../include/m_pd.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef NT
#pragma warning( disable : 4244 )
#define inline
#endif

#define WAVETABLESIZE 65536 //2^16
#define DEFAULT_NPARTIALS 100
#define DEFAULT_interp_incr 0.0045 //per sample, this is 20 ms @ 44k sr 

static t_class *ats_sinnoi_class;

//RANDI UGen data structure by Pablo Di Liscia
typedef struct { 
  float size; 		//size of the frame in samples this should be sr/freq.
  float a1;   		//current amplitude value
  float a2;   		//next    amplitude value
  float cnt;  		//sample  position counter
  float dif;       	//to save time computing a2-a1
  float dif_d_siz;      //to save time computing a2-a1 / size
  float sr; 		//sampling rate  
  float freq;		//frequency
} RANDI;

//t_partial represents one partial member in the bank
typedef struct _partial
{
  int   index;
  float fCurr;
  float freq;
  float fIncr;
  float aCurr;
  float amp;
  float aIncr;
  float rCurr;
  float res;
  float rIncr;
  float phase;
  unsigned long  nInterp;
} t_partial;

typedef struct _ats_sinnoi
{
  t_object x_obj; 	
  float    *wavetable;	
  int      wavetablesize;
  int      got_a_table;
  t_partial *pBank;
  RANDI	 **randi;
  float    infreq;	
  float    inamp;
  float    inres;
  float    sampleRate;
  float    sampleperiod; 
  float    interp_incr;
  long     interpSamples;
  int      sp;
  int      nPartials;
  float    bwscal;
} t_ats_sinnoi;

/*RANDI Audio Functions Prototypes by Pablo Di Liscia*/
RANDI *set_randi(int sr, float freq);
float randi(RANDI *radat);
///////////////////////////////////////////////////////////////////////////////////////
/*----- Interpolation Time -----
  milleseconds to interpolate over; so samples = (n*SR)/1000
  divide only when converting the interp time to samples(here),since it 
  is only used as a denominator to find the increment proportion:
  SP= 1/SR, 1/(n*SR/1000) = (1000*SP)/n
*/
///////////////////////////////////////////////////////////////////////////////////////
static void ats_sinnoi_bwscal(t_ats_sinnoi *x, t_floatarg n)
{
    
  if(n <= 0. || n >= 1.) {
    post("ATS_SINNOI~: wrong value for bandwidth scaling, using default value (0.1)");
    n=0.1;
    return;
  }
  x->bwscal = n;
  post("ATS_SINNOI~: bandwidth scaler changed (%f)", x->bwscal);
}
///////////////////////////////////////////////////////////////////////////////////////
static void ats_sinnoi_interpMs(t_ats_sinnoi *x, t_floatarg n)
{
    
  if(n > 0) x->interp_incr =(1000* x->sampleperiod)/ n ;
  else x->interp_incr = x->sampleperiod; 
    
  x->interpSamples = (unsigned long)((n *.001) * x->sampleRate);
}
///////////////////////////////////////////////////////////////////////////////////////
static void ats_sinnoi_nPartials(t_ats_sinnoi *x, t_floatarg n)
{
  int i;
  x->pBank = (t_partial *)resizebytes( x->pBank, x->nPartials * sizeof(t_partial), n * sizeof(t_partial));
	
  for(i =0; i < x->nPartials; i++) {
    freebytes(x->randi[i], sizeof(RANDI));
  }
  freebytes(x->randi, sizeof(RANDI*)*x->nPartials);
  x->nPartials = n;
   
  x->randi=(RANDI**)getbytes(x->nPartials * sizeof(RANDI*));
  for(i = 0; i < x->nPartials; i++) {
    x->randi[i]=set_randi(x->sampleRate, 10.);
    x->pBank[i].res = 0.;
    x->pBank[i].rCurr = 0.;
  }
	
  post("ATS_SINNOI~: max partials: %d", x->nPartials);
}
///////////////////////////////////////////////////////////////////////////////////////
static void ats_sinnoi_index(t_ats_sinnoi *x, t_floatarg in)
{
  int i, iindex;
  iindex = (int)in;
  t_partial *bank = x->pBank; 
	    
  if( iindex < 0) {
    error("negative index rejected");
    return;
  }

  for(i =0; i < x->nPartials; i++)
    {
      if( bank[i].index == iindex)
	{//recalculate increment slope from current interpolated positions and update goal
	  if(bank[i].aCurr == 0) bank[i].aCurr = 0.0000001;
	  if(bank[i].rCurr == 0) bank[i].rCurr = 0.0000001;
	  bank[i].fIncr = (x->infreq - bank[i].fCurr) * x->interp_incr;
	  bank[i].aIncr = (x->inamp - bank[i].aCurr) * x->interp_incr;
	  bank[i].rIncr = x->inres==0.0 ? 0.0 :(x->inres - bank[i].rCurr) * x->interp_incr;
	  
	  bank[i].freq = x->infreq;
	  bank[i].amp  = x->inamp;
	  bank[i].res  = x->inres;
	  bank[i].nInterp = x->interpSamples;
	  return;
	}
    } //end continuing partial

  //new partial, see if there is an empty slot for the new partial
  for(i =0; i < x->nPartials; i++)
    {
      if(bank[i].aCurr == 0)
	{ //new partial, only ramp amp from zero, 
	  bank[i].index = iindex;
	  bank[i].fCurr = x->infreq;
	  bank[i].fIncr = 0;
	  bank[i].freq = x->infreq;
	  bank[i].amp = x->inamp;
	  bank[i].res = x->inres;
	  bank[i].nInterp = x->interpSamples;
	  bank[i].aCurr = 0.0000001; 
	  bank[i].rCurr = 0.0000001;
	  bank[i].aIncr = x->inamp * x->interp_incr;
	  bank[i].rIncr = x->inres==0.0 ? 0.0 : x->inres * x->interp_incr;
	  return;
	}
    } //end new partial for
 
  //ats_sinnoi is full, steal oldest partial (creates a pop) and  ramp amp from zero
  bank[x->sp].index = iindex;
  bank[x->sp].fCurr = x->infreq;
  bank[x->sp].fIncr = 0;
  bank[x->sp].freq  = x->infreq;
  bank[x->sp].amp   = x->inamp;
  bank[x->sp].nInterp = x->interpSamples;
  bank[x->sp].aCurr = 0.0000001;
  bank[x->sp].rCurr = 0.0000001;
  bank[x->sp].aIncr = x->inamp * x->interp_incr;
  bank[x->sp].rIncr = x->inres==0.0 ? 0.0 : x->inres * x->interp_incr;  
  x->sp++;
  if(x->sp == x->nPartials) x->sp = 0;
}
///////////////////////////////////////////////////////////////////////////////////////
static void ats_sinnoi_table(t_ats_sinnoi *x, t_symbol *tablename)
{
  if(!x->got_a_table) free(x->wavetable);

  t_garray *a;

  if (!(a = (t_garray *)pd_findbyclass(tablename, garray_class)))
    pd_error(x, "%s: no such array", tablename->s_name);
  else if (!garray_getfloatarray(a, &x->wavetablesize, &x->wavetable))
    pd_error(x, "%s: bad template for tabread", tablename->s_name);
  else //table exists
    {
      post("wavetablesize: %d", x->wavetablesize );
    }
  x->got_a_table = 1;
}
///////////////////////////////////////////////////////////////////////////////////////
static void ats_sinnoi_print(t_ats_sinnoi *x)
{
  t_partial *bank = x->pBank; 
    
  post("#.  Index,  Freq,  Amp");
  int i;
  //for every partial
  for(i=0; i < x->nPartials; i++) {
    if(bank[i].aCurr) {
      post("%d. index: %d,freq: %f,amp: %f", i, bank[i].index, 
	   bank[i].freq,  bank[i].amp );
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////
static void ats_sinnoi_reset(t_ats_sinnoi *x)
{
  int i;
  /*just mute the output by setting all amplitudes to zero*/
 post("ats_sinnoi~: reseting...");
  for(i=0; i < x->nPartials; i++) {
    x->pBank[i].amp = 0.;
    x->pBank[i].res = 0.;
    x->pBank[i].aCurr = 0.;
    x->pBank[i].rCurr = 0.;
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////
static t_int *ats_sinnoi_perform(t_int *w)
{
  t_ats_sinnoi *x = (t_ats_sinnoi *)(w[1]);
  t_float *outd = (t_float *)(w[2]);
  t_float *outr = (t_float *)(w[3]);
  t_int n = (t_int)(w[4]);
  t_int i, sample;    
  t_float phaseincrement;
  t_float sample_sum, freq, amp, osc_sig;
  t_int	lookup;
  t_partial *bank = x->pBank; 
  float lowscal= x->bwscal / 4.;

  //clear output buffers
  memset( outd , 0., n *sizeof( t_float ));
  memset( outr , 0., n *sizeof( t_float ));
	
  for(i=0; i < x->nPartials; i++) { //for every partial 
   if(bank[i].aCurr != 0.0) {
    for(sample = 0; sample < n; sample++) {//and every sample..
      if(bank[i].nInterp > 0) {
	bank[i].fCurr += bank[i].fIncr;
	bank[i].aCurr += bank[i].aIncr;
	bank[i].rCurr += bank[i].rIncr;
	--bank[i].nInterp;
      }
      else {
	bank[i].fCurr = bank[i].freq;
	bank[i].aCurr = bank[i].amp;
	bank[i].rCurr = bank[i].res;
      }
      // get the phase increment freq = cyc/sec,
      //sr = samp/sec, phaseinc = cyc/samp = freq/sr = freq * sampleperiod
      phaseincrement = bank[i].fCurr * x->sampleperiod;
      bank[i].phase += phaseincrement;
      while(bank[i].phase >= 1.0f) //..and wrap
	bank[i].phase -= 1.0f;
      while(bank[i].phase < 0.0f)
	bank[i].phase += 1.0f;
		
      lookup  = (int)(x->wavetablesize * bank[i].phase);
      x->randi[i]->freq = (bank[i].fCurr < 500.? bank[i].fCurr * lowscal : bank[i].fCurr * x->bwscal);
      osc_sig = *(x->wavetable + lookup);
      *(outd+sample) += osc_sig * bank[i].aCurr; 
      *(outr+sample) += osc_sig * bank[i].rCurr * randi(x->randi[i]); 
    }//end for samples
   } //end if x->index
  }//end for partials
  return (w+5);
}
///////////////////////////////////////////////////////////////////////////////////////
static void ats_sinnoi_dsp(t_ats_sinnoi *x, t_signal **sp)
{
  x->sampleRate =  sp[0]->s_sr;
  x->sampleperiod = 1 / x->sampleRate;
  dsp_add(ats_sinnoi_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}
///////////////////////////////////////////////////////////////////////////////////////
static void *ats_sinnoi_new(void)
{
  t_ats_sinnoi *x = (t_ats_sinnoi *)pd_new(ats_sinnoi_class);
    
  float twopi, size;
  int i;

  outlet_new(&x->x_obj, gensym("signal"));
  outlet_new(&x->x_obj, gensym("signal"));
  floatinlet_new(&x->x_obj, &x->infreq);
  floatinlet_new(&x->x_obj, &x->inamp);
  floatinlet_new(&x->x_obj, &x->inres);
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("interp"));

  //hardcoded because dsp hasn't been turned on yet
  //prevents divide by zero in ats_sinnoi_index()
  x->sampleRate = sys_getsr(); //48000;
  x->sampleperiod = 1 / x->sampleRate;
  ats_sinnoi_interpMs( x, 5.0); 

  x->got_a_table = 0;
  x->sp = 0;
  x->nPartials = DEFAULT_NPARTIALS;
  x->pBank = (t_partial *)getbytes( x->nPartials * sizeof(t_partial));
  memset(x->pBank, 0, x->nPartials * sizeof(t_partial));
	
  x->randi=(RANDI**)getbytes(x->nPartials * sizeof(RANDI*));
  for(i = 0; i < x->nPartials; i++) {
    x->randi[i]=set_randi(x->sampleRate, 10.);
  }
	   
  twopi = 8.0f * atan(1.0f);
  x->wavetablesize = WAVETABLESIZE;
  float *sinewave;
  sinewave = (t_float *)malloc(x->wavetablesize * sizeof(t_float));
  for(i = 0; i < x->wavetablesize; i++) {
    sinewave[i] = sin(twopi * (float)i/ x->wavetablesize);
  }
  x->wavetable = &sinewave[0];    
  x->bwscal=0.1;
  x->inres=0.0;

	
  return (x);
}
///////////////////////////////////////////////////////////////////////////////////////
static void ats_sinnoi_free(t_ats_sinnoi *x)
{
  free(x->pBank);
  free(x->randi);
  if(!x->got_a_table)
    free(x->wavetable);
}
///////////////////////////////////////////////////////////////////////////////////////
void ats_sinnoi_tilde_setup(void)
{
  ats_sinnoi_class = class_new(gensym("ats_sinnoi~"),(t_newmethod)ats_sinnoi_new,\
			       (t_method)ats_sinnoi_free,sizeof(t_ats_sinnoi), 0, A_DEFFLOAT, 0);
  class_addfloat(ats_sinnoi_class, ats_sinnoi_index);
  class_addmethod(ats_sinnoi_class, (t_method)ats_sinnoi_table, gensym("table"), A_SYMBOL);
  class_addmethod(ats_sinnoi_class, (t_method)ats_sinnoi_interpMs, gensym("interp"), A_FLOAT, 0);
  class_addmethod(ats_sinnoi_class, (t_method)ats_sinnoi_dsp, gensym("dsp"), (t_atomtype)0);
  class_addmethod(ats_sinnoi_class, (t_method)ats_sinnoi_print, gensym("print"), 0);
  class_addmethod(ats_sinnoi_class, (t_method)ats_sinnoi_reset, gensym("reset"), 0);
  class_addmethod(ats_sinnoi_class, (t_method)ats_sinnoi_nPartials, gensym("partials"), A_FLOAT, 0);
  class_addmethod(ats_sinnoi_class, (t_method)ats_sinnoi_bwscal, gensym("bwscal"), A_FLOAT, 0);
  post(" ATS_SINNOI~ external, V. 0.1, by Richie Eakin and Pablo Di Liscia, 2012");
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
/*RANDI Audio functions definitions by Pablo Di Liscia*/
/********************************************************/
RANDI *set_randi(int sr, float freq)
{
  RANDI *randi;

  randi=(RANDI*)malloc(sizeof(RANDI));

  randi->size= (float)((float)sr / freq);
  randi->a1  = (float)rand();
  randi->a2  = (float)rand();
  randi->cnt =  0.;
  randi->freq=  freq;
  randi->sr  =  sr;
  /*precompute first segment data to save audio rate computation time*/
  randi->dif =  randi->a2 - randi->a1; 
  randi->dif_d_siz= randi->dif / randi->size;

  return(randi);
}
///////////////////////////////////////////////////////////////////////////////////////
/*
  randi outputs random numbers in the range of 1 and -1 
  getting a new number at frequency freq and linearly 
  interpolating the intermediate values.
*/
float randi(RANDI *radat)
{
  float samp;
  
  if(radat->freq == (float)radat->sr) { //just white noise
    return( (float)rand()/ (float)RAND_MAX);
  }

  if(radat->cnt++ >= radat->size) { /*new segment*/
    radat->a1  = radat->a2;     /*target is now start*/
    radat->a2  = (float)rand(); /*get a new target value randomly*/
    radat->cnt = 0.;
    /*precompute segment data to save audio rate computation time*/
    radat->dif =  radat->a2 - radat->a1;  
    radat->size= (float) ((float)radat->sr / radat->freq);
    radat->dif_d_siz= radat->dif / radat->size;
  }
  samp=(float)(radat->dif_d_siz * radat->cnt) + radat->a1;
  
  return(1. - ((float)(samp /RAND_MAX) * 2.)); //normalize output
}
/********************************************************************/
