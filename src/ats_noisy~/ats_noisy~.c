/*ats_noisy external by Oscar Pablo Di Liscia*/
/*
odiliscia@unq.edu.ar pdiliscia@iuna.edu.ar
Programa de Investigacion "Sistemas Temporales y Sintesis Espacial en el Arte Sonoro"
http://stseas.web.unq.edu.ar/
odiliscia@unq.edu.ar

Escuela Universitaria de Artes, Universidad Nacional de Quilmes, Argentina
*/
/*odiliscia@unq.edu.ar pdiliscia@iuna.edu.ar*/
/*ats synthesis external for pd*/

#include "../../include/m_pd.h"
#include "ats_noisy.h"
#include <math.h>
/******************************************************************/
/******************************************************************/
static t_class *ats_noisy_tilde_class;
/******************************************************************/
/******************************************************************/
double band_edges[NB_RES+1] = ATSA_CRITICAL_BAND_EDGES;
				
typedef struct _ats_noisy_tilde {
  t_object x_obj;
  t_sample f;
  float    amp;
  INTERP   ra[NB_RES];  
  int      active;
  int	   nb_res;
  RANDI    **randi;
  int	   nb_det;
  OSCIL    **oscil;
  int	   tlen;
  float	   *table;
  float    band_centers[NB_RES];
  float    tstep;
  int      sstep;
  int      scnt;
  float    sr;
  float    sampleperiod; 
  float    amp_max;
} t_ats_noisy_tilde;

/******************************************************************/
void ats_noisy_tilde_destroy(t_ats_noisy_tilde *x);
/******************************************************************/
/******************************************************************/
/******************************************************************/
void ats_noisy_reinit(t_ats_noisy_tilde *x) {
  int i; 

  x->tstep = INTERP_TIME;                /*default interp period*/     
  x->sstep =(int)(x->tstep * x->sr);     /*calculate interp period in samples*/
  x->scnt=0; 
  /*initialize data for amplitude interpolation*/
  memset(&x->ra, 0, x->nb_res * sizeof(INTERP));

  return;
}
/******************************************************************/
void ats_noisy_tilde_any(t_ats_noisy_tilde *x, t_symbol *s, t_int argc, t_atom *argv)
{
  int i;
  float f;
  t_symbol *temp;

  /*stop*/
  if(strcmp(s->s_name, "stop") == 0) {
    x->active=0;
    post(" ats_noisy~: stopped");
    return;
  }
  /*start*/
  if(strcmp(s->s_name, "start") == 0) {
    x->active=1;
    post(" ats_noisy~: active");
    return;
  }
  
  /*reinit*/
  if(strcmp(s->s_name, "reinit") == 0) {
    post(" ats_noisy~: reinitializing");
    ats_noisy_reinit(x);
    return;
  }
  
  /*max amp*/
  if(strcmp(s->s_name, "amp") == 0) {
    post(" ats_noisy~: max amp value %f", x->amp_max);
    return;
  }
  
  return;
}
/*******************************************************************/
/*******************************************************************/
t_int *ats_noisy_tilde_perform(t_int *w)
{
  t_ats_noisy_tilde *x = (t_ats_noisy_tilde *)(w[1]);
  t_sample  *out =     (t_sample *)(w[2]);
  int          n  =    (int) (w[3]);
  int i;
  float output;
  
  //clear output buffer 
  //memset( out , 0., n *sizeof( t_float ));
	
  if(x->active) {
    while(n--) {
      output=0.;
      if (x->scnt) {
	for(i=0; i<x->nb_res; ++i) {
	  if(x->ra[i].curr !=0.)
	    output+= x->ra[i].curr*oscil(x->oscil[i])*randi(x->randi[i]);
	  if(x->ra[i].step != 0.)
	    x->ra[i].curr += x->ra[i].step;
	}
	--x->scnt;
      }
      else {
	for(i=0; i<x->nb_res; ++i) {
	  if(x->ra[i].curr !=0.)
	    output+=x->ra[i].curr*oscil(x->oscil[i])*randi(x->randi[i]);
	}
      }
      *out++ =output;
      x->amp_max=(fabs(output) > x->amp_max ? fabs(output) : x->amp_max);	
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  else{
    while(n--) {
      *out++ =0.;
    }
  }  

  return (w+4);
}
/******************************************************************/
void ats_noisy_tilde_dsp(t_ats_noisy_tilde *x, t_signal **sp)
{ 
  dsp_add(ats_noisy_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}
/******************************************************************/
void noise_amp(t_ats_noisy_tilde *x, t_symbol *s, int argc, t_atom *argv){
	
  t_symbol *temp;
  int i;
	
  for(i = 0; i < argc; i++){
    temp = atom_getsymbol(&argv[i]);
    if(strcmp(temp->s_name, "float") == 0){
      x->ra[i].val  = atom_getfloat(&argv[i]);
      x->ra[i].next = x->ra[i].val;
      x->ra[i].step = (x->ra[i].next - x->ra[i].curr) / (float)x->sstep; /*compute step*/
    }
  }
  x->scnt=x->sstep;
	
  return;
}
/******************************************************************/
static void interp(t_ats_noisy_tilde *x, t_floatarg n)
{
    
  x->tstep =( n > 0.0 ? n*0.001 : x->sampleperiod);  	/*convert msecs. to secs.*/     
  x->sstep =( n > 0.0 ? (int)(x->tstep * x->sr) : 1); /*calculate time step in samples*/
  post(" ats_noisy: interpolation time changed to %f msecs., %d samples ", n, x->sstep);
	
  return;
}
/******************************************************************/
void *ats_noisy_tilde_new(void)
{ 
  int i;
  
  t_ats_noisy_tilde *x = (t_ats_noisy_tilde *)pd_new(ats_noisy_tilde_class);
  x->active=0;
  x->randi=NULL;
  x->oscil=NULL;
  x->amp_max=0.;
    
  /*initialize values for interpolation*/
  x->tlen=TLEN;
  x->sr  =sys_getsr();
  post("Actual Sampling Rate is: %f", x->sr);
  x->sampleperiod= 1./x->sr;
  x->tstep = INTERP_TIME;  				 /*default interp period*/     
  x->sstep =(int)(x->tstep* x->sr);      /*calculate interp period in samples*/
  x->scnt=0; /*very important in the beggining!*/
  
  x->nb_res=x->nb_det=NB_RES;	/*for noise ony*/
  /*compute sine table for lookup oscilators*/
  x->table=make_sine_table(1., 0., x->tlen);
  /*allocate memory for each RANDI and OSCIL unit generators*/
  x->randi=(RANDI**)getbytes(sizeof(RANDI*)*x->nb_res);
  x->oscil=(OSCIL**)getbytes(sizeof(OSCIL*)*x->nb_res);

  for(i=0; i < x->nb_res; ++i) {
    /*use the same loop to compute the critical band frequency centers*/
    x->band_centers[i]= band_edges[i] + ((band_edges[i+1]-band_edges[i])*0.5);
    /*allocate and set randi UGs*/
    x->randi[i]=(RANDI*)set_randi((int)x->sr, (band_edges[i+1]-band_edges[i]));
    /*allocate and set oscil UGs*/
    x->oscil[i]=(OSCIL*)set_oscil(x->tlen, 0.0, (int)x->sr, x->table);
    /*as for noise only the frequency of these will not change, we set their increment here*/
    x->oscil[i]->incr=F2INCR(x->band_centers[i], x->oscil[i]->lsr);
  }
  /*initialize data for amplitude interpolation*/
  memset(&x->ra, 0, x->nb_res * sizeof(INTERP));
  /*control inlet for a list of the noise bands amplitudes*/
  inlet_new(&x->x_obj,&x->x_obj.ob_pd,gensym("list"),gensym("noise_amp"));
  /*control inlet for the interpolation time in msecs.*/
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("interp"));
  /* AUDIO OUTLETS*/
  outlet_new(&x->x_obj, &s_signal);
    
  return (void *)x;
}
/******************************************************************/
/******************************************************************/
void ats_noisy_tilde_setup(void) {

  ats_noisy_tilde_class = class_new(gensym("ats_noisy~"),
				    (t_newmethod)ats_noisy_tilde_new,
				    (t_method)ats_noisy_tilde_destroy, sizeof(t_ats_noisy_tilde),
				    CLASS_DEFAULT, 0);	

  class_addanything(ats_noisy_tilde_class, ats_noisy_tilde_any);
  class_addmethod(ats_noisy_tilde_class,(t_method)ats_noisy_tilde_dsp, gensym("dsp"),0);
  class_addmethod(ats_noisy_tilde_class,(t_method)noise_amp, gensym("noise_amp"),A_GIMME,0);
  class_addmethod(ats_noisy_tilde_class,(t_method)interp, gensym("interp"), A_FLOAT, 0);
  
  post("\n ats_noisy~ external, Pablo Di Liscia, \n UNQ, Argentina, 2013");

  return;
}
/******************************************************************/
void ats_noisy_tilde_destroy(t_ats_noisy_tilde *x) 
{
  int i;

  for(i=0; i < x->nb_res; ++i) {
    freebytes(x->randi[i], sizeof(RANDI));
    freebytes(x->oscil[i], sizeof(OSCIL));
  }
  freebytes(x->randi, sizeof(RANDI*)*x->nb_res);
  freebytes(x->oscil, sizeof(OSCIL*)*x->nb_det);
  freebytes(x->table, sizeof(float)*(x->tlen+1));
  
  return;
}
/******************************************************************/
/******************************************************************/
/**************DEFINITIONS OF AUDIO FUNCTIONS *********************/
/********************************************************/
float *make_sine_table(float amp, float phase, int length)
{
  float *tp=NULL;
  int i;
  double theta=0.;
  double incr = (double)PI2 / ((double)(length));
  
  tp= (float *) malloc((length+1) * sizeof(float));
  
  theta+= (double)phase * PI2; /*phase value from 0 to 1*/
  
  for(i=0; i < length; i++) {
    tp[i] = (float)sin(theta)*(double)amp;
    theta +=incr;
  }
  tp[length]=tp[0]; /*extended guard point*/
  
  return(tp);
}
/******************************************************************/
OSCIL *set_oscil(int L, float phase, int sr, float *table)
{
  OSCIL *oscil;

  oscil=(OSCIL*) malloc(sizeof(OSCIL));

  oscil->L=(double)L;
  oscil->incr =0.;
  oscil->phase=(double)phase;
  oscil->index=(double)phase*oscil->L; /*initial phase*/
  oscil->lsr  =oscil->L/(double)sr;    /*to compute increment from frequency*/ 
  oscil->table=table;
  return(oscil);
}
/********************************************************/
/*non interpolating table lookup oscilator for static frequency*/
/*increment must be set externally in the beggining to save computation time */
inline float oscil(OSCIL *oscil)
{
  float sample=0.;
  while ( oscil->index >= oscil->L ) oscil->index -= oscil->L ;
  while ( oscil->index < 0. )        oscil->index += oscil->L ;

  sample=oscil->table[(int)oscil->index];
  oscil->index+=oscil->incr;

  return(sample);
}
/********************************************************************/
/********************************************************************/
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
  randi->siz_d_sr= randi->size / (float)randi->sr;
  randi->incr    =(float)randi->freq * randi->siz_d_sr;

  return(randi);
}
/********************************************************************/
/******************************************************************/
/*
  randi outputs random numbers in the range of 1, -1 
  getting a new number at frequency freq and linearly 
  interpolating the intermediate values.
  In this version, frequency changes are NOT acknowledged to save computation time 
*/
inline float randi(RANDI *radat)
{
  float samp;
  
  if(radat->freq == (float)radat->sr) { //just white noise
    return( (float)rand()/ (float)RAND_MAX);
  }

  if(radat->cnt >= radat->size) { /*new segment*/
    radat->a1  = radat->a2;     /*target is now start*/
    radat->a2  = (float)rand(); /*get a new target value randomly*/
    radat->cnt = 0.;
    /*precompute segment data to save audio rate computation time*/
    radat->dif =  radat->a2 - radat->a1;  
    radat->dif_d_siz= radat->dif / radat->size;
    radat->siz_d_sr = radat->size / (float)radat->sr;
  }
  samp=(float)(radat->dif_d_siz * radat->cnt) + radat->a1;
  
  radat->cnt+=radat->incr;
  
  return(1. - ((float)(samp /RAND_MAX) * 2.)); 
}
