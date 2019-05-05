/* ats-pd version 0.2
 * December, 2014
 * atsread by Alex Norman with modifications by Pablo Di Liscia.

odiliscia@unq.edu.ar pdiliscia@iuna.edu.ar
Programa de Investigacion "Sistemas Temporales y Sintesis Espacial en el Arte Sonoro"
http://stseas.web.unq.edu.ar/
odiliscia@unq.edu.ar
Escuela Universitaria de Artes, Universidad Nacional de Quilmes, Argentina

 * read the README file for info about the opcode. 
 The modifications  to atsread that were included by Pablo Di Liscia are meant basically to better connection
 with the synthesis externals (ats_noisy~, ats_sinnoi~ and oscbank~) and are listed below:
 1-	Inclusion of an extra outlet for the output of the header data of the ATS file. 
 At least knowledge of the type (basically to know whether residual data is present or not), 
 the duration and the number of partials of the opened file is needed in order to properly set 
 the corresponding synthesis units and to send to them the time data. The header data is sent in 
 the form of a list of floats once a file is opened.
 2-	Inclusion of a function to compute the residual noise values corresponding to each partial, 
 for the case in which residual information is present and both, deterministic and residual synthesis, are required . 
 If the file has residual data, then this function is call once it is opened, and the resulting values 
 for each partial at each frame are stored in memory.
 Such C function (Nband_energy_to_res) was taken from the synthesis engine code of the 
 program WATSH (by J. Pampin, O. P. Di Liscia and P. Moss) with improvements made by
 Alex Norman (in ugnorman.c from the Csound program source Code). 
 3-	Inclusion of an extra outlet for the output of the noise values for each partial computed in 2.
 4-	Inclusion of an extra outlet for the output of the index (i.e., the partial number). 
 The partial number if sent as a float value before its amplitude, frequency and noise (if any) data are sent. 
 5-	Modification of the output format of the amplitude, frequency, phase and noise values. 
 These are now sent as independent float values (with their partial index preceding them), 
 instead of as lists of float values. 
*/
#include "../../include/m_pd.h"
#include <stdio.h>
#include <math.h> 
#include <malloc.h>
#include <string.h>

#define ATSA_NOISE_VARIANCE 0.04
#define ATSA_CRITICAL_BAND_EDGES {0.0, 100.0, 200.0, 300.0, 400.0, 510.0, \
      630.0, 770.0, 920.0, 1080.0, 1270.0,				\
      1480.0, 1720.0, 2000.0, 2320.0, 2700.0,				\
      3150.0, 3700.0, 4400.0, 5300.0, 6400.0,				\
      7700.0, 9500.0, 12000.0, 15500.0, 20000.0}
#define ATSA_CRITICAL_BANDS 25
#define ENG_RMS(val, ws) sqrt(val /(ws * (float)ATSA_NOISE_VARIANCE))

typedef struct _atsdataloc {
  double amp;
  double freq;
}	t_atsdataloc;

typedef struct _atshead {
  double  magic;  		/* ats magic number */
  double  sampr;  		/* sampling rate */
  double  frmsz;  		/* frame size in samples */
  double  winsz;  		/* window size in samples */
  double  npartials;    /* number of partials */
  double  nfrms;  		/* number of frames */
  double  ampmax; 		/* max amplitude */
  double  freqmax;  	/*max frequency*/
  double  dur;    		/* duration seconds */
  double  type;   		/* Ats Frame type 1-4 */
}	t_atshead;

static t_class *atsread_class;

typedef struct _atsread {
  t_object  x_obj;
  t_float timepnt;
  int berrflg;	//time pointer bounderror flag
  t_atsdataloc ** data;
  double ** nzdata;
  double ** nzpdata;
  t_atshead atshead;
  double framemult;
  int * outpartials;	//partials to output
  int outnzbands[ATSA_CRITICAL_BANDS];
  int outnz;	//should we output the noise?
  int outsines;	//should we output the sines?
  	
  t_atom * amplist;
  t_atom * freqlist;
  t_atom * reslist;
  t_atom nzlist[ATSA_CRITICAL_BANDS];
  t_atom dlist[3];	
  t_outlet *freq_out, *amp_out, *nz_out, *res_out, *data_out, *ind_out;
  
} t_atsread;


#ifdef BYTESWAP
//byte swaps a double
double bswap(double * swap_me){
  unsigned char *sw_in, sw_out[8];	//for swaping
  sw_in = (unsigned char *)swap_me;
  sw_out[0] = sw_in[7];
  sw_out[1] = sw_in[6];
  sw_out[2] = sw_in[5];
  sw_out[3] = sw_in[4];
  sw_out[4] = sw_in[3];
  sw_out[5] = sw_in[2];
  sw_out[6] = sw_in[1];
  sw_out[7] = sw_in[0];
  return *((double *)sw_out);
}
#endif

//function prototypes
void readdatafile(t_atsread *, const char *);
void Nband_energy_to_res(t_atsread *p, float *maxval);
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
void atsread_timepnt(t_atsread *x, t_floatarg time){
  double mytime = (double)time;
  int i, frame, outpnz;	//outp is the number of partials to output
  double frac, tempamp, tempfreq, tempnz, tempres, framedur= x->atshead.dur / x->atshead.nfrms;
  t_atom *amplist, *freqlist, *nzlist, *reslist;
  amplist = x->amplist;
  freqlist = x->freqlist;
  nzlist = x->nzlist;
  reslist = x->reslist;

  if(x->atshead.magic != 123){
    post("ATSREAD: you must first open an Atsfile");
    return;
  }
  outpnz = 0;	//set the output noiseband count to zero	
  frame = (int)floor(mytime * x->framemult);

  /////////////////////////////////////////////////////////////////////////
  if(frame < 0){
    if(x->berrflg){
      post("ATSREAD: negative time pointers not allowed, setting to zero");
      x->berrflg = 0;
    }
    if(x->outsines == 1){
      for(i = 0; i < x->atshead.npartials; i++){
  			if(x->outpartials[i] == 1){
  	  		outlet_float(x->ind_out,  (float)(i+1));
  	  		outlet_float(x->freq_out, (float)(x->data[0][i]).freq);
  	  		outlet_float(x->amp_out,  (float)(x->data[0][i]).amp);
  	  		if(x->outnz==1 && (int)x->atshead.type > 2) { /*deterministic and noise*/
  	    		outlet_float(x->res_out,(float)x->nzpdata[0][i]);
  	  		}
  			}
      }
    }
    if(x->outnz == 1 && (int)x->atshead.type > 2){
      //we can and should output noise
      	//if(!x->outsines) {  /*only residual*/
  				for(i = 0; i < ATSA_CRITICAL_BANDS; i++){
  	  			if(x->outnzbands[i] == 1){
  	    			tempnz = ENG_RMS((float)x->nzdata[0][i], x->atshead.winsz);
  	    			SETFLOAT(&(nzlist[outpnz]),(float)tempnz);
  	    			outpnz++;
  	  			}
  				}
      	//}
  	 	}	
  	}
  /////////////////////////////////////////////////////////////////////////
  else if(frame >= ((int)(x->atshead.nfrms) - 1)){
    frame = (int)(x->atshead.nfrms) - 1;
    if(x->berrflg){
      post("ATSREAD: timepointer exceeds last frame, setting to last frame");
      x->berrflg = 0;
    }
    if(x->outsines == 1){
      for(i = 0; i < x->atshead.npartials; i++){
  			if(x->outpartials[i] == 1){
  	  		outlet_float(x->ind_out,  (float)(i+1));
  	  		outlet_float(x->freq_out, (float)(x->data[frame][i]).freq);
  	  		outlet_float(x->amp_out,  (float)(x->data[frame][i]).amp);
  	  		if(x->outnz==1 && (int)x->atshead.type > 2) { /*deterministic and noise*/
  	    		outlet_float(x->res_out,(float)x->nzpdata[frame][i]);
  	  		}
  			}
      }
    }
    if(x->outnz == 1 && (int)x->atshead.type > 2){
      //if(!x->outsines) {  /*only residual*/
  			for(i = 0; i < ATSA_CRITICAL_BANDS; i++){
  	  		if(x->outnzbands[i] == 1){
  	    		tempnz = ENG_RMS((float)x->nzdata[frame][i], x->atshead.winsz);
  	    		SETFLOAT(&(nzlist[outpnz]),(float)tempnz);
  	    		outpnz++;
  	  		}
  			}
      //}
    }
  }
  /////////////////////////////////////////////////////////////////////////
	else{
    x->berrflg = 1; 
    frac = (mytime - ((double)(frame))/x->framemult) / framedur;		
    if(x->outsines == 1){
      for(i = 0; i < x->atshead.npartials; i++){
				if(x->outpartials[i] == 1){
	  			tempamp = (x->data[frame][i]).amp + frac * ((x->data[1 + frame][i]).amp - (x->data[frame][i]).amp);
	  			tempfreq = (x->data[frame][i]).freq + frac * ((x->data[1 + frame][i]).freq - (x->data[frame][i]).freq);
	  			outlet_float(x->ind_out,  (float)(i+1));
	  			outlet_float(x->freq_out, (float)tempfreq);
	  			outlet_float(x->amp_out,  (float)tempamp);
	  			if(x->outnz==1 && (int)x->atshead.type > 2) { /*deterministic and noise*/
	    			tempres = x->nzpdata[frame][i] + frac * ((x->nzpdata[1 + frame][i]) - (x->nzpdata[frame][i]));
	    			outlet_float(x->res_out,  (float)tempres);
	  			}
				}
      }
    }
    if(x->outnz == 1 && (int)x->atshead.type > 2){
      //we can and should output noise
      //if(!x->outsines) {  /*only residual*/
				for(i = 0; i < ATSA_CRITICAL_BANDS; i++){
	  			if(x->outnzbands[i] == 1){
	    			tempnz = (x->nzdata[frame][i]) + frac * ((x->nzdata[frame + 1][i]) - (x->nzdata[frame][i]));
	    			tempnz = ENG_RMS(tempnz, x->atshead.winsz);
	    			SETFLOAT(&(nzlist[outpnz]),(float)tempnz);
	    			outpnz++;
	  			}			
				}
      //}
    }
   }
  /////////////////////////////////////////////////////////////////////////
  /*this is the only list output case*/
  if(x->outnz == 1 && /*x->outsines == 0*/(int)x->atshead.type > 2)
    outlet_list(x->nz_out, gensym("list"), outpnz, nzlist);
}
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
void atsread_any(t_atsread *x, t_symbol *s, int argc, t_atom *argv){
  char * filename;
  int i, j, from, to;
  t_symbol * temp;

  if(strcmp(s->s_name, "open") == 0){
    //open data file
    if(argc < 1){
      post("ATSREAD: you need to specify a file name");
      return;
    }
		
    temp = atom_getsymbol(&argv[0]);
		
    if(temp == NULL){
      post("ATSREAD: not a valid filename");
      return;
    }
    readdatafile(x, temp->s_name);
    x->berrflg = 1;
    return;
  }
  if(strcmp(s->s_name, "noise") == 0){
    x->outnz = 1;
    post("ATSREAD: outputing noise, if there is any");
    return;
  }
  if(strcmp(s->s_name, "nonoise") == 0){
    x->outnz = 0;
    post("ATSREAD: not outputing noise");
    return;
  }
  if(strcmp(s->s_name, "sines") == 0){
    x->outsines = 1;
    post("ATSREAD: outputing sines");
    return;
  }
  if(strcmp(s->s_name, "nosines") == 0){
    x->outsines = 0;
    post("ATSREAD: not outputing sines");
    return;
  }
	
  /* BELOW HERE WE DEAL WITH PARTIAL and BAND ADDING/SETTING/REMOVING */
	
  if(x->atshead.magic != 123){
    post("ATSREAD: you must open a file before setting up the partials");
    return;
  }
  if(argc < 1){
    post("ATSREAD: you need to specify at least 1 partial/band to set/add/remove");
    return;
  }
	
  if(strcmp(s->s_name, "set") == 0 || strcmp(s->s_name, "add") == 0){
    //set or add partials to output
    if(strcmp(s->s_name, "set") == 0){
      //clear the data table
      for(i = 0; i < (int)x->atshead.npartials; i++)
	x->outpartials[i] = 0;
    }
		
    for(i = 0; i < argc; i++){
      temp = atom_getsymbol(&argv[i]);
      if(strcmp(temp->s_name, "float") == 0){
	j = (int)atom_getfloat(&argv[i]);
	if(j < 1)
	  j = 1;
	else if(j >= (int)x->atshead.npartials)
	  j = (int)x->atshead.npartials;
	x->outpartials[j - 1] = 1;
      }
      else{
	if(sscanf(temp->s_name, "%i..%i", &from, &to) == 2){
	  if(from > to){
	    j = from;
	    from = to;
	    to = j;
	  }
	  if(from < 1)
	    from = 1;
	  if(to > (int)x->atshead.npartials)
	    to = (int)x->atshead.npartials;
	  for(j = from; j <= to; j++)
	    x->outpartials[j-1] = 1;
	}
	else
	  post("ATSREAD: %s is not a valid range (ie 1..42) or partial number (ie 12)", temp->s_name);
      }
    }
    return;
  }
  if(strcmp(s->s_name, "remove") == 0){
    //remove partials
    for(i = 0; i < argc; i++){
      temp = atom_getsymbol(&argv[i]);
      if(strcmp(temp->s_name, "float") == 0){
	j = (int)atom_getfloat(&argv[i]);
	if(j < 1)
	  j = 1;
	else if(j >= (int)x->atshead.npartials)
	  j = (int)x->atshead.npartials;
	x->outpartials[j - 1] = 0;
      }
      else{
	if(sscanf(temp->s_name, "%i..%i", &from, &to) == 2){
	  if(from > to){
	    j = from;
	    from = to;
	    to = j;
	  }
	  if(from < 1)
	    from = 1;
	  if(to > (int)x->atshead.npartials)
	    to = (int)x->atshead.npartials;
	  for(j = from; j <= to; j++)
	    x->outpartials[j-1] = 0;
	}
	else
	  post("ATSREAD: %s is not a valid range (ie 1..42) or partial number (ie 12)", temp->s_name);
      }
    }
    return;
  }
  if(strcmp(s->s_name, "setnz") == 0 || strcmp(s->s_name, "addnz") == 0){
    //set or add bands to output
    if(strcmp(s->s_name, "setnz") == 0){
      //clear the data table
      for(i = 0; i < ATSA_CRITICAL_BANDS; i++)
	x->outnzbands[i] = 0;
    }
		
    for(i = 0; i < argc; i++){
      temp = atom_getsymbol(&argv[i]);
      if(strcmp(temp->s_name, "float") == 0){
				j = (int)atom_getfloat(&argv[i]);
				if(j < 1)
	  			j = 1;
				else if(j >= ATSA_CRITICAL_BANDS)
	  			j = ATSA_CRITICAL_BANDS;
				x->outnzbands[j - 1] = 1;
      }
      else{
				if(sscanf(temp->s_name, "%i..%i", &from, &to) == 2){
	  			if(from > to){
	    			j = from;
	    			from = to;
	    			to = j;
	  			}
	  		if(from < 1)
	    		from = 1;
	  		if(to > ATSA_CRITICAL_BANDS)
	    		to = ATSA_CRITICAL_BANDS;
	  		for(j = from; j <= to; j++)
	    		x->outnzbands[j-1] = 1;
				}
			else
	  		post("ATSREAD: %s is not a valid range (ie 1..25) or band number (ie 2)", temp->s_name);
      }
    }
    return;
  }
  if(strcmp(s->s_name, "removenz") == 0){
    //remove bands
    for(i = 0; i < argc; i++){
      temp = atom_getsymbol(&argv[i]);
      if(strcmp(temp->s_name, "float") == 0){
	j = (int)atom_getfloat(&argv[i]);
	if(j < 1)
	  j = 1;
	else if(j >= ATSA_CRITICAL_BANDS)
	  j = ATSA_CRITICAL_BANDS;
	x->outnzbands[j - 1] = 0;
      }
      else{
	if(sscanf(temp->s_name, "%i..%i", &from, &to) == 2){
	  if(from > to){
	    j = from;
	    from = to;
	    to = j;
	  }
	  if(from < 1)
	    from = 1;
	  if(to > ATSA_CRITICAL_BANDS)
	    to = ATSA_CRITICAL_BANDS;
	  for(j = from; j <= to; j++)
	    x->outnzbands[j-1] = 0;
	}
	else
	  post("ATSREAD: %s is not a valid range (ie 1..25) or band number (ie 12)", temp->s_name);
      }
    }
    return;
  }
}
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
void *atsread_new(t_symbol *s, int argc, t_atom *argv){
  t_symbol * temp;
  int i;
	
  t_atsread *x = (t_atsread *)pd_new(atsread_class);
  	
  // initalize
  x->berrflg = 1;
  x->data = NULL;
  x->nzdata = NULL;
  x->nzpdata = NULL;
  x->atshead.magic = 0;
  x->outpartials = NULL;
  x->amplist = NULL;
  x->freqlist = NULL;
  x->reslist = NULL;
  x->outnz = 1;
  x->outsines = 1;
	
  //outlets
  x->ind_out = outlet_new(&x->x_obj, &s_float);
  x->freq_out = outlet_new(&x->x_obj, &s_float);
  x->amp_out = outlet_new(&x->x_obj, &s_float);
  x->nz_out = outlet_new(&x->x_obj, &s_float);
  x->res_out = outlet_new(&x->x_obj, &s_float);  /*Added by Oscar Pablo Di Liscia*/
  x->data_out = outlet_new(&x->x_obj, &s_float); /*Added by Oscar Pablo Di Liscia*/

  /* deal with creation time arguments */
	/*
  for(i = 0; i < argc; i++){
    temp = atom_getsymbol(&argv[i]);
    if(strcmp(temp->s_name, "noise") == 0){
      x->outnz = 1;
      post("ATSREAD: outputing noise, if there is any");
    }
    else if(strcmp(temp->s_name, "nonoise") == 0){
      x->outnz = 0;
      post("ATSREAD: not outputing noise");
    }
    else if(strcmp(temp->s_name, "sines") == 0){
      x->outsines = 1;
      post("ATSREAD: outputing sines");
    }
    else if(strcmp(temp->s_name, "nosines") == 0){
      x->outsines = 0;
      post("ATSREAD: not outputing sines");
    }
    else{
      //open data file
      if(temp == NULL){
	post("ATSREAD: no filename requested to open");
      }
      else
	readdatafile(x, temp->s_name);
    }
  }*/
  return (void *)x;
}
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
void atsread_destroy(t_atsread *x){
  int i;
  //clean up allocated data
  if(x->data != NULL){
    for(i = 0; i < (int)(x->atshead).nfrms; i++)
      free(x->data[i]);
    free(x->data);
    x->data = NULL;
  }
  if(x->nzdata != NULL){
    for(i = 0; i < (int)(x->atshead).nfrms; i++)
      free(x->nzdata[i]);
    free(x->nzdata);
    x->nzdata = NULL;
  }
  
  if(x->nzpdata != NULL){
    for(i = 0; i < (int)(x->atshead).nfrms; i++)
      free(x->nzpdata[i]);
    free(x->nzpdata);
    x->nzpdata = NULL;
  }
  
  if(x->outpartials != NULL)
    free(x->outpartials);
  if(x->amplist != NULL)
    free(x->amplist);
  if(x->freqlist != NULL)
    free(x->freqlist);
  if(x->reslist != NULL)
    free(x->reslist);
}
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
void atsread_setup(void) {
  atsread_class = class_new(gensym("atsread"),
			    (t_newmethod)atsread_new, (t_method)atsread_destroy, sizeof(t_atsread), CLASS_DEFAULT, A_GIMME, 0);
	
  class_addfloat(atsread_class, atsread_timepnt);
  class_addanything(atsread_class, atsread_any);
  class_sethelpsymbol(atsread_class, gensym("help-atsread"));
  post("\n ATSREAD external, version 2012 by Alex Norman \n with improvements by Oscar Pablo Di Liscia\n");
}
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
// set up data file to be read
void readdatafile(t_atsread *x, const char * filename){
  FILE *fin;
  double temp;
  float maxval=0.;
  int i,j;
  t_atom *dlist;
  dlist = x->dlist;

#ifdef BYTESWAP
  int swapped = 0;
  int k;
#endif
	
  if((fin=fopen(filename,"rb"))==0){
    post("ATSREAD: cannot open file %s", filename);
    return;
  }
  //delete any previous data
  if(x->data != NULL){
    for(i = 0; i < (int)(x->atshead).nfrms; i++)
      free(x->data[i]);
    free(x->data);
    x->data = NULL;
  }
  if(x->nzdata != NULL){
    for(i = 0; i < (int)(x->atshead).nfrms; i++)
      free(x->nzdata[i]);
    free(x->nzdata);
    x->nzdata = NULL;
  }
  if(x->nzpdata != NULL){
    for(i = 0; i < (int)(x->atshead).nfrms; i++)
      free(x->nzpdata[i]);
    free(x->nzpdata);
    x->nzpdata = NULL;
  }       
  // read the header
  fread(&(x->atshead),sizeof(t_atshead),1,fin);

  //make sure this is an ats file
  if((int)x->atshead.magic != 123){

#ifndef BYTESWAP
    post("ATSREAD - Magic Number not correct in %s", filename);
    fclose(fin);
    return;
#else
    //attempt to byte swap
    if((int)bswap(&x->atshead.magic) != 123){	//not an ats file
      post("ATSREAD - Magic Number not correct in %s", filename);
      fclose(fin);
      return;
    } else {
      post("ATSREAD - %s is of the wrong endianness, byteswapping.", filename);
      swapped = 1;	//true (for later use)
			//byte swap
      for(i = 0; i < (sizeof(t_atshead) / sizeof(double)); i++)
	((double *)(&x->atshead))[i]  = bswap(&(((double *)(&x->atshead))[i]));
    }
#endif
  }
  post("\t %s stats: \n \t\tType: %i, Frames: %i, Partials: %i, Sampling rate: %i Hz, Duration %f sec, Max Freq: %f Hz, Max Amp: %f",
       filename, (int)x->atshead.type, (int)x->atshead.nfrms, (int)x->atshead.npartials, (int)x->atshead.sampr, (float)x->atshead.dur, 
       (float)x->atshead.freqmax, (float)x->atshead.ampmax);

  if((int)x->atshead.type > 4 || (int)x->atshead.type < 1){
    post("ATSREAD - %i type Ats files not supported yet",(int)x->atshead.type);
    x->atshead.magic = 0;
    fclose(fin);
    return;
  }
        
  //allocate memory for a new data table
  x->data = (t_atsdataloc **)malloc(sizeof(t_atsdataloc *) * (int)x->atshead.nfrms);
  for(i = 0; i < (int)x->atshead.nfrms; i++)
    x->data[i] = (t_atsdataloc *)malloc(sizeof(t_atsdataloc) * (int)x->atshead.npartials);
                
  // store the data into table(s)
        
  switch ((int)x->atshead.type){
  case 1 :
    x->outnz=0;
    for(i = 0; i < (int)x->atshead.nfrms; i++){
      //eat time data
      fread(&temp,sizeof(double),1,fin);
      //get a whole frame of partial data
      fread(&x->data[i][0], sizeof(t_atsdataloc),(int)x->atshead.npartials,fin);
#ifdef BYTESWAP
      if(swapped) {	//byte swap if neccesary
	for(k = 0; k < (int)x->atshead.npartials; k++){
	  x->data[i][k].amp  = bswap(&x->data[i][k].amp);
	  x->data[i][k].freq  = bswap(&x->data[i][k].freq);
	}
      }
#endif
    }
    break;
  case 2 :
    x->outnz=0;
    for(i = 0; i < (int)x->atshead.nfrms; i++){
      //eat time data
      fread(&temp,sizeof(double),1,fin);
      //get partial data
      for(j = 0; j < (int)x->atshead.npartials; j++){
	fread(&x->data[i][j], sizeof(t_atsdataloc),1,fin);
	//eat phase info
	fread(&temp, sizeof(double),1,fin);
#ifdef BYTESWAP
	if(swapped) {	//byte swap if nessisary
	  x->data[i][j].amp = bswap(&x->data[i][j].amp);
	  x->data[i][j].freq = bswap(&x->data[i][j].freq);
	}
#endif
      }
    }
    break;
  case 3 :
    //allocate mem for the noise data
    x->nzdata  = (double **)malloc(sizeof(double *) * (int)x->atshead.nfrms);
    x->nzpdata = (double **)malloc(sizeof(double *) * (int)x->atshead.nfrms);
	
    for(i = 0; i < (int)x->atshead.nfrms; i++){
      //allocate this row
      x->nzdata[i]  = (double *)malloc(sizeof(double) * ATSA_CRITICAL_BANDS);
      x->nzpdata[i] = (double *)malloc(sizeof(double) * (int)x->atshead.npartials);
      //eat time data
      fread(&temp,sizeof(double),1,fin);
      //get a whole frame of partial data
      fread(&x->data[i][0], sizeof(t_atsdataloc),(int)x->atshead.npartials,fin);
      //get the noise
      fread(&x->nzdata[i][0], sizeof(double),ATSA_CRITICAL_BANDS,fin);
#ifdef BYTESWAP
      if(swapped) {	//byte swap if nessisary
	for(k = 0; k < (int)x->atshead.npartials; k++){
	  x->data[i][k].freq = bswap(&x->data[i][k].freq);
	  x->data[i][k].amp = bswap(&x->data[i][k].amp);
	}
	for(k = 0; k < ATSA_CRITICAL_BANDS; k++)
	  x->nzdata[i][k] = bswap(&x->nzdata[i][k]);
      }
#endif
    }
    break;
  case 4 :
    //allocate mem for the noise data
    x->nzdata = (double **)malloc(sizeof(double *) * (int)x->atshead.nfrms);
    x->nzpdata = (double **)malloc(sizeof(double *) * (int)x->atshead.nfrms);
    for(i = 0; i < (int)x->atshead.nfrms; i++){
      //allocate noise for this row
      x->nzdata[i] = (double *)malloc(sizeof(double) * ATSA_CRITICAL_BANDS);
      x->nzpdata[i] = (double *)malloc(sizeof(double) * (int)x->atshead.npartials);
      //eat time data
      fread(&temp,sizeof(double),1,fin);
      //get partial data
      for(j = 0; j < (int)x->atshead.npartials; j++){
	fread(&x->data[i][j], sizeof(t_atsdataloc),1,fin);
	//eat phase info
	fread(&temp, sizeof(double),1,fin);
#ifdef BYTESWAP
	if(swapped) {	//byte swap if nessisary
	  x->data[i][j].freq = bswap(&x->data[i][j].freq);
	  x->data[i][j].amp = bswap(&x->data[i][j].amp);
	}
#endif
      }
      //get the noise
      fread(&x->nzdata[i][0], sizeof(double),ATSA_CRITICAL_BANDS,fin);
#ifdef BYTESWAP
      if(swapped) {	//byte swap if nessisary
	for(k = 0; k < ATSA_CRITICAL_BANDS; k++)
	  x->nzdata[i][k] = bswap(&x->nzdata[i][k]);
      }
#endif
    }
    break;
  }
  /*if the file has a residual part, compute the energy of the noise for each partial*/
  /*added by Pablo Di Liscia, 2012*/
  if((int)x->atshead.type == 3 || (int)x->atshead.type == 4) {
    x->outnz=1;
    post("ATSREAD: computing noise energy for each partial");
	Nband_energy_to_res(x, &maxval);	
    post("ATSREAD: max. RMS value found %f", maxval);
  }
 
  // set up the frame multiplier
  x->framemult = (x->atshead.nfrms - (double)1.)/ x->atshead.dur;

  if(x->outpartials != NULL)
    free(x->outpartials);
  x->outpartials = (int *)malloc(sizeof(int) * (int)(x->atshead.npartials));
  for(i = 0; i < (int)x->atshead.npartials; i++)
    x->outpartials[i] = 1;
  for(i = 0; i < ATSA_CRITICAL_BANDS; i++)
    x->outnzbands[i] = 1;

  if(x->amplist != NULL)
    free(x->amplist);
  if(x->freqlist != NULL)
    free(x->freqlist);
  if(x->reslist != NULL)
    free(x->reslist);
  x->amplist  = (t_atom *)malloc(sizeof(t_atom) * (int)(x->atshead.npartials));
  x->freqlist = (t_atom *)malloc(sizeof(t_atom) * (int)(x->atshead.npartials));
  x->reslist  = (t_atom *)malloc(sizeof(t_atom) * (int)(x->atshead.npartials));
  
  /*output duration, number of partials and file type, added by Oscar Pablo Di Liscia, 2012*/
  SETFLOAT(&(dlist[0]),(float)x->atshead.dur);
  SETFLOAT(&(dlist[1]),(float)x->atshead.npartials);
  SETFLOAT(&(dlist[2]),(float)x->atshead.type);
  outlet_list(x->data_out, gensym("list"), 3, dlist);
  
  //close the input file
  fclose(fin);
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
void Nband_energy_to_res(t_atsread *p, float *maxval)
{
    int     i, j, k;
    float   edges[] = ATSA_CRITICAL_BAND_EDGES;
    double  bandsum[ATSA_CRITICAL_BANDS];
    double  partialfreq;
    double  partialamp;
    double  **partialband;
    int     *bandnum;

	*maxval=0.0;
    partialband = (double **) malloc(sizeof(double*)* (int) p->atshead.npartials);
    bandnum = (int *) malloc(sizeof(int)* (int) p->atshead.npartials);

    for (i = 0; i < (int) p->atshead.nfrms; i++) {
      /* init sums */
      memset(bandsum, 0.0, ATSA_CRITICAL_BANDS*sizeof(double));
      /* find sums per band */
      for (j = 0; j < (int) p->atshead.npartials; j++) {
        partialfreq = p->data[i][j].freq; 
        partialamp =  p->data[i][j].amp;   
        for (k = 0; k < ATSA_CRITICAL_BANDS; k++) {
          if ((partialfreq < edges[k + 1]) && (partialfreq >= edges[k])) {
            bandsum[k] += partialamp;
            bandnum[j] = k;
            partialband[j] = &p->nzdata[i][k];
            break;
          }
        }
      }
      /* compute energy per partial */
      for (j = 0; j < (int) p->atshead.npartials; j++) {
        if (bandsum[bandnum[j]] > 0.0) {
			p->nzpdata[i][j]= ENG_RMS(p->data[i][j].amp * (*partialband[j]) / bandsum[bandnum[j]],p->atshead.winsz);
			*maxval=((float)p->nzpdata[i][j] > (*maxval) ? p->nzpdata[i][j] : (*maxval));
		}
        else
			p->nzpdata[i][j]=0.0;
      }
    }
    free(partialband);
    free(bandnum);
	return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
