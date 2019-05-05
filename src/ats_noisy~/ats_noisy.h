/*ats_noisy.h by Pablo Di Liscia*/

#ifndef ATS_NOISY_H 
#define ATS_NOISY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/********************************************************/
/*SYMBOLIC CONSTANTS*/
/********************************************************/
#define PI    	3.14159265358979323846264338327
#define PI2   	6.28318530717958647692528676654
#define TLEN	65536
#define INTERP_TIME 0.0075
#define NB_RES 	25 /*number of critical bands*/
/* array of critical band frequency edges base on data from: 
 * Zwicker, Fastl (1990) "Psychoacoustics Facts and Models", 
 * Berlin ; New York : Springer-Verlag
 */
#define ATSA_CRITICAL_BAND_EDGES {0.0, 100.0, 200.0, 300.0, 400.0, 510.0, \
      630.0, 770.0, 920.0, 1080.0, 1270.0,				\
      1480.0, 1720.0, 2000.0, 2320.0, 2700.0,				\
      3150.0, 3700.0, 4400.0, 5300.0, 6400.0,				\
      7700.0, 9500.0, 12000.0, 15500.0, 20000.0}
#define  ATS_NOISE_THRESHOLD -120
#define  ATS_CRITICAL_BANDS   25
#define  ATS_NOISE_VARIANCE  0.04 

/********************************************************/
/*MACROS*/
/********************************************************/
/*conversions*/
/*dB value from linear value*/
#define LIN2dB(x) (double)(20. * log10(x))
/*compute oscilator increment from frequency*/ 
#define ITOF(frec, L, sr)  frec * ((float)L/(float)sr)
/*compute an increment of s semitones*/
#define S2INCR(s) (pow(2.,(double)s/12.))
/*compute increment from freq.*/
#define F2INCR(f, lsr) lsr*f 
/*convert angles*/
#define DEG2RAD(deg) (deg/360.)* PI2
#define RAD2DEG(rad) (rad/PI2) * 360.
/*get a random value between 0 and "range"*/
#define RAND_FLOAT(range)  range*((float)rand()/ (float)RAND_MAX)
#define BIRAND_FLOAT(range)  (range*2.*((float)rand()/ (float)RAND_MAX)) - range
/*get an int random value between 0 and "range-1"*/
#define RAND_INT(range) rand()%range 
/********************************************************/
/*DATA STRUCTURES for each UG*/
/********************************************************/
typedef struct {
  double  index;
  double  incr;
  double  phase;
  double  L;
  double  lsr;
  float   *table;
}OSCIL;
typedef struct { 
  float size; //size of the frame in samples this should be sr/freq.
  float a1;   //current amplitude value
  float a2;   //next    amplitude value
  float cnt;  //sample  position counter
  float dif;         //to save time computing a2-a1
  float dif_d_siz;   //to save time computing a2-a1 / size
  float siz_d_sr;    // size/sr  
  float sr;   
  float freq;
  float incr;
} RANDI;
/*this is used to interpolate amplitude values*/
typedef struct {
  float val;
  float curr;
  float next;
  float step;
}INTERP;
/***************************************************************************/
/**************PROTOTYPES OF FUNCTIONS FOR AUDIO SYNTHESIS******************/
/***************************************************************************/
float *make_sine_table(float amp, float phase, int length);
/*OSCIL*/
OSCIL   		*set_oscil(int L, float phase, int sr, float *table);
inline float    oscil(OSCIL *oscil);
/*RANDI*/
RANDI 			*set_randi(int sr, float freq);
inline float 	randi(RANDI *radat);
#endif
