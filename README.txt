ATS-PD by Pablo Di Liscia
odiliscia@unq.edu.ar 
Programa de Investigacion "Sistemas Temporales y Sintesis Espacial en el Arte Sonoro"
http://stseas.web.unq.edu.ar/
odiliscia@unq.edu.ar

This is a software package for ATS Synthesys (Juan Pampin) using 
Pure Data (Miller Puckette).
The externals included are:

1-atsread, by Alex Norman, with improvements by Pablo Di Liscia
2-ats_noisy~ by Pablo Di Liscia.
3-ats_sinnoi~ by Richie Eakin and Pablo Di Liscia.
4-oscbank~ by Richie Eakin.

Make sure you have all the externals properly compiled
and instaled in your system before running the examples.
The oscbank~ external by Richie Eakin is also needed:
https://github.com/pd-l2ork/pd/tree/master/externals/oscbank~

Download the 'ats-files' package in case you want have some ATS analysis
files to use with the example patches.

For a full study of the ATS technique, you may consult:
https://dxarts.washington.edu/sites/dxarts/files/documents/wiki/ats_theory.pdf
http://lac.linuxaudio.org/2013/papers/26.pdf

Compiling the externals:
(note that you may need to edit the Makefiles to meet your system settings and
the compiler you will use)

a-Linux: 
	-Copy all the  ats-pd folder in the path you want.
	-cd your_folder/ats-pd/
	-make -f Makefile.linux

b-Windows (you will need CygWin properly installed in your system):
	-Copy all the  ats-pd folder in the path you want.
	-cd your_folder/ats-pd/
	-make -f Makefile.win

Please, email me if you have questions or suggestions.

Pablo Di Liscia
odiliscia@unq.edu.ar
Universidad Nacional de Quilmes
Argentina

Information Updated at 01/03/2019

  
