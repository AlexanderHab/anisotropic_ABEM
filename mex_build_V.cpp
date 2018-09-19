/***************************************************************************/
/* mex_build_v.cpp                                                         */
/*                                                                         */
/*   Definiert die MEX Funktion mex_build_V(coo,ele,zeta,typ)              */
/*     coo  : Koordinaten Matrix                                           */
/*     ele  : Element Matrix                                               */
/*     zeta : Steuerungsvariable der Zulaessigeitsbedingungen              */
/*     typ  : Bestimmt die Art der Berechnung                              */
/*                                                                         */
/*     TYPES----                                                           */
/*     0 : nur die Abstaende der Elemente werden berechnet                 */
/*       zeta wird ignoriert                                               */
/*                                                                         */
/*     1 : Voll Analytisch                                                 */
/*       zeta wird ignoriert                                               */
/*                                                                         */
/*     2 : Semianalytisch im Fernfeld Quadratur ueber beide Elemente       */
/*       Quadratur bei : max{dia(T1),dia(T2)} < zeta[0] * dist(T1,T2)      */
/*                                                                         */
/*     3 : Semianalytisch wie (2) und Quadratur ueber eine Achse           */
/*       Quadratur bei : min{x[1],y[1]} < zeta[0] * dist(x[1],y[1])        */
/*                                                                         */
/*                                                                         */
/***************************************************************************/
/* Author: Peter Schaefer                             schaeferpm@gmail.com */
/* Version: 1.0  (2012)                                                    */
/***************************************************************************/

/*
 * wenn gesetzt wird, paralleles Rechnen benutzt
 */

//Alexx
// #define PARALLEL
#include <cmath>
#include <mex.h>
// 
// #ifdef PARALLEL
// #include <omp.h>
// #endif

#include "slpRectangle_pone.hpp"

/*
 * Gibt vom Vektor vec die Dimension != 0 zurueck
 */
int inline dimOfVec(double* vec) {
	if (vec[2] != 0)
		return 2;
	else if (vec[1] != 0)
		return 1;
	else if (vec[0] != 0)
		return 0;
	else {
		mexErrMsgTxt("length of Site is zero");
		return -1;
	}
}

/*
 * Gibt von [0 1 2] die fehlende Zahl zurueck
 */
int inline dimOfThird(int a, int b) {
	switch (a + b) {
	case 1:
		return 2;
	case 2:
		return 1;
	case 3:
		return 0;
	default:
		return -1;
	}
}

/*
 * addiert Vektor b zu Vektor a
 */
void inline add(double* a, double* b) {
	for (int i = 0; i < 3; ++i)
		a[i] += b[i];
}

/*
 * vergleicht zwei Vektoren a, b
 */
bool inline iseq(double* a, double* b) {
	for (int i = 0; i < 3; ++i)
		if (a[i] != b[i])
			return false;
	return true;
}

/*
 * vergleicht Vektor a mit Konstante b
 */
bool inline iseq(double* a, double b) {
	for (int i = 0; i < 3; ++i)
		if (a[i] != b)
			return false;
	return true;
}

/*
 * speichert die Summe der Vektoren a,b in Vektor c
 */
void inline add(double* a, double* b, double* c) {
	for (int i = 0; i < 3; ++i)
		a[i] = b[i] + c[i];
}

/*
 * subtrahiert Vektor b von Vektor a
 */
void inline sub(double* a, double* b) {
	for (int i = 0; i < 3; ++i)
		a[i] -= b[i];
}

/*
 * speichert die Differenz der Vektoren a,b in Vektor c
 */
void inline sub(double* a, double* b, double* c) {
	for (int i = 0; i < 3; ++i)
		a[i] = b[i] - c[i];
}

/*
 * speichert Vektor b in Vektor a
 */
void inline set(double* a, double* b) {
	for (int i = 0; i < 3; ++i)
		a[i] = b[i];
}

/*
 * setzt alle Werte in Vektor a auf Konstante b
 */
void inline set(double* a, double b) {
	for (int i = 0; i < 3; ++i)
		a[i] = b;
}

/*
 * speichert die kleinste Zeile von Matrix x in Vektor d
 */
void inline getSCorner(double* d, double x[4][3]) {
	set(d, x[0]);
	for (int i = 1; i < 4; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (d[j] > x[i][j])
				set(d, x[i]);
			else if (d[j] < x[i][j])
				break;
		}
	}
}

/*
 * berechnet den Abstand von zwei parallelen Elementen * 4Pi
 */
double inline distT_par(double b, double d, double t, double v, double d1,
		double d2, double d3, double * dummy) {
	double dis1 = 0, dis2 = 0;

	if (d1 < 0)
		dis1 = -d1 - t;
	else if (d1 > 0)
		dis1 = d1 - b;
	if (dis1 < 0)
		dis1 = 0;

	if (d2 < 0)
		dis2 = -d2 - v;
	else if (d2 > 0)
		dis2 = d2 - d;
	if (dis2 < 0)
		dis2 = 0;

	return sqrt(dis1 * dis1 + dis2 * dis2 + d3 * d3) * (4 * M_PI);
}

/*
 * berechnet den Abstand von zwei orthogonalen Elementen * 4Pi
 */
double inline distT_ort(double b, double d, double t, double v, double d1,
		double d2, double d3, double * dummy) {
	double dis1 = 0, dis2 = 0, dis3 = 0;
	if (d1 < 0)
		dis1 = -d1;
	else
		dis1 = d1 - b;
	if (dis1 < 0)
		dis1 = 0;

	if (d2 < 0)
		dis2 = -d2 - t;
	else
		dis2 = d2 - d;
	if (dis2 < 0)
		dis2 = 0;

	if (d3 < 0)
		dis3 = -d3 - v;
	else
		dis3 = d3;
	if (dis3 < 0)
		dis3 = 0;

	return sqrt(dis1 * dis1 + dis2 * dis2 + dis3 * dis3) * (4 * M_PI);
}

/*
 * mexFunktion zum Aufbauen der Galerkin Matrix
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, k; //Schleifenindizes
	bool count = false;

	/*
	 * Sicherheitsueberpruefung der Parameter
	 */
	if ((nrhs != 4))
		mexErrMsgTxt(
				"erwarte (coordinates(Nx3), elements(Mx4), zeta(double), type(int))");
	if (nlhs > 2 || nlhs == 0)
		mexErrMsgTxt("Funktion hat maximal 2 Rueckgabewerte [A counts]");

	int cm = mxGetM(prhs[0]);
	int cn = mxGetN(prhs[0]);
	if (cn != 3)
		mexErrMsgTxt("erwarte Koordinatenmatrix im (Nx3) Format");

	int em = mxGetM(prhs[1]);
	int en = mxGetN(prhs[1]);
	if (en != 4)
		mexErrMsgTxt("erwarte Elementmatrix im (Mx4) Format");

// Alexx
// 	/*
// 	 * Parameter fuer parallele Berechnung setzen
// 	 */
// #ifdef PARALLEL
// #define MINSIZE_PER_WORKER 3
// #define MAX_WORKER 10
// 
// 	int actualNumberOfThreads = omp_get_max_threads();
// 
// 	if (MAX_WORKER < actualNumberOfThreads)
// 	actualNumberOfThreads = MAX_WORKER;
// 
// 	int firstRow = 0, lastRow = -1;
// 	int targetSize = em * (em + 1) / (2 * actualNumberOfThreads);
// 
// 	if (targetSize < MINSIZE_PER_WORKER) {
// 		actualNumberOfThreads = (int) em * (em + 1) / (2 * MINSIZE_PER_WORKER);
// 		if (actualNumberOfThreads < 1) {
// 			actualNumberOfThreads = 1;
// 			targetSize = em;
// 		} else {
// 			targetSize = MINSIZE_PER_WORKER;
// 		}
// 	}
// #endif

	/*
	 * Variablen initialisieren
	 */
    
    
    
	plhs[0] = mxCreateDoubleMatrix(em, em, mxREAL);

	double * A = mxGetPr(plhs[0]);
	double * C = mxGetPr(prhs[0]);      // get Coordinates
	double * E = mxGetPr(prhs[1]);      // get Elements
	double * COUNT = NULL;

	if(nlhs==2){                       
		count = true;
		plhs[1] = mxCreateDoubleMatrix(2,2,mxREAL);
		COUNT = mxGetPr(plhs[1]);
	}

	int type = (int) *(mxGetPr(prhs[3]));
	double * zeta = { 0 };

	if (mxGetN(prhs[2]) > 1)
		zeta = mxGetPr(prhs[2]) + 1;

	//setzen der Quadraturpunkte
	int quad = (int) *mxGetPr(prhs[2]);
	quad = setQuad(quad);
    

	//Art der Berechnung bestimmen
	double (*ctypeP)(double, double, double, double, double, double, double,
			double*);
	double (*ctypeO)(double, double, double, double, double, double, double,
			double*);

	switch (type) {
	default: // voll analytisch
		ctypeP = cParO1;
		ctypeO = cOrtO1;
		break;
	case 2: // analytisch oder Qadratur
		ctypeP = cParO2;
		ctypeO = cOrtO2;
		break;
	case 3: // analytisch oder Qadratur oder Qy[1]x2
		ctypeP = cParO3;
		ctypeO = cOrtO3;
		break;
	case 4: // analytisch oder Qy[1]x2
		ctypeP = cParO4;
		ctypeO = cOrtO4;
//		ctypeO = cOrtO2;
		break;
	case 0: // Distanz der Elemente
		ctypeP = distT_par;
		ctypeO = distT_ort;
		break;
	}

//      printf("test c-file 3 \n");
     
     
	//Lageinformationen
	int rx, rxa, rxb, ry, rya, ryb;
	double tmp;
	double x[4][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
	double y[4][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
	double xa[3] = { 0, 0, 0 };
	double xb[3] = { 0, 0, 0 };
	double ya[3] = { 0, 0, 0 };
	double yb[3] = { 0, 0, 0 };
	double d[3] = { 0, 0, 0 };

	/*
	 * Schleife ueber Elemente initialisieren
	 */

// Alexx
// #ifdef PARALLEL
// 	{
// #pragma omp parallel for schedule(static,1) private(i,j,k,tmp,x,y,xa,xb,ya,yb,d,rx,rxa,rxb,ry,rya,ryb,lastRow,firstRow) shared(C,E,ctypeO,ctypeP,em,targetSize,actualNumberOfThreads)
// 
// 		for (i = 0; i < actualNumberOfThreads; ++i) {
// 			firstRow = (int) sqrt(2 * i * targetSize);
// 			if (i == actualNumberOfThreads - 1) {
// 				lastRow = em - 1;
// 			} else {
// 				lastRow = (int) sqrt(2 * (i + 1) * targetSize) - 1;
// 			}
// 
// 			for (j = firstRow; j <= lastRow; ++j) {
// #else
// 	{
// 		{
 			for (j = 0; j < em; ++j) {
// 
// #endif

            
				// Eckdaten zwischenspeichern
				x[0][0] = C[(int) E[j] - 1];        // beruecksichtige spaltenweise speicherung
				x[0][1] = C[cm + (int) E[j] - 1];
				x[0][2] = C[2 * cm + (int) E[j] - 1];

				x[1][0] = C[(int) E[em + j] - 1];
				x[1][1] = C[cm + (int) E[em + j] - 1];
				x[1][2] = C[2 * cm + (int) E[em + j] - 1];

				x[2][0] = C[(int) E[2 * em + j] - 1];
				x[2][1] = C[cm + (int) E[2 * em + j] - 1];
				x[2][2] = C[2 * cm + (int) E[2 * em + j] - 1];

				x[3][0] = C[(int) E[3 * em + j] - 1];
				x[3][1] = C[cm + (int) E[3 * em + j] - 1];
				x[3][2] = C[2 * cm + (int) E[3 * em + j] - 1];

				// Seitenvektoren aufbauen
				sub(xa, x[1], x[0]);    // erg is stored in xa
				sub(xb, x[3], x[0]);

				// Lageeigenschaften des Flaechenstuecks
				rxa = dimOfVec(xa);                 // index der nichtnull-koordinate
				rxb = dimOfVec(xb);                 // index der nichtnull-koordinate
				rx = dimOfThird(rxa, rxb);
                
//                 printf("test c-file 5 \n");

				if (xa[rxa] <= 0 || xb[rxb] <= 0)
					mexWarnMsgTxt("X Element hat Seitenlaengen <=0!");

				if (rxa == rxb)
					mexWarnMsgTxt("X Laenge ist nur in eine Richtung");

				for (k = 0; k <= j; ++k) {

					//Eckdaten zwischenspeichern
					y[0][0] = C[(int) E[k] - 1];
					y[0][1] = C[cm + (int) E[k] - 1];
					y[0][2] = C[2 * cm + (int) E[k] - 1];

					y[1][0] = C[(int) E[em + k] - 1];
					y[1][1] = C[cm + (int) E[em + k] - 1];
					y[1][2] = C[2 * cm + (int) E[em + k] - 1];

					y[2][0] = C[(int) E[2 * em + k] - 1];
					y[2][1] = C[cm + (int) E[2 * em + k] - 1];
					y[2][2] = C[2 * cm + (int) E[2 * em + k] - 1];

					y[3][0] = C[(int) E[3 * em + k] - 1];
					y[3][1] = C[cm + (int) E[3 * em + k] - 1];
					y[3][2] = C[2 * cm + (int) E[3 * em + k] - 1];

					// Seiten Vektoren aufbauen
					sub(ya, y[1], y[0]);
					sub(yb, y[3], y[0]);

					// Lageeigenschaften des Flaechenstuecks
					rya = dimOfVec(ya);                     // kann ev. zu problemen führen
					ryb = dimOfVec(yb);                     // wenn eintrag = 0 aber trotzen dim 
					ry = dimOfThird(rya, ryb);

					if (ya[rya] <= 0 || yb[ryb] <= 0)
						mexWarnMsgTxt("Y Element hat Seitenlaengen <=0!");

					if (rya == ryb)
						mexWarnMsgTxt("Y Laenge ist nur in einer Richtung");

					//delta berechnen
					sub(d, y[0], x[0]);                     // is sicher nicht das  minimum
                    
//                     printf("test c-file 6 \n");

					/*
					 * Matrix Eintrag berechnen
					 */
					if (rx == ry) {
						if (rxa == rya) {
							tmp = ctypeP((xa[rxa]), (xb[rxb]), (ya[rxa]),
									(yb[rxb]), d[rxa], d[rxb], d[rx], zeta);

						} else {
							tmp = ctypeP((xa[rxa]), (xb[rxb]), (yb[rxa]),
									(ya[rxb]), d[rxa], d[rxb], d[rx], zeta);
						}

					} else {

						if (rxa == rya) {
							tmp = ctypeO((xb[rxb]), (xa[rxa]), (ya[rya]),
									(yb[ryb]), d[rxb], d[rxa], d[rx], zeta);
						} else if (rxa == ryb) {
							tmp = ctypeO((xb[rxb]), (xa[rxa]), (yb[ryb]),
									(ya[rya]), d[rxb], d[rxa], d[rx], zeta);
						} else if (rxb == rya) {
							tmp = ctypeO((xa[rxa]), (xb[rxb]), (ya[rya]),
									(yb[ryb]), d[rxa], d[rxb], d[rx], zeta);
						} else {
							tmp = ctypeO((xa[rxa]), (xb[rxb]), (yb[ryb]),
									(ya[rya]), d[rxa], d[rxb], d[rx], zeta);
						}

					}
// Alexx
// #ifdef PARALLEL
// #pragma omp critical
// #endif
// 
 					A[(k * em) + j] = 1. / (4 * M_PI) * tmp;
 					if (k != j)
 						A[(j * em) + k] = 1. / (4 * M_PI) * tmp;
// 
// #ifdef PARALLEL
// #pragma omp end critical
// #endif
				}
//			}
//		}
	}
//     printf("test c-file 7 \n");
	return;
}
