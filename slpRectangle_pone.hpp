/***************************************************************************/
/* slpRectangle_pone.hpp                                                   */
/*   Implementiert alle Funktionen zur Berechnung der Galerkin Matrix      */
/*   for discontinouos P1 elements                                         */
/*                                                                         */
/*                                                                         */
/*                                                                         */
/***************************************************************************/
/* Author: Peter Schaefer                             schaeferpm@gmail.com */
/* Alexx 2018                                                              */
/***************************************************************************/


#ifndef HILBERT3D_LAPLACE_SLPRECTANGLE_HPP_GUARD_
#define HILBERT3D_LAPLACE_SLPRECTANGLE_HPP_GUARD_

/*
 * voll analytisch rechnen
 */
double cParO1(double, double, double, double, double, double, double, double*, double, double, double, double);
double cOrtO1(double, double, double, double, double, double, double, double*, double, double, double, double);

/*
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cParO2(double, double, double, double, double, double, double, double*, double, double, double, double);
double cOrtO2(double, double, double, double, double, double, double, double*, double, double, double, double);

/*
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cParO3(double, double, double, double, double, double, double, double*, double, double, double, double);
double cOrtO3(double, double, double, double, double, double, double, double*, double, double, double, double);

/*
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * analytisch sonst
 */
double cParO4(double, double, double, double, double, double, double, double*, double, double, double, double);
double cOrtO4(double, double, double, double, double, double, double, double*, double, double, double, double);


/*
 * setzt die 2^Anzahl der Auswertungsstellen
 */
int setQuad(int);




/*
 *
 *Now the old stuff (for P0)
 *
 *
 */

/*
 * voll analytisch rechnen
 */
double cParO1(double, double, double, double, double, double, double, double*);
double cOrtO1(double, double, double, double, double, double, double, double*);

/*
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cParO2(double, double, double, double, double, double, double, double*);
double cOrtO2(double, double, double, double, double, double, double, double*);

/*
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cParO3(double, double, double, double, double, double, double, double*);
double cOrtO3(double, double, double, double, double, double, double, double*);

/*
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * analytisch sonst
 */
double cParO4(double, double, double, double, double, double, double, double*);
double cOrtO4(double, double, double, double, double, double, double, double*);




#endif
