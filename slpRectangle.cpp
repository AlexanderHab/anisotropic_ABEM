/***************************************************************************/
/* slpRectangle.cpp                                                        */
/*   Implementiert alle Funktionen zum Berechnen der GalerkinMatrix        */
/*                                                                         */
/*                                                                         */
/*                                                                         */
/*                                                                         */
/***************************************************************************/
/* Author: Peter Schaefer                             schaeferpm@gmail.com */
/* Version: 1.0  (2012)                                                    */
/***************************************************************************/

#include <math.h>

#include <mex.h>

#include "slpRectangle.hpp"
#include "gauss.hpp"



int GAUSS_size = GAUSS_SIZE[2];
gauss * GAUSS_nodes = GAUSS_NODES[2];



/*
 * setzt den Quadraturgrad, also 2^q Auswertungsstellen
 */
int setQuad(int q) {
	GAUSS_size = GAUSS_SIZE[q];
	GAUSS_nodes = GAUSS_NODES[q];
	return GAUSS_size;
}

/*
 * gibt das Signum von x zurueck
 */
int inline sign(double x) {
	return x > 0 ? 1 : (x < 0 ? -1 : 0);
}

/*
 * gibt das Maximum der beiden Zahlen x,y zurueck
 */
double inline max(double x, double y) {
	return x < y ? y : x;
}

/*
 * gibt das Minimum der beiden Zahlen x,y zurueck
 */
double inline min(double x, double y) {
	return x > y ? y : x;
}

/*
 * kleinere Seiten nach vorn und dadurch kleinerer Durchmesser nach vorn
 */
void inline switch_site_par(double& b, double& d, double& t, double& v,
		double& d1, double& d2) {
	double tmp = 0;

	// kleine Seite nach vorn (x1 y1)
	if (b > t) {
		tmp = b;
		b = t;
		t = tmp;
		d1 = -d1;
	}

	// kleine Seite nach vorn (x2 y2)
	if (d > v) {
		tmp = d;
		d = v;
		v = tmp;
		d2 = -d2;
	}
}

/*
 * kleinere Seiten nach vorn und dadurch kleinerer Durchmesser nach vorn
 */
void inline switch_site_ort(double& b, double& d, double& t, double& v,
		double& d1, double& d2, double& d3) {
	double tmp = 0;

	// kleine Seite nach vorn (x2 y2)
	if (d > t) {
		tmp = d;
		d = t;
		t = tmp;
		d2 = -d2;
	}

	// kleine Seite nach vorn (x1 y3)
	if (b > v) {
		tmp = b;
		b = v;
		v = tmp;
		tmp = -d1;
		d1 = -d3;
		d3 = tmp;
	}
}

/*
 * kleinere Achse nach vorn
 */
void inline switch_dim_par(double& b, double& d, double& t, double& v,
		double& d1, double& d2) {
	//
	if (b * b + t * t > d * d + v * v) {
		double tmp = 0;

		tmp = b;
		b = d;
		d = tmp;
		tmp = t;
		t = v;
		v = tmp;

		tmp = d1;
		d1 = d2;
		d2 = tmp;
	}
}

/*
 * gibt den Abstand zum Quadrat zurueck
 */
double inline dist2_par(double b, double d, double t, double v, double d1,
		double d2, double d3) {
	double dis1 = 0, dis2 = 0;
	if (d1 < 0)
		dis1 = -d1 - t;
	else
		dis1 = d1 - b;
	if (dis1 < 0)
		dis1 = 0;

	if (d2 < 0)
		dis2 = -d2 - v;
	else
		dis2 = d2 - d;
	if (dis2 < 0)
		dis2 = 0;

	return dis1 * dis1 + dis2 * dis2 + d3 * d3;
}

/*
 * gibt den Abstand zum Quadrat zurueck
 */
double inline dist2_ort(double b, double d, double t, double v, double d1,
		double d2, double d3) {
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

	return dis1 * dis1 + dis2 * dis2 + dis3 * dis3;
}

/*
 * gibt den Abstand in einer Richtung zum Quadrat zurueck
 */
double inline dist_s2(double b, double t, double d1) {
	double dis1 = 0;
	if (d1 < 0)
		dis1 = -d1 - t;
	else
		dis1 = d1 - b;
	if (dis1 < 0)
		dis1 = 0;

	return dis1 * dis1;
}

/*
 * Kernfunktion 1 /|x-y|
 */
double inline f_par(double x1, double x2, double y1, double y2, double d1,
		double d2, double d3) {
	return 1.
			/ sqrt(
					(x1 - y1 - d1) * (x1 - y1 - d1)
							+ (x2 - y2 - d2) * (x2 - y2 - d2) + d3 * d3);
}
/*
 * Kernfunktion 1 /|x-y|
 */
double inline f_ort(double x1, double x2, double y2, double y3, double d1,
		double d2, double d3) {
	return 1.
			/ sqrt(
					(x1 - d1) * (x1 - d1) + (x2 - y2 - d2) * (x2 - y2 - d2)
							+ (y3 + d3) * (y3 + d3));
}

/*
 * 1 mal integrierte Kernfunktion int 1 /|x-y| dy (QUADRATUR)
 */
double inline g_QY(double p, double y, double x, double l) {
	double sol = 0;
	for (int i = 0; i < GAUSS_size; ++i) {
		sol += pow(
				(x - GAUSS_nodes[i].n * y) * (x - GAUSS_nodes[i].n * y) + l * l,
				p) * GAUSS_nodes[i].w * y;
	}

	return sol;
}

/*
 * Stammfunktion g
 * 1 mal integrierte Kernfunktion int 1 /|x-y| dy (ANALYTISCH)
 */
double inline g_AY(double p, double y, double x, double l) {
	double sol = 0;

	if (l != 0) {
		if (p == 0.5) {
			sol = (y - x) / 2 * sqrt((y - x) * (y - x) + l * l)
					+ l * l / 2 * asinh((y - x) / fabs(l));
		} else if (p == 0)
			sol = y - x;
		else if (p == -0.5)
			sol = asinh((y - x) / fabs(l));
		else if (p == -1)
			sol = atan((y - x) / fabs(l)) / fabs(l);
		else if (p == -1.5)
			sol = (y - x) / ((l * l) * sqrt((y - x) * (y - x) + l * l));
		else
			sol = (y - x) * pow((y - x) * (y - x) + l * l, p)
					+ 2 * p * l * l * g_AY(p - 1, y, x, l) / (2 * p + 1);
	} else {
		if (p == -0.5)
			sol = sign(y - x) * log(fabs(y - x));
		else
			sol = (y - x) * pow(fabs(y - x), 2 * p) / (2 * p + 1);
	}

	return sol;
}

/*
 * 2 mal integrierte Kernfunktion int int 1 /|x-y| dy_2 dy_1  (QUADRATUR)
 */
double inline G_QY1Y2(double p, double y1, double y2, double x1, double x2,
		double l) {

	double sol = 0;
	for (int i = 0; i < GAUSS_size; ++i) {
		for (int j = 0; j < GAUSS_size; ++j) {
			sol += pow(
					(x1 - y1 * GAUSS_nodes[i].n) * (x1 - y1 * GAUSS_nodes[i].n)
							+ (x2 - y2 * GAUSS_nodes[j].n)
									* (x2 - y2 * GAUSS_nodes[j].n) + l * l, p)
					* y1 * GAUSS_nodes[i].w * y2 * GAUSS_nodes[j].w;
		}
	}
	return sol;
}

/*
 * 2 mal integrierte Kernfunktion int int 1 /|x-y| dy_2 dx_2  (ANALYTISCH)
 */
double inline G_AY2X2(double y1, double y2, double x1, double x2, double d3) {
	double hlp = ((x1 - y1) * (x1 - y1) + (x2 - y2) * (x2 - y2) + d3 * d3);
	double sol = sqrt(hlp);

	if ((x2 - y2) != 0)
		sol += (x2 - y2) * log(sqrt(hlp) - (x2 - y2));

	return sol;
}

/*
 * Stammfunktion G
 * 2 mal integrierte Kernfunktion int int 1 /|x-y| dy_2 dy_1  (ANALYTISCH)
 */
double inline G_AY1Y2(double p, double y1, double y2, double x1, double x2,
		double l) {
	double pt = p + 1.5;
	double sol = 0;
	if (pt == 0) {
		if (l == 0) {
			sol = -sqrt((y1 - x1) * (y1 - x1) + (y2 - x2) * (y2 - x2))
					/ ((y1 - x1) * (y2 - x2));
		} else {
			sol = sign((y1 - x1) * (y2 - x2)) / (2 * fabs(l))
					* acos(
							-2 * (y1 - x1) * (y1 - x1) * (y2 - x2) * (y2 - x2)
									/ (((y1 - x1) * (y1 - x1) + l * l)
											* ((y2 - x2) * (y2 - x2) + l * l))
									+ 1);
		}
	} else if ((pt > 0) && ((int) pt == pt)) {
		if (l != 0)
			sol = 2 * p * l * l * G_AY1Y2(p - 1, y1, y2, x1, x2, l);
		if ((y1 - x1) != 0)
			sol += (y1 - x1)
					* g_AY(p, y2, x2, sqrt((y1 - x1) * (y1 - x1) + l * l));
		if ((y2 - x2) != 0)
			sol += (y2 - x2)
					* g_AY(p, y1, x1, sqrt((y2 - x2) * (y2 - x2) + l * l));
		sol /= 2 * p + 2;
	} else {
		//Fall wird nicht benoetigt
		sol = NAN;
		mexWarnMsgTxt("G_AY1Y2 fuer Parameter nicht implementiert!");
	}

	return sol;
}

/*
 * Stammfunktion F
 * 4 mal integrierte Kernfunktion int int 1 /|x-y| dy dx  (ANALYTISCH)
 */
double inline F_par(double x1, double x2, double y1, double y2, double d1,
		double d2, double d3) {
	double sol = (x1 - y1 - d1) * (x2 - y2 - d2);

	if (sol != 0)
		sol *= G_AY1Y2(-0.5, x1, x2, y1 + d1, y2 + d2, d3);

	if ((x1 - y1 - d1) != 0)
		sol -= (x1 - y1 - d1)
				* g_AY(0.5, x1, y1 + d1,
						sqrt((x2 - y2 - d2) * (x2 - y2 - d2) + d3 * d3));

	if ((x2 - y2 - d2) != 0)
		sol -= (x2 - y2 - d2)
				* g_AY(0.5, x2, y2 + d2,
						sqrt((x1 - y1 - d1) * (x1 - y1 - d1) + d3 * d3));

	double hlp = ((x1 - y1 - d1) * (x1 - y1 - d1)
			+ (x2 - y2 - d2) * (x2 - y2 - d2) + d3 * d3);
	sol += 1. / 3 * hlp * sqrt(hlp);
	return sol;
}

/*
 * Stammfunktion F
 * 4 mal integrierte Kernfunktion int int 1 /|x-y| dy dx  (ANALYTISCH)
 */
double inline F_ort(double x1, double x2, double y2, double y3, double d1,
		double d2, double d3) {
	double sol = -G_AY1Y2(0.5, y3, x1, -d3, d1, x2 - y2 - d2);

	if ((x1 - d1) * (x2 - y2 - d2) != 0)
		sol -= (x1 - d1) * (x2 - y2 - d2)
				* G_AY1Y2(-0.5, x2, y3, y2 + d2, -d3, x1 - d1);

	if ((x1 - d1) != 0)
		sol += (x1 - d1)
				* g_AY(0.5, y3, -d3,
						sqrt(
								(x1 - d1) * (x1 - d1)
										+ (x2 - y2 - d2) * (x2 - y2 - d2)));

	if ((y3 + d3) * (x2 - y2 - d2) != 0)
		sol -= (y3 + d3) * (x2 - y2 - d2)
				* G_AY1Y2(-0.5, x1, x2, d1, y2 + d2, -y3 - d3);

	if ((y3 + d3) != 0)
		sol += (y3 + d3)
				* g_AY(0.5, x1, d1,
						sqrt(
								(x2 - y2 - d2) * (x2 - y2 - d2)
										+ (y3 + d3) * (y3 + d3)));

	return sol / 2.;
}

/*
 * Stammfunktion (fuer Quadratur)
 * 2 mal integrierte Kernfunktion int int 1 /|x-y| dy_2 dy_1  (ANALYTISCH)
 */
double inline FY1Y2_par(double x1, double x2, double y1, double y2, double d1,
		double d2, double d3) {
	return G_AY1Y2(-0.5, y1, y2, x1 - d1, x2 - d2, d3);
}

/*
 * Stammfunktion (fuer Quadratur)
 * 2 mal integrierte Kernfunktion int int 1 /|x-y| dy_2 dy_3  (ANALYTISCH)
 */
double inline FY2Y3_ort(double x1, double x2, double y2, double y3, double d1,
		double d2, double d3) {
	return G_AY1Y2(-0.5, y2, y3, x2 - d2, -d3, x1 - d1);
}

/*
 * Quadratur ueber zwei Integrale und Grenzen in Stammfunktion ueber zwei Integrale
 * b,d - Grenzen fuer Quadratur
 * t,v - Grenzen in Stammfunktion
 */
template<double (F)(double, double, double, double, double, double, double)>
double intQ2A2(double b, double d, double t, double v, double d1, double d2,
		double d3) {
	double sol = 0;
	double w1, x1, x2;

	for (int i = 0; i < GAUSS_size; ++i) {
		x1 = b * GAUSS_nodes[i].n;
		w1 = GAUSS_nodes[i].w;
		for (int j = 0; j < GAUSS_size; ++j) {
			x2 = d * GAUSS_nodes[j].n;
			sol += (F(x1, x2, t, v, d1, d2, d3) - F(x1, x2, 0, v, d1, d2, d3)
					- F(x1, x2, t, 0, d1, d2, d3) + F(x1, x2, 0, 0, d1, d2, d3))
					* w1 * GAUSS_nodes[j].w;
		}
	}

	return sol * b * d;
}

/*
 * Quadratur ueber vier Integrale
 * b,d,t,v - Grenzen fuer Quadratur
 */
template<double (F)(double, double, double, double, double, double, double)>
double intQ4(double b, double d, double t, double v, double d1, double d2,
		double d3) {

	double sol = 0;
	int i, j, k, l;
	double x1, x2, y1;
	double w1, w2, w3;

	for (i = 0; i < GAUSS_size; ++i) {
		x1 = b * GAUSS_nodes[i].n;
		w1 = GAUSS_nodes[i].w;
		for (j = 0; j < GAUSS_size; ++j) {
			x2 = d * GAUSS_nodes[j].n;
			w2 = w1 * GAUSS_nodes[j].w;
			for (k = 0; k < GAUSS_size; ++k) {
				y1 = t * GAUSS_nodes[k].n;
				w3 = w2 * GAUSS_nodes[k].w;
				for (l = 0; l < GAUSS_size; ++l) {
					sol += F(x1, x2, y1, v * GAUSS_nodes[l].n, d1, d2, d3) * w3
							* GAUSS_nodes[l].w;
				}
			}
		}
	}

	return sol * b * d * t * v;
}

/*
 * Grenzen in Stammfunktion ueber vier Integrale
 * b,d,t,v - Grenzen fuer Stammfunktion
 */
template<double (F)(double, double, double, double, double, double, double)>
double inline intA4(double b, double d, double t, double v, double d1,
		double d2, double d3) {

	return F(b, d, t, v, d1, d2, d3) - F(b, d, t, 0, d1, d2, d3)
			- F(b, d, 0, v, d1, d2, d3) + F(b, d, 0, 0, d1, d2, d3)
			- F(b, 0, t, v, d1, d2, d3) + F(b, 0, t, 0, d1, d2, d3)
			+ F(b, 0, 0, v, d1, d2, d3) - F(b, 0, 0, 0, d1, d2, d3)
			- F(0, d, t, v, d1, d2, d3) + F(0, d, t, 0, d1, d2, d3)
			+ F(0, d, 0, v, d1, d2, d3) - F(0, d, 0, 0, d1, d2, d3)
			+ F(0, 0, t, v, d1, d2, d3) - F(0, 0, t, 0, d1, d2, d3)
			- F(0, 0, 0, v, d1, d2, d3) + F(0, 0, 0, 0, d1, d2, d3);

}

/*
 * fuer parallele Elemente
 * voll analytisch
 */
double cParO1(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta) {
	return intA4<F_par>(b, d, t, v, d1, d2, d3);
}
/*
 * fuer orthogonale Elemente
 * voll analytisch
 */
double cOrtO1(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta) {
	return intA4<F_ort>(b, d, t, v, d1, d2, d3);
}

/*
 * fuer parallele Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cParO2(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta) {
	double tmp = 0;

	//kurze Seite nach vorn
	switch_site_par(b, d, t, v, d1, d2);

	if (zeta[0] * zeta[0] * (t * t + v * v) < dist2_par(b, d, t, v, d1, d2, d3))
		return intQ4<f_par>(b, d, t, v, d1, d2, d3);

	return intA4<F_par>(b, d, t, v, d1, d2, d3);
}

/*
 * fuer orthogonale Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cOrtO2(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta) {
	double tmp = 0, dis2 = 0;

	//kurze Seite nach vorn
	switch_site_ort(b, d, t, v, d1, d2, d3);

	if (zeta[0] * zeta[0] * (t * t + v * v) < dist2_ort(b, d, t, v, d1, d2, d3))
		return intQ4<f_ort>(b, d, t, v, d1, d2, d3);

	return intA4<F_ort>(b, d, t, v, d1, d2, d3);
}

/*
 * fuer parallele Elemente
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cParO3(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta) {

	//kurze Seite nach vorn
	switch_site_par(b, d, t, v, d1, d2);

	if (zeta[1] * zeta[1] * (b * b + d * d) < dist2_par(b, d, t, v, d1, d2, d3))
		return intQ2A2<FY1Y2_par>(b, d, t, v, d1, d2, d3);

	if (zeta[0] * zeta[0] * (t * t + v * v) < dist2_par(b, d, t, v, d1, d2, d3))
		return intQ4<f_par>(b, d, t, v, d1, d2, d3);

	return intA4<F_par>(b, d, t, v, d1, d2, d3);
}

/*
 * fuer orthogonale Elemente
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cOrtO3(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta) {
	double tmp = 0, dis2 = 0;

	//kurze Seite nach vorn
	switch_site_ort(b, d, t, v, d1, d2, d3);

	if (zeta[1] * zeta[1] * (b * b + d * d) < dist2_ort(b, d, t, v, d1, d2, d3))
		return intQ2A2<FY2Y3_ort>(b, d, t, v, d1, d2, d3);

	if (zeta[0] * zeta[0] * (t * t + v * v) < dist2_ort(b, d, t, v, d1, d2, d3))
		return intQ4<f_ort>(b, d, t, v, d1, d2, d3);

	return intA4<F_ort>(b, d, t, v, d1, d2, d3);
}

/*
 * fuer parallele Elemente
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * analytisch sonst
 */
double cParO4(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta) {

	//kurze Seite nach vorn
	switch_site_par(b, d, t, v, d1, d2);

	if (zeta[1] * zeta[1] * (b * b + d * d) < dist2_par(b, d, t, v, d1, d2, d3))
		return intQ2A2<FY1Y2_par>(b, d, t, v, d1, d2, d3);

	return intA4<F_par>(b, d, t, v, d1, d2, d3);
}
/*
 * fuer orthogonale Elemente
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * analytisch sonst
 */
double cOrtO4(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta) {
	double tmp = 0, dis2 = 0;

	//kurze Seite nach vorn
	switch_site_ort(b, d, t, v, d1, d2, d3);

	if (zeta[1] * zeta[1] * (b * b + d * d) < dist2_ort(b, d, t, v, d1, d2, d3))
		return intQ2A2<FY2Y3_ort>(b, d, t, v, d1, d2, d3);

	return intA4<F_ort>(b, d, t, v, d1, d2, d3);
}
