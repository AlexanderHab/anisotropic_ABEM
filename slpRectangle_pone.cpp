/***************************************************************************/
/* slpRectangle.cpp                                                        */
/*   Implementiert alle Funktionen zum Berechnen der GalerkinMatrix        */
/*                                                                         */
/*                                                                         */
/*                                                                         */
/*                                                                         */
/***************************************************************************/
/* Code is based on the bachelor thesis of Peter Schäfer                   */
/* adaptiert von Alexx 2018                                                */
/***************************************************************************/


#include <math.h>
#include <mex.h>

#include "slpRectangle_pone.hpp"
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
 * Stammfunktion g
 * 1 mal integrierte Kernfunktion int 1 /|x-y| dy (ANALYTISCH)
 * 1 parameter k
 * Formula on Maischak p3. (bottom (10))
 */
double inline g_neu_pone(double p, double y, double x, double a, double k) {

    double sol = 0;
// printf("g_neu_pone, p= %f, y=%f, x=%f, a=%f, k=%f \n", p,y,x,a, k);

    // if k=0 solve directly formula (13)
    if (k==0){
        if(a!=0){
            if ((y-x) != 0) {
            
                if (p == 0)
                    sol = y - x;
                else if (p == -0.5)
                    sol = asinh((y - x) / fabs(a));

                else if (p == -1)
                    sol = atan((y - x) / fabs(a));

                else if (p == -1.5)
                    sol = (y - x) / ((a * a) * sqrt((y - x) * (y - x) + a * a));
                
                else {
                    sol = (y - x) * pow((y - x) * (y - x) + a * a, p)
					+ 2 * p * a * a * g_neu_pone(p - 1, y, x, a,k);
                    sol = sol / (2 * p + 1);
                }
                    
            }
            else{
                sol = 0;
            }
        }      
        else{
           printf("a=0 in g_neu_pone mit k=0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        }
        
    }   
     //otherwise (k!=0) start recursion (formula (10))
     else{
         if(k+2*p+1 !=0){
//               printf("standard case \n");
         
             sol = pow(y,k-1)*pow((y-x)*(y-x)+a*a,p+1);

             if(k+p !=0 && x!=0){
                sol += 2*x*(k+p)*g_neu_pone(p,y,x,a,k-1);
             }

             if(k-1 >0){
                 sol += (-1)*(k-1)*(x*x+a*a)* g_neu_pone(p,y,x,a,k-2);
             }

             sol = sol/(k+2*p+1);
         }
         else if(k+2*p+1==0){
             printf("new case \n");
             
             if(2*p+2 != 0){
             
                sol = pow(y,k-1)*pow((y-x)*(y-x)+a*a,p+1) ;        
             
                 if(k>1){
                     sol += (k-1)* g_neu_pone(p+1,y,x,a,k-2);
                 }
                sol = sol/(2*p+2);
             
             }            
             else if( 2*p+2 ==0 && k==1){
                 sol = log( sqrt((y-x)*(y-x)+a*a)) +x*g_neu_pone(-1,y,x,a,0);
                 
             }
             
         }
         else{
             printf("something went wrong in g_neu \n");
         
         }


     }
    
//	printf("g_neu_pone sol: %f \n",sol); 
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
 *lt. Maischak p. 4 (14)
 * resp. p5
 */
    


double inline G_2_neu(double p, double y1, double y2, double x1, double x2, 
        double a, double k, double l) {
      
    double sol = 0;
    double zz =0; 

// printf("G_2_neu, p= %f, y1=%f, y2=%f, x1=%f, x2=%f, a=%f, k=%f, l=%f \n", p,y1,y2,x1,x2,a,k,l);    
    if(p != -3./2){
        
        if( fabs(y1-x1) !=0){
            sol = pow(y1,k)*(y1-x1)*g_neu_pone(p,y2,x2,sqrt( (y1-x1)*(y1-x1)+a*a ),l);
//             if (isnan(g_neu_pone(p,y2,x2,sqrt( (y1-x1)*(y1-x1)+a*a ),l)))
//                 printf("g 1: %f, %f %f %f %f, \n", p,y2,x2,sqrt( (y1-x1)*(y1-x1)+a*a ),l); 
        }

        if(k != 0 && x1!=0){
            sol += k*x1*G_2_neu(p,y1,y2,x1,x2,a,k-1,l);           
//             if( isnan(G_2_neu(p,y1,y2,x1,x2,a,k-1,l)))
//                 printf("G_2_neu(p,y1,y2,x1,x2,a,k-1,l): %f %f %f %f %f %f %f %f  \n", p,y1,y2,x1,x2,a,k-1,l); 

        }

        if( fabs(y2-x2) !=0 ){
            sol += pow(y2,l)*(y2-x2)*g_neu_pone(p,y1,x1,sqrt( (y2-x2)*(y2-x2)+a*a ),k);
//             if (isnan(g_neu_pone(p,y1,x1,sqrt( (y2-x2)*(y2-x2)+a*a ),k)))
//                 printf("g 1: %f, %f %f %f %f, \n", p,y1,x1,sqrt( (y2-x2)*(y2-x2)+a*a ),k); 
        }


        if(l != 0 && x2 !=0 ){
            sol += l*x2*G_2_neu(p,y1,y2,x1,x2,a,k,l-1);  
//             if( isnan(G_2_neu(p,y1,y2,x1,x2,a,k,l-1)))
//                 printf("G2(p,y1,y2,x1,x2,a,k,l-1) sol normal 4:  %f mit %f %f %f %f  %f %f %f %f \n", sol,p,y1,y2,x1,x2,a,k,l-1); 

        }
        if( (p != 0) && (a !=0)  ){
            sol += 2*p*a*a*G_2_neu(p-1,y1,y2,x1,x2,a,k,l);
//             if( isnan(G_2_neu(p-1,y1,y2,x1,x2,a,k,l)))
//                 printf("G2(p-1,y1,y2,x1,x2,a,k,l))) sol normal 5: %f mit %f %f %f %f  %f %f %f %f \n", sol,p-1,y1,y2,x1,x2,a,k,l ); 
         }

        sol = sol/(k+l+2*p+2); 
    }
    else if(p == -3./2){
        if(k==0 && l==0){
    //         printf("a= %f, x1 = %f, x2 = %f, y1= %f, y2=%f \n ",a,x1,x2,y1,y2); 
		sol = sign((y1 - x1) * (y2 - x2)*a) / (2 *a)  * acos( (-2)*(y1-x1)*(y1-x1)*(y2-x2)*(y2-x2)/ ( ((y1-x1)*(y1-x1)+a*a)*((y2-x2)*(y2-x2)+a*a)) +1);
//             if( isnan(sol))
//                 printf("G2 sol k=0 l=0: %f \n", sol); 
             }
        else if(k==0 && l==1){
            sol = (-1)*asinh( (y1-x1) /  sqrt((y2-x2)*(y2-x2) +a*a) ) + x2*G_2_neu(p,y1,y2,x1,x2,a,0,0);
//             if( isnan(G_2_neu(p,y1,y2,x1,x2,a,0,0)))
//                 printf("G2(p,y1,y2,x1,x2,a,0,0) sol k=0 l=1: %f mit %f %f %f %f  %f %f %f %f \n", sol,p,y1,y2,x1,x2,a,zz,zz); 
//             if( isnan( asinh( (y1-x1) / ( sqrt((y2-x2)*(y2-x2) +a*a))) ) ){
//                 printf("asin is nan mit x1= %f, y1=%f, x2=%f y2=%f, a=%f \n",x1,y1,x2,y2,a);
//                 printf( "sqrt((y2-x2)*(y2-x2) +a*a)= %f", sqrt((y2-x2)*(y2-x2) +a*a));
//            }
        }
        else if(k==0 && l>1){
            sol = (y1-x1)*g_neu_pone(-1./2,y2,x2,sqrt((y1-x1)*(y1-x1)+a*a),l-2);
            sol += 2*x2*G_2_neu(p,y1,y2,x1,x2,a,0,l-1);
            sol += (-1)*(x2*x2 + a*a) * G_2_neu(p,y1,y2,x1,x2,a,0,l-2);    
//             if( isnan(G_2_neu(p,y1,y2,x1,x2,a,0,l-1)))
//                 printf("G_2_neu(p,y1,y2,x1,x2,a,0,l-1) sol k=0 l>1: %f mit %f %f %f %f  %f %f %f %f\n", sol,p,y1,y2,x1,x2,a,zz,l-1);
//              if( isnan(g_neu_pone(-1./2,y2,x2,sqrt((y1-x1)*(y1-x1)+a*a),l-2)))
//                 printf("g_neu_pone(-1./2,y2,x2,sqrt((y1-x1)*(y1-x1)+a*a),l-2) sol k=0 l>1: %f mit %f %f %f %f  %f \n", sol,-1./2,y2,x2,sqrt((y1-x1)*(y1-x1)+a*a),l-2);
//              if( isnan(G_2_neu(p,y1,y2,x1,x2,a,0,l-2)))
//                 printf("G_2_neu(p,y1,y2,x1,x2,a,0,l-2) sol k=0 l>1: %f mit %f %f %f %f  %f %f %f %f\n", sol,p,y1,y2,x1,x2,a,zz,l-2);
        }
        else{
            sol = (-1)*pow(y1,k-1)*g_neu_pone(-1./2,y2,x2,sqrt((y1-x1)*(y1-x1) + a*a),l);
            
            if(x1 != 0){
                sol += x1*G_2_neu(p,y1,y2,x1,x2,a,k-1,l);
            }
            if( isnan(sol)){
                printf("G2 sol else k=0: %f \n", sol);
//             if(isnan(g_neu_pone(-1./2,y2,x2,sqrt((y1-x1)*(y1-x1) + a*a),l)))
//                 printf("g_neu_pone(-1./2,y2,x2,sqrt((y1-x1)*(y1-x1) + a*a),l) sol k=0 l>1: %f mit %f %f %f %f  %f \n", sol,1./2,y2,x2,sqrt((y1-x1)*(y1-x1) + a*a),l);
            }        
            if (k>1){
                sol += (k-1)*G_2_neu(-1./2,y1,y2,x1,x2,a,k-2,l);                               
//              if( isnan(sol))
//                  printf("G2 sol else k>1: %f \n", sol) ;
            }
        }
    
        
    }
    else{
        printf("something went wrong with the p in G_2_neu \n");
    }

        
   
    
    
//    printf("sol: %f, \n", sol); 
    return sol;   
    
}  


/*
 * Stammfunktion G
 * 3 mal integrierte Kernfunktion int int 1 /|x-y| dy_3 dy_2 dy_1
 * 3 parameter l,k,n
 *  Formula p. 8 (25)
 */

double inline G_3_neu(double p, double x2, double y1, double y2, double d1,
        double d2,double a, double l, double m, double n) {
    
    //double x1_d1 = x1-d1   
    double sol = 0;
//     printf("G_3_neu: p= %f, l=%f, m=%f, n=%f \n", p,l,m,n);
    
	if(y2 !=0 )    
		sol += pow(y2,n+1)*G_2_neu(p,y1, x2, d1,y2+d2,a,l,m);


    
//     if(sol>0)
//         printf("!!!!!! first part in G_3: %f \n", G_2_neu(p,y1, x2, d1,y2+d2,a,l,m));

    sol += pow(x2,l)*G_2_neu(p,y1,y2,d1,x2-d2,a,m,n+1); 

    
//     if(sol>0)
//         printf("!!!!!! second part in G_3: %f mit  p= %f, y1=%f, y2=%f, d1 = %f , d2=%f, a=%f, m=%f, n=%f \n", G_2_neu(p,y1,y2,d1,x2-d2,a,m,n+1),p,y1,y2,d1,x2-d2,a,m,n+1);
            
    if(l!=0){
        sol += (-1)*l*G_3_neu(p,x2,y1,y2,d1,d2,a,l-1,m,n+1);
//         
//         if(sol>0)
//             printf("!!!!!! third part in G_3: %f \n", G_3_neu(p,x2,y1,y2,d1,d2,a,l-1,m,n+1));
    }
    
    
    sol = 1./(n+1)* sol;
    
    return sol;    
    
    
    
    
}


/*
 * Stammfunktion F
 * 4 mal integrierte Kernfunktion int int 1 /|x-y| dy dx  (ANALYTISCH)
 */
// the Formula is given in the Maischak Paper p.9 (26)
//
//
double inline F_par_pone(double x1, double x2, double y1, double y2, double d1,
		double d2, double a, double k, double l, double m, double n) {
    
    double p = -1./2;
    double sol = 0;

      if(k+l+m+n == 0)
  //       printf("F_par_pone: p= %f, x1=%f, x2=%f, y1=%f, y2=%f, d1=%f, d2=%f, a=%f, k=%f,l=%f, m=%f, n=%f \n", p,x1,x2,y1,y2,d1,d2,a,k,l,m,n);

      sol += pow(y1,m+1)*G_3_neu(p,x2,x1,y2,y1+d1,d2,a,l,k,n);

      if(isnan(sol))
        printf(" G_3_neu(p,x2, x1,y2,y1+d1,d2,a,l,k,n): %f %f %f %f %f %f %f %f %f %f \n", p,x2, x1,y2,y1+d1,d2,a,l,k,n);

      sol += pow(x1,k)*G_3_neu(p,x2,y1,y2,x1-d1,d2,a,l,m+1,n);    //x1 is not sure!! x2 should be y2, typo on in the paper
      if(isnan(sol))
        printf(" G_3_neu(p,y2,y1,y2,x1-d1,d2,a,l,m+1,n): %f %f %f %f %f %f %f %f %f %f \n", p,y2,y1,y2,x1-d1,d2,a,l,m+1,n);

      if(k!=0){
            sol += (-1)* k*F_par_pone(x1,x2,y1,y2,d1,d2,a,k-1,l,m+1,n);
            if(isnan(sol))
                printf(" F_par_pone(x1,x2,y1,y2,d1,d2,a,k-1,l,m+1,n): %f  %f %f %f %f %f %f %f %f %f %f \n",x1,x2,y1,y2,d1,d2,a,k-1,l,m+1,n);
        }


    sol = sol/(m+1);
    
    return sol;
    
}

/*
 * Stammfunktion F
 * 4 mal integrierte Kernfunktion int int 1 /|x-y| dy dx  (ANALYTISCH)
 */
double inline F_ort_pone(double x1, double x2, double y2, double y3, double y1_d1,
		double d2, double x3_d3, double k, double l, double m, double n) {
    //y1_d1 = y1 + d1
    //x3_d3 = x3 - d3    
    double p = -1./2;    
    double sol = 0;
    
//    printf("F_ort_pone: p= %f, k=%f,l=%f, m=%f, n=%f \n", p,k,l,m,n);
    
    
      sol += pow(y2,m+1)*G_3_neu(p,x1,x2,y3,y1_d1,y2+d2,-x3_d3,k,l,n); 
    
      sol += pow(x2,l)*G_3_neu(p,x1,y2,y3,y1_d1,x2-d2,-x3_d3,k,m+1,n);    //x1 is not sure!! typo on in the paper
      
      if(l!=0){    
            sol += (-1)* l*F_ort_pone(x1,x2,y2,y3,y1_d1,d2,-x3_d3,k,l-1,m+1,n);
        }
    
    
    sol = sol/(m+1);
    
    return sol;
    
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
 *für p1 
 */
template<double (F)(double, double, double, double, double, double, double, double, double, double, double)>
double inline intA4(double b, double d, double t, double v, double d1,
		double d2, double d3, double k, double l, double m, double n) {
    
   // printf("test f %f \n",		F(b, d, t, v, d1, d2, d3, k, l, m, n)	 );
    
	return    F(b, d, t, v, d1, d2, d3, k, l, m, n) - F(b, d, t, 0, d1, d2, d3, k, l, m, n)
			- F(b, d, 0, v, d1, d2, d3, k, l, m, n) + F(b, d, 0, 0, d1, d2, d3, k, l, m, n)
			- F(b, 0, t, v, d1, d2, d3, k, l, m, n) + F(b, 0, t, 0, d1, d2, d3, k, l, m, n)
			+ F(b, 0, 0, v, d1, d2, d3, k, l, m, n) - F(b, 0, 0, 0, d1, d2, d3, k, l, m, n)
			- F(0, d, t, v, d1, d2, d3, k, l, m, n) + F(0, d, t, 0, d1, d2, d3, k, l, m, n)
			+ F(0, d, 0, v, d1, d2, d3, k, l, m, n) - F(0, d, 0, 0, d1, d2, d3, k, l, m, n)
			+ F(0, 0, t, v, d1, d2, d3, k, l, m, n) - F(0, 0, t, 0, d1, d2, d3, k, l, m, n)
			- F(0, 0, 0, v, d1, d2, d3, k, l, m, n) + F(0, 0, 0, 0, d1, d2, d3, k, l, m, n);

}



/*
 * fuer parallele Elemente
 * voll analytisch
 */
double cParO1(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta, double k, double l, double m, double n) {
	return intA4<F_par_pone>(b, d, t, v, d1, d2, d3, k, l, m, n);
}
/*
 * fuer orthogonale Elemente
 * voll analytisch
 */
double cOrtO1(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta, double k, double l, double m, double n) {
	return intA4<F_ort_pone>(b, d, t, v, d1, d2, d3, k, l, m, n);
}

/*
 * fuer parallele Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cParO2(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta, double k, double l, double m, double n) {
	double tmp = 0;

	//kurze Seite nach vorn
	switch_site_par(b, d, t, v, d1, d2);

	if (zeta[0] * zeta[0] * (t * t + v * v) < dist2_par(b, d, t, v, d1, d2, d3))
		return intQ4<f_par>(b, d, t, v, d1, d2, d3);

	return intA4<F_par_pone>(b, d, t, v, d1, d2, d3,k,l,m,n);
}

/*
 * fuer orthogonale Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cOrtO2(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta, double k, double l, double m, double n) {
	double tmp = 0, dis2 = 0;

	//kurze Seite nach vorn
	switch_site_ort(b, d, t, v, d1, d2, d3);

	if (zeta[0] * zeta[0] * (t * t + v * v) < dist2_ort(b, d, t, v, d1, d2, d3))
		return intQ4<f_ort>(b, d, t, v, d1, d2, d3);

	return intA4<F_ort_pone>(b, d, t, v, d1, d2, d3,k,l,m,n);
}

/*
 * fuer parallele Elemente
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cParO3(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta, double k, double l, double m, double n) {

	//kurze Seite nach vorn
	switch_site_par(b, d, t, v, d1, d2);

	if (zeta[1] * zeta[1] * (b * b + d * d) < dist2_par(b, d, t, v, d1, d2, d3))
		return intQ2A2<FY1Y2_par>(b, d, t, v, d1, d2, d3);

	if (zeta[0] * zeta[0] * (t * t + v * v) < dist2_par(b, d, t, v, d1, d2, d3))
		return intQ4<f_par>(b, d, t, v, d1, d2, d3);

	return intA4<F_par_pone>(b, d, t, v, d1, d2, d3,k,l,m,n);
}

/*
 * fuer orthogonale Elemente
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * volle Quadratur fuer zeta_Q zulaessige Elemente
 * analytisch sonst
 */
double cOrtO3(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta, double k, double l, double m, double n) {
	double tmp = 0, dis2 = 0;

	//kurze Seite nach vorn
	switch_site_ort(b, d, t, v, d1, d2, d3);

	if (zeta[1] * zeta[1] * (b * b + d * d) < dist2_ort(b, d, t, v, d1, d2, d3))
		return intQ2A2<FY2Y3_ort>(b, d, t, v, d1, d2, d3);

	if (zeta[0] * zeta[0] * (t * t + v * v) < dist2_ort(b, d, t, v, d1, d2, d3))
		return intQ4<f_ort>(b, d, t, v, d1, d2, d3);

	return intA4<F_ort_pone>(b, d, t, v, d1, d2, d3,k,l,m,n);
}

/*
 * fuer parallele Elemente
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * analytisch sonst
 */
double cParO4(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta, double k, double l, double m, double n) {

	//kurze Seite nach vorn
	switch_site_par(b, d, t, v, d1, d2);

	if (zeta[1] * zeta[1] * (b * b + d * d) < dist2_par(b, d, t, v, d1, d2, d3))
		return intQ2A2<FY1Y2_par>(b, d, t, v, d1, d2, d3);

	return intA4<F_par_pone>(b, d, t, v, d1, d2, d3,k,l,m,n);
}
/*
 * fuer orthogonale Elemente
 * Quadratur ueber ein Element fuer zeta_E zulaessige Elemente
 * analytisch sonst
 */
double cOrtO4(double b, double d, double t, double v, double d1, double d2,
		double d3, double* zeta, double k, double l, double m, double n) {
	double tmp = 0, dis2 = 0;

	//kurze Seite nach vorn
	switch_site_ort(b, d, t, v, d1, d2, d3);

	if (zeta[1] * zeta[1] * (b * b + d * d) < dist2_ort(b, d, t, v, d1, d2, d3))
		return intQ2A2<FY2Y3_ort>(b, d, t, v, d1, d2, d3);

	return intA4<F_ort_pone>(b, d, t, v, d1, d2, d3,k,l,m,n);
}





/*
 *
 *
 *
 *Functions for P0 elements
 *
 *
 */

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


