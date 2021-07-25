#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

void _rinvGauss(double *rig, double *runif, double* rnorm, int n, double mu, double lambda){
    int i;
	double v, x, b = 0.5*mu/lambda, a = mu*b, c = 4.0*mu*lambda, d = mu*mu;

	for (i=0; i<n; i++) {
		if (mu<=0 || lambda<=0)
			rig[i] = 0.0;
		else {
			v   = rnorm[i]*rnorm[i];
			x   = mu + a*v - b*sqrt(c*v + d*v*v);
			rig[i]=(runif[i]<(mu/(mu+x)))?x:d/x;
			if (rig[i]<0.0)
				v = x;
		}
	}
}

double mean(double *x, int n)
{
	int i;
	double avg = 0.0;
	for(i=0;i<n;i++)
	{
		avg += x[i];
	}
	avg /= n;
	return avg;
}

void _muhat_Poisson(double *x, int n, double *mu)
{

    double mu0;

  	mu0     = mean(x, n);
    mu[0]   = mu0;
  	mu[1]   = mu0*(mu0+1);
  	mu[2]   = mu0*(pow(mu0,2) + 3*mu0+1);
  	mu[3]   = mu0*(pow(mu0,3) + 6*pow(mu0,2) + 7*mu0 + 1);
  	mu[4]   = mu0*(pow(mu0,4) + 10*pow(mu0,3) + 25*pow(mu0,2) + 15*mu0+1);
  	mu[5]   = mu0*(pow(mu0,5) + 15*pow(mu0,4) + 65*pow(mu0,3) + 90*pow(mu0,2) + 31*mu0+1);
}

void _muhat_IG(double *x, int n, double *mu)
{
	int k;
	double lambda;
	double *y;
	y = (double*)malloc(sizeof(double)*n);
	for(k=0;k<n;k++)
	{
		y[k] = 1.0/(x[k]);
	}
  	lambda = 1.0/(mean(y, n)-1/mean(x, n));
  	mu[0] = mean(x, n);
  	mu[1] = (mu[0]+lambda)*pow(mu[0],2)/lambda;
  	for(k=2;k<6;k++)
  	{
        mu[k] = pow(mu[0],2)*mu[k-2]+pow(mu[0],2)*(2*k-1)*mu[k-1]/lambda;
  	}
    free(y);
}

double _sigma_hat(double *L, double gamma)
{
	double sigma_hat;
	sigma_hat = pow(L[0],2)*L[5]
               + (4*gamma-6)*pow(L[0],3)*L[4]
               - 4*gamma*L[0]*L[1]*L[4]
               + (2*gamma-5)*(2*gamma-5)*pow(L[0],4)*L[3]
               - (4*gamma-4)*(2*gamma-3)*pow(L[0],2)*L[1]*L[3]
               + 4*gamma*gamma*pow(L[1],2)*L[3]
               + 2*L[0]*L[2]*L[3]
               + (1-4*gamma)*L[1]*pow(L[2],2)
               - (16*gamma*gamma-40*gamma+12)*L[0]*pow(L[1],2)*L[2]
               + 32*(gamma-1)*(gamma-2)*pow(L[0],3)*L[1]*L[2]
               + (4*gamma-10)*pow(L[0],2)*pow(L[2],2)
               + (16-8*gamma)*(2*gamma-3)*pow(L[0],5)*L[2]
               + 12*(gamma-1)*(2*gamma-3)*pow(L[0],2)*pow(L[1],3)
               + (35-18*gamma)*(2*gamma-3)*pow(L[0],4)*pow(L[1],2)
               - 4*gamma*gamma*pow(L[1],4)
               + (8-4*gamma)*(8-4*gamma)*pow(L[0],6)*L[1];
	return sigma_hat;
}

void cal_L(double *x, double *L, int n)
{
    int i;
    double tmp;

    for(i=0; i<4; i++) L[i] = 0.0;

    for (i = 0; i < n; i++)
    {
        tmp     = x[i];
        L[0]    += tmp;
        tmp     *= x[i];
        L[1]    += tmp;
        tmp     *= x[i];
        L[2]    += tmp;
        tmp     *= x[i];
        L[3]    += tmp;
    }
}

double leave_one_S_hat(double xleave, double *L, int n, double *coe)
{
    double S_hat, L1, L2, L3, L4, l1, l2, l3, l4;

    S_hat   = xleave;
    L1      = L[0] - S_hat;
    S_hat   *= xleave;
    L2      = L[1] - S_hat;
    S_hat   *= xleave;
    L3      = L[2] - S_hat;
    S_hat   *= xleave;
    L4      = L[3] - S_hat;
    l1      = L1/(n-1.0);
    l2      = L2/(n-1.0);
    l3      = L3/(n-1.0);
    l4      = L4/(n-1.0);

    S_hat   = coe[0]*pow(l1, 4);
    S_hat   += coe[1]*l2*l2;
    S_hat   += coe[2]*l1*l1*l2;
    S_hat   += coe[3]*l4;
    S_hat   += coe[4]*l1*l3;

    return S_hat;
}

double S_minus(double *x, int n, double gamma, double *sminus)
{
    int i;
    double *coe, *L, L0;

    n--;
    L       = (double*)malloc(sizeof(double)*4);
    coe     = (double*)malloc(sizeof(double)*5);
    coe[0]  = (2-gamma)*(n/(n-1.0))*(n/(n-2.0))*(n/(n-3.0));
    coe[1]  = 3*(n/(n-2.0))*(1/(n-3.0))-gamma*(n/(n-1.0))*((n*n-3*n+3)/((n-2.0)*(n-3.0)));
    coe[2]  = -3*((n+1.0)/(n-1.0))*(n/(n-2.0))*(n/(n-3.0))+2*gamma*(n/(n-1.0))*(n/(n-2.0))*(n/(n-3.0));
    coe[3]  = -(n/(n-1.0))*((n+1.0)/(n-2.0))*(1/(n-3.0))+gamma*(n/(n-2.0))*(1/(n-3.0));
    coe[4]  = (n/(n-1.0))*((n*n+n+4)/((n-2.0)*(n-3.0)))-4*gamma*(n/(n-2.0))*(1/(n-3.0));
    n++;


    cal_L(x, L, n);


    for (i = 0; i < n; i++)
    {
        sminus[i] = leave_one_S_hat(x[i], L, n, coe);
    }

    L0 = L[0];
    free(coe);
    free(L);

    return L0;
}

void leave_one_x_bar(double *x, double *xbar, int n, double L1)
{
    int i, j;

    j = n-1;
    for (i = 0; i < n; i++) {
        xbar[i] = (L1 - x[i])/j;
    }
}

double exp_approx(double x){
    x = 1.0 + x/256;
    for (int i = 0; i < 4; i++){
        x *= x;
    }
    return x;
}

double T_n_w_1(double *x, int n, double gamma, double a) //nT
{
    double tmp, L0, T = 0.0;
    int i, j;
    double *xbar, *Sminus;

    xbar    = (double*)malloc(sizeof(double)*n);
    Sminus  = (double*)malloc(sizeof(double)*n);


    L0      = S_minus(x, n, gamma, Sminus);
    leave_one_x_bar(x, xbar, n, L0);


    for (i = 0; i < n; i++)
        T += Sminus[i]*Sminus[i];

    for (i = 0; i < n-1; i++) {
        for (j = i+1; j < n; j++) {
            tmp = (xbar[i]-xbar[j])/2.0;
            tmp = exp(-tmp*tmp/a);
            T   += Sminus[i]*Sminus[j]*tmp*2;
        }
    }

    T /= n;

    free(xbar);
    free(Sminus);

    return T;
}

double T_n_w_1B(double *x, int n, double gamma, double a) //nT
{
    double tmp, L0, T = 0.0;
    int i, j;
    double *xbar, *Sminus;

    xbar    = (double*)malloc(sizeof(double)*n);
    Sminus  = (double*)malloc(sizeof(double)*n);


    L0      = S_minus(x, n, gamma, Sminus);
    leave_one_x_bar(x, xbar, n, L0);


    for (i = 0; i < n; i++)
        T += Sminus[i]*Sminus[i];
    if(n>1000){
        for (i = 0; i < n-1; i++) {
            for (j = i+1; j < n; j++) {
                tmp = (xbar[i]-xbar[j])/2.0;
                tmp = exp_approx(-tmp*tmp/a);
                T   += Sminus[i]*Sminus[j]*tmp*2;
            }
        }
    }
    else{
        for (i = 0; i < n-1; i++) {
            for (j = i+1; j < n; j++) {
                tmp = (xbar[i]-xbar[j])/2.0;
                tmp = exp(-tmp*tmp/a);
                T   += Sminus[i]*Sminus[j]*tmp*2;
            }
        }
    }

    T /= n;

    free(xbar);
    free(Sminus);

    return T;
}

double T_n_w_2(double *x, int n, double gamma, double a)
{
    double tmp, L0, T = 0.0;
    int i, j;
    double *xbar, *L, *Sminus;

    xbar    = (double*)malloc(sizeof(double)*n);
    Sminus  = (double*)malloc(sizeof(double)*n);


    L0      = S_minus(x, n, gamma, Sminus);
    leave_one_x_bar(x, xbar, n, L0);



    for (i = 0; i < n; i++)
        T += Sminus[i]*Sminus[i];

    for (i = 0; i < n-1; i++) {
        for (j = i+1; j < n; j++) {
            tmp = (xbar[i]-xbar[j]);
            tmp = a*a/(a*a + tmp*tmp);
            T   += Sminus[i]*Sminus[j]*tmp*2;
        }
    }

    T /= n;

    free(xbar);
    free(Sminus);

    return T;
}

SEXP RINV_GAUSS(SEXP RUNIF, SEXP RNORM, SEXP N, SEXP MU, SEXP LAM)
{
    int n = INTEGER(N)[0];
	SEXP rIG;
	PROTECT(rIG = allocVector(REALSXP, n));
    _rinvGauss(REAL(rIG), REAL(RUNIF), REAL(RNORM), n, REAL(MU)[0], REAL(LAM)[0]);

    UNPROTECT(1);
	return rIG;
}

SEXP sigma_hat(SEXP L, SEXP GAMMA)
{
	SEXP Sigma;
	PROTECT(Sigma = allocVector(REALSXP, 1));
	REAL(Sigma)[0] = _sigma_hat(REAL(L), REAL(GAMMA)[0]);
    UNPROTECT(1);
	return Sigma;
}

SEXP muhat_Ig(SEXP X, SEXP N)
{
    SEXP Mu;
    PROTECT(Mu = allocVector(REALSXP, 6));

    _muhat_IG(REAL(X), INTEGER(N)[0], REAL(Mu));

    UNPROTECT(1);
    return Mu;
}

SEXP muhat_poisson(SEXP X, SEXP N)
{
    SEXP Mu;
    PROTECT(Mu = allocVector(REALSXP, 6));

    _muhat_Poisson(REAL(X), INTEGER(N)[0], REAL(Mu));

    UNPROTECT(1);
    return Mu;
}

SEXP _Tnw(SEXP X, SEXP N, SEXP GAMMA, SEXP A, SEXP WEIGHT)
{
    SEXP Test;
    PROTECT(Test = allocVector(REALSXP, 1));

    if (INTEGER(WEIGHT)[0] == 1)
        REAL(Test)[0] = T_n_w_1(REAL(X), INTEGER(N)[0], REAL(GAMMA)[0], REAL(A)[0]);
    else
        REAL(Test)[0] = T_n_w_2(REAL(X), INTEGER(N)[0], REAL(GAMMA)[0], REAL(A)[0]);

    UNPROTECT(1);
    return Test;
}

SEXP _TnwB(SEXP X, SEXP N, SEXP GAMMA, SEXP A, SEXP WEIGHT)
{
    SEXP Test;
    PROTECT(Test = allocVector(REALSXP, 1));

    if (INTEGER(WEIGHT)[0] == 1)
        REAL(Test)[0] = T_n_w_1B(REAL(X), INTEGER(N)[0], REAL(GAMMA)[0], REAL(A)[0]);
    else
        REAL(Test)[0] = T_n_w_2(REAL(X), INTEGER(N)[0], REAL(GAMMA)[0], REAL(A)[0]);

    UNPROTECT(1);
    return Test;
}