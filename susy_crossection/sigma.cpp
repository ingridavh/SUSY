#include "sigma.h"
#include <armadillo>
#include <math.h>

using namespace arma;
using namespace std;

sigma::sigma(double s) :
    m_s(s),
    m_mg(400.), //Lower limit
    m_mq(400.), //Lower limit
    m_pi(arma::datum::pi),
    m_alpha(0.1184), //Value for aplha(m_Z)
    m_i(1)
{
    m_m2 = m_mg*m_mg-m_mq*m_mq;
}

double sigma::sigma_diff(double x1, double x2, double s_hat){
    //Differential cross section for quark quark - squark squark
    //To be integrated over convolution integral

    //Get rid of negative arguments
    if ((1-4.*m_mq*m_mq/s_hat) < 0)
    {
        //cout << "The argument to beta and L1 is negative, so it gives no contribution." << endl;
        //cout << "Error for: x1= " << x1 << ", x2 = " << x2 << endl;
        return 0;
    }

    //Set variables that depend on s
    double betaq = sqrt(1-4.*m_mq*m_mq/s_hat);

    double L1 = log((s_hat + 2.*m_m2 - s_hat*betaq)/(s_hat + 2.*m_m2 + s_hat*betaq));

    //Set last term, only valid for i=j (t-channel)
    double A = 1;
    if (m_i == 1){
        A = (m_pi*m_alpha*m_alpha/s_hat)*(8.*m_mg*m_mg/(27.*(s_hat+2*m_m2)))*L1;
    }    else {
        A = 0;
    }

    double dsigma = (m_pi*m_alpha*m_alpha/s_hat)*(betaq*(-4./9.-4.*m_m2*m_m2/(9.*(m_mg*m_mg*s_hat+m_m2*m_m2)))+ (-4./9.-8.*m_m2/(9.*s_hat))*L1)+A;

    return dsigma;
}


double sigma::f(double x)
    //Parton distribution function (PDF)
{
    //just a sample, change this later
    return pow((1-x),4);
}

double sigma::integrate_simpson(double x_max, double x_min, double N){
    //Use Simpson's method to integrate over f(x1) and f(x2)

    double h = (x_max - x_min)/N;

    cout << "You are now inside the Simpson integrator, with step size " << h << endl;
    if (m_i == 1) {
        cout << "This cross section is for same flavor initial quarks." << endl;
    } else {
        cout << "This cross section is for different flavor intial quarks." << endl;
    }

    double fa = f(x_min)*f(x_min)*sigma_diff(x_min, x_min, m_s*x_min*x_min);
    double fc = f(x_max)*f(x_max)*sigma_diff(x_max, x_max, m_s*x_max*x_max);

    //integration variables
    double x1; double x2; double s_hat;
    double F; double coeff;

    for (int i = 1; i < N; i++ ) {
        x1 = x_min + i*h;

        for (int j=0; j < N+1 ; j++ ) {
            x2 = x_min + j*h;

            s_hat = m_s*x1*x2;

            if (j==0 && i%2 == 1) {coeff = 4;} else if (j==0 && i%2 == 0){ coeff = 2;}
            else if (j==N && i%2 == 1) {coeff = 4;} else if (j==N && i%2==0) { coeff = 2;}
            else if (i%2 == 1 && j%2 == 1) { coeff = 16;}
            else if (i%2 == 1 && j%2 == 0) { coeff = 8;}
            else if (i%2 == 0 && j%2 == 1) { coeff = 8;}
            else if (i%2 == 0 && i%2 == 0) {coeff = 4;}

            F += coeff*f(x1)*f(x2)*sigma_diff(x1, x2, s_hat);
        }
    }

    F += fa + fc;
    F /= h*h;
    F *= 1/(2.5819)*pow(10,12); //Convert to femtobarn^(-1)
    return F;
}

void sigma::same_quarks()
{
    m_i = 1;
}

void sigma::different_quarks()
{
    m_i = 0;
}

