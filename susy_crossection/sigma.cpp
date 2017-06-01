#include "LHAPDF/LHAPDF.h"
#include "sigma.h"
#include <armadillo>
#include <math.h>

using namespace LHAPDF;
using namespace arma;
using namespace std;

sigma::sigma(double s) :
    m_s(s),
    m_mg(800.), //From article-022
    m_mq(1000.), //Lower limit
    m_pi(arma::datum::pi),
    m_alpha(0.1148), //Value for alpha(m_Z)
    m_error_count(0),
    m_no_prob(0),
    m_pid1(1),
    m_pid2(2)
{
    m_m2 = m_mg*m_mg-m_mq*m_mq;
}

double sigma::sigma_diff(double x1, double x2, double s_hat){
    //Differential cross section for quark quark - squark squark
    //To be integrated over convolution integral

    double betaq; double L1;
    //Get rid of negative arguments
    if ((1-4.*m_mq*m_mq/s_hat) <= 0)
    {
        //cout << "The argument to beta and L1 is negative, so it gives no contribution." << endl;
        //cout << "Error for: x1= " << x1 << ", x2 = " << x2 << endl;
        betaq = 0;
        L1 = 0;

        m_error_count += 1;
    } else {
        betaq = sqrt(1-4.*m_mq*m_mq/s_hat);
        L1 = log((s_hat + 2.*m_m2 - s_hat*betaq)/(s_hat + 2.*m_m2 + s_hat*betaq));

        m_no_prob += 1;

    }


    //double L1 = log((s_hat + 2.*m_m2 - s_hat*betaq)/(s_hat + 2.*m_m2 + s_hat*betaq));

    //Set last term, only valid for i=j (t-channel)
    double A = 1;
    if (m_pid1 == m_pid2){
        A = (m_pi*m_alpha*m_alpha/s_hat)*(8.*m_mg*m_mg/(27.*(s_hat+2*m_m2)))*L1;
    }    else {
        A = 0;
    }



    double dsigma = (m_pi*m_alpha*m_alpha/s_hat)*(betaq*(-4./9.-4.*m_m2*m_m2/(9.*(m_mg*m_mg*s_hat+m_m2*m_m2)))+ (-4./9.-8.*m_m2/(9.*s_hat))*L1)+A;

    if (dsigma < 0) {
        cout << "A : " << A << endl;
        cout << "betaq : " << betaq << endl;
    };


    return dsigma;
}

double sigma::f(double x)
    //Parton distribution function (PDF)
{
    //just a sample, change this later
    return pow((1-x),4);
}

double sigma::integrate_simpson(double x_max, double x_min, double N){
    //Use Rectangle method to integrate over f(x1) and f(x2)

    double h = abs(x_max - x_min)/N;

    //Set right parameters for LHAPDF
    const string setname = "cteq66";
    const PDF* pdf = mkPDF(setname, 0);
    double q2 = m_mq*m_mq; //GeV - from article

    cout << "You are now inside the Simpson integrator, with step size " << h << endl;

    //integration variables
    double x1; double x2; double s_hat;
    double F = 0;

    const double xf_test = pdf->xfxQ2(m_pid1, x_min, q2);

    cout << "The minimum xf-value is " << xf_test << endl;
    cout << "The minimum sigma-value is " << sigma_diff(x_min, x_min, x_min*x_min*m_s) << endl;

    cout << "m_q= " << m_mq << ", m_g= " << m_mg << ", s= " << m_s << endl;

    //Coefficient to account for not counting ingoing quarks twice
    double coeff;

    for (int i = 0; i < N; i++ ) {
        x1 = (i+0.5)*h;
        for (int j=0; j < N ; j++ ) {
            x2 = (j+0.5)*h;
            s_hat = m_s*x1*x2;

            //sum over initial quark flavors, pid=1,2,3,4
            for (int pid1 = 1; pid1 < 5; pid1++){
                for (int pid2 = 1; pid2 < 5; pid2++) {
                    if (pid1 == pid2) { coeff = 1.;} else {coeff = 0.5;}

                    m_pid1 = pid1; m_pid2 = pid2;

                    double xf1 = pdf->xfxQ2(m_pid1, x1, q2);
                    double xf2 = pdf->xfxQ2(m_pid2, x2, q2);

                    if (coeff*xf1*xf2*sigma_diff(x1, x2, s_hat) < 0) {
                        cout << "sigma : " << sigma_diff(x1, x2, s_hat) << endl;
                        cout << "x1 : " << x1 << ", x2 : " << x2 << endl;
                    };

                    F += coeff*xf1*xf2*sigma_diff(x1, x2, s_hat)/(x1*x2);
                }
            }
        }
    }

    F = F*h*h;
    F *= 1/(2.5819)*pow(10,9); //Convert to picobarn
    return F;
}
