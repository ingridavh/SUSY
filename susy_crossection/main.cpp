#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include "sigma.h"
#include <string>
#include <fstream>

//To compile and run this program in the terminal
//>> g++ -c  sigma.cpp sigma.h `lhapdf-config --cflags --ldflags`
//>> g++ main.cpp sigma.o -o kjor `lhapdf-config --cflags --ldflags`
//>> ./kjor

using namespace LHAPDF;
using namespace std;

//A program for calculating the cross section for squark pair production at the LHC.
//Energy units are [GeV], and the final cross section is given in [fb^(-1)].

int main()
{
    cout << "Hello Subatomic World!" << endl;

    double s = 7000*7000;

    //Number of integration points N
    double N = 100;

    //Minimum x-value
    double x_min = 1e-5;

    //Maximum x_value
    double x_max = 1.0;

    sigma my_sigma(s);

//    cout << "Well, it looks like we have an instance called my_sigma" << endl;

//    my_sigma.m_mq = 1000;
//    double crossection = my_sigma.integrate_simpson(x_max, x_min, N);
//    cout << "The final cross section for mg = 800 GeV, mq = 1000 GeV and sqrt(s)= 7000: " << crossection << endl;

    //Test for single cross section with mq = 1000 GeV and sqrt(s) = 7000

    my_sigma.m_mq = 1000;
    for (int s_ = 1000; s_ < 13000; s_ += 200){
        my_sigma.m_s = s_*s_;
        cout << my_sigma.sigma_diff(1,1,s_*s_) << ", ";

    }
    cout << "end" << endl;

    //Plot xf function
    const string setname = "cteq66";
    const PDF* pdf = mkPDF(setname, 0);
    vector<int> pids = pdf->flavors();
    double q2 = my_sigma.m_mq*my_sigma.m_mq; //GeV - from article
    int pid = pids[1];

    for (double x_ = 0; x_ < 1; x_+= 0.01){
        const double xf_test = pdf->xfxQ2(pid, x_, q2);
        cout << xf_test << ", ";

    }

    cout << "end 2" << endl;


    //Write cross sections for several s to file

    //Open the output file
//    string filename = "sd_vary_mq_mg800.txt";
//    ofstream ofile;
//    if (!ofile.is_open()){
//        ofile.open(filename.c_str(), ofstream::out);
//        if (!ofile.good()){
//            cout << "Error opening file " << filename << ", aborting!" << endl;
//            return 0;
//        }
//    }

//    ofile << "Integrated using Simpsons method using N=100 points." << endl;
//    ofile << "s Gluino mass set to 800 GeV" << endl;
//    ofile << "m_q " << "s " << " sigma " << endl;

//    for (double mq = 600; mq < 2100; mq += 200) {
//        my_sigma.m_mq = mq;

//        for (double sqrt_s = 5000; sqrt_s <= 13000; sqrt_s += 200){
//            double my_s = sqrt_s*sqrt_s;
//            my_sigma.m_s = my_s;
//            cout << "Center of mass energy " << sqrt_s << endl;
//            double my_cross = my_sigma.integrate_simpson(x_max, x_min, N);
//            ofile << mq << " " << sqrt_s << " " << my_cross << endl;
//        }
//    }

//    ofile.close();

//    cout << "Error count " << my_sigma.m_error_count << endl;

    return 0;
}

