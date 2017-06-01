#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include "sigma.h"
#include <string>
#include <fstream>

//To compile and run this program in the terminal
//>> g++ -c  sigma.cpp sigma.h `lhapdf-config --cflags --ldflags`
//>> g++ main.cpp sigma.o -o kjor `lhapdf-config --cflags --ldflags`
//>> ./kjor

using namespace std;

//A program for calculating the cross section for squark pair production at the LHC.
//Energy units are [GeV], and the final cross section is given in [fb^(-1)].

int main()
{
    //Center of mass energy squared, s
    double s = 7000*7000;

    //Number of integration points N
    double N = 1000;

    //Minimum x-value
    double x_min = 0.00001;

    cout << "x_min : " << x_min << endl;

    //Maximum x_value
    double x_max = 1.0;

    sigma my_sigma(s);

    //Testing sigma for 7000 GeV
    //-----------------------------------------------
//    my_sigma.m_mq = 1000;
//    double crossection = my_sigma.integrate_simpson(x_max, x_min, N);
//    cout << "The final cross section for mg = 800 GeV, mq = 1000 GeV and sqrt(s)= 7000: " << crossection << endl;
    //-----------------------------------------------


    //Plot xf function
    //-----------------------------------------------
    //Open the output file
//    string filename = "test_fx_pid1.txt";
//    ofstream ofile;
//    if (!ofile.is_open()){
//        ofile.open(filename.c_str(), ofstream::out);
//        if (!ofile.good()){
//            cout << "Error opening file " << filename << ", aborting!" << endl;
//            return 0;
//        }
//    }

//    const string setname = "cteq66";
//    const PDF* pdf = mkPDF(setname, 0);
//    double q2 = my_sigma.m_mq*my_sigma.m_mq; //GeV - from article
//    int my_pid = 2;

//    ofile << "x" << " " << "fx" << endl;
//    for (double x_ = 0; x_ < 1; x_+= 0.00001){
//        const double xf_test = pdf->xfxQ2(my_pid, x_, q2);
//        ofile << x_ << " " << xf_test << endl;

//    }

//    ofile.close();


//    //-----------------------------------------------

//    //Write cross sections for several s to file

//    //Open the output file
//    string filename = "sd_vary_mq_mg2000.txt";
//    ofstream ofile;
//    if (!ofile.is_open()){
//        ofile.open(filename.c_str(), ofstream::out);
//        if (!ofile.good()){
//            cout << "Error opening file " << filename << ", aborting!" << endl;
//            return 0;
//        }
//    }

//    my_sigma.m_mg = 2000; //GeV
//    ofile << "Integrated using Simpsons method using N=1000 points." << endl;
//    ofile << "s Gluino mass set to 1200 GeV" << endl;
//    ofile << "m_q " << "s " << " sigma " << endl;

//    for (double mq = 800; mq < 2100; mq += 200) {
//        my_sigma.m_mq = mq;

//        for (double sqrt_s = 3000; sqrt_s <= 13000; sqrt_s += 500){
//            double my_s = sqrt_s*sqrt_s;
//            my_sigma.m_s = my_s;
//            cout << "Center of mass energy " << sqrt_s << endl;
//            double my_cross = my_sigma.integrate_simpson(x_max, x_min, N);
//            ofile << mq << " " << sqrt_s << " " << my_cross << endl;
//        }
//    }

//    ofile.close();
//    //--------------------------------------------


    //-----------------------------------------------

    //Write cross sections for several s to file

    //Open the output file
    string filename = "sd_sigma_low.txt";
    ofstream ofile;
    if (!ofile.is_open()){
        ofile.open(filename.c_str(), ofstream::out);
        if (!ofile.good()){
            cout << "Error opening file " << filename << ", aborting!" << endl;
            return 0;
        }
    }

    my_sigma.m_mg = 2000; //GeV
    my_sigma.m_mq = 800; //GeV
    ofile << "Integrated using Simpsons method using N=1000 points." << endl;
    ofile << "s Gluino mass set to 2000 GeV, squark mass: 800 GeV" << endl;
    ofile << "s " << " sigma " << endl;

    for (double sqrt_s = 2000; sqrt_s <= 10000; sqrt_s += 100){
        double my_s = sqrt_s*sqrt_s;
        my_sigma.m_s = my_s;
        cout << "Center of mass energy " << sqrt_s << endl;
        double my_cross = my_sigma.integrate_simpson(x_max, x_min, N);
        //double my_cross = my_sigma.sigma_diff(1,1,my_s)*1/(2.5819)*pow(10,9);
        ofile << sqrt_s << " " << my_cross << endl;
    }

    ofile.close();
    //--------------------------------------------




    return 0;
}

