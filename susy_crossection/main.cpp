#include <iostream>
#include "sigma.h"
#include <string>
#include <fstream>

using namespace std;

//A program for calculating the cross section for squark pair production at the LHC.
//Energy units are [GeV], and the final cross section is given in [fb^(-1)].

int main()
{
    cout << "Hello Subatomic World!" << endl;

    //Center of mass energy squared
    double s = 1000*1000;

    //Number of integration points N
    double N = 10;

    //Minimum x-value
    double x_min = 0.05;

    //Maximum x_value
    double x_max = 1.0;

    sigma my_sigma(s);

    cout << "Well, it looks like we have an instance called my_sigma" << endl;

    //integrate_simpson(x_max, x_min, N)

    my_sigma.different_quarks();

    double cross_section = my_sigma.integrate_simpson(x_max, x_min, N);

    cout << "And the cross section for center of mass energy s = " << s << " sigma = " << cross_section << endl;

    //Write cross sections for several s to file

    //Open the output file
    string filename = "susy_data1.txt";
    ofstream ofile;
    if (!ofile.is_open()){
        ofile.open(filename.c_str(), ofstream::out);
        if (!ofile.good()){
            cout << "Error opening file " << filename << ", aborting!" << endl;
            return 0;
        }
    }

    ofile << "s " << " sigma " << endl;

    for (double sqrt_s = 500; sqrt_s <= 3000; sqrt_s += 500){
        double my_s = sqrt_s*sqrt_s;
        my_sigma.m_s = my_s;
        cout << "Center of mass energy " << sqrt_s << endl;
        double my_cross = my_sigma.integrate_simpson(x_max, x_min, N);
        ofile << sqrt_s << " " << my_cross << endl;
    }

    ofile.close();

    return 0;
}

