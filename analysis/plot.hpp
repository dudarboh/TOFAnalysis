#include "TF1.h"
#include <string>

// in mm/ns
#define SPEED_OF_LIGHT 299.792458
// const double rInner = 1804.8;


double calibrateTOF(double tof, int nHits,string method){
    if ( method == "Closest"){
        TF1 calib = TF1("closest", "[0]+[1]/(x-[2])");
        calib.SetParameters(1.00005e+00, 8.27272e-04, -3.37065e+00);
        return tof/calib.Eval(nHits);
    }
    else if (method == "Fastest"){
        TF1 calib = TF1("fastest", "[0]+[1]/(x-[2])");
        calib.SetParameters(1.00002e+00, 6.16835e-04, -1.15474e+00);
        return tof/calib.Eval(nHits);
    }
    else if (method == "Frank"){
        TF1 calib = TF1("frank", "pol2");
        if(nHits < 45) calib.SetParameters(9.99898e-01, -5.00921e-06, 6.92759e-08);
        else calib.SetParameters(9.99753e-01, 1.45337e-06, -2.28934e-09);
        return tof/calib.Eval(nHits);
    }
    else if (method == "Cyl"){
        TF1 calib = TF1("cyl", "[0] + [1]*sqrt(x)");
        calib.SetParameters(9.99883e-01, 2.16310e-05);
        return tof/calib.Eval(nHits);
    }
    else throw std::string("Wrong method passed to the calibrateTOF");
}
