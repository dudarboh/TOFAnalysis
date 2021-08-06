#ifndef TofAnaUtils_h
#define TofAnaUtils_h 1

#include "EVENT/ReconstructedParticle.h"
#include <utility>
#include "TString.h"
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"

class TofMethod{
    public:
    TofMethod(double smearing){
        _smearing = smearing;
    }

    virtual double calculate( EVENT::ReconstructedParticle* ) {return 0.;};

    TString _branch_name;
    double _tof, _smearing;
};

class TofClosest : public TofMethod{
    public:
    TofClosest(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_closest_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo );
};

class TofFastest : public TofMethod{
    public:
    TofFastest(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_fastest_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo );
};

class TofFrankFit : public TofMethod{
    public:
    TofFrankFit(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_frank_fit_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo );
};

class TofFrankAvg : public TofMethod{
    public:
    TofFrankAvg(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_frank_avg_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo );
};

class TofSet : public TofMethod{
    public:
    TofSet(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_set_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo );
};


std::pair<double, double> getTpcR(const dd4hep::Detector& detector);
dd4hep::rec::Vector3D calculateMomentum(EVENT::ReconstructedParticle* pfo, int location, double bField);
double calculateTrackLength(EVENT::ReconstructedParticle* pfo, int location);
double calculateTrackLengthIntegral(EVENT::ReconstructedParticle* pfo);

double calculateTrackLengthSET(EVENT::ReconstructedParticle* pfo, int location);
double calculateTrackLengthIntegralSET(EVENT::ReconstructedParticle* pfo);

int findShowerStart(EVENT::ReconstructedParticle* pfo);

#endif
