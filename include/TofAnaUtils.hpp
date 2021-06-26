#ifndef TofAnaUtils_h
#define TofAnaUtils_h 1

#include "EVENT/ReconstructedParticle.h"
#include <utility>
#include "TString.h"
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"

struct TofMethod{
    TofMethod(double smearing){
        _smearing = smearing;
    }

    virtual double calculate( EVENT::ReconstructedParticle* pfo ) = 0;

    TString _branch_name;
    double _tof, _smearing;
};

struct TofClosest : public TofMethod{
    TofClosest(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_closest_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo );
};

struct TofFastest : public TofMethod{
    TofFastest(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_fastest_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo );
};

struct TofFrankFit : public TofMethod{
    TofFrankFit(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_frank_fit_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo );
};

struct TofFrankAvg : public TofMethod{
    TofFrankAvg(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_frank_avg_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo );
};

struct TofSet : public TofMethod{
    TofSet(double smearing) : TofMethod(smearing) {
        _branch_name = Form( "tof_set_%d", int(_smearing) );
    };
    double calculate( EVENT::ReconstructedParticle* pfo, double rTpcOuter );
};


std::pair<double, double> getTpcR(const dd4hep::Detector& detector);
dd4hep::rec::Vector3D calculateMomentum(EVENT::ReconstructedParticle* pfo, int location, double bField);
double calculateTrackLength(EVENT::ReconstructedParticle* pfo, int location);
double calculateTrackLengthIntegral(EVENT::ReconstructedParticle* pfo);
void findShowerStart(EVENT::ReconstructedParticle* pfo);

#endif
