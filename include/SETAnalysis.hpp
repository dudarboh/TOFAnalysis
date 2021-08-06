#ifndef SETAnalysis_h
#define SETAnalysis_h 1

#include <memory>
#include <string>
#include <vector>

#include "marlin/Processor.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"

#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "UTIL/LCRelationNavigator.h"

using marlin::Processor;
using std::string, std::cout, std::endl, std::unique_ptr, std::vector, std::pair;
using ROOT::Math::XYZVector;

class SETAnalysis : public Processor {
    public:
        SETAnalysis();
        Processor* newProcessor() {return new SETAnalysis;}
        void init();
        void processEvent(LCEvent* evt);
        void end();

        MCParticle* getMcMaxWeight(LCRelationNavigator pfoToMc, ReconstructedParticle* pfo);
        pair<double, double> getTpcR();
        vector<TrackerHit*> getSetHits(Track* track, double tpcROuter);
        double getTrackLength(Track* track, int from=TrackState::AtIP, int to=TrackState::AtCalorimeter);
        double getTofClosest( Cluster* cluster, XYZVector posTrackAtCalo, double smearing=0. );
        double getTofFastest( Cluster* cluster, XYZVector posTrackAtCalo, double smearing=0. );
        double getTofFrankFit( Cluster* cluster, XYZVector posTrackAtCalo, XYZVector momTrackAtCalo, double smearing=0., unsigned int nLayers=10 );
        double getTofFrankAvg( Cluster* cluster, XYZVector posTrackAtCalo, XYZVector momTrackAtCalo, double smearing=0., unsigned int nLayers=10 );


        unique_ptr<TFile> _file;
        unique_ptr<TTree> _tree;
        string _outputFileName;

        XYZVector _tsLastPos;
        XYZVector _tsLastMom;
        double _tsLastOmega;
        double _tsLastTanL;
        double _tsLastPhi;
        double _tsLastD0;
        double _tsLastZ0;

        XYZVector _tsCaloPos;
        XYZVector _tsCaloMom;
        double _tsCaloOmega;
        double _tsCaloTanL;
        double _tsCaloPhi;
        double _tsCaloD0;
        double _tsCaloZ0;

        int _nSetHits;
        XYZVector _setHitPos;
        double _setHitTime;
        XYZVector _setPosTrue;

        XYZVector _caloPosTrue;
        double _caloTofTrue;
        double _caloTofClosest;
        double _caloTofFastest;
        double _caloTofFrankFit;
        double _caloTofFrankAvg;

        double _trackLengthSet;
        double _trackLengthCalo;

        double _bField;
        double _tpcROuter;
};


#endif
