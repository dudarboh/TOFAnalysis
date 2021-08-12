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
#include "EVENT/SimCalorimeterHit.h"
#include "UTIL/LCRelationNavigator.h"

using marlin::Processor;
using std::string, std::cout, std::endl, std::unique_ptr, std::vector, std::pair;
using ROOT::Math::XYZVectorF;
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
        TrackerHit* getSetHit(Track* track, double tpcROuter);
        double getTrackLength(Track* track, int from=TrackState::AtIP, int to=TrackState::AtCalorimeter);
        double getTrackLengthIntegral(Track* track);

        CalorimeterHit* getClosestHit( Cluster* cluster, XYZVectorF posTrackAtCalo);
        pair<XYZVectorF, double> getFastestHit( Cluster* cluster, double smearing=0.);
        double getTofFrankFit( Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing=0., unsigned int nLayers=10 );
        double getTofFrankAvg( Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing=0., unsigned int nLayers=10 );

        int getNEcalHits(Cluster* cluster);

        unique_ptr<TFile> _file;
        unique_ptr<TTree> _tree;
        string _outputFileName;

        int _nEvent;
        int _pdg;
        XYZVectorF _tsLastPos{};
        XYZVectorF _tsLastMom{};
        double _tsLastOmega;
        double _tsLastTanL;
        double _tsLastPhi;
        double _tsLastD0;
        double _tsLastZ0;

        XYZVectorF _tsCaloPos{};
        XYZVectorF _tsCaloMom{};
        double _tsCaloOmega;
        double _tsCaloTanL;
        double _tsCaloPhi;
        double _tsCaloD0;
        double _tsCaloZ0;

        bool _hasSetHit;
        int _nEcalHits;
        XYZVector _posSetHit{};

        XYZVectorF _posClosest{};
        XYZVectorF _posFastest[5];

        double _trackLengthSet;
        double _trackLengthCalo;
        double _trackLengthIntegral;

        double _smearings[5] = {0., 10., 30., 50., 100.};
        double _tofSetFront[5];
        double _tofSetBack[5];
        double _tofClosest[5];
        double _tofFastest[5];
        double _tofFrankFit[5];
        double _tofFrankAvg[5];

        double _bField;
        double _tpcROuter;
};


#endif
