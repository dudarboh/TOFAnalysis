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
        double getTrackLengthIntegral(Track* track);

        CalorimeterHit* getFastestHit( Cluster* cluster);
        CalorimeterHit* getClosestHit( Cluster* cluster, XYZVector posTrackAtCalo);
        double getTofFrankFit( Cluster* cluster, XYZVector posTrackAtCalo, XYZVector momTrackAtCalo, double smearing=0., unsigned int nLayers=10 );
        double getTofFrankAvg( Cluster* cluster, XYZVector posTrackAtCalo, XYZVector momTrackAtCalo, double smearing=0., unsigned int nLayers=10 );

        XYZVector getFastestContPos(SimCalorimeterHit* hit);
        int getNEcalHits(Cluster* cluster);

        unique_ptr<TFile> _file;
        unique_ptr<TTree> _tree;
        string _outputFileName;

        int _nEvent;
        int _pdg;
        XYZVector _tsLastPos{};
        XYZVector _tsLastMom{};
        double _tsLastOmega;
        double _tsLastTanL;
        double _tsLastPhi;
        double _tsLastD0;
        double _tsLastZ0;

        XYZVector _tsCaloPos{};
        XYZVector _tsCaloMom{};
        double _tsCaloOmega;
        double _tsCaloTanL;
        double _tsCaloPhi;
        double _tsCaloD0;
        double _tsCaloZ0;

        int _nSetHits;
        XYZVector _setHitPos{};
        double _setHitTime;
        XYZVector _setPosTrue{};

        int _nEcalHits;
        XYZVector _posClosest{};
        XYZVector _posClosestSim{};
        double _tofClosest;
        double _tofClosestSim;
        XYZVector _posFastest{};
        XYZVector _posFastestSim{};
        double _tofFastest;
        double _tofFastestSim;
        double _tofFrankFit;
        double _tofFrankAvg;

        double _trackLengthSet;
        double _trackLengthCalo;
        double _trackLengthIntegral;

        double _bField;
        double _tpcROuter;
        double _smearing;
};


#endif
