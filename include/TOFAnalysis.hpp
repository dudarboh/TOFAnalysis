#ifndef TOFAnalysis_h
#define TOFAnalysis_h 1

#include <memory>
#include <string>
#include <vector>
#include <numeric>

#include "marlin/Processor.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/Vector3D.h"
#include "TH1F.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/SimCalorimeterHit.h"
#include "UTIL/LCRelationNavigator.h"
#include <DD4hep/Detector.h>


namespace MarlinTrk{
    class IMarlinTrkSystem;
}

using marlin::Processor;
using std::string, std::cout, std::endl, std::unique_ptr, std::vector, std::pair;
using ROOT::Math::XYZVectorF;
using ROOT::Math::XYZVector;

class TOFAnalysis : public Processor {
    public:
        TOFAnalysis();
        Processor* newProcessor() {return new TOFAnalysis;}
        void init();
        void processEvent(LCEvent* evt);
        void end();

        MCParticle* getMcMaxWeight(LCRelationNavigator pfoToMc, ReconstructedParticle* pfo);
        pair<double, double> getTpcR();
        TrackerHit* getSetHit(Track* track, double tpcROuter);

        CalorimeterHit* getClosestHit( Cluster* cluster, XYZVectorF posTrackAtCalo);
        pair<XYZVectorF, double> getFastestHit( Cluster* cluster, double smearing=0.);
        double getTofFrankFit( Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing=0., unsigned int nLayers=10 );
        double getTofFrankAvg( Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing=0., unsigned int nLayers=10 );

        int getNEcalHits(Cluster* cluster);

        void drawPfo(Track* track, Cluster* cluster);


        unique_ptr<TFile> _file;
        unique_ptr<TTree> _tree;
        string _outputFileName;

        int _nEvent;
        int _pdg;

        double dummy = std::numeric_limits<double>::min();
        XYZVectorF dummyVec = XYZVectorF(dummy, dummy, dummy);
        std::map< std::string, XYZVectorF > _tsPos = { {"ip", dummyVec}, {"first", dummyVec}, {"last", dummyVec}, {"ecal", dummyVec} };
        std::map< std::string, XYZVectorF > _tsMom = { {"ip", dummyVec}, {"first", dummyVec}, {"last", dummyVec}, {"ecal", dummyVec} };
        std::map< std::string, double > _tsOmega = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };
        std::map< std::string, double > _tsTanL = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };
        std::map< std::string, double > _tsPhi = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };
        std::map< std::string, double > _tsD0 = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };
        std::map< std::string, double > _tsZ0 = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };

        bool _hasSetHit;
        int _nEcalHits;
        int _nFitHits;
        XYZVector _posSetHit{};

        XYZVectorF _posClosest{};
        XYZVectorF _posFastest[5];

        std::map< std::string, double > _mom = { {"hmSet", 0.}, {"hmEcal", 0.} };

        std::map< std::string, double > _trackLength = { {"set", dummy}, {"ecal", dummy}  };
        std::map< std::string, double > _phiCurl = { {"set", dummy}, {"ecal", dummy}  };

        double _smearings[5] = {0., 10., 30., 50., 100.};
        double _tofSetFront[5];
        double _tofSetBack[5];
        double _tofClosest[5];
        double _tofFastest[5];
        double _tofFrankFit[5];
        double _tofFrankAvg[5];

        // MarlinTrk v02-00 release notes
        // USERS SHOULD NO LONGER DELETE THE IMarlinTrkSystem POINTER IN THEIR CODE (Marlin processor)
        MarlinTrk::IMarlinTrkSystem* _trkSystem;
        double _bField;
        double _tpcROuter;
        dd4hep::Detector& _theDetector = dd4hep::Detector::getInstance();

};


#endif
