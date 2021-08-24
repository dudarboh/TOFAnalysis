#ifndef SETAnalysis_h
#define SETAnalysis_h 1

#include <memory>
#include <string>
#include <vector>
#include <numeric>

#include "marlin/Processor.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"
#include "TH1F.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/SimCalorimeterHit.h"
#include "UTIL/LCRelationNavigator.h"

namespace MarlinTrk{
    class IMarlinTrkSystem;
}

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

        double estimateTrackLengthTanL(double phi1, double phi2, double omega, double tanL);
        double estimateTrackLengthZ(double phi1, double phi2, double omega, double z1, double z2);

        double getTrackLength(Track* track, bool extrapolateToIp=true, bool extrapolateToEcal=true, std::string method="dz");
        double getMom2Harmonic(Track* track, bool extrapolateToIp, bool extrapolateToEcal, std::string method);

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
        XYZVector _posSetHit{};

        XYZVectorF _posClosest{};
        XYZVectorF _posFastest[5];

        std::map< std::string, double > _mom = { {"sqrHmTanL", 0.}, {"sqrHmTanLSet", 0.}, {"sqrHmDz", 0.}, {"sqrHmDzSet", 0.} };

        std::map< std::string, double > _trackLength = { {"ip", dummy}, {"set", dummy}, {"calo", dummy}, {"refitTanL", dummy}, {"refitZ", dummy}, {"setRefitZ", dummy}  };

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
};


#endif
