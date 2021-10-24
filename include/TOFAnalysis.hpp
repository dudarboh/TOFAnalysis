#ifndef TOFAnalysis_h
#define TOFAnalysis_h 1

#include <memory>
#include <map>

#include "marlin/Processor.h"
#include "TFile.h"
#include "TTree.h"

#include "DD4hep/Detector.h"
#include "Math/Vector3D.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

class TOFAnalysis : public marlin::Processor {
    public:
        TOFAnalysis();
        marlin::Processor* newProcessor() {return new TOFAnalysis;}
        void init();
        void processEvent(LCEvent* evt);
        void end();

        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;
        std::string _outputFileName;

        int _nEvent;
        int _pdg;

        double dummy = std::numeric_limits<double>::min();
        ROOT::Math::XYZVectorF dummyVec = ROOT::Math::XYZVectorF(dummy, dummy, dummy);
        std::map< std::string, ROOT::Math::XYZVectorF > _tsPos = { {"ip", dummyVec}, {"first", dummyVec}, {"last", dummyVec}, {"ecal", dummyVec} };
        std::map< std::string, ROOT::Math::XYZVectorF > _tsMom = { {"ip", dummyVec}, {"first", dummyVec}, {"last", dummyVec}, {"ecal", dummyVec} };
        std::map< std::string, double > _tsOmega = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };
        std::map< std::string, double > _tsTanL = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };
        std::map< std::string, double > _tsPhi = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };
        std::map< std::string, double > _tsD0 = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };
        std::map< std::string, double > _tsZ0 = { {"ip", dummy}, {"first", dummy}, {"last", dummy}, {"ecal", dummy} };

        bool _hasSetHit;
        int _nEcalHits;
        int _nFitHits;
        ROOT::Math::XYZVector _posSetHit{};

        ROOT::Math::XYZVectorF _posClosest{};
        ROOT::Math::XYZVectorF _posFastest[5];

        std::map< std::string, double > _mom = { {"hmSet", 0.}, {"hmEcal", 0.} };

        std::map< std::string, double > _trackLength = { {"set", dummy}, {"ecal", dummy}  };
        std::map< std::string, double > _nCurls = { {"set", dummy}, {"ecal", dummy}  };

        double _smearings[5] = {0., 10., 30., 50., 100.};
        double _tofSetFront[5];
        double _tofSetBack[5];
        double _tofClosest[5];
        double _tofFastest[5];
        double _tofFrankFit[5];
        double _tofFrankAvg[5];

        // MarlinTrk v02-00 release notes:
        // USERS SHOULD NO LONGER DELETE THE IMarlinTrkSystem POINTER IN THEIR CODE (Marlin processor)
        MarlinTrk::IMarlinTrkSystem* _trkSystem;
        dd4hep::Detector& _theDetector = dd4hep::Detector::getInstance();
        double _bField;
        double _tpcROuter;
};


#endif
