#ifndef TOFAnalysis_h
#define TOFAnalysis_h 1

#include <string>
#include <vector>
#include "marlin/Processor.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"

class TOFAnalysis : public marlin::Processor {
    public:

        TOFAnalysis(const TOFAnalysis&) = delete;
        TOFAnalysis& operator=(const TOFAnalysis&) = delete;
        marlin::Processor* newProcessor() { return new TOFAnalysis; }

        TOFAnalysis();
        void init();
        void processEvent(EVENT::LCEvent* evt);
        void end();

    private:
        int _nEvent;

        std::unique_ptr<TFile> _file;
        std::unique_ptr<TTree> _tree;

        // Variables for tree branches
        int _pdg;
        bool _hasSetHit;
        int _nEcalHits;
        ROOT::Math::XYZVectorF _posECAL;
        ROOT::Math::XYZVectorF _momECAL;
        double _d0ECAL;
        double _z0ECAL;

        double _trackLengthSET;
        double _momentumSET;
        double _trackLengthECAL;
        double _momentumECAL;
        std::vector<double> _tof_closest = std::vector<double>(11, 0.);
        std::vector<double> _tof_set = std::vector<double>(11, 0.);
        std::vector< std::vector<double> > _tof_avg = std::vector< std::vector<double> >(11, std::vector<double>(30, 0.));
        std::vector< std::vector<double> > _tof_fit = std::vector< std::vector<double> >(11, std::vector<double>(30, 0.));

        MarlinTrk::IMarlinTrkSystem* _trkSystem = nullptr;

        double _bField{};
        double _tpcOuterR{};
};

#endif
