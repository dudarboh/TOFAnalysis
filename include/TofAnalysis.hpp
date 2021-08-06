/**
    @file TOFAnalysis.h
    @author Bohdan Dudar
    @brief TOFAnalysis class for extracting data from slcio into root file
*/

#ifndef TofAnalysis_h
#define TofAnalysis_h 1

#include <string>
#include <vector>
#include <memory>
#include <utility>

#include "marlin/Processor.h"
// #include "Math/Vector3D.h"
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"
#include "EVENT/ReconstructedParticle.h"

#include "TFile.h"
#include "TTree.h"
// #include "TString.h"




class TofAnalysis : public marlin::Processor {
    public:
        TofAnalysis();
        Processor* newProcessor() {return new TofAnalysis;}
        void init();
        void processEvent(LCEvent* evt);
        void end();

        //Utility functions
        void writeTpcHits(EVENT::ReconstructedParticle* pfo);
        void writeEcalHits(EVENT::ReconstructedParticle* pfo);
        void writeSetHits(EVENT::ReconstructedParticle* pfo);
        double calcTofClosest(EVENT::ReconstructedParticle* pfo, double smearing);
        double calcTofFastest(EVENT::ReconstructedParticle* pfo, double smearing);
        double calcTofFrankFit(EVENT::ReconstructedParticle* pfo, double smearing);
        double calcTofFrankAvg(EVENT::ReconstructedParticle* pfo, double smearing);
        double calcTofSet(EVENT::ReconstructedParticle* pfo, double smearing);

        std::pair<double, double> getTpcR(const dd4hep::Detector& detector);
        dd4hep::rec::Vector3D calcMomentum(EVENT::ReconstructedParticle* pfo, int location, double bField);
        double calcTrackLength(EVENT::ReconstructedParticle* pfo, int location);
        double calcTrackLengthIntegral(EVENT::ReconstructedParticle* pfo);

        double calcTrackLengthSET(EVENT::ReconstructedParticle* pfo, int location);
        double calcTrackLengthIntegralSET(EVENT::ReconstructedParticle* pfo);

        int findShowerStart(EVENT::ReconstructedParticle* pfo);


    private:
        //Steering variables///
        std::string _outputFileName;
        bool _writeTpcHits;
        bool _writeEcalHits;
        bool _writeSetHits;

        int _nEvt;
        std::pair <double, double> _tpcR;
        double _bField;

        std::unique_ptr <TFile> _file;
        std::unique_ptr <TTree> _tree;

        /////ROOT branches/////
        int _pdg;

        int _nTPCHits;
        std::vector <dd4hep::rec::Vector3D> _posTPCHit;
        std::vector <float> _tTPCHit;

        int _nECALHits;
        std::vector <dd4hep::rec::Vector3D> _posECALHit;
        std::vector <float> _tECALHit;
        std::vector <int> _layerECALHit;
        std::vector <float> _eECALHit;

        int _nSETHits;
        std::vector <dd4hep::rec::Vector3D> _posSETHit;
        std::vector <float> _tSETHit;

        dd4hep::rec::Vector3D _momIP;
        dd4hep::rec::Vector3D _momECAL;



        double _trackLengthIP;
        double _trackLengthECAL;
        double _trackLengthIntegral;

        double _trackLengthIPSET;
        double _trackLengthECALSET;
        double _trackLengthIntegralSET;

        double _tofClosest;
        double _tofFastest;
        double _tofFrankFit;
        double _tofFrankAvg;
        double _tofSet;

};


#endif
