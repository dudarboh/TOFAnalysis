/**
    @file TOFAnalysis.h
    @author Bohdan Dudar
    @date June 2020
    @brief TOFAnalysis class for extracting data from slcio into root file
*/

#ifndef TofAnalysis_h
#define TofAnalysis_h 1

#include <string>
#include <vector>
#include <memory>
#include <string>
#include <utility>

#include "marlin/Processor.h"
#include "DDRec/Vector3D.h"
#include "EVENT/ReconstructedParticle.h"

#include "TFile.h"
#include "TTree.h"

using marlin::Processor;


class TofAnalysis : public Processor {
    public:
        TofAnalysis();
        Processor* newProcessor() {return new TofAnalysis;}
        void init();
        void processEvent(LCEvent* evt);
        void end();

        void writeTPCHits(EVENT::ReconstructedParticle* pfo);
        void writeECALHits(EVENT::ReconstructedParticle* pfo);
        void writeSETHits(EVENT::ReconstructedParticle* pfo);
    private:
        ///////////////////////
        //Steering variables///
        ///////////////////////
        std::string _outputFileName;
        bool _writeTPCHits;
        bool _writeECALHits;
        bool _writeSETHits;

        int _nEvt;
        std::pair <double, double> _tpcR;
        double _bField;

        std::unique_ptr <TFile> _file;
        std::unique_ptr <TTree> _tree;

        ///////////////////////
        /////ROOT branches/////
        ///////////////////////
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

        // vector<tof_method> _tofs = {TofClosest(50.),
        //                             tof_method("fastest", 0),
        //                             tof_method("frankFit", 0),
        //                             tof_method("frankAvg", 0)};

};


#endif
