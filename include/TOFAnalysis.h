/**
    @file TOFAnalysis.h
    @author Bohdan Dudar
    @date September 2020
    @brief TOFAnalysis class for extracting data from slcio into root file for further python script TOFAnalysis
*/

#ifndef TOFAnalysis_h
#define TOFAnalysis_h 1

#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"
#include <EVENT/ReconstructedParticle.h>
#include "EVENT/CalorimeterHit.h"
#include "marlinutil/CalorimeterHitType.h"

#include "UTIL/ILDConf.h"
#include <UTIL/PIDHandler.h>
#include <chrono>
#include <iostream>
#include "lcio.h"
#include <TFile.h>
#include <TTree.h>
#include "DDRec/Vector3D.h"

#include "TInterpreter.h"

using namespace dd4hep::rec;
using namespace std;
using namespace lcio;
using namespace marlin;

class TOFAnalysis : public Processor {
public:
    /**
    @brief Constructor of the TOFAnalysis processor

    Register input LCCollection type and name from slcio file using parent class method
    <a href="https://ilcsoft.desy.de/Marlin/current/doc/html/classmarlin_1_1Processor.html#a60e489b83117049d67ada1d0b2558199">registerInputCollection()</a>.
    Write a name into @p _PFOs variable.

    Constructor starts twice. Caused by a parent class behaviour.

    */
    TOFAnalysis();
    ~TOFAnalysis();

    Processor* newProcessor() {return new TOFAnalysis;}
    void init();
    void processEvent(LCEvent* evt);
    void end();

    /**
    @brief Very brief documentation example

    */
    Double_t getAlgorithmParameter(ReconstructedParticle* pfo, int algorithmID, int parameterID, PIDHandler& handler);

protected:
    //Time and status progress
	int _nEvt;
    chrono::time_point<chrono::system_clock> _start;

    string _outputFileName;

    //root file data
    double _p;
    float _charge;

    float _d0;
    float _phi;
    float _omega;
    float _z0;
    float _tanL;
    float _chi2;
    int _ndf;
    float _dEdX;

    float _d0FirstState;
    float _phiFirstState;
    float _omegaFirstState;
    float _z0FirstState;
    float _tanLFirstState;
    float _xRefFirstState;
    float _yRefFirstState;
    float _zRefFirstState;

    float _d0LastState;
    float _phiLastState;
    float _omegaLastState;
    float _z0LastState;
    float _tanLLastState;
    float _xRefLastState;
    float _yRefLastState;
    float _zRefLastState;

    float _d0CalState;
    float _phiCalState;
    float _omegaCalState;
    float _z0CalState;
    float _tanLCalState;
    float _xRefCalState;
    float _yRefCalState;
    float _zRefCalState;

    int _nTrHits;
    int _nTPCHits;
    vector <double> _xTrHit;
    vector <double> _yTrHit;
    vector <double> _zTrHit;
    vector <float> _tTrHit;

    // Only hits with some time info and from ECAL
    int _nCalHits;
    vector <float> _xCalHit;
    vector <float> _yCalHit;
    vector <float> _zCalHit;
    vector <float> _tCalHit;
    vector <int> _layerCalHit;
    vector <double> _dToLineCalHit;
    vector <double> _dToRefPointCalHit;

    // lcio//DD4HEP specifics
    const int _TPCindex;

    TFile* _file;
    TTree* _tree;
};


#endif
