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
    int _nPFOs;
    int _nGoodPFOs;
    vector <float> _p;
    vector <float> _charge;

    vector <float> _d0;
    vector <float> _phi;
    vector <float> _omega;
    vector <float> _z0;
    vector <float> _tanL;
    vector <float> _chi2;
    vector <int> _ndf;
    vector <float> _dEdX;
    vector <float> _length;

    vector <float> _d0Last;
    vector <float> _phiLast;
    vector <float> _omegaLast;
    vector <float> _z0Last;
    vector <float> _tanLLast;
    vector <float> _xRefLast;
    vector <float> _yRefLast;
    vector <float> _zRefLast;

    vector <float> _d0Calo;
    vector <float> _phiCalo;
    vector <float> _omegaCalo;
    vector <float> _z0Calo;
    vector <float> _tanLCalo;
    vector <float> _xRefCalo;
    vector <float> _yRefCalo;
    vector <float> _zRefCalo;

    vector <int> _nHitsTrack;
    vector <int> _nHitsTPC;
    vector < vector<double> > _xHit;
    vector < vector<double> > _yHit;
    vector < vector<double> > _zHit;
    vector < vector<float> > _tHit;

    vector <float> _xCluster;
    vector <float> _yCluster;
    vector <float> _zCluster;
    vector <float> _phiCluster;
    vector <float> _thetaCluster;

    // Only hits with some time info and from ECAL
    vector <int> _nHitsCluster;
    vector < vector<float> > _xHitCluster;
    vector < vector<float> > _yHitCluster;
    vector < vector<float> > _zHitCluster;
    vector < vector<float> > _tHitCluster;
    vector < vector<int> > _layerHitCluster;
    vector < vector<float> > _dToLineHitCluster;
    vector < vector<float> > _dToRefPointHitCluster;


    // lcio//DD4HEP specifics
    const int _TPCindex;

    TFile* _file;
    TTree* _tree;
};


#endif
