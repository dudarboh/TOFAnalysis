/**
    @file ExtractPFO.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractPFO class for extracting data from slcio into root file
*/

#ifndef ExtractPFO_hpp
#define ExtractPFO_hpp 1

#include "marlin/Processor.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <chrono>
#include <string>
#include <vector>
using marlin::Processor;
using namespace ROOT::Math;
using namespace std::chrono;
using std::vector, std::string, std::unique_ptr;


class ExtractPFO : public Processor {
public:
    ExtractPFO();

    Processor* newProcessor() {return new ExtractPFO;}
    void init();
    void processEvent(LCEvent* evt);
    void end();

private:
    unique_ptr <TFile> _file;
    unique_ptr <TTree> _tree;
	int _nEvt;
    time_point <system_clock> _start;

    string _outputFileName;

    PxPyPzEVector _pReco;
    float _charge;
    int _nTracks;
    int _nMC;
    vector <float> _weightMC;
    vector <int> _PDG;
    vector <float> _chargeMC;
    vector <float> _mMC;
    vector <XYZTVector> _vtxMC;
    vector <PxPyPzEVector> _pMC;
    vector <int> _isCreatedInSimulation;
    vector <int> _isBackscatter;
    vector <int> _vertexIsNotEndpointOfParent;
    vector <int> _isDecayedInTracker;
    vector <int> _isDecayedInCalorimeter;
    vector <int> _hasLeftDetector;
    vector <int> _isStopped;
    vector <int> _isOverlay;
};


#endif
