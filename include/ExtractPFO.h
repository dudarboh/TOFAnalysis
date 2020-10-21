/**
    @file ExtractPFO.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractPFO class for extracting data from slcio into root file
*/

#ifndef ExtractPFO_h
#define ExtractPFO_h 1

#include "TFile.h"
#include "TTree.h"
#include "marlin/Processor.h"
#include <chrono>
#include <vector>
#include <string>
using marlin::Processor;
using namespace std::chrono;
using std::vector, std::string;

class ExtractPFO : public Processor {
public:
    ExtractPFO();
    ~ExtractPFO();

    Processor* newProcessor() {return new ExtractPFO;}
    void init();
    void processEvent(LCEvent* evt);
    void end();

protected:
    //Time and status progress
	int _nEvt;
    time_point <system_clock> _start;

    string _outputFileName;

    TFile* _file;
    TTree* _tree;

    float _charge;
    double _p;
    double _pt;
    int _nTracks;

    int _nMC;

    vector<float> _weightMC;
    vector<float> _energyMC;
    vector<int> _PDG;
    vector<double> _xMC;
    vector<double> _yMC;
    vector<double> _zMC;
    vector<float> _tMC;
    vector<double> _pMC;
    vector<double> _ptMC;
    vector<float> _massMC;
    vector<float> _chargeMC;

    vector<int> _isCreatedInSimulation;
    vector<int> _isBackscatter;
    vector<int> _vertexIsNotEndpointOfParent;
    vector<int> _isDecayedInTracker;
    vector<int> _isDecayedInCalorimeter;
    vector<int> _hasLeftDetector;
    vector<int> _isStopped;
    vector<int> _isOverlay;
};


#endif
