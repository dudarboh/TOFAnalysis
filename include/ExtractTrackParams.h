/**
    @file ExtractTrackParams.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractTrackParams class for extracting data from slcio into root file
*/

#ifndef ExtractTrackParams_h
#define ExtractTrackParams_h 1

#include "TFile.h"
#include "TTree.h"
#include "marlin/Processor.h"
#include <chrono>
#include <vector>
#include <string>
using marlin::Processor;
using namespace std::chrono;
using std::vector, std::string;

#include <EVENT/TrackState.h>
// using EVENT::TrackState;
// #include "lcio.h"

class ExtractTrackParams : public Processor {
public:
    ExtractTrackParams();
    ~ExtractTrackParams();

    Processor* newProcessor() {return new ExtractTrackParams;}
    void init();
    void processEvent(LCEvent* evt);
    void end();

protected:
    //Time and status progress
	int _nEvt;
    time_point <system_clock> _start;

    string _outputFileName;
    string _massAssumption;

    TFile* _file;
    TTree* _tree;

    double _bField;

    float _chi2;
    int _ndf;
    float _dEdX;

    //TrackState parameters for each state
    static const int _nTrackStates = 4;
    int _trackStates[_nTrackStates] = {TrackState::AtIP, TrackState::AtFirstHit, TrackState::AtLastHit, TrackState::AtCalorimeter};

    float _p[_nTrackStates];
    float _pt[_nTrackStates];
    float _d0[_nTrackStates];
    float _phi[_nTrackStates];
    float _omega[_nTrackStates];
    float _z0[_nTrackStates];
    float _tanL[_nTrackStates];
    float _xRef[_nTrackStates];
    float _yRef[_nTrackStates];
    float _zRef[_nTrackStates];
};


#endif
