/**
    @file ExtractTrackerHits.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractTrackerHits class for extracting tracker hits from slcio to root file
*/

#ifndef ExtractTrackerHits_h
#define ExtractTrackerHits_h 1

#include "TFile.h"
#include "TTree.h"
#include "marlin/Processor.h"
#include <chrono>
#include <vector>
#include <string>
using marlin::Processor;
using namespace std::chrono;
using std::vector, std::string;

class ExtractTrackerHits : public Processor {
public:
    ExtractTrackerHits();
    ~ExtractTrackerHits();

    Processor* newProcessor() {return new ExtractTrackerHits;}
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

    double _rTPCInner;
    double _rTPCOuter;

    static const int _nTrackerRegions = 3;
    int _nHits[_nTrackerRegions];
    vector <double> _x[_nTrackerRegions];
    vector <double> _y[_nTrackerRegions];
    vector <double> _z[_nTrackerRegions];
    vector <float> _t[_nTrackerRegions];
};


#endif
