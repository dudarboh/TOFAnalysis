/**
    @file ExtractTrackerHits.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractTrackerHits class for extracting tracker hits from slcio to root file
*/

#ifndef ExtractTrackerHits_hpp
#define ExtractTrackerHits_hpp 1

#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "marlin/Processor.h"
#include <chrono>
#include <vector>
#include <string>
using marlin::Processor;
using namespace ROOT::Math;
using namespace std::chrono;
using std::vector, std::string, std::unique_ptr;

class ExtractTrackerHits : public Processor {
public:
    ExtractTrackerHits();

    Processor* newProcessor() {return new ExtractTrackerHits;}
    void init();
    void processEvent(LCEvent* evt);
    void end();

protected:
    unique_ptr <TFile> _file;
    unique_ptr <TTree> _tree;
	int _nEvt;
    time_point <system_clock> _start;

    string _outputFileName;

    double _rTPCInner;
    double _rTPCOuter;

    static const int _nTrackerRegions = 3;
    int _nHits[_nTrackerRegions];
    vector <XYZTVector> _pos[_nTrackerRegions];
    vector <int> _nMC[_nTrackerRegions];
    vector <XYZTVector> _posMC[_nTrackerRegions];
    vector <PxPyPzEVector> _pMC[_nTrackerRegions];
    vector <float> _pathLengthMC[_nTrackerRegions];
    vector <float> _isProducedBySecondary[_nTrackerRegions];
};


#endif
