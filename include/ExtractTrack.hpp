/**
    @file ExtractTrack.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractTrack class for extracting data from slcio into root file
*/

#ifndef ExtractTrack_hpp
#define ExtractTrack_hpp 1

#include "marlin/Processor.h"
#include "EVENT/TrackState.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include <chrono>
#include <vector>
#include <string>
using marlin::Processor;
using namespace ROOT::Math;
using namespace std::chrono;
using std::vector, std::string, std::unique_ptr;


class ExtractTrack : public Processor {
public:
    ExtractTrack();

    Processor* newProcessor() {return new ExtractTrack;}
    void init();
    void processEvent(LCEvent* evt);
    void end();

protected:
    unique_ptr <TFile> _file;
    unique_ptr <TTree> _tree;
	int _nEvt;
    time_point <system_clock> _start;

    double _bField;
    string _outputFileName;
    string _massAssumption;

    float _chi2;
    int _ndf;
    float _dEdX;
    static const int _nTrackStates = 4;
    int _trackStates[_nTrackStates] = {TrackState::AtIP, TrackState::AtFirstHit, TrackState::AtLastHit, TrackState::AtCalorimeter};
    XYZVector _p[_nTrackStates];
    XYZPoint _ref[_nTrackStates];
    float _d0[_nTrackStates];
    float _phi[_nTrackStates];
    float _omega[_nTrackStates];
    float _z0[_nTrackStates];
    float _tanL[_nTrackStates];
};


#endif
