/**
    @file TOFAnalysis.h
    @author Bohdan Dudar
    @date October 2020
    @brief TOFAnalysis class for extracting data from slcio into root file
*/

#ifndef TOFAnalysis_hpp
#define TOFAnalysis_hpp 1

#include "marlin/Processor.h"
#include "TFile.h"
#include "TTree.h"
#include <chrono>
#include <string>
#include <vector>
#include "DDRec/Vector3D.h"

using marlin::Processor;
using namespace std::chrono;
using std::vector, std::string, std::unique_ptr;
using dd4hep::rec::Vector3D;

class TOFAnalysis : public Processor {
public:
    TOFAnalysis();
    Processor* newProcessor() {return new TOFAnalysis;}
    void init();
    void processEvent(LCEvent* evt);
    void end();

private:
    unique_ptr <TFile> _file;
    unique_ptr <TTree> _tree;
	int _nEvt;
    time_point <system_clock> _start;

    string _outputFileName;

    // SET hits
    int _nSETHits;
    vector <Vector3D> _xyzSETHit;
    vector <float> _tSETHit;
    // ECAL hits
    int _nECALHits;
    vector <Vector3D> _xyzECALHit;
    vector <float> _tECALHit;
    vector <int> _layerECALHit;
    vector <float> _eECALHit;
    // Track
    float _chi2Track;
    int _ndfTrack;
    float _dEdXTrack;
    float _lengthTrackIP;
    float _lengthTrackCalo;
    float _lengthTrackIntegral;
    // Track States
    Vector3D _pTrackAtIP;
    Vector3D _pTrackAtCalo;
    Vector3D _xyzTrackAtCalo;
    float _d0TrackAtCalo;
    float _z0TrackAtCalo;
    // Cluster
    Vector3D _xyzCluster;
    //MCParticle from PFO
    float _weightMC;
    Vector3D _xyzVtxMC;
    float _tVtxMC;
    Vector3D _pMC;
    int _PDG;
    int _isBackscatter;
    int _isDecayedInTracker;

    float _rTPCOuter;
    double _bField[3];
};


#endif
