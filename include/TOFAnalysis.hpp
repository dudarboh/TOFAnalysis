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
#include "Math/Vector3D.h"
#include "UTIL/ILDConf.h"

using marlin::Processor;
using namespace std::chrono;
using std::vector, std::string, std::unique_ptr;
using namespace ROOT::Math;

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
    int _nSETHits, _nTPCHits;
    vector <XYZVector> _xyzSETHit;
    vector <double> _tSETHit;
    // ECAL hits
    int _nECALHits;
    vector <XYZVector> _xyzECALHit;
    vector <double> _tECALHit;
    vector <int> _layerECALHit;
    vector <double> _eECALHit;
    // Track
    double _chi2Track;
    int _ndfTrack;
    double _dEdXTrack;
    double _lengthTrackIP;
    double _lengthTrackCalo;
    double _lengthTrackIntegral;
    // Track States
    XYZVector _pTrackAtIP;
    XYZVector _pTrackAtCalo;
    XYZVector _xyzTrackAtCalo;
    double _d0TrackAtIP;
    double _z0TrackAtIP;
    double _d0TrackAtCalo;
    double _z0TrackAtCalo;
    // Cluster
    XYZVector _xyzCluster;
    //MCParticle from PFO
    double _weightMC;
    XYZVector _xyzVtxMC;
    double _tVtxMC;
    XYZVector _pMC;
    int _PDG;

    double _rTPCOuter;
    double _bField[3];
    // lcio//DD4HEP specifics
    const int _TPCID = ILDDetID::TPC;
    const int _TPCindex = _TPCID*2-2;
};


#endif
