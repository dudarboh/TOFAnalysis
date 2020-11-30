/**
    @file ExtractCalorimeterHits.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractCalorimeterHits class for extracting tracker hits from slcio to root file
*/

#ifndef ExtractCalorimeterHits_h
#define ExtractCalorimeterHits_h 1

#include "TFile.h"
#include "TTree.h"
#include "marlin/Processor.h"
#include <chrono>
#include <vector>
#include <string>
using marlin::Processor;
using namespace std::chrono;
using std::vector, std::string;

class ExtractCalorimeterHits : public Processor {
public:
    ExtractCalorimeterHits();
    ~ExtractCalorimeterHits();

    Processor* newProcessor() {return new ExtractCalorimeterHits;}
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

    // static const int _nCalorimeterRegions = 1;
    int _nHits;
    vector <float> _x;
    vector <float> _y;
    vector <float> _z;
    vector <float> _t;
    vector <float> _tMCFastest;
    vector <int> _layer;

};


#endif
