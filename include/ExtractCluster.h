/**
    @file ExtractCluster.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractCluster class for extracting tracker hits from slcio to root file
*/

#ifndef ExtractCluster_h
#define ExtractCluster_h 1

#include "TFile.h"
#include "TTree.h"
#include "marlin/Processor.h"
#include <chrono>
#include <vector>
#include <string>
using marlin::Processor;
using namespace std::chrono;
using std::vector, std::string;

class ExtractCluster : public Processor {
public:
    ExtractCluster();
    ~ExtractCluster();

    Processor* newProcessor() {return new ExtractCluster;}
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

    float _energy;
    float _x;
    float _y;
    float _z;
    // These phi and theta define direction of the shower development.
    // Not its position!
    float _phi;
    float _theta;
};


#endif
