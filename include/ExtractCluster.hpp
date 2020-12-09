/**
    @file ExtractCluster.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractCluster class for extracting tracker hits from slcio to root file
*/

#ifndef ExtractCluster_hpp
#define ExtractCluster_hpp 1

#include "TFile.h"
#include "TTree.h"
#include "Math/Point3D.h"
#include "marlin/Processor.h"
#include <chrono>
#include <vector>
#include <string>
using marlin::Processor;
using namespace ROOT::Math;
using namespace std::chrono;
using std::vector, std::string, std::unique_ptr;

class ExtractCluster : public Processor {
public:
    ExtractCluster();

    Processor* newProcessor() {return new ExtractCluster;}
    void init();
    void processEvent(LCEvent* evt);
    void end();

protected:
    //Time and status progress
	int _nEvt;
    time_point <system_clock> _start;

    string _outputFileName;
    unique_ptr<TFile> _file;
    unique_ptr<TTree> _tree;

    float _energy;
    XYZPoint _pos;
    // These phi and theta define direction of the shower development.
    // Not its position!
    float _phi;
    float _theta;
};


#endif
