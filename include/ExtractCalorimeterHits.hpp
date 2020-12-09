/**
    @file ExtractCalorimeterHits.h
    @author Bohdan Dudar
    @date October 2020
    @brief ExtractCalorimeterHits class for extracting tracker hits from slcio to root file
*/

#ifndef ExtractCalorimeterHits_hpp
#define ExtractCalorimeterHits_hpp 1

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

class ExtractCalorimeterHits : public Processor {
public:
    ExtractCalorimeterHits();

    Processor* newProcessor() {return new ExtractCalorimeterHits;}
    void init();
    void processEvent(LCEvent* evt);
    void end();

protected:
    unique_ptr <TFile> _file;
    unique_ptr <TTree> _tree;
	int _nEvt;
    time_point <system_clock> _start;

    string _outputFileName;

    int _nHits;
    vector <XYZTVector> _pos;
    vector <int> _layer;
    vector <float> _energy;
};


#endif
