
#include "ExtractTrackerHits.h"

// #include <iostream>
using std::cout, std::endl;

//TPC radii
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
using dd4hep::Detector, dd4hep::DetElement, dd4hep::mm;

#include "DDRec/DetectorData.h"
using dd4hep::rec::FixedPadSizeTPCData;

#include "EVENT/LCCollection.h"
using EVENT::LCCollection;
#include "EVENT/ReconstructedParticle.h"
using EVENT::ReconstructedParticle;

#include <EVENT/SimTrackerHit.h>
using EVENT::SimTrackerHit;

ExtractTrackerHits aExtractTrackerHits;

ExtractTrackerHits::ExtractTrackerHits() : Processor("ExtractTrackerHits"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("TrackerHits.root"));
}

ExtractTrackerHits::~ExtractTrackerHits(){
    delete _tree;
    delete _file;
}

void ExtractTrackerHits::init(){
    _nEvt = 0;
    _start = system_clock::now();

    _file = new TFile(_outputFileName.c_str(), "RECREATE");
    _tree = new TTree("TrackerHits", "Tree with tracker hits");

    _tree->Branch("nInnerHits", &_nHits[0], "nInnerHits/I");
    _tree->Branch("xInner", &_x[0]);
    _tree->Branch("yInner", &_y[0]);
    _tree->Branch("zInner", &_z[0]);
    _tree->Branch("tInner", &_t[0]);

    _tree->Branch("nTPCHits", &_nHits[1], "nTPCHits/I");
    _tree->Branch("xTPC", &_x[1]);
    _tree->Branch("yTPC", &_y[1]);
    _tree->Branch("zTPC", &_z[1]);

    _tree->Branch("nSETHits", &_nHits[2], "nSETHits/I");
    _tree->Branch("xSET", &_x[2]);
    _tree->Branch("ySET", &_y[2]);
    _tree->Branch("zSET", &_z[2]);
    _tree->Branch("tSET", &_t[2]);

    //get VXD, SIT, TPC, SET radii to write hits based on their subdetector
    const Detector& detector = Detector::getInstance();
    const DetElement tpcDet = detector.detector("TPC");
    const FixedPadSizeTPCData* tpc = tpcDet.extension <FixedPadSizeTPCData>() ;
    _rTPCInner = tpc->rMinReadout/mm;
    _rTPCOuter = tpc->rMaxReadout/mm;
}

void ExtractTrackerHits::processEvent(LCEvent* evt){
    ++_nEvt;
    if(_nEvt%10 == 0){
        double elapsedTime = duration<double>(system_clock::now() - _start).count();
        cout << "Event: "<<_nEvt<<"   Elapsed Time: "<<elapsedTime<<" sec     Avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    }

    LCCollection* colPFO = evt->getCollection("PandoraPFOs");

    for (int i=0; i<colPFO->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        // Only simple cases of PFOs
        if( nClusters != 1 || nTracks > 1) continue;

        if (nTracks == 0){
            for (int j = 0; j < _nTrackerRegions; ++j){
                _nHits[j] = 0;
            }
            _tree->Fill();
            continue;
        }

        const Track* track = pfo->getTracks()[0];

        for (auto&& hit : track->getTrackerHits() ){
            // const vector<LCObject*>& objVec = hit->getRawHits();
            // cout<<objVec.size()<<endl;

            const double* pos = hit->getPosition();
            const double rho = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
            // 0 - Inner, 1 - TPC, 2 - SET
            const int trackerIdx = (_rTPCInner <= rho && rho <= _rTPCOuter) + 2*(rho > _rTPCOuter);
            _x[trackerIdx].push_back(pos[0]);
            _y[trackerIdx].push_back(pos[1]);
            _z[trackerIdx].push_back(pos[2]);
            _t[trackerIdx].push_back(hit->getTime());
        }
        for (int j = 0; j < _nTrackerRegions; ++j) _nHits[j] = _x[j].size();
        _tree->Fill();

        //Clear all the vectors before next PFO
        for (int j = 0; j < _nTrackerRegions; ++j){
            _x[j].clear();
            _y[j].clear();
            _z[j].clear();
            _t[j].clear();
        }
    } // end of PFOs loop
}

void ExtractTrackerHits::end(){
    double elapsedTime = duration<double>(system_clock::now() - _start).count();
    cout<<"Finished writing Tracker Hits"<<endl;
    cout<<"Total number of events: "<<_nEvt<<endl;
    cout<<"Total elapsed time: "<<elapsedTime<<" sec"<<endl;
    cout<<"Total avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    _file->Write();
    cout<<_file->GetName()<<"   file is written in the current directory"<<endl;
}
