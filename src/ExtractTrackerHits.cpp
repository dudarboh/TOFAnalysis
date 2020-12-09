
#include "ExtractTrackerHits.hpp"


#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include <EVENT/SimTrackerHit.h>
#include <UTIL/LCRelationNavigator.h>
using dd4hep::Detector, dd4hep::DetElement, dd4hep::mm, dd4hep::rec::FixedPadSizeTPCData;
using EVENT::LCCollection, EVENT::ReconstructedParticle, EVENT::SimTrackerHit;
using std::cout, std::endl;


ExtractTrackerHits aExtractTrackerHits;

ExtractTrackerHits::ExtractTrackerHits() : Processor("ExtractTrackerHits"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("TrackerHits.root"));
}

void ExtractTrackerHits::init(){
    _nEvt = 0;
    _start = system_clock::now();

    _file.reset( new TFile(_outputFileName.c_str(), "RECREATE") );
    _tree.reset( new TTree("TrackerHits", "Tree with tracker hits") );

    _tree->Branch("nInnerHits", &_nHits[0]);
    _tree->Branch("posInner", &_pos[0]);

    _tree->Branch("nTPCHits", &_nHits[1]);
    _tree->Branch("posTPC", &_pos[1]);

    _tree->Branch("nTPCMC", &_nMC[1]);
    _tree->Branch("posTPCMC", &_posMC[1]);
    _tree->Branch("pTPCMC", &_pMC[1]);
    _tree->Branch("pathLengthTPCMC", &_pathLengthMC[1]);
    _tree->Branch("isProducedBySecondaryTPCMC", &_isProducedBySecondary[1]);

    _tree->Branch("nSETHits", &_nHits[2]);
    _tree->Branch("posSET", &_pos[2]);

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

    LCCollection* colRelationTPC = evt->getCollection("TPCTrackerHitRelations");
    LCRelationNavigator relationTPC(colRelationTPC);
    LCCollection* colRelationSET = evt->getCollection("SETTrackerHitRelations");
    LCRelationNavigator relationSET(colRelationSET);

    for (int i=0; i<colPFO->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        // Only simple cases of PFOs
        if( nClusters != 1 || nTracks > 1) continue;

        if (nTracks == 0){
            for (int j = 0; j < _nTrackerRegions; ++j) _nHits[j] = 0;
            _tree->Fill();
            continue;
        }

        const Track* track = pfo->getTracks()[0];

        for (auto&& hit : track->getTrackerHits() ){
            const double* xyz = hit->getPosition();
            XYZTVector pos(xyz[0], xyz[1], xyz[2], hit->getTime());
            double rho = pos.rho();
            // 0 - Inner, 1 - TPC, 2 - SET
            const int trackerIdx = (_rTPCInner <= rho && rho <= _rTPCOuter) + 2*(rho > _rTPCOuter);
            _pos[trackerIdx].push_back(pos);

            if (trackerIdx != 1) continue;

            vector <LCObject*> relationObjects = relationTPC.getRelatedToObjects(hit);
            _nMC[trackerIdx].push_back(relationObjects.size());
            if (relationObjects.size() == 0) continue;
            // Push only 1st hit info. Managing 2d array in root files is pain
            SimTrackerHit* mcHit = dynamic_cast<SimTrackerHit*>(relationObjects[0]);
            const double* xyzMC = mcHit->getPosition();
            _posMC[trackerIdx].push_back( XYZTVector(xyzMC[0], xyzMC[1], xyzMC[2], mcHit->getTime()) );
            const float* momMC = mcHit->getMomentum();
            _pMC[trackerIdx].push_back( PxPyPzEVector(momMC[0], momMC[1], momMC[2], mcHit->getEDep()) );
            _pathLengthMC[trackerIdx].push_back(mcHit->getPathLength());
            _isProducedBySecondary[trackerIdx].push_back(mcHit->isProducedBySecondary());
        }
        for (int j = 0; j < _nTrackerRegions; ++j) _nHits[j] = _pos[j].size();
        _tree->Fill();

        for (int j = 0; j < _nTrackerRegions; ++j){
            _pos[j].clear();
            _nMC[j].clear();
            _posMC[j].clear();
            _pMC[j].clear();
            _pathLengthMC[j].clear();
            _isProducedBySecondary[j].clear();
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
