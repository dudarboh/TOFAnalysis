
#include "ExtractCalorimeterHits.hpp"

#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include <EVENT/SimCalorimeterHit.h>
#include <UTIL/LCRelationNavigator.h>
#include "marlinutil/CalorimeterHitType.h"
using EVENT::LCCollection, EVENT::ReconstructedParticle, EVENT::SimCalorimeterHit;
using std::cout, std::endl;


ExtractCalorimeterHits aExtractCalorimeterHits;

ExtractCalorimeterHits::ExtractCalorimeterHits() : Processor("ExtractCalorimeterHits"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("CalorimeterHits.root"));
}


void ExtractCalorimeterHits::init(){
    _nEvt = 0;
    _start = system_clock::now();

    _file.reset( new TFile(_outputFileName.c_str(), "RECREATE") );
    _tree.reset( new TTree("ECALHits", "Tree with ECAL hits") );

    _tree->Branch("nHits", &_nHits);
    _tree->Branch("posECALHit", &_pos);
    _tree->Branch("layer", &_layer);
    _tree->Branch("energy", &_energy);
}

void ExtractCalorimeterHits::processEvent(LCEvent* evt){
    ++_nEvt;
    LCCollection* colPFO = evt->getCollection("PandoraPFOs");
    LCCollection* colRelationBarrel = evt->getCollection("EcalBarrelRelationsSimRec");
    LCRelationNavigator relationBarrel(colRelationBarrel);

    for (int i=0; i<colPFO->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        // Only simple cases of PFOs
        if( nClusters != 1 || nTracks > 1) continue;

        const Cluster* cluster = pfo->getClusters()[0];
        //Count only ECAl
        for (auto&& hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isEcal = (hitType.caloID() == CHT::ecal);
            if (!isEcal) continue;
            const float* pos = hit->getPosition();
            _pos.push_back( XYZTVector(pos[0], pos[1], pos[2], hit->getTime()) );
            _layer.push_back(hitType.layer());
            _energy.push_back(hit->getEnergy());
        }
        _nHits = _pos.size();
        _tree->Fill();

        _pos.clear();
        _layer.clear();
        _energy.clear();
    } // end of PFOs loop
}

void ExtractCalorimeterHits::end(){
    double elapsedTime = duration<double>(system_clock::now() - _start).count();
    cout<<"Finished writing Calorimeter Hits"<<endl;
    cout<<"Total number of events: "<<_nEvt<<endl;
    cout<<"Total elapsed time: "<<elapsedTime<<" sec"<<endl;
    cout<<"Total avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    _file->Write();
    cout<<_file->GetName()<<"   file is written in the current directory"<<endl;
}
