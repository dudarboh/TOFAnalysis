
#include "ExtractCalorimeterHits.h"

#include <iostream>
using std::cout, std::endl;

#include "EVENT/LCCollection.h"
using EVENT::LCCollection;
#include "EVENT/ReconstructedParticle.h"
using EVENT::ReconstructedParticle;

// CHT class
#include "marlinutil/CalorimeterHitType.h"

#include "marlin/VerbosityLevels.h"

ExtractCalorimeterHits aExtractCalorimeterHits;

ExtractCalorimeterHits::ExtractCalorimeterHits() : Processor("ExtractCalorimeterHits"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("CalorimeterHits.root"));
}

ExtractCalorimeterHits::~ExtractCalorimeterHits(){
    delete _tree;
    delete _file;
}

void ExtractCalorimeterHits::init(){
    _nEvt = 0;
    _start = system_clock::now();

    _file = new TFile(_outputFileName.c_str(), "RECREATE");
    _tree = new TTree("ECALHits", "Tree with ECAL hits");

    _tree->Branch("nHits", &_nHits, "nHits/I");
    _tree->Branch("x", &_x);
    _tree->Branch("y", &_y);
    _tree->Branch("z", &_z);
    _tree->Branch("t", &_t);
    _tree->Branch("layer", &_layer);

}

void ExtractCalorimeterHits::processEvent(LCEvent* evt){
    //Pring status
    ++_nEvt;
    if(_nEvt%10 == 0){
        double elapsedTime = duration<double>(system_clock::now() - _start).count();
        cout<<"Event: "<<_nEvt<<"   Elapsed Time: "<<elapsedTime<<" sec     Avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    }

    //Get collection of PFOs for this event
    LCCollection* col = nullptr;
    try{
        col = evt->getCollection("PandoraPFOs");
    }
    catch (...){
        cout<<"Event "<<_nEvt<<" has no PandoarPFOs collection. Skip event"<<endl;
        return;
    }

    for (int p=0; p<col->getNumberOfElements(); ++p){
        const ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col->getElementAt(p));
        // Look only at PFOs with 1 cluster and 1 track
        if(pfo->getClusters().size() != 1 || pfo->getTracks().size() != 1) continue;

        const Cluster* cluster = pfo->getClusters()[0];
        //Count only ECAl
        for (auto&& hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isEcal = (hitType.caloID() == CHT::ecal);
            if (!isEcal) continue;

            const float* pos = hit->getPosition();
            _x.push_back(pos[0]);
            _y.push_back(pos[1]);
            _z.push_back(pos[2]);
            _t.push_back(hit->getTime());
            _layer.push_back(hitType.layer());
        }
        _nHits = _x.size();
        _tree->Fill();

        //Clear all the vectors before next PFO
        _x.clear();
        _y.clear();
        _z.clear();
        _t.clear();
        _layer.clear();
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