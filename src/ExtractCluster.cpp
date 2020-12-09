#include "ExtractCluster.hpp"

#include <iostream>
using std::cout, std::endl;

#include "EVENT/LCCollection.h"
using EVENT::LCCollection;
#include "EVENT/ReconstructedParticle.h"
using EVENT::ReconstructedParticle;

ExtractCluster aExtractCluster;

ExtractCluster::ExtractCluster() : Processor("ExtractCluster"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("Cluster.root"));
}

void ExtractCluster::init(){
    _nEvt = 0;
    _start = system_clock::now();

    _file.reset( new TFile(_outputFileName.c_str(), "RECREATE") );
    _tree.reset( new TTree("Cluster", "Tree with cluster info") );

    _tree->Branch("energy", &_energy);
    _tree->Branch("pos", &_pos);
    _tree->Branch("phi", &_phi);
    _tree->Branch("theta", &_theta);

}

void ExtractCluster::processEvent(LCEvent* evt){
    ++_nEvt;
    LCCollection* colPFO = evt->getCollection("PandoraPFOs");

    for (int i=0; i<colPFO->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        // Only simple cases of PFOs
        if( nClusters != 1 || nTracks > 1) continue;

        const Cluster* cluster = pfo->getClusters()[0];

        _energy = cluster->getEnergy();
        const float* pos = cluster->getPosition();
        _pos = XYZPoint(pos[0], pos[1], pos[2]);
        _phi = cluster->getIPhi();
        _theta = cluster->getITheta();
        _tree->Fill();
    } // end of PFOs loop
}

void ExtractCluster::end(){
    double elapsedTime = duration<double>(system_clock::now() - _start).count();
    cout<<"Finished writing Calorimeter Hits"<<endl;
    cout<<"Total number of events: "<<_nEvt<<endl;
    cout<<"Total elapsed time: "<<elapsedTime<<" sec"<<endl;
    cout<<"Total avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    _file->Write();
    cout<<_file->GetName()<<"   file is written in the current directory"<<endl;
}
