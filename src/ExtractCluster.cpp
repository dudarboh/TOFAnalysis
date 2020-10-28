#include "ExtractCluster.h"

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

ExtractCluster::~ExtractCluster(){
    delete _tree;
    delete _file;
}

void ExtractCluster::init(){
    _nEvt = 0;
    _start = system_clock::now();

    _file = new TFile(_outputFileName.c_str(), "RECREATE");
    _tree = new TTree("Cluster", "Tree with cluster info");

    _tree->Branch("energy", &_energy, "energy/F");
    _tree->Branch("x", &_x, "x/F");
    _tree->Branch("y", &_y, "y/F");
    _tree->Branch("z", &_z, "z/F");
    _tree->Branch("phi", &_phi, "phi/F");
    _tree->Branch("theta", &_theta, "theta/F");

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
        _x = pos[0];
        _y = pos[1];
        _z = pos[2];
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
