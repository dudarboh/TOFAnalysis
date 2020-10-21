#include "ExtractPFO.h"

// #include <iostream>
using std::cout, std::endl;

#include "EVENT/LCCollection.h"
using EVENT::LCCollection;
#include "EVENT/ReconstructedParticle.h"
using EVENT::ReconstructedParticle;

#include <EVENT/MCParticle.h>
using EVENT::MCParticle;

#include <UTIL/LCRelationNavigator.h>

ExtractPFO aExtractPFO;

ExtractPFO::ExtractPFO() : Processor("ExtractPFO"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("PFO.root"));
}

ExtractPFO::~ExtractPFO(){
    delete _tree;
    delete _file;
}

void ExtractPFO::init(){
    _nEvt = 0;
    _start = system_clock::now();

    _file = new TFile(_outputFileName.c_str(), "RECREATE");

    _tree = new TTree("PFO", "PFO parameters");

    _tree->Branch("charge", &_charge, "charge/F");
    _tree->Branch("p", &_p, "p/D");
    _tree->Branch("pt", &_pt, "pt/D");
    _tree->Branch("nTracks", &_nTracks, "nTracks/I");
    _tree->Branch("nMC", &_nMC, "nMC/I");
    _tree->Branch("weightMC", &_weightMC);
    _tree->Branch("energyMC", &_energyMC);
    _tree->Branch("PDG", &_PDG);
    _tree->Branch("xMC", &_xMC);
    _tree->Branch("yMC", &_yMC);
    _tree->Branch("zMC", &_zMC);
    _tree->Branch("tMC", &_tMC);
    _tree->Branch("pMC", &_pMC);
    _tree->Branch("ptMC", &_ptMC);
    _tree->Branch("massMC", &_massMC);
    _tree->Branch("chargeMC", &_chargeMC);

    _tree->Branch("isCreatedInSimulation", &_isCreatedInSimulation);
    _tree->Branch("isBackscatter", &_isBackscatter);
    _tree->Branch("vertexIsNotEndpointOfParent", &_vertexIsNotEndpointOfParent);
    _tree->Branch("isDecayedInTracker", &_isDecayedInTracker);
    _tree->Branch("isDecayedInCalorimeter", &_isDecayedInCalorimeter);
    _tree->Branch("hasLeftDetector", &_hasLeftDetector);
    _tree->Branch("isStopped", &_isStopped);
    _tree->Branch("isOverlay", &_isOverlay);
}

void ExtractPFO::processEvent(LCEvent* evt){
    ++_nEvt;
    if(_nEvt%10 == 0){
        double elapsedTime = duration<double>(system_clock::now() - _start).count();
        cout << "Event: "<<_nEvt<<"   Elapsed Time: "<<elapsedTime<<" sec     Avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    }

    LCCollection* colPFO = evt->getCollection("PandoraPFOs");

    LCCollection* colRelation = evt->getCollection("RecoMCTruthLink");
    LCRelationNavigator relation(colRelation);

    for (int i=0; i<colPFO->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
        int nClusters = pfo->getClusters().size();
        _nTracks = pfo->getTracks().size();

        // Only simple cases of PFOs
        if( nClusters != 1 || _nTracks > 1) continue;

        //Fill branch variables
        _charge = pfo->getCharge();
        const double* mom = pfo->getMomentum();
        _p = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
        _pt = sqrt(mom[0]*mom[0] + mom[1]*mom[1]);

        const vector <LCObject*>& relationObjects = relation.getRelatedToObjects(pfo);
        const vector <float>& relationWeights = relation.getRelatedToWeights(pfo);
        _nMC = relationObjects.size();
        for(int j=0; j<_nMC; ++j){
            MCParticle* mcPFO = dynamic_cast<MCParticle*> (relationObjects[j]);
            _weightMC.push_back(relationWeights[j]);
            _energyMC.push_back(mcPFO->getEnergy());
            _PDG.push_back(mcPFO->getPDG());
            const double* pos = mcPFO->getVertex();
            _xMC.push_back(pos[0]);
            _yMC.push_back(pos[1]);
            _zMC.push_back(pos[2]);
            _tMC.push_back(mcPFO->getTime());
            const double* momMC = mcPFO->getMomentum();
            _pMC.push_back(sqrt(momMC[0]*momMC[0] + momMC[1]*momMC[1] + momMC[2]*momMC[2]));
            _ptMC.push_back(sqrt(momMC[0]*momMC[0] + momMC[1]*momMC[1]));
            _massMC.push_back(mcPFO->getMass());
            _chargeMC.push_back(mcPFO->getCharge());

            _isCreatedInSimulation.push_back(mcPFO->isCreatedInSimulation());
            _isBackscatter.push_back(mcPFO->isBackscatter());
            _vertexIsNotEndpointOfParent.push_back(mcPFO->vertexIsNotEndpointOfParent());
            _isDecayedInTracker.push_back(mcPFO->isDecayedInTracker());
            _isDecayedInCalorimeter.push_back(mcPFO->isDecayedInCalorimeter());
            _hasLeftDetector.push_back(mcPFO->hasLeftDetector());
            _isStopped.push_back(mcPFO->isStopped());
            _isOverlay.push_back(mcPFO->isOverlay());
        }
        _tree->Fill();

        _weightMC.clear();
        _energyMC.clear();
        _PDG.clear();
        _xMC.clear();
        _yMC.clear();
        _zMC.clear();
        _tMC.clear();
        _pMC.clear();
        _ptMC.clear();
        _massMC.clear();
        _chargeMC.clear();

        _isCreatedInSimulation.clear();
        _isBackscatter.clear();
        _vertexIsNotEndpointOfParent.clear();
        _isDecayedInTracker.clear();
        _isDecayedInCalorimeter.clear();
        _hasLeftDetector.clear();
        _isStopped.clear();
        _isOverlay.clear();
    } // end of PFOs loop
}

void ExtractPFO::end(){
    double elapsedTime = duration<double>(system_clock::now() - _start).count();
    cout<<"Finished writing Track parameters"<<endl;
    cout<<"Total number of events: "<<_nEvt<<endl;
    cout<<"Total elapsed time: "<<elapsedTime<<" sec"<<endl;
    cout<<"Total avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    _file->Write();
    cout<<_file->GetName()<<"   file is written in the current directory"<<endl;
}
