#include "ExtractPFO.h"

// #include <iostream>
using std::cout, std::endl, std::stringstream, std::runtime_error;

#include "EVENT/LCCollection.h"
using EVENT::LCCollection;
#include "EVENT/ReconstructedParticle.h"
using EVENT::ReconstructedParticle;

#include <EVENT/MCParticle.h>
using EVENT::MCParticle;

#include <UTIL/LCRelationNavigator.h>


//This part to get inner radius of ECAL
#include <DDRec/DetectorData.h>
using dd4hep::rec::LayeredCalorimeterData;

#include "DD4hep/Detector.h"
#include <DD4hep/DetType.h>
using dd4hep::Detector, dd4hep::DetType, dd4hep::DetElement;
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DD4hepUnits.h"
using dd4hep::mm;

//This function is only to check rInner of ECAL barrel
LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag=0) {
    LayeredCalorimeterData * theExtension = 0;
    Detector & mainDetector = Detector::getInstance();
    const vector<DetElement>& theDetectors = dd4hep::DetectorSelector(mainDetector).detectors(  includeFlag, excludeFlag );

    if( theDetectors.size()  != 1 ){
        stringstream es ;
        es << " getExtension: selection is not unique (or empty)  includeFlag: " << DetType( includeFlag ) << " excludeFlag: " << DetType( excludeFlag )
        << " --- found detectors : " ;
        for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ) es << theDetectors.at(i).name() << ", " ;
        throw runtime_error( es.str() ) ;
    }

    theExtension = theDetectors.at(0).extension<LayeredCalorimeterData>();
    return theExtension;
}

ExtractPFO aExtractPFO;

ExtractPFO::ExtractPFO() : Processor("ExtractPFO"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("PFO.root"));
}

ExtractPFO::~ExtractPFO(){
    delete _tree;
    delete _file;
}

void ExtractPFO::init(){
    //This is only to check rInner of ECAL barrel
    const LayeredCalorimeterData* eCalBarrelExtension = getExtension( (DetType::CALORIMETER|DetType::ELECTROMAGNETIC|DetType::BARREL), (DetType::AUXILIARY|DetType::FORWARD) );
    const double rInner = eCalBarrelExtension->extent[0]/dd4hep::mm;
    cout<<"Inner radius: "<<rInner<<endl;

    _nEvt = 0;
    _start = system_clock::now();

    _file = new TFile(_outputFileName.c_str(), "RECREATE");

    _tree = new TTree("PFO", "PFO parameters");

    _tree->Branch("charge", &_charge, "charge/F");
    _tree->Branch("px", &_px, "px/D");
    _tree->Branch("py", &_py, "py/D");
    _tree->Branch("pz", &_pz, "pz/D");
    _tree->Branch("nTracks", &_nTracks, "nTracks/I");
    _tree->Branch("nMC", &_nMC, "nMC/I");
    _tree->Branch("weightMC", &_weightMC);
    _tree->Branch("energyMC", &_energyMC);
    _tree->Branch("PDG", &_PDG);
    _tree->Branch("xMC", &_xMC);
    _tree->Branch("yMC", &_yMC);
    _tree->Branch("zMC", &_zMC);
    _tree->Branch("tMC", &_tMC);
    _tree->Branch("pxMC", &_pxMC);
    _tree->Branch("pyMC", &_pyMC);
    _tree->Branch("pzMC", &_pzMC);
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
        _px = mom[0];
        _py = mom[1];
        _pz = mom[2];

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
            _pxMC.push_back(momMC[0]);
            _pyMC.push_back(momMC[1]);
            _pzMC.push_back(momMC[2]);
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
        _pxMC.clear();
        _pyMC.clear();
        _pzMC.clear();
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
