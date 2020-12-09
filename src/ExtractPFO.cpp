#include "ExtractPFO.hpp"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "UTIL/LCRelationNavigator.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DD4hepUnits.h"
using EVENT::LCCollection, EVENT::ReconstructedParticle, EVENT::MCParticle;
using dd4hep::rec::LayeredCalorimeterData;
using dd4hep::Detector, dd4hep::DetType, dd4hep::DetElement, dd4hep::mm;
using std::cout, std::endl, std::stringstream, std::runtime_error;


//This function is only to check rInner of ECAL barrel
LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag=0) {
    Detector& mainDetector = Detector::getInstance();
    const vector<DetElement>& theDetectors = dd4hep::DetectorSelector(mainDetector).detectors(  includeFlag, excludeFlag );

    if( theDetectors.size()  != 1 ){
        stringstream es ;
        es << " getExtension: selection is not unique (or empty)  includeFlag: " << DetType( includeFlag ) << " excludeFlag: " << DetType( excludeFlag )
        << " --- found detectors : " ;
        for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ) es << theDetectors.at(i).name() << ", " ;
        throw runtime_error( es.str() ) ;
    }
    return theDetectors.at(0).extension<LayeredCalorimeterData>();
}

ExtractPFO aExtractPFO;

ExtractPFO::ExtractPFO() : Processor("ExtractPFO"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("PFO.root"));
}

void ExtractPFO::init(){
    //This is only to check rInner of ECAL barrel
    const LayeredCalorimeterData* eCalBarrelExtension = getExtension( (DetType::CALORIMETER|DetType::ELECTROMAGNETIC|DetType::BARREL), (DetType::AUXILIARY|DetType::FORWARD) );
    const double rInner = eCalBarrelExtension->extent[0]/dd4hep::mm;
    cout<<"Inner radius: "<<rInner<<endl;

    _nEvt = 0;
    _start = system_clock::now();

    _file.reset(new TFile(_outputFileName.c_str(), "RECREATE"));
    _tree.reset(new TTree("PFO", "Pandora PFO properties"));

    _tree->Branch("pReco", &_pReco);
    _tree->Branch("charge", &_charge);
    _tree->Branch("nTracks", &_nTracks);

    _tree->Branch("nMC", &_nMC);
    _tree->Branch("weightMC", &_weightMC);
    _tree->Branch("PDG", &_PDG);
    _tree->Branch("chargeMC", &_chargeMC);
    _tree->Branch("mMC", &_mMC);
    _tree->Branch("vtxMC", &_vtxMC);
    _tree->Branch("pMC", &_pMC);
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
        ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> ( colPFO->getElementAt(i) );
        int nClusters = pfo->getClusters().size();
        _nTracks = pfo->getTracks().size();

        // Only simple cases of PFOs
        if( nClusters != 1 || _nTracks > 1) continue;

        const double* mom = pfo->getMomentum();
        _pReco = PxPyPzEVector(mom[0], mom[1], mom[2], pfo->getEnergy());
        _charge = pfo->getCharge();

        const vector <LCObject*>& relationObjects = relation.getRelatedToObjects(pfo);
        const vector <float>& relationWeights = relation.getRelatedToWeights(pfo);

        _nMC = relationObjects.size();
        for(int j=0; j<_nMC; ++j){
            _weightMC.push_back( relationWeights[j] );
            MCParticle* mcPFO = dynamic_cast <MCParticle*> ( relationObjects[j] );
            const double* vtx = mcPFO->getVertex();
            const double* momMC = mcPFO->getMomentum();

            _vtxMC.push_back( XYZTVector(vtx[0], vtx[1], vtx[2], mcPFO->getTime()) );
            _pMC.push_back( PxPyPzEVector(momMC[0], momMC[1], momMC[2], mcPFO->getEnergy()) );
            _PDG.push_back( mcPFO->getPDG() );
            _mMC.push_back(mcPFO->getMass());
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
        _vtxMC.clear();
        _pMC.clear();
        _PDG.clear();
        _mMC.clear();
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
