#include "TofAnalysis.hpp"
#include "TofAnaUtils.hpp"

#include <marlin/Global.h>
#include "DD4hep/Detector.h"
#include "EVENT/LCCollection.h"
#include <marlinutil/GeometryUtil.h>
#include <marlin/ProcessorEventSeeder.h>
#include "CLHEP/Random/Randomize.h"
#include "marlinutil/CalorimeterHitType.h"

// #include "LCIO/EVENT/MCParticle.h"
// #include "LCIO/EVENT/TrackerHit.h"
// #include "UTIL/LCRelationNavigator.h"
// #include "DDRec/DetectorData.h"
// #include "DD4hep/DetType.h"
// #include "DD4hep/DetectorSelector.h"
// #include "DD4hep/DD4hepUnits.h"
// #include "HelixClass.h"
//
//
// EVENT::MCParticle, EVENT::TrackerHit;
// using dd4hep::rec::LayeredCalorimeterData;
// using dd4hep::Detector, dd4hep::DetType, dd4hep::DetElement, dd4hep::rec::FixedPadSizeTPCData;

using std::string;
using std::vector;
using dd4hep::Detector;
using dd4hep::rec::Vector3D;
using EVENT::LCCollection;
using EVENT::ReconstructedParticle;


TofAnalysis aTofAnalysis;

TofAnalysis::TofAnalysis() : Processor("TofAnalysis"){
    registerProcessorParameter(string("outputFile"),
                               string("Name of the output root file"),
                               _outputFileName,
                               string("TofAnalysis_RENAME.root"));

    registerProcessorParameter(string("writeTPCHits"),
                              string("Write TPC hits to the output root file"),
                              _writeTPCHits,
                              bool(false));

    registerProcessorParameter(string("writeECALHits"),
                              string("Write ECAL hits to the output root file"),
                              _writeECALHits,
                              bool(false));

    registerProcessorParameter(string("writeSETHits"),
                              string("Write SET hits to the output root file"),
                              _writeSETHits,
                              bool(false));

}

void TofAnalysis::init(){
    printParameters();

    _nEvt = 0;
    const Detector& detector = Detector::getInstance();
    _tpcR = getTpcR(detector);
    _bField = MarlinUtil::getBzAtOrigin();

    _file.reset( new TFile(_outputFileName.c_str(), "RECREATE") );
    _tree.reset( new TTree("TOFAnalysis", "TOFAnalysis") );

    _tree->Branch("pdg", &_pdg);
    if(_writeTPCHits){
        _tree->Branch("nTPCHits", &_nTPCHits);
        _tree->Branch("posTPCHit", &_posTPCHit);
        _tree->Branch("tTPCHit", &_tTPCHit);
    }
    if(_writeECALHits){
        _tree->Branch("nECALHits", &_nECALHits);
        _tree->Branch("posECALHit", &_posECALHit);
        _tree->Branch("tECALHit", &_tECALHit);
        _tree->Branch("layerECALHit", &_layerECALHit);
        _tree->Branch("eECALHit", &_eECALHit);
    }
    if(_writeSETHits){
        _tree->Branch("nSETHits", &_nSETHits);
        _tree->Branch("posSETHit", &_posSETHit);
        _tree->Branch("tSETHit", &_tSETHit);
    }


    marlin::Global::EVENTSEEDER->registerProcessor(this);
    return;
}

void TofAnalysis::processEvent(LCEvent* evt){
    CLHEP::HepRandom::setTheSeed( marlin::Global::EVENTSEEDER->getSeed(this) );
    ++_nEvt;

    LCCollection* colPFO = evt->getCollection("PandoraPFOs");
    for (int i=0; i<colPFO->getNumberOfElements(); ++i){

        ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> ( colPFO->getElementAt(i) );
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();
        // Only simple cases of PFOs
        if( nClusters != 1 || nTracks > 1) continue;

        if(_writeTPCHits) writeTPCHits(pfo);
        if(_writeECALHits) writeECALHits(pfo);
        if(_writeSETHits) writeSETHits(pfo);

        // for(auto& m:_tofs) m.var = calculateTOF(pfo, m.name, m.smearing);

        _tree->Fill();
    }
}

void TofAnalysis::end(){
    _file->Write();
    std::cout<<_file->GetName()<<"   file is written in the current directory"<<std::endl;
    std::cout<<"Processed number of events: "<<_nEvt<<std::endl;
    // std::cout<<"Total elapsed time: "<<elapsedTime<<" sec"<<std::endl;
    // std::cout<<"Total avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<std::endl;
}


void TofAnalysis::writeTPCHits(ReconstructedParticle* pfo){
    _posTPCHit.clear();
    _tTPCHit.clear();

    if (pfo->getTracks().size() == 0){
        _nTPCHits = 0;
        return;
    }

    const Track* track = pfo->getTracks()[0];
    for ( const auto& hit : track->getTrackerHits() ){
        Vector3D pos( hit->getPosition() );
        bool isTPCHit = _tpcR.first <= pos.rho() && pos.rho() <= _tpcR.second;

        if (isTPCHit){
            float time = hit->getTime();

            _posTPCHit.push_back(pos);
            _tTPCHit.push_back(time);
        }
    }
    _nTPCHits = _posTPCHit.size();
    return;
}


void TofAnalysis::writeECALHits(ReconstructedParticle* pfo){
    _posECALHit.clear();
    _tECALHit.clear();
    _layerECALHit.clear();
    _eECALHit.clear();

    const Cluster* cluster = pfo->getClusters()[0];
    for ( const auto& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );

        if (isECALHit){
            Vector3D pos( hit->getPosition() );
            float time = hit->getTime();
            int layer = hitType.layer();
            float energy = hit->getEnergy();

            _posECALHit.push_back(pos);
            _tECALHit.push_back(time);
            _layerECALHit.push_back(layer);
            _eECALHit.push_back(energy);
        }
    }
    _nECALHits = _posECALHit.size();
    return;
}


void TofAnalysis::writeSETHits(ReconstructedParticle* pfo){
    _posSETHit.clear();
    _tSETHit.clear();

    if (pfo->getTracks().size() == 0){
        _nSETHits = 0;
        return;
    }

    const Track* track = pfo->getTracks()[0];
    for ( const auto& hit : track->getTrackerHits() ){
        Vector3D pos( hit->getPosition() );
        bool isSETHit = pos.rho() > _tpcR.second;

        if (isSETHit){
            float time = hit->getTime();

            _posSETHit.push_back(pos);
            _tSETHit.push_back(time);
        }
    }
    _nSETHits = _posSETHit.size();
    return;
}
