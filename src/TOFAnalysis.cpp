#include "TOFAnalysis.hpp"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "EVENT/TrackerHit.h"
#include "UTIL/LCRelationNavigator.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DD4hepUnits.h"
#include "marlinutil/CalorimeterHitType.h"
#include "HelixClass.h"
using EVENT::LCCollection, EVENT::ReconstructedParticle, EVENT::MCParticle, EVENT::TrackerHit;
using dd4hep::rec::LayeredCalorimeterData;
using dd4hep::Detector, dd4hep::DetType, dd4hep::DetElement, dd4hep::rec::FixedPadSizeTPCData;
using std::cout, std::endl;

bool sortbyr(TrackerHit* a, TrackerHit* b) {
    const double* posA = a->getPosition();
    const double* posB = b->getPosition();
    return XYZVector(posA[0], posA[1], posA[2]).r() < XYZVector(posB[0], posB[1], posB[2]).r();}


//This function is only to check rInner of ECAL barrel
LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag=0) {
    Detector& mainDetector = Detector::getInstance();
    const vector<DetElement>& theDetectors = dd4hep::DetectorSelector(mainDetector).detectors(  includeFlag, excludeFlag );

    if( theDetectors.size()  != 1 ){
        cout << " getExtension: selection is not unique (or empty)  includeFlag: " << dd4hep::DetType( includeFlag ) << " excludeFlag: " << dd4hep::DetType( excludeFlag )
 		<< " --- found detectors : " ;
        for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ) cout << theDetectors.at(i).name() << ", " ;
        assert(0);
    }

    return theDetectors.at(0).extension<LayeredCalorimeterData>();
}

TOFAnalysis aTOFAnalysis;

TOFAnalysis::TOFAnalysis() : Processor("TOFAnalysis"){
    registerProcessorParameter(string("outputFile"),
                               string("Name of the output root file"),
                               _outputFileName,
                               string("TOFAnalysis_RENAME.root"));
}

void TOFAnalysis::init(){
    //Check outer TPC radius to collect only SET hits
    const Detector& detector = Detector::getInstance();
    const DetElement tpcDet = detector.detector("TPC");
    const FixedPadSizeTPCData* tpc = tpcDet.extension <FixedPadSizeTPCData>();
    _rTPCOuter = tpc->rMaxReadout/dd4hep::mm;
    detector.field().magneticField({0., 0., 0.}, _bField);

    //This is only to check rInner of ECAL barrel
    const LayeredCalorimeterData* eCalBarrelExtension = getExtension( (DetType::CALORIMETER|DetType::ELECTROMAGNETIC|DetType::BARREL), (DetType::AUXILIARY|DetType::FORWARD) );
    const double rInner = eCalBarrelExtension->extent[0]/dd4hep::mm;
    cout<<"Inner radius: "<<rInner<<endl;

    _nEvt = 0;
    _start = system_clock::now();

    _file.reset(new TFile(_outputFileName.c_str(), "RECREATE"));
    _tree.reset(new TTree("TOFAnalysis", "ROOT file for TOF analysis"));

    //SET
    _tree->Branch("nSETHits", &_nSETHits);
    _tree->Branch("xyzSETHit", &_xyzSETHit);
    _tree->Branch("tSETHit", &_tSETHit);
    //ECAL hits
    _tree->Branch("nECALHits", &_nECALHits);
    _tree->Branch("xyzECALHit", &_xyzECALHit);
    _tree->Branch("tECALHit", &_tECALHit);
    _tree->Branch("layerECALHit", &_layerECALHit);
    _tree->Branch("eECALHit", &_eECALHit);
    //Track
    _tree->Branch("chi2Track", &_chi2Track);
    _tree->Branch("ndfTrack", &_ndfTrack);
    _tree->Branch("dEdXTrack", &_dEdXTrack);
    _tree->Branch("lengthTrackIP", &_lengthTrackIP);
    _tree->Branch("lengthTrackCalo", &_lengthTrackCalo);
    _tree->Branch("lengthTrackIntegral", &_lengthTrackIntegral);
    //Track States
    _tree->Branch("pTrackAtIP", &_pTrackAtIP);
    _tree->Branch("pTrackAtCalo", &_pTrackAtCalo);
    _tree->Branch("xyzTrackAtCalo", &_xyzTrackAtCalo);
    _tree->Branch("d0TrackAtCalo", &_d0TrackAtCalo);
    _tree->Branch("z0TrackAtCalo", &_z0TrackAtCalo);
    //Cluster
    _tree->Branch("xyzCluster", &_xyzCluster);
    //MCParticle
    _tree->Branch("weightMC", &_weightMC);
    _tree->Branch("PDG", &_PDG);
    _tree->Branch("xyzVtxMC", &_xyzVtxMC);
    _tree->Branch("tVtxMC", &_tVtxMC);
    _tree->Branch("pMC", &_pMC);
}

void TOFAnalysis::processEvent(LCEvent* evt){
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
        int nTracks = pfo->getTracks().size();

        // Only simple cases of PFOs
        if( nClusters != 1 || nTracks > 1) continue;

        //SET
        if (nTracks == 0) _nSETHits = 0;
        else{
            const Track* track = pfo->getTracks()[0];

            for (const auto& hit:track->getTrackerHits()){
                XYZVector xyzHit = XYZVector( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
                if ( xyzHit.rho() > _rTPCOuter ){
                    _xyzSETHit.push_back(xyzHit);
                    _tSETHit.push_back(hit->getTime());
                }
            }
            _nSETHits = _xyzSETHit.size();
        }

        //Cluster
        const Cluster* cluster = pfo->getClusters()[0];
        _xyzCluster = XYZVector( cluster->getPosition()[0], cluster->getPosition()[1], cluster->getPosition()[2] );
        // ECAL hits
        for (const auto& hit:cluster->getCalorimeterHits()){
            //Count only ECAL hits
            CHT hitType( hit->getType() );
            bool isEcal = (hitType.caloID() == CHT::ecal);
            if (!isEcal) continue;
            XYZVector xyzHit = XYZVector( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
            _xyzECALHit.push_back(xyzHit);
            _tECALHit.push_back( hit->getTime() );
            _layerECALHit.push_back( hitType.layer() );
            _eECALHit.push_back( hit->getEnergy() );
        }
        _nECALHits = _xyzECALHit.size();
        // Track
        if (nTracks == 0){
            _chi2Track = 0.;
            _ndfTrack = 0;
            _dEdXTrack = 0.;
            _lengthTrackIP = 0.;
            _lengthTrackCalo = 0.;
            _lengthTrackIntegral = 0.;
            _pTrackAtIP = XYZVector();
            _pTrackAtCalo = XYZVector();
            _xyzTrackAtCalo = XYZVector();
            _d0TrackAtCalo = 0.;
            _z0TrackAtCalo = 0.;
        }
        else{
            const Track* track = pfo->getTracks()[0];

            _chi2Track = track->getChi2();
            _ndfTrack = track->getNdf();
            _dEdXTrack = track->getdEdx();

            const TrackState* tsIP = track->getTrackState(TrackState::AtIP);
            double phiIP = tsIP->getPhi();
            double omegaIP = tsIP->getOmega();
            double tanLIP = tsIP->getTanLambda();
            double d0IP = tsIP->getD0();
            double z0IP = tsIP->getZ0();
            const TrackState* tsCalo = track->getTrackState(TrackState::AtCalorimeter);
            double phiCalo = tsCalo->getPhi();
            double omegaCalo = tsCalo->getOmega();
            double tanLCalo = tsCalo->getTanLambda();
            _d0TrackAtCalo = tsCalo->getD0();
            _z0TrackAtCalo = tsCalo->getZ0();

            _lengthTrackIP = abs( (phiIP - phiCalo)/omegaIP )*sqrt(1. + tanLIP*tanLIP);
            _lengthTrackCalo = abs( (phiIP - phiCalo)/omegaCalo )*sqrt(1. + tanLCalo*tanLCalo);

            const TrackState* tsFirst = track->getTrackState(TrackState::AtFirstHit);
            double phiFirst = tsFirst->getPhi();
            const TrackState* tsLast = track->getTrackState(TrackState::AtLastHit);
            double phiLast = tsLast->getPhi();
            _lengthTrackIntegral = abs( (phiIP - phiFirst)/omegaIP )*sqrt(1. + tanLIP*tanLIP) + abs( (phiLast - phiCalo)/omegaCalo )*sqrt(1. + tanLCalo*tanLCalo);

            vector <TrackerHit*> trackHits = track->getTrackerHits();
            sort(trackHits.begin(), trackHits.end(), sortbyr);
            for (size_t j=1; j < trackHits.size(); ++j){
                _lengthTrackIntegral += (XYZVector(trackHits[j]->getPosition()[0], trackHits[j]->getPosition()[1], trackHits[j]->getPosition()[2]) - XYZVector(trackHits[j-1]->getPosition()[0], trackHits[j-1]->getPosition()[1], trackHits[j-1]->getPosition()[2])).r();
            }

            HelixClass helixIP;
            helixIP.Initialize_Canonical(phiIP, d0IP, z0IP, omegaIP, tanLIP, _bField[2]/dd4hep::tesla);
            _pTrackAtIP = XYZVector( helixIP.getMomentum()[0], helixIP.getMomentum()[1], helixIP.getMomentum()[2] );

            HelixClass helixCalo;
            helixCalo.Initialize_Canonical(phiCalo, _d0TrackAtCalo, _z0TrackAtCalo, omegaCalo, tanLCalo, _bField[2]/dd4hep::tesla);
            _pTrackAtCalo = XYZVector( helixCalo.getMomentum()[0], helixCalo.getMomentum()[1], helixCalo.getMomentum()[2] );

            _xyzTrackAtCalo = XYZVector( tsCalo->getReferencePoint()[0], tsCalo->getReferencePoint()[1], tsCalo->getReferencePoint()[2] );
        }
        //PFOs MCParticle
        const vector <LCObject*>& relationObjects = relation.getRelatedToObjects(pfo);
        const vector <float>& relationWeights = relation.getRelatedToWeights(pfo);
        int idxMaxWeight = std::distance(relationWeights.begin(), std::max_element(relationWeights.begin(), relationWeights.end()));
        MCParticle* mcPFO = dynamic_cast <MCParticle*> ( relationObjects[idxMaxWeight] );

        _weightMC = relationWeights[idxMaxWeight];
        _PDG = mcPFO->getPDG();
        _xyzVtxMC = XYZVector( mcPFO->getVertex()[0], mcPFO->getVertex()[1], mcPFO->getVertex()[2] );
        _tVtxMC = mcPFO->getTime();
        _pMC = XYZVector( mcPFO->getMomentum()[0], mcPFO->getMomentum()[1], mcPFO->getMomentum()[2] );

        // Fill tree after everything is done
        _tree->Fill();
        //clear vectors area
        _xyzSETHit.clear();
        _tSETHit.clear();
        _xyzECALHit.clear();
        _tECALHit.clear();
        _layerECALHit.clear();
        _eECALHit.clear();
    } //end of PFO loop
}

void TOFAnalysis::end(){
    double elapsedTime = duration<double>(system_clock::now() - _start).count();
    cout<<"Total number of events: "<<_nEvt<<endl;
    cout<<"Total elapsed time: "<<elapsedTime<<" sec"<<endl;
    cout<<"Total avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    _file->Write();
    cout<<_file->GetName()<<"   file is written in the current directory"<<endl;
}
