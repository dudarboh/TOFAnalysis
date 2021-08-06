#include "TofAnalysis.hpp"

#include <limits>
#include <algorithm>

#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "EVENT/TrackState.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/MCParticle.h"

#include "HelixClass.h"

#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <marlin/Global.h>
#include <marlin/ProcessorEventSeeder.h>
#include "marlinutil/CalorimeterHitType.h"
#include "UTIL/LCRelationNavigator.h"
#include "marlinutil/GeometryUtil.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using EVENT::LCCollection;
using EVENT::ReconstructedParticle;
using EVENT::Track;
using EVENT::TrackState;
using EVENT::TrackerHit;
using EVENT::CalorimeterHit;
using EVENT::Cluster;

using marlin::Processor;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::FixedPadSizeTPCData;
using dd4hep::rec::Vector3D;
// using ROOT::Math::XYZVector;




TofAnalysis aTofAnalysis;

TofAnalysis::TofAnalysis() : Processor("TofAnalysis"){
    registerProcessorParameter(string("outputFile"),
                               string("Name of the output root file"),
                               _outputFileName,
                               string("TofAnalysis_RENAME.root"));

    registerProcessorParameter(string("writeTpcHits"),
                              string("Write TPC hits to the output root file"),
                              _writeTpcHits,
                              bool(false));

    registerProcessorParameter(string("writeEcalHits"),
                              string("Write ECAL hits to the output root file"),
                              _writeEcalHits,
                              bool(false));

    registerProcessorParameter(string("writeSetHits"),
                              string("Write SET hits to the output root file"),
                              _writeSetHits,
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
    if(_writeTpcHits){
        _tree->Branch("nTPCHits", &_nTPCHits);
        _tree->Branch("posTPCHit", &_posTPCHit);
        _tree->Branch("tTPCHit", &_tTPCHit);
    }
    if(_writeEcalHits){
        _tree->Branch("nECALHits", &_nECALHits);
        _tree->Branch("posECALHit", &_posECALHit);
        _tree->Branch("tECALHit", &_tECALHit);
        _tree->Branch("layerECALHit", &_layerECALHit);
        _tree->Branch("eECALHit", &_eECALHit);
    }
    if(_writeSetHits){
        _tree->Branch("nSETHits", &_nSETHits);
        _tree->Branch("posSETHit", &_posSETHit);
        _tree->Branch("tSETHit", &_tSETHit);
    }
    _tree->Branch("momIP", &_momIP);
    _tree->Branch("momECAL", &_momECAL);

    _tree->Branch("trackLengthIP", &_trackLengthIP);
    _tree->Branch("trackLengthECAL", &_trackLengthECAL);
    _tree->Branch("trackLengthIntegral", &_trackLengthIntegral);

    _tree->Branch("trackLengthIPSET", &_trackLengthIPSET);
    _tree->Branch("trackLengthECALSET", &_trackLengthECALSET);
    _tree->Branch("trackLengthIntegralSET", &_trackLengthIntegralSET);

    _tree->Branch("tofClosest", &_tofClosest);
    _tree->Branch("tofFastest", &_tofFastest);
    _tree->Branch("tofFrankFit", &_tofFrankFit);
    _tree->Branch("tofFrankAvg", &_tofFrankAvg);
    _tree->Branch("tofSet", &_tofSet);


    marlin::Global::EVENTSEEDER->registerProcessor(this);
    return;
}

void TofAnalysis::processEvent(LCEvent* evt){
    CLHEP::HepRandom::setTheSeed( marlin::Global::EVENTSEEDER->getSeed(this) );
    ++_nEvt;

    LCCollection* colPFO = evt->getCollection("PandoraPFOs");
    LCCollection* colRelation = evt->getCollection("RecoMCTruthLink");
    LCRelationNavigator relation(colRelation);

    for (int i=0; i<colPFO->getNumberOfElements(); ++i){

        ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> ( colPFO->getElementAt(i) );

        const vector <LCObject*>& pfoMCs = relation.getRelatedToObjects(pfo);
        const vector <float>& pfoWeights = relation.getRelatedToWeights(pfo);
        int maxWeight = std::max_element(pfoWeights.begin(), pfoWeights.end()) - pfoWeights.begin();
        MCParticle* mcPFO = dynamic_cast <MCParticle*> ( pfoMCs[maxWeight] );


        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();
        // Only simple cases of PFOs
        if( nClusters != 1 || nTracks != 1) continue;


        if(_writeTpcHits) writeTpcHits(pfo);
        if(_writeEcalHits) writeEcalHits(pfo);
        if(_writeSetHits) writeSetHits(pfo);

        _pdg = mcPFO->getPDG();
        _momIP = calcMomentum(pfo, TrackState::AtIP, _bField);
        _momECAL = calcMomentum(pfo, TrackState::AtCalorimeter, _bField);
        _trackLengthIP = calcTrackLength(pfo, TrackState::AtIP);
        _trackLengthECAL = calcTrackLength(pfo, TrackState::AtCalorimeter);
        _trackLengthIntegral = calcTrackLengthIntegral(pfo);
        _trackLengthIPSET = calcTrackLengthSET(pfo, TrackState::AtIP);
        _trackLengthECALSET = calcTrackLengthSET(pfo, TrackState::AtLastHit);
        _trackLengthIntegralSET = calcTrackLengthIntegralSET(pfo);

        _tofClosest = calcTofClosest(pfo, 0.);
        _tofFastest = calcTofFastest(pfo, 0.);
        // _tofFrankFit = calcTofFrankFit(pfo, 0.);
        _tofFrankAvg = calcTofFrankAvg(pfo, 0.);
        _tofSet = calcTofSet(pfo, 0.);


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


void TofAnalysis::writeTpcHits(ReconstructedParticle* pfo){
    _posTPCHit.clear();
    _tTPCHit.clear();

    if (pfo->getTracks().size() == 0){
        _nTPCHits = 0;
        return;
    }

    const Track* track = pfo->getTracks()[0];
    for ( const auto& hit : track->getTrackerHits() ){
        Vector3D pos( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
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


void TofAnalysis::writeEcalHits(ReconstructedParticle* pfo){
    _posECALHit.clear();
    _tECALHit.clear();
    _layerECALHit.clear();
    _eECALHit.clear();

    const Cluster* cluster = pfo->getClusters()[0];
    for ( const auto& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );

        if (isECALHit){
            Vector3D pos( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
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


void TofAnalysis::writeSetHits(ReconstructedParticle* pfo){
    _posSETHit.clear();
    _tSETHit.clear();

    if (pfo->getTracks().size() == 0){
        _nSETHits = 0;
        return;
    }

    const Track* track = pfo->getTracks()[0];
    for ( const auto& hit : track->getTrackerHits() ){
        Vector3D pos( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
        bool isSETHit = pos.rho() > _tpcR.second;
        if (isSETHit){
            cout<<"hit type:   "<<hit->getType()<<endl;
            for ( const auto& raw_hit : hit->getRawHits() ){
                // cout<<raw_hit->getType()<<endl;

            }

            float time = hit->getTime();

            _posSETHit.push_back(pos);
            _tSETHit.push_back(time);
        }
    }
    _nSETHits = _posSETHit.size();
    return;
}


std::pair<double, double> TofAnalysis::getTpcR(const Detector& detector){
    // const Detector& detector = Detector::getInstance();
    const DetElement tpcDet = detector.detector("TPC");
    const FixedPadSizeTPCData* tpc = tpcDet.extension <FixedPadSizeTPCData>();
    double rInner = tpc->rMinReadout/dd4hep::mm;
    double rOuter = tpc->rMaxReadout/dd4hep::mm;
    return std::make_pair(rInner, rOuter);
}


Vector3D TofAnalysis::calcMomentum(ReconstructedParticle* pfo, int location, double bField){
    if ( pfo->getTracks().size() == 0 ) return Vector3D();

    const Track* track = pfo->getTracks()[0];

    const TrackState* ts = track->getTrackState(location);
    double phi = ts->getPhi();
    double omega = ts->getOmega();
    double tanL = ts->getTanLambda();
    double d0 = ts->getD0();
    double z0 = ts->getZ0();

    HelixClass helix;
    helix.Initialize_Canonical(phi, d0, z0, omega, tanL, bField);

    return Vector3D( helix.getMomentum() );
}


double TofAnalysis::calcTrackLength(ReconstructedParticle* pfo, int location){
    if ( pfo->getTracks().size() == 0 ) return 0.;

    const Track* track = pfo->getTracks()[0];
    const TrackState* tsIp = track->getTrackState(TrackState::AtIP);
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);

    double phiIp = tsIp->getPhi();
    double phiEcal = tsEcal->getPhi();

    double omega, tanL;
    if (location == TrackState::AtIP){
        omega = tsIp->getOmega();
        tanL = tsIp->getTanLambda();
    }
    else if (location  == TrackState::AtCalorimeter){
        omega = tsEcal->getOmega();
        tanL = tsEcal->getTanLambda();
    }
    else return 0.;

    return std::abs( (phiIp - phiEcal)/omega )*std::sqrt(1. + tanL*tanL);
}


double TofAnalysis::calcTrackLengthIntegral(ReconstructedParticle* pfo){
    if ( pfo->getTracks().size() == 0 ) return 0.;

    const Track* track = pfo->getTracks()[0];
    const TrackState* tsIP = track->getTrackState(TrackState::AtIP);
    const TrackState* tsFirst = track->getTrackState(TrackState::AtFirstHit);
    const TrackState* tsLast = track->getTrackState(TrackState::AtLastHit);
    const TrackState* tsECAL = track->getTrackState(TrackState::AtCalorimeter);

    double phiIP = tsIP->getPhi();
    double phiFirst = tsFirst->getPhi();
    double phiLast = tsLast->getPhi();
    double phiECAL = tsECAL->getPhi();

    double omegaIP = tsIP->getOmega();
    double omegaECAL = tsECAL->getOmega();
    double tanLIP = tsIP->getTanLambda();
    double tanLECAL = tsECAL->getTanLambda();


    double tr_len = 0;
    tr_len += std::abs( (phiIP - phiFirst)/omegaIP )*std::sqrt(1. + tanLIP*tanLIP);
    tr_len += std::abs( (phiLast - phiECAL)/omegaECAL )*std::sqrt(1. + tanLECAL*tanLECAL);

    std::vector <TrackerHit*> trackHits = track->getTrackerHits();

    auto sortbyr = [](const TrackerHit* a, const TrackerHit* b){
        Vector3D posA( a->getPosition() );
        Vector3D posB( b->getPosition() );
        return posA.r() < posB.r();
    };

    std::sort(trackHits.begin(), trackHits.end(), sortbyr);

    for (size_t j=1; j < trackHits.size(); ++j){
        Vector3D pos1( trackHits[j-1]->getPosition() );
        Vector3D pos2( trackHits[j]->getPosition() );
        tr_len += (pos2-pos1).r();
    }
    return tr_len;
}


double TofAnalysis::calcTrackLengthSET(ReconstructedParticle* pfo, int location){
    if ( pfo->getTracks().size() == 0 ) return 0.;

    const Track* track = pfo->getTracks()[0];
    const TrackState* tsIp = track->getTrackState(TrackState::AtIP);
    const TrackState* tsLast = track->getTrackState(TrackState::AtFirstHit);

    double phiIp = tsIp->getPhi();
    double phiLast = tsLast->getPhi();

    double omega, tanL;
    if (location == TrackState::AtIP){
        omega = tsIp->getOmega();
        tanL = tsIp->getTanLambda();
    }
    else if (location  == TrackState::AtFirstHit){
        omega = tsLast->getOmega();
        tanL = tsLast->getTanLambda();
    }
    else return 0.;

    return std::abs( (phiIp - phiLast)/omega )*std::sqrt(1. + tanL*tanL);
}


double TofAnalysis::calcTrackLengthIntegralSET(ReconstructedParticle* pfo){
    if ( pfo->getTracks().size() == 0 ) return 0.;

    const Track* track = pfo->getTracks()[0];
    const TrackState* tsIP = track->getTrackState(TrackState::AtIP);
    const TrackState* tsFirst = track->getTrackState(TrackState::AtFirstHit);

    double phiIP = tsIP->getPhi();
    double phiFirst = tsFirst->getPhi();

    double omegaIP = tsIP->getOmega();
    double tanLIP = tsIP->getTanLambda();


    double tr_len = 0;
    tr_len += std::abs( (phiIP - phiFirst)/omegaIP )*std::sqrt(1. + tanLIP*tanLIP);

    std::vector <TrackerHit*> trackHits = track->getTrackerHits();

    auto sortbyr = [](const TrackerHit* a, const TrackerHit* b){
        Vector3D posA( a->getPosition() );
        Vector3D posB( b->getPosition() );
        return posA.r() < posB.r();
    };

    std::sort(trackHits.begin(), trackHits.end(), sortbyr);

    for (size_t j=1; j < trackHits.size(); ++j){
        Vector3D pos1( trackHits[j-1]->getPosition() );
        Vector3D pos2( trackHits[j]->getPosition() );
        tr_len += (pos2-pos1).r();
    }
    return tr_len;
}


double TofAnalysis::calcTofClosest( ReconstructedParticle* pfo, double smearing){
    //No track no "closest to track" time
    if ( pfo->getTracks().size() == 0 ) return 0.;

    // Get Track position at ECAL
    const Track* track = pfo->getTracks()[0];
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D posTrackAtCalo( ts->getReferencePoint() );

    double closestDistance = std::numeric_limits<double>::max();
    double time = std::numeric_limits<double>::max();

    //Loop over ECAL hits
    const Cluster* cluster = pfo->getClusters()[0];
    for ( const auto& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        Vector3D pos( hit->getPosition() );
        if( (pos - posTrackAtCalo).r() < closestDistance ){
            closestDistance = (pos - posTrackAtCalo).r();
            time = CLHEP::RandGauss::shoot(hit->getTime(), smearing);
        }
    }
    return time - closestDistance/CLHEP::c_light;
}


double TofAnalysis::calcTofFastest( ReconstructedParticle* pfo, double smearing){
    //No track --- no correction for distance is possible time
    if ( pfo->getTracks().size() == 0 ) return 0.;

    // Get Track position at ECAL
    const Track* track = pfo->getTracks()[0];
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D posTrackAtCalo( ts->getReferencePoint() );

    double distance = std::numeric_limits<double>::max();
    double time = std::numeric_limits<double>::max();

    //Loop over ECAL hits
    const Cluster* cluster = pfo->getClusters()[0];
    for ( const auto& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        Vector3D pos( hit->getPosition() );
        if( hit->getTime() < time ){
            time = CLHEP::RandGauss::shoot(hit->getTime(), smearing);
            distance = (pos - posTrackAtCalo).r();
        }
    }

    return time - distance/CLHEP::c_light;
}


double TofAnalysis::calcTofFrankFit( ReconstructedParticle* pfo, double smearing){
    //No track --- no correction for distance is possible time
    if ( pfo->getTracks().size() == 0 ) return 0;

    // Get Track position at ECAL
    const Track* track = pfo->getTracks()[0];
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D posTrackAtCalo( ts->getReferencePoint() );


    double bField = MarlinUtil::getBzAtOrigin();
    Vector3D momAtECAL = calcMomentum(pfo, TrackState::AtCalorimeter, bField);

    const Cluster* cluster = pfo->getClusters()[0];

    unsigned int nLayers = 10;
    std::vector <double> d(nLayers);
    std::vector <double> time(nLayers);
    std::vector <double> d_err(nLayers);
    std::vector <double> time_err(nLayers);
    std::vector <double> closestDistanceToLine(nLayers);

    for (unsigned int l=0; l < nLayers; ++l){
        closestDistanceToLine[l] = std::numeric_limits<double>::max();
        d_err[l] = 0.;
        time_err[l] = 0.1;
    }

    for (unsigned int l=0; l < nLayers; ++l){
        for ( const auto& hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if (! isECALHit) continue;
            if (hitType.layer() != l) continue;

            Vector3D pos( hit->getPosition() );
            double dToLine = (pos - posTrackAtCalo).cross(momAtECAL.unit()).r();
            if (dToLine < closestDistanceToLine[l]){
                closestDistanceToLine[l] = dToLine;
                time[l] = CLHEP::RandGauss::shoot(hit->getTime(), smearing);
                d[l] = (pos - posTrackAtCalo).r();
            }
        }
    }

    TGraphErrors gr(nLayers, &d[0], &time[0], &d_err[0], &time_err[0]);
    TF1* fit = new TF1("fit", "pol1");
    gr.Fit(fit, "Q");
    fit = gr.GetFunction("fit");
    return fit->GetParameter(0);
}


double TofAnalysis::calcTofFrankAvg( ReconstructedParticle* pfo, double smearing){
    //No track --- no correction for distance is possible time
    if ( pfo->getTracks().size() == 0 ) return 0;

    // Get Track position at ECAL
    const Track* track = pfo->getTracks()[0];
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D posTrackAtCalo( ts->getReferencePoint() );


    double bField = MarlinUtil::getBzAtOrigin();
    Vector3D momAtECAL = calcMomentum(pfo, TrackState::AtCalorimeter, bField);

    const Cluster* cluster = pfo->getClusters()[0];

    unsigned int nLayers = 10;
    std::vector <double> d(nLayers);
    std::vector <double> time(nLayers);
    std::vector <double> closestDistanceToLine(nLayers);

    double tof = 0.;

    for (unsigned int l=0; l < nLayers; ++l){
        closestDistanceToLine[l] = std::numeric_limits<double>::max();
    }

    for (unsigned int l=0; l < nLayers; ++l){
        for ( const auto& hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if (! isECALHit) continue;
            if (hitType.layer() != l) continue;

            Vector3D pos( hit->getPosition() );
            double dToLine = (pos - posTrackAtCalo).cross(momAtECAL.unit()).r();
            if (dToLine < closestDistanceToLine[l]){
                closestDistanceToLine[l] = dToLine;
                time[l] = CLHEP::RandGauss::shoot(hit->getTime(), smearing);
                d[l] = (pos - posTrackAtCalo).r();
            }
        }
        tof += (time[l] - d[l]/CLHEP::c_light) / nLayers;
    }
    return tof;
}

double TofAnalysis::calcTofSet( ReconstructedParticle* pfo, double smearing){
    if (pfo->getTracks().size() == 0){
        return 0.;
    }

    const Detector& detector = Detector::getInstance();
    const DetElement tpcDet = detector.detector("TPC");
    const FixedPadSizeTPCData* tpc = tpcDet.extension <FixedPadSizeTPCData>();
    double rOuter = tpc->rMaxReadout/dd4hep::mm;

    const Track* track = pfo->getTracks()[0];
    double tof = std::numeric_limits<double>::max();
    for ( const auto& hit : track->getTrackerHits() ){
        Vector3D pos( hit->getPosition() );
        bool isSETHit = pos.rho() > rOuter;

        if (isSETHit){
            if (hit->getTime() < tof) tof = hit->getTime();
        }
    }
    return tof;
}


int TofAnalysis::findShowerStart(ReconstructedParticle* pfo){
    const Cluster* cluster = pfo->getClusters()[0];

    std::cout<<"Maybe one day . . ."<<std::endl;
    // auto findIter = integrations.find( _integration_method ) ;
    return 0.;
}
