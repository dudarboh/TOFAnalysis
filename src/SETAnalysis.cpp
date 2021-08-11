#include "SETAnalysis.hpp"

#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/SimTrackerHit.h"
#include "HelixClass.h"

#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"

#include "marlinutil/GeometryUtil.h"
#include "marlinutil/CalorimeterHitType.h"

#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "Math/Vector3D.h"
#include "TGraphErrors.h"
#include "TF1.h"

using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::FixedPadSizeTPCData;

SETAnalysis aSETAnalysis;

SETAnalysis::SETAnalysis() : Processor("SETAnalysis"){
    registerProcessorParameter(string("outputFile"),
                               string("Name of the output root file"),
                               _outputFileName,
                               string("SETAnalysis_RENAME.root"));

    registerProcessorParameter(string("smearing"),
                               string("Time resolution in ***ps***"),
                               _smearing,
                               double(0.));

}

void SETAnalysis::init(){
    _nEvent = 0;
    _file.reset( new TFile(_outputFileName.c_str(), "RECREATE") );
    _tree.reset( new TTree("SETAnalysis", "SETAnalysis") );

    _tree->Branch("pdg", &_pdg);
    _tree->Branch("ts_last_pos", &_tsLastPos);
    _tree->Branch("ts_last_mom", &_tsLastMom);
    _tree->Branch("ts_last_omega", &_tsLastOmega);
    _tree->Branch("ts_last_tanL", &_tsLastTanL);
    _tree->Branch("ts_last_phi", &_tsLastPhi);
    _tree->Branch("ts_last_d0", &_tsLastD0);
    _tree->Branch("ts_last_z0", &_tsLastZ0);

    _tree->Branch("ts_calo_pos", &_tsCaloPos);
    _tree->Branch("ts_calo_mom", &_tsCaloMom);
    _tree->Branch("ts_calo_omega", &_tsCaloOmega);
    _tree->Branch("ts_calo_tanL", &_tsCaloTanL);
    _tree->Branch("ts_calo_phi", &_tsCaloPhi);
    _tree->Branch("ts_calo_d0", &_tsCaloD0);
    _tree->Branch("ts_calo_z0", &_tsCaloZ0);

    _tree->Branch("n_ecal_hits", &_nEcalHits);
    _tree->Branch("pos_closest", &_posClosest);
    _tree->Branch("tof_closest", &_tofClosest);

    _tree->Branch("pos_fastest", &_posFastest);
    _tree->Branch("tof_fastest", &_tofFastest);

    _tree->Branch("tof_frank_fit", &_tofFrankFit);
    _tree->Branch("tof_frank_avg", &_tofFrankAvg);

    _tree->Branch("n_set_hits", &_nSetHits);
    _tree->Branch("set_hit_pos", &_setHitPos);
    _tree->Branch("set_hit_time", &_setHitTime);
    _tree->Branch("set_pos_true", &_setPosTrue);

    _tree->Branch("track_length_set", &_trackLengthSet);
    _tree->Branch("track_length_calo", &_trackLengthCalo);
    // _tree->Branch("track_length_integral", &_trackLengthIntegral);

    _bField = MarlinUtil::getBzAtOrigin();
    _tpcROuter = getTpcR().second;
}


void SETAnalysis::processEvent(LCEvent* evt){
    //profiler notes:
    // LCRelationNavigators - are heavy. Remove if possible
    //
    ++_nEvent;
    cout<<"****** Event: "<<_nEvent<<endl;
    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator pfoToMc( evt->getCollection("RecoMCTruthLink") );

    //Don't need for now.. Might need in the future
    LCRelationNavigator spPointToSimHit( evt->getCollection("SETSpacePointRelations") );
    // LCRelationNavigator caloHitToSimHit( evt->getCollection("EcalBarrelRelationsSimRec") );

    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        MCParticle* mc = getMcMaxWeight(pfoToMc, pfo);
        const vector<Track*>& tracks = pfo->getTracks();
        const vector<Cluster*>& clusters = pfo->getClusters();

        if(mc == nullptr || tracks.size() != 1 || clusters.size() != 1) continue;

        Track* track = tracks[0];
        Cluster* cluster = clusters[0];
        // const TrackState* tsIp = track->getTrackState(TrackState::AtIP);
        // const TrackState* tsFirst = track->getTrackState(TrackState::AtFirstHit);
        const TrackState* tsLast = track->getTrackState(TrackState::AtLastHit);
        const TrackState* tsCalo = track->getTrackState(TrackState::AtCalorimeter);

        vector<TrackerHit*> setHits = getSetHits(track, _tpcROuter);
        _nSetHits = setHits.size();
        _nEcalHits = getNEcalHits(cluster);
        if (_nSetHits == 0 || _nEcalHits == 0) continue;

        _pdg = mc->getPDG();
        _setHitPos.SetCoordinates( setHits[0]->getPosition() );
        _setHitTime = CLHEP::RandGauss::shoot( setHits[0]->getTime(), _smearing / 1000. ) ;

        SimTrackerHit* setSimHit = dynamic_cast <SimTrackerHit*> (spPointToSimHit.getRelatedToObjects( setHits[0] )[0] );
        _setPosTrue.SetCoordinates( setSimHit->getPosition() );

        _tsLastOmega = tsLast->getOmega();
        _tsLastTanL = tsLast->getTanLambda();
        _tsLastPhi = tsLast->getPhi();
        _tsLastD0 = tsLast->getD0();
        _tsLastZ0 = tsLast->getZ0();
        _tsLastPos.SetCoordinates( tsLast->getReferencePoint() );
        HelixClass helixLast;
        helixLast.Initialize_Canonical(_tsLastPhi, _tsLastD0, _tsLastZ0, _tsLastOmega, _tsLastTanL, _bField);
        _tsLastMom.SetCoordinates( helixLast.getMomentum() );

        _tsCaloOmega = tsCalo->getOmega();
        _tsCaloTanL = tsCalo->getTanLambda();
        _tsCaloPhi = tsCalo->getPhi();
        _tsCaloD0 = tsCalo->getD0();
        _tsCaloZ0 = tsCalo->getZ0();
        _tsCaloPos.SetCoordinates( tsCalo->getReferencePoint() );
        HelixClass helixCalo;
        helixCalo.Initialize_Canonical(_tsCaloPhi, _tsCaloD0, _tsCaloZ0, _tsCaloOmega, _tsCaloTanL, _bField);
        _tsCaloMom.SetCoordinates( helixCalo.getMomentum() );

        _trackLengthSet = getTrackLength(track, TrackState::AtIP, TrackState::AtLastHit);
        _trackLengthCalo = getTrackLength(track, TrackState::AtIP, TrackState::AtCalorimeter);
        // _trackLengthIntegral = getTrackLengthIntegral(track);
        _tofFrankFit = getTofFrankFit( cluster, _tsCaloPos, _tsCaloMom, _smearing / 1000. );
        _tofFrankAvg = getTofFrankAvg( cluster, _tsCaloPos, _tsCaloMom, _smearing / 1000. );


        CalorimeterHit* closestHit =  getClosestHit( cluster, _tsCaloPos );
        _posClosest .SetCoordinates( closestHit->getPosition() );
        _tofClosest = CLHEP::RandGauss::shoot( closestHit->getTime(), _smearing / 1000. ) - ( _posClosest - _tsCaloPos ).r()/CLHEP::c_light;

        CalorimeterHit* fastestHit =  getFastestHit( cluster );
        _posFastest.SetCoordinates( fastestHit->getPosition() );
        _tofFastest = CLHEP::RandGauss::shoot( fastestHit->getTime(), _smearing / 1000. ) - ( _posFastest - _tsCaloPos ).r()/CLHEP::c_light;

        _tree->Fill();
    }

}

void SETAnalysis::end(){
    _file->Write();
}


MCParticle* SETAnalysis::getMcMaxWeight(LCRelationNavigator pfoToMc, ReconstructedParticle* pfo){
    MCParticle* mc = nullptr;
    const vector <LCObject*>& mcs = pfoToMc.getRelatedToObjects(pfo);
    const vector <float>& mcWeights = pfoToMc.getRelatedToWeights(pfo);
    if (mcs.size() == 0) return mc;
    int maxW = std::max_element(mcWeights.begin(), mcWeights.end()) - mcWeights.begin();
    mc = dynamic_cast <MCParticle*> ( mcs[maxW] );
    return mc;
}


pair<double, double> SETAnalysis::getTpcR(){
    const Detector& detector = Detector::getInstance();
    const DetElement tpcDet = detector.detector("TPC");
    const FixedPadSizeTPCData* tpc = tpcDet.extension <FixedPadSizeTPCData>();
    double rInner = tpc->rMinReadout/dd4hep::mm;
    double rOuter = tpc->rMaxReadout/dd4hep::mm;
    return std::make_pair(rInner, rOuter);
}

vector<TrackerHit*> SETAnalysis::getSetHits(Track* track, double tpcROuter){
    const vector<TrackerHit*>& hits = track->getTrackerHits();
    vector<TrackerHit*> setHits;
    XYZVector pos;
    for(const auto hit: hits){
        pos.SetCoordinates( hit->getPosition() );
        if ( pos.rho() > tpcROuter ) setHits.push_back(hit);
    }
    return setHits;
}


double SETAnalysis::getTrackLength(Track* track, int from, int to){
    // AtLastHit should give the identical omega and tanL as at the AtCalorimeter.
    // I did check that.
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    double omega = tsEcal->getOmega();
    double tanL = tsEcal->getTanLambda();

    const TrackState* ts1 = track->getTrackState(from);
    const TrackState* ts2 = track->getTrackState(to);
    double phi1 = ts1->getPhi();
    double phi2 = ts2->getPhi();

    return std::abs( (phi1 - phi2)/omega )*std::sqrt(1. + tanL*tanL);
}


double SETAnalysis::getTrackLengthIntegral(Track* track){
    double trackLength = 0.;

    //*** track length from IP state to 1st tracker hit state***
    const TrackState* tsIp = track->getTrackState(TrackState::AtIP);
    double phiIP = tsIp->getPhi();
    double omegaIP = tsIp->getOmega();
    double tanLIP = tsIp->getTanLambda();
    const TrackState* tsFirst = track->getTrackState(TrackState::AtFirstHit);
    double phiFirst = tsFirst->getPhi();
    trackLength += std::abs( (phiIP - phiFirst)/omegaIP )*std::sqrt(1. + tanLIP*tanLIP);

    //*** in the tracker region sum hit-by-hit assuming straight line path between hits***
    vector <TrackerHit*> trackHits = track->getTrackerHits();

    auto sortByR = [](TrackerHit* a, TrackerHit* b) {
        XYZVector posA, posB;
        posA.SetCoordinates( a->getPosition() );
        posB.SetCoordinates( b->getPosition() );
        return posA.r() < posB.r();
    };
    sort(trackHits.begin(), trackHits.end(), sortByR);

    for (size_t j=1; j < trackHits.size(); ++j){
        XYZVector pos1, pos2;
        pos1.SetCoordinates( trackHits[j-1]->getPosition() );
        pos2.SetCoordinates( trackHits[j]->getPosition() );
        trackLength += (pos2-pos1).r();
    }


    // track length from last tracker hit to the ECAL surface
    const TrackState* tsLast = track->getTrackState(TrackState::AtLastHit);
    double phiLast = tsLast->getPhi();
    const TrackState* tsCalo = track->getTrackState(TrackState::AtCalorimeter);
    double phiCalo = tsCalo->getPhi();
    double omegaCalo = tsCalo->getOmega();
    double tanLCalo = tsCalo->getTanLambda();

    trackLength += std::abs( (phiLast - phiCalo)/omegaCalo )*std::sqrt(1. + tanLCalo*tanLCalo);
    return trackLength;
}


CalorimeterHit* SETAnalysis::getClosestHit( Cluster* cluster, XYZVectorF posTrackAtCalo){
    CalorimeterHit* closestHit = nullptr;

    double closestDistance = std::numeric_limits<double>::max();
    for ( const auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        XYZVectorF hitPos;
        hitPos.SetCoordinates( hit->getPosition() );
        double dToEntry = (hitPos - posTrackAtCalo).r();
        if( dToEntry < closestDistance ){
            closestDistance = dToEntry;
            closestHit = hit;
        }
    }
    return closestHit;
}


CalorimeterHit* SETAnalysis::getFastestHit( Cluster* cluster ){
    CalorimeterHit* earliestHit = nullptr;

    double earliestTime = std::numeric_limits<double>::max();
    for ( const auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        double timeHit = hit->getTime();
        if( timeHit < earliestTime ){
            earliestTime = timeHit;
            earliestHit = hit;
        }
    }
    return earliestHit;
}


double SETAnalysis::getTofFrankFit( Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing, unsigned int nLayers ){
    vector <double> d;
    vector <double> time;
    vector <double> d_err;
    vector <double> time_err;

    for (unsigned int l=0; l < nLayers; ++l){
        double closestDistance = std::numeric_limits<double>::max();
        double closestTime = std::numeric_limits<double>::max();
        double closestDistanceToLine = std::numeric_limits<double>::max();

        for ( const auto& hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if ( (! isECALHit) || (hitType.layer() != l) ) continue;

            XYZVectorF pos;
            pos.SetCoordinates( hit->getPosition() );
            double dToLine = (pos - posTrackAtCalo).Cross(momTrackAtCalo.unit()).r();
            if (dToLine < closestDistanceToLine){
                closestDistance = (pos - posTrackAtCalo).r();
                closestTime = CLHEP::RandGauss::shoot(hit->getTime(), smearing);
                closestDistanceToLine = dToLine;
            }
        }
        if ( closestDistanceToLine == std::numeric_limits<double>::max() ) continue;
        d.push_back(closestDistance);
        time.push_back(closestTime);
        d_err.push_back(0.);
        time_err.push_back(0.1);
    }
    // Can't fit 0 or 1 point. Must return something meaningfull
    if ( d.size() == 0 ) return 0.;
    else if ( d.size() == 1 ) return time[0] - d[0]/CLHEP::c_light;

    TGraphErrors gr(d.size(), &d[0], &time[0], &d_err[0], &time_err[0]);
    gr.Fit("pol1", "Q");
    return gr.GetFunction("pol1")->GetParameter(0);
}


double SETAnalysis::getTofFrankAvg( Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing, unsigned int nLayers ){
    double tof = 0.;
    int nHits = 0;

    for (unsigned int l=0; l < nLayers; ++l){
        double closestDistance = std::numeric_limits<double>::max();
        double closestTime = std::numeric_limits<double>::max();
        double closestDistanceToLine = std::numeric_limits<double>::max();

        for ( const auto& hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if ( (! isECALHit) || (hitType.layer() != l) ) continue;

            XYZVectorF pos;
            pos.SetCoordinates( hit->getPosition() );
            double dToLine = (pos - posTrackAtCalo).Cross(momTrackAtCalo.unit()).r();
            if (dToLine < closestDistanceToLine){
                closestDistance = (pos - posTrackAtCalo).r();
                closestTime = CLHEP::RandGauss::shoot(hit->getTime(), smearing);
                closestDistanceToLine = dToLine;
            }
        }
        if ( closestDistanceToLine == std::numeric_limits<double>::max() ) continue;

        tof += closestTime - closestDistance/CLHEP::c_light;
        ++nHits;
    }
    if (nHits == 0) return 0.;

    return tof / nHits;
}



int SETAnalysis::getNEcalHits(Cluster* cluster){
    int nHits = 0;
    for ( const auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (isECALHit) ++nHits;
    }
    return nHits;
}
