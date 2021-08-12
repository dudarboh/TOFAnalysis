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
#include "TString.h"

using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::FixedPadSizeTPCData;

SETAnalysis aSETAnalysis;

SETAnalysis::SETAnalysis() : Processor("SETAnalysis"){
    registerProcessorParameter(string("outputFile"),
                               string("Name of the output root file"),
                               _outputFileName,
                               string("SETAnalysis_RENAME.root"));

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

    _tree->Branch("has_set_hit", &_hasSetHit);
    _tree->Branch("n_ecal_hits", &_nEcalHits);

    _tree->Branch("pos_set_hit", &_posSetHit);
    _tree->Branch("pos_closest", &_posClosest);

    _tree->Branch("track_length_set", &_trackLengthSet);
    _tree->Branch("track_length_calo", &_trackLengthCalo);
    // _tree->Branch("track_length_integral", &_trackLengthIntegral);

    for(unsigned int i=0; i < std::size(_smearings); ++i ){
        _tree->Branch(Form("pos_fastest_%d", int(_smearings[i]) ), &( _posFastest[i] ) );
        _tree->Branch(Form("tof_fastest_%d", int(_smearings[i]) ), &( _tofFastest[i] ) );
        _tree->Branch(Form("tof_closest_%d", int(_smearings[i]) ), &( _tofClosest[i] ) );
        _tree->Branch(Form("tof_set_front_%d", int(_smearings[i]) ), &( _tofSetFront[i] ) );
        _tree->Branch(Form("tof_set_back_%d", int(_smearings[i]) ), &( _tofSetBack[i] ) );
        _tree->Branch(Form("tof_frank_fit_%d", int(_smearings[i]) ), &( _tofFrankFit[i] ) );
        _tree->Branch(Form("tof_frank_avg_%d", int(_smearings[i]) ), &( _tofFrankAvg[i] ) );
    }


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
        /////////////////////////////////// Load all info from PFO
        ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

        MCParticle* mc = getMcMaxWeight(pfoToMc, pfo);
        if(mc == nullptr) continue;

        const vector<Cluster*>& clusters = pfo->getClusters();
        if(clusters.size() != 1) continue;
        Cluster* cluster = clusters.at(0);
        _nEcalHits = getNEcalHits(cluster);
        if (_nEcalHits == 0) continue;

        const vector<Track*>& tracks = pfo->getTracks();
        if(tracks.size() != 1) continue;
        Track* track = tracks.at(0);
        TrackerHit* setHit = getSetHit(track, _tpcROuter);
        _hasSetHit = (setHit != nullptr);
        if (!_hasSetHit) continue;

        ///////////////////////////////////
        _pdg = mc->getPDG();
        _posSetHit.SetCoordinates( setHit->getPosition() );

        const LCObjectVec& setSimHits = spPointToSimHit.getRelatedToObjects( setHit );
        if (setSimHits.size() != 2) cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAA: "<<setSimHits.size()<<endl;



        /////////////////WRITE TRACK STATES/////////////////////////
        const TrackState* tsLast = track->getTrackState(TrackState::AtLastHit);
        _tsLastOmega = tsLast->getOmega();
        _tsLastTanL = tsLast->getTanLambda();
        _tsLastPhi = tsLast->getPhi();
        _tsLastD0 = tsLast->getD0();
        _tsLastZ0 = tsLast->getZ0();
        _tsLastPos.SetCoordinates( tsLast->getReferencePoint() );
        HelixClass helixLast;
        helixLast.Initialize_Canonical(_tsLastPhi, _tsLastD0, _tsLastZ0, _tsLastOmega, _tsLastTanL, _bField);
        _tsLastMom.SetCoordinates( helixLast.getMomentum() );

        const TrackState* tsCalo = track->getTrackState(TrackState::AtCalorimeter);
        _tsCaloOmega = tsCalo->getOmega();
        _tsCaloTanL = tsCalo->getTanLambda();
        _tsCaloPhi = tsCalo->getPhi();
        _tsCaloD0 = tsCalo->getD0();
        _tsCaloZ0 = tsCalo->getZ0();
        _tsCaloPos.SetCoordinates( tsCalo->getReferencePoint() );
        HelixClass helixCalo;
        helixCalo.Initialize_Canonical(_tsCaloPhi, _tsCaloD0, _tsCaloZ0, _tsCaloOmega, _tsCaloTanL, _bField);
        _tsCaloMom.SetCoordinates( helixCalo.getMomentum() );

        ///////////////////////WRITE TRACK LENGTHS/////////////////////
        _trackLengthSet = getTrackLength(track, TrackState::AtIP, TrackState::AtLastHit);
        _trackLengthCalo = getTrackLength(track, TrackState::AtIP, TrackState::AtCalorimeter);
        // _trackLengthIntegral = getTrackLengthIntegral(track);


        SimTrackerHit* setSimHitFront = dynamic_cast <SimTrackerHit*> (setSimHits.at(0) );
        SimTrackerHit* setSimHitBack = dynamic_cast <SimTrackerHit*> (setSimHits.at(1) );


        CalorimeterHit* closestHit = getClosestHit( cluster, _tsCaloPos );
        _posClosest.SetCoordinates( closestHit->getPosition() );

        for(unsigned int j=0; j < std::size(_smearings); ++j ){
            pair<XYZVectorF, double> fastestHit = getFastestHit( cluster, _smearings[j] / 1000. );
            _posFastest[j] = fastestHit.first;
            _tofFastest[j] = fastestHit.second - ( _posFastest[j] - _tsCaloPos ).r()/CLHEP::c_light;

            _tofClosest[j] = CLHEP::RandGauss::shoot( closestHit->getTime(), _smearings[j] / 1000. ) - ( _posClosest - _tsCaloPos ).r()/CLHEP::c_light;
            _tofSetFront[j] = CLHEP::RandGauss::shoot( setSimHitFront->getTime(), _smearings[j] / 1000. );
            _tofSetBack[j] =  CLHEP::RandGauss::shoot( setSimHitBack->getTime(), _smearings[j] / 1000. );
            _tofFrankFit[j] = getTofFrankFit( cluster, _tsCaloPos, _tsCaloMom, _smearings[j] / 1000. );
            _tofFrankAvg[j] = getTofFrankAvg( cluster, _tsCaloPos, _tsCaloMom, _smearings[j] / 1000. );
        }
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

TrackerHit* SETAnalysis::getSetHit(Track* track, double tpcROuter){
    const vector<TrackerHit*>& hits = track->getTrackerHits();
    XYZVector pos;
    //Performance: loop from the end. SET hits at the end!
    for (int i=hits.size()-1; i>=0; --i){
        pos.SetCoordinates( hits[i]->getPosition() );
        if ( pos.rho() > tpcROuter ) return hits[i];
    }
    return nullptr;
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


pair<XYZVectorF, double> SETAnalysis::getFastestHit( Cluster* cluster, double smearing ){
    pair<XYZVectorF, double> earliestHit{};

    earliestHit.second = std::numeric_limits<double>::max();
    for ( const auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        double timeHit = CLHEP::RandGauss::shoot( hit->getTime(), smearing );
        if( timeHit < earliestHit.second ){
            earliestHit.first.SetCoordinates( hit->getPosition() );
            earliestHit.second = timeHit;
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
                closestTime = CLHEP::RandGauss::shoot( hit->getTime() , smearing );
                closestDistanceToLine = dToLine;
            }
        }
        if ( closestDistanceToLine == std::numeric_limits<double>::max() ) continue;
        d.push_back(closestDistance);
        time.push_back(closestTime);
        d_err.push_back(0.);
        time_err.push_back(0.3);
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
