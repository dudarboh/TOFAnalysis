#include "SETAnalysis.hpp"

#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/SimTrackerHit.h"
#include "HelixClass.h"

#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Vector3D.h"

#include "marlinutil/GeometryUtil.h"
#include "marlinutil/CalorimeterHitType.h"

#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "Math/Vector3D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TString.h"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackImpl.h>

#include "streamlog/streamlog.h"


using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::FixedPadSizeTPCData;

using namespace MarlinTrk;

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
    _tree->Branch("track_length_refit", &_trackLengthRefit);
    _tree->Branch("track_length_adrian", &_trackLengthAdrian);

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
    _trkSystem = Factory::createMarlinTrkSystem( "DDKalTest", nullptr, "" ) ;

    _trkSystem->setOption( IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init() ;

}


void SETAnalysis::processEvent(LCEvent* evt){
    ++_nEvent;
    cout<<"****** Event: "<<_nEvent<<endl;

    // set the correct configuration for the tracking system for this event
    TrkSysConfig< IMarlinTrkSystem::CFG::useQMS> mson( _trkSystem, true );
    TrkSysConfig< IMarlinTrkSystem::CFG::usedEdx> elosson( _trkSystem, true);
    TrkSysConfig< IMarlinTrkSystem::CFG::useSmoothing> smoothon( _trkSystem, true);

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
        _trackLengthAdrian = getTrackLengthAdrian(track);
        _trackLengthRefit = getTrackLengthRefit(track);
        cout<<"Calo track length: "<<_trackLengthCalo<<endl;
        cout<<"Adrian track length: "<<_trackLengthAdrian<<endl;
        cout<<"Refit track length: "<<_trackLengthRefit<<endl;

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


double SETAnalysis::getTrackLengthAdrian(Track* track){
    // version of Mitaroff arxiv:2107.02031v1
    // approach 1, using only the track states
    // equation 20 with phi*=phi
    double trackLength = 0.;

    //*** track length from IP state to 1st tracker hit state***
    const TrackState* tsIp = track->getTrackState(TrackState::AtIP);
    double phiIp = tsIp->getPhi();
    double d0Ip = tsIp->getD0();
    double z0Ip = tsIp->getZ0();
    double omegaIp = tsIp->getOmega();
    double tanLIp = tsIp->getTanLambda();
    HelixClass helixIp;
    helixIp.Initialize_Canonical(phiIp, d0Ip, z0Ip, omegaIp, tanLIp, _bField);
    double qIp = helixIp.getCharge();
    double rhIp = helixIp.getRadius();
    dd4hep::rec::Vector3D pIp= helixIp.getMomentum();
    double khIp = 1./rhIp;

    XYZVector posIp = XYZVector(tsIp->getReferencePoint()[0], tsIp->getReferencePoint()[1], tsIp->getReferencePoint()[2]);
    // streamlog_out(DEBUG) << "AtIP q="<<qIp<<" Bz="<<B<<" Ku="<<CLHEP::c_light<<" z0="<<z0Ip<<" p="<<pIp.r()<<" rh="<<rhIp<<" posR:"<<posIp.r()<<std::endl;

    std::vector <TrackState*> trackSt = track->getTrackStates();
    double phi[1000], k[1000], z[1000];

    phi[0]=phiIp;
    k[0]=khIp;
    z[0]=tsIp->getReferencePoint()[2];

    auto sortByR = [&](TrackState* ts_a, TrackState* ts_b) {
        XYZVector posA = XYZVector(ts_a->getReferencePoint()[0], ts_a->getReferencePoint()[1], ts_a->getReferencePoint()[2]);
        XYZVector posB = XYZVector(ts_b->getReferencePoint()[0], ts_b->getReferencePoint()[1], ts_b->getReferencePoint()[2]);
        return posA.r() < posB.r();
    };
    sort(trackSt.begin(), trackSt.end(), sortByR);

    for (size_t j=1; j < trackSt.size(); ++j){
        HelixClass helix;
        double phi_helix = trackSt.at(j)->getPhi();
        double d0 = trackSt.at(j)->getD0();
        double z0 = trackSt.at(j)->getZ0();
        double omega = trackSt.at(j)->getOmega();
        double tanL = trackSt.at(j)->getTanLambda();
        helix.Initialize_Canonical(phi_helix, d0, z0, omega, tanL, _bField);

        double q = helix.getCharge();
        double rh = helix.getRadius();
        dd4hep::rec::Vector3D p = helix.getMomentum();
        double kh = 1./rh;
        XYZVector pos = XYZVector(trackSt.at(j)->getReferencePoint()[0], trackSt.at(j)->getReferencePoint()[1], trackSt.at(j)->getReferencePoint()[2]);

        // streamlog_out(DEBUG) << "TrackState at:"<<trackSt.at(j)->getLocation()<<" q:"<<q<<" Bz="<<B<<" Ku="<<CLHEP::c_light<<" z0="<<z0<<" p="<<p.r()<<" rh="<<rh<<" posR="<<pos.r()<<std::endl;
        phi[j]=phi_helix;
        k[j]=kh;
        z[j]=trackSt.at(j)->getReferencePoint()[2];

        trackLength += sqrt ( pow( (phi[j]-phi[j-1])/k[j-1],2) + pow(z[j]-z[j-1],2) ) ;
    }

    const TrackState* tsCalo= track->getTrackState(TrackState::AtCalorimeter);
    double phiCalo = tsCalo->getPhi();
    double omegaCalo = tsCalo->getOmega();
    double tanLCalo = tsCalo->getTanLambda();

    // streamlog_out(DEBUG)<< "Simplest one:"<<((phiIp - phiCalo)/omegaCalo)*sqrt(1. + tanLCalo*tanLCalo) << " Mitaroff: "<<trackLength<<std::endl;

    return trackLength;
}


double SETAnalysis::getTrackLengthRefit(Track* track){

    vector <TrackerHit*> trackHits = track->getTrackerHits();

    auto sortByRho = [](TrackerHit* a, TrackerHit* b) {
        XYZVector posA, posB;
        posA.SetCoordinates( a->getPosition() );
        posB.SetCoordinates( b->getPosition() );
        return posA.rho() < posB.rho();
    };
    sort(trackHits.begin(), trackHits.end(), sortByRho);

    // setup initial dummy covariance matrix
    vector<float> covMatrix(15);
    // initialize variances
    covMatrix[0]  = ( 1e+06 ); //sigma_d0^2
    covMatrix[2]  = ( 100. ); //sigma_phi0^2
    covMatrix[5]  = ( 0.00001 ); //sigma_omega^2
    covMatrix[9]  = ( 1e+06 ); //sigma_z0^2
    covMatrix[14] = ( 100. ); //sigma_tanl^2

    double maxChi2PerHit = 100.;
    IMarlinTrack* marlinTrk = _trkSystem->createTrack();
    TrackImpl refittedTrack;


    //Unfortunately this option doesn't work!!!
    // createFinalisedLCIOTrack(marlinTrk, trackHits, &refittedTrack, IMarlinTrack::backward, 0 , _bField, maxChi2PerHit);
    //Need to initialize manually at the last hit
    TrackState* tsInit = const_cast<EVENT::TrackState* > ( track->getTrackState( TrackState::AtLastHit ) ) ;
    TrackStateImpl preFit( *tsInit ) ;
    preFit.setCovMatrix( covMatrix )  ;
    createFinalisedLCIOTrack(marlinTrk, trackHits, &refittedTrack, IMarlinTrack::backward, &preFit , _bField, maxChi2PerHit);

    vector<pair<TrackerHit*, double> > hitsInFit;
    vector<pair<TrackerHit*, double> > outliers;

    marlinTrk->getHitsInFit(hitsInFit);
    marlinTrk->getOutliers(outliers);
    if(outliers.size() == 0) cout<< "No outliers... Good"<<endl;

    double dummyChi2;
    int dummyNdf;
    TrackStateImpl* ts1 = new TrackStateImpl();
    // ts1 = track->getTrackState(TrackState::AtIP);
    TrackerHit* firstHit = hitsInFit.back().first;
    marlinTrk->propagate(Vector3D(), firstHit, *ts1, dummyChi2, dummyNdf);
    double trackLength = 0.;
    for(int i=0; i < trackHits.size() ; i+=10 ){
        TrackStateImpl* ts2 = new TrackStateImpl();
        marlinTrk->getTrackState(trackHits[i], *ts2, dummyChi2, dummyNdf);

        double phi = ts1->getPhi();
        double phi2 = ts2->getPhi();
        double omega = ts1->getOmega();
        double tanL = ts1->getTanLambda();
        trackLength += std::abs( (phi2 - phi)/omega )*std::sqrt(1. + tanL*tanL);
        ts1 = ts2;
    }

    double phi = ts1->getPhi();
    double phi2 = track->getTrackState(TrackState::AtCalorimeter)->getPhi();
    double omega = ts1->getOmega();
    double tanL = ts1->getTanLambda();
    trackLength += std::abs( (phi2 - phi)/omega )* sqrt(1. + tanL*tanL);

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
