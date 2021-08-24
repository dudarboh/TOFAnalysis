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
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackImpl.h>

#include "streamlog/streamlog.h"

#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>


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
    std::cout.precision(7);
    _nEvent = 0;
    _file.reset( new TFile(_outputFileName.c_str(), "RECREATE") );
    _tree.reset( new TTree("SETAnalysis", "SETAnalysis") );

    _tree->Branch("pdg", &_pdg);
    _tree->Branch("has_set_hit", &_hasSetHit);
    _tree->Branch("n_ecal_hits", &_nEcalHits);

    std::vector<std::string> tsNames{"ip", "first", "last", "ecal"};
    for (auto ts : tsNames){
        _tree->Branch( Form("ts_%s_pos", ts.c_str()), &_tsPos[ts] );
        _tree->Branch( Form("ts_%s_mom", ts.c_str()), &_tsMom[ts] );
        _tree->Branch( Form("ts_%s_omega", ts.c_str()), &_tsOmega[ts] );
        _tree->Branch( Form("ts_%s_tanL", ts.c_str()), &_tsTanL[ts] );
        _tree->Branch( Form("ts_%s_phi", ts.c_str()), &_tsPhi[ts] );
        _tree->Branch( Form("ts_%s_d0", ts.c_str()), &_tsD0[ts] );
        _tree->Branch( Form("ts_%s_z0", ts.c_str()), &_tsZ0[ts] );
    }

    _tree->Branch("track_length_set", &_trackLength["set"]);
    _tree->Branch("track_length_refit_set_z", &_trackLength["setRefitZ"]);
    _tree->Branch("track_length_ip", &_trackLength["ip"]);
    _tree->Branch("track_length_calo", &_trackLength["calo"]);
    _tree->Branch("track_length_refit_tanL", &_trackLength["refitTanL"]);
    _tree->Branch("track_length_refit_z", &_trackLength["refitZ"]);

    _tree->Branch("pos_set_hit", &_posSetHit);
    _tree->Branch("pos_closest", &_posClosest);

    _tree->Branch("mom2_hm_TanL", &_mom["sqrHmTanL"]);
    _tree->Branch("mom2_hm_tanL_set", &_mom["sqrHmTanLSet"]);
    _tree->Branch("mom2_hm_dz", &_mom["sqrHmDz"]);
    _tree->Branch("mom2_hm_dz_set", &_mom["sqrHmDzSet"]);

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
    streamlog_out(MESSAGE)<<"****** Event: "<<_nEvent<<endl;

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
        // if (!_hasSetHit) continue;

        ///////////////////////////////////
        _pdg = mc->getPDG();

        if (_hasSetHit) _posSetHit.SetCoordinates( setHit->getPosition() );


        // if (_hasSetHit)
        const LCObjectVec& setSimHits = spPointToSimHit.getRelatedToObjects( setHit );

        /////////////////WRITE TRACK STATES/////////////////////////
        // should return at Ip, at First, at Last, at Calo
        const vector<TrackState*> trackStates = track->getTrackStates();
        std::vector<std::string> tsNames{"ip", "first", "last", "ecal"};

        for(unsigned int j=0; j < trackStates.size(); ++j){
            _tsOmega[ tsNames[j] ] = trackStates[j]->getOmega();
            _tsTanL[ tsNames[j] ] = trackStates[j]->getTanLambda();
            _tsPhi[ tsNames[j] ] = trackStates[j]->getPhi();
            _tsD0[ tsNames[j] ] = trackStates[j]->getD0();
            _tsZ0[ tsNames[j] ] = trackStates[j]->getZ0();
            _tsPos[ tsNames[j] ].SetCoordinates( trackStates[j]->getReferencePoint() );
            HelixClass helix;
            helix.Initialize_Canonical(_tsPhi[tsNames[j]], _tsD0[tsNames[j]], _tsZ0[tsNames[j]], _tsOmega[tsNames[j]], _tsTanL[tsNames[j]], _bField);
            _tsMom[ tsNames[j] ].SetCoordinates( helix.getMomentum() );
        }

        ///////////////////////WRITE TRACK LENGTHS/////////////////////
        _trackLength["ip"] = getTrackLength(track, true, true, "TDR");
        _trackLength["calo"] = getTrackLength(track, true, true, "helixEcal");
        _trackLength["set"] = getTrackLength(track, true, false, "helixEcal");
        _trackLength["setRefitZ"] = getTrackLength(track, true, false, "dz");
        _trackLength["refitTanL"] = getTrackLength(track, true, true, "tanL");
        _trackLength["refitZ"] = getTrackLength(track, true, true, "dz");

        _mom["sqrHmTanL"] = getMom2Harmonic( track, true, true, "tanL" );
        _mom["sqrHmTanLSet"] = getMom2Harmonic( track, true, false, "tanL" );
        _mom["sqrHmDz"] = getMom2Harmonic( track, true, true, "dz" );
        _mom["sqrHmDzSet"] = getMom2Harmonic( track, true, false, "dz" );


        streamlog_out(DEBUG)<<"Refit: "<<_trackLength["refitTanL"]<<endl;
        streamlog_out(DEBUG)<<"Winni: "<<_trackLength["refitZ"]<<endl;
        streamlog_out(DEBUG)<<"Diff: "<<_trackLength["refitZ"]-_trackLength["refitTanL"]<<endl;

        SimTrackerHit* setSimHitFront = nullptr;
        SimTrackerHit* setSimHitBack = nullptr;
        if (_hasSetHit){
            setSimHitFront = dynamic_cast <SimTrackerHit*> (setSimHits.at(0) );
            setSimHitBack = dynamic_cast <SimTrackerHit*> (setSimHits.at(1) );
        }

        CalorimeterHit* closestHit = getClosestHit( cluster, _tsPos["ecal"] );
        _posClosest.SetCoordinates( closestHit->getPosition() );

        for(unsigned int j=0; j < std::size(_smearings); ++j ){
            pair<XYZVectorF, double> fastestHit = getFastestHit( cluster, _smearings[j] / 1000. );
            _posFastest[j] = fastestHit.first;
            _tofFastest[j] = fastestHit.second - ( _posFastest[j] - _tsPos["ecal"] ).r()/CLHEP::c_light;

            _tofClosest[j] = CLHEP::RandGauss::shoot( closestHit->getTime(), _smearings[j] / 1000. ) - ( _posClosest - _tsPos["ecal"] ).r()/CLHEP::c_light;
            if (_hasSetHit){
                _tofSetFront[j] = CLHEP::RandGauss::shoot( setSimHitFront->getTime(), _smearings[j] / 1000. );
                _tofSetBack[j] =  CLHEP::RandGauss::shoot( setSimHitBack->getTime(), _smearings[j] / 1000. );
            }
            else{
                _tofSetFront[j] = 0.;
                _tofSetBack[j] = 0.;
            }
            _tofFrankFit[j] = getTofFrankFit( cluster, _tsPos["ecal"], _tsMom["ecal"], _smearings[j] / 1000. );
            _tofFrankAvg[j] = getTofFrankAvg( cluster, _tsPos["ecal"], _tsMom["ecal"], _smearings[j] / 1000. );
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


double SETAnalysis::estimateTrackLengthTanL(double phi1, double phi2, double omega, double tanL) {
    bool phiFlipped = std::abs(phi2-phi1) > M_PI;
    double dPhi;
    if (phiFlipped) dPhi = 2*M_PI - std::abs(phi2 - phi1);
    else dPhi = std::abs(phi2 - phi1);

    return dPhi/std::abs(omega) * std::sqrt(1.+tanL*tanL);
}


double SETAnalysis::estimateTrackLengthZ(double phi1, double phi2, double omega, double z1, double z2) {
    bool phiFlipped = std::abs(phi2-phi1) > M_PI;
    double dPhi;
    if (phiFlipped) dPhi = 2*M_PI - std::abs(phi2 - phi1);
    else dPhi = std::abs(phi2 - phi1);

    return std::sqrt( std::pow(dPhi/omega, 2) + std::pow(z2-z1, 2) );
}


double SETAnalysis::getTrackLength(Track* track, bool extrapolateToIp, bool extrapolateToEcal, std::string method){
    streamlog_out(DEBUG)<<"____________________________________________________________________________________________________________________________________"<<endl;

    if (method == "TDR"){
        //OBSOLETE, don't use
        // This is TDR production release from 2018 for improvement comparisons since I joined...
        // I literally copied from TofUtils.. Just renamed variables
        const TrackState* tsIp = track->getTrackState( TrackState::AtIP ) ;
        float phiIp = tsIp->getPhi() ;
        float omega = tsIp->getOmega()  ;
        float tanL = tsIp->getTanLambda() ;

        float phiEcal;
        if (extrapolateToEcal) {
            const TrackState* tsEcal = track->getTrackState( TrackState::AtCalorimeter ) ;
            phiEcal = tsEcal->getPhi();
        }
        else{
            const TrackState* tsLast= track->getTrackState(TrackState::AtLastHit);
            phiEcal = tsLast->getPhi();
        }

        float length = (phiIp-phiEcal)*(1/omega) * sqrt( 1 + tanL * tanL ) ;
        return length;
    }
    else if (method == "helixEcal"){
        //don't refit hit-by-hit. Used before. Use trackState at Calo for curvature
        const TrackState* ts1 = track->getTrackState(TrackState::AtIP);
        double phi1 = ts1->getPhi();
        const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
        double omega = tsEcal->getOmega();
        double tanL = tsEcal->getTanLambda();

        double phi2;
        if (extrapolateToEcal) phi2 = tsEcal->getPhi();
        else{
            const TrackState* tsLast= track->getTrackState(TrackState::AtLastHit);
            phi2 = tsLast->getPhi();
        }
        return std::abs( (phi1 - phi2)/omega )*std::sqrt(1. + tanL*tanL);
    }

    std::vector <TrackerHit*> trackHits = track->getTrackerHits();

    //It is probably required by the fitter...
    auto sortByR = [](TrackerHit* a, TrackerHit* b) {
        XYZVector posA, posB;
        posA.SetCoordinates( a->getPosition() );
        posB.SetCoordinates( b->getPosition() );
        return posA.rho() < posB.rho();
    };
    std::sort(trackHits.begin(), trackHits.end(), sortByR);

    // setup initial dummy covariance matrix
    vector<float> covMatrix(15);
    // initialize variances
    covMatrix[0]  = 1e+06; //sigma_d0^2
    covMatrix[2]  = 100.; //sigma_phi0^2
    covMatrix[5]  = 0.00001; //sigma_omega^2
    covMatrix[9]  = 1e+06; //sigma_z0^2
    covMatrix[14] = 100.; //sigma_tanl^2
    double maxChi2PerHit = 100.;
    IMarlinTrack* marlinTrk = _trkSystem->createTrack();
    TrackImpl refittedTrack;

    //Need to initialize trackState at last hit
    TrackStateImpl preFit = *(track->getTrackState( TrackState::AtLastHit ));
    preFit.setCovMatrix( covMatrix );
    createFinalisedLCIOTrack(marlinTrk, trackHits, &refittedTrack, IMarlinTrack::backward, &preFit , _bField, maxChi2PerHit);
    //fit is finished, collect the hits

    vector<pair<TrackerHit*, double> > hitsInFit;
    vector<pair<TrackerHit*, double> > outliers;
    marlinTrk->getHitsInFit(hitsInFit);
    marlinTrk->getOutliers(outliers);

    vector<TrackStateImpl> trackStates;

    if (extrapolateToEcal) trackStates.push_back( *(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter))));
    //add track states at SET hit. If it doesn't exist, then it copies the most outer tpc hit (shouldn't affect track length)
    trackStates.push_back( *(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );


    // loop over all hits in the order of how they were fitted (fitted backwards and hits sorted by rho).
    int nHits = hitsInFit.size();
    for(int i=0; i < nHits; ++i ){
        TrackStateImpl ts;
        double chi2Tmp;
        int ndfTmp;
        marlinTrk->getTrackState(hitsInFit[i].first, ts, chi2Tmp, ndfTmp);
        trackStates.push_back(ts);
    }
    if (extrapolateToIp) trackStates.push_back(*(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));

    //track parameters
    vector<double> phi, d0, z0, omega, tanL, refx, refy, refz, p, pt, pz, z, trackLength;

    for(auto ts:trackStates){
        phi.push_back( ts.getPhi() );
        d0.push_back( ts.getD0() );
        z0.push_back( ts.getZ0() );
        omega.push_back( ts.getOmega() );
        tanL.push_back( ts.getTanLambda() );
        refx.push_back(ts.getReferencePoint()[0]);
        refy.push_back(ts.getReferencePoint()[1]);
        refz.push_back(ts.getReferencePoint()[2]);
        z.push_back(ts.getReferencePoint()[2] + ts.getZ0() );
        HelixClass helix;
        helix.Initialize_Canonical(phi.back(), d0.back(), z0.back(), omega.back(), tanL.back(), _bField);
        p.push_back( dd4hep::rec::Vector3D( helix.getMomentum() ).r() );
        pt.push_back( dd4hep::rec::Vector3D( helix.getMomentum() ).trans() );
        pz.push_back( dd4hep::rec::Vector3D( helix.getMomentum() ).z() );

        if (phi.size() == 1) continue; //can't calculate track length from 1 element
        if (method == "tanL") trackLength.push_back( estimateTrackLengthTanL( phi.rbegin()[1], phi.back(), omega.rbegin()[1], tanL.rbegin()[1] ) );
        else if (method == "dz") trackLength.push_back( estimateTrackLengthZ( phi.rbegin()[1], phi.back(), omega.rbegin()[1], z.rbegin()[1], z.back() ) );
        else throw std::string("Uncknown method, available are: TDR, helixEcal, tanL, dz");
    }
    double totalTrackLength = std::accumulate(trackLength.begin(), trackLength.end(), 0.);

    streamlog_out(DEBUG)<<"nHits="<<nHits<<"     "<<endl;
    for(auto it= trackStates.rbegin(); it != trackStates.rend(); ++it){
        bool idxIp = it == trackStates.rbegin();
        bool idxLast = it == trackStates.rend() - 2;
        bool idxEcal = it == trackStates.rend() - 1;
        int idx = std::distance( trackStates.begin(), it.base() ) - 1;
        if ( !(idxIp || idxLast || idxEcal ||idx % 20 == 0) || idx == 0 ) continue;
        streamlog_out(DEBUG)<<"****************************"<<endl;
        if ( idxIp ) {streamlog_out(DEBUG)<<"atIp:   "<<endl;}
        else if ( idxLast ) {streamlog_out(DEBUG)<<"atLast:   "<<endl;}
        else if ( idxEcal ) {streamlog_out(DEBUG)<<"atEcal:   "<<endl;}
        for(int i=0; i<2; ++i){
            idx -= i;
            streamlog_out(DEBUG)<<"phi[deg]="<< phi.at(idx)*180./M_PI<<"    " ;
            streamlog_out(DEBUG)<<"d0[mm]="<< d0.at(idx) <<"    " ;
            streamlog_out(DEBUG)<<"z0[mm]="<< z0.at(idx) <<"    " ;
            streamlog_out(DEBUG)<<"R[mm]="<< 1. / omega.at(idx) <<"    " ;
            streamlog_out(DEBUG)<<"theta[deg]="<< std::atan( 1./tanL.at(idx) )*180./M_PI <<"    " <<endl;
            streamlog_out(DEBUG)<<"ref[mm]=("<<refx.at(idx)<<"  "<<refy.at(idx)<<"  "<<refz.at(idx)<<")    " ;
            streamlog_out(DEBUG)<<"z[mm]="<< z.at(idx) <<"    " ;
            streamlog_out(DEBUG)<<"p[GeV]="<< p.at(idx) <<"    " ;
            streamlog_out(DEBUG)<<"pt[GeV]="<< pt.at(idx) <<"    " ;
            streamlog_out(DEBUG)<<"pz[GeV]="<< pz.at(idx) <<"    "<<endl;
        }
        streamlog_out(DEBUG)<<"Track length="<<trackLength[idx - 1]<<endl;
        streamlog_out(DEBUG)<<endl;
    }

    delete marlinTrk;

    return totalTrackLength;
}


double SETAnalysis::getMom2Harmonic(Track* track, bool extrapolateToIp, bool extrapolateToEcal, std::string method){
    //Return harmonic mean of squared momentum from Winni's paper

    std::vector <TrackerHit*> trackHits = track->getTrackerHits();

    //It is probably required by the fitter...
    auto sortByR = [](TrackerHit* a, TrackerHit* b) {
        XYZVector posA, posB;
        posA.SetCoordinates( a->getPosition() );
        posB.SetCoordinates( b->getPosition() );
        return posA.rho() < posB.rho();
    };
    std::sort(trackHits.begin(), trackHits.end(), sortByR);

    // setup initial dummy covariance matrix
    vector<float> covMatrix(15);
    // initialize variances
    covMatrix[0]  = 1e+06; //sigma_d0^2
    covMatrix[2]  = 100.; //sigma_phi0^2
    covMatrix[5]  = 0.00001; //sigma_omega^2
    covMatrix[9]  = 1e+06; //sigma_z0^2
    covMatrix[14] = 100.; //sigma_tanl^2
    double maxChi2PerHit = 100.;
    IMarlinTrack* marlinTrk = _trkSystem->createTrack();
    TrackImpl refittedTrack;

    //Need to initialize trackState at last hit
    TrackStateImpl preFit = *(track->getTrackState( TrackState::AtLastHit ));
    preFit.setCovMatrix( covMatrix );
    createFinalisedLCIOTrack(marlinTrk, trackHits, &refittedTrack, IMarlinTrack::backward, &preFit , _bField, maxChi2PerHit);
    //fit is finished, collect the hits

    vector<pair<TrackerHit*, double> > hitsInFit;
    marlinTrk->getHitsInFit(hitsInFit);

    vector<TrackStateImpl> trackStates;

    if (extrapolateToEcal) trackStates.push_back( *(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter))));
    //add track states at SET hit. If it doesn't exist, then it copies the most outer tpc hit (shouldn't affect track length)
    trackStates.push_back( *(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );


    // loop over all hits in the order of how they were fitted (fitted backwards and hits sorted by rho).
    int nHits = hitsInFit.size();
    for(int i=0; i < nHits; ++i ){
        TrackStateImpl ts;
        double chi2Tmp;
        int ndfTmp;
        marlinTrk->getTrackState(hitsInFit[i].first, ts, chi2Tmp, ndfTmp);
        trackStates.push_back(ts);
    }
    if (extrapolateToIp) trackStates.push_back(*(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));

    //track parameters
    std::vector<double> phi, d0, z0, omega, tanL, refx, refy, refz, z, p, pt, pz, trackLength, pWeighted;

    for(auto ts:trackStates){
        phi.push_back( ts.getPhi() );
        d0.push_back( ts.getD0() );
        z0.push_back( ts.getZ0() );
        omega.push_back( ts.getOmega() );
        tanL.push_back( ts.getTanLambda() );
        refx.push_back(ts.getReferencePoint()[0]);
        refy.push_back(ts.getReferencePoint()[1]);
        refz.push_back(ts.getReferencePoint()[2]);
        z.push_back(ts.getReferencePoint()[2] + ts.getZ0() );
        HelixClass helix;
        helix.Initialize_Canonical(phi.back(), d0.back(), z0.back(), omega.back(), tanL.back(), _bField);
        p.push_back( dd4hep::rec::Vector3D( helix.getMomentum() ).r() );
        pt.push_back( dd4hep::rec::Vector3D( helix.getMomentum() ).trans() );
        pz.push_back( dd4hep::rec::Vector3D( helix.getMomentum() ).z() );


        if (phi.size() == 1) continue; //can't calculate track length from 1 element
        if (method == "tanL") trackLength.push_back( estimateTrackLengthTanL( phi.rbegin()[1], phi.back(), omega.rbegin()[1], tanL.rbegin()[1] ) );
        else if (method == "dz") trackLength.push_back( estimateTrackLengthZ( phi.rbegin()[1], phi.back(), omega.rbegin()[1], z.rbegin()[1], z.back() ) );
        else throw std::string("Uncknown method, available are: tanL, dz");
        pWeighted.push_back( trackLength.back() /(p.back()*p.back()) );
    }
    double totalTrackLength = std::accumulate(trackLength.begin(), trackLength.end(), 0.);
    double totalWeight = std::accumulate(pWeighted.begin(), pWeighted.end(), 0.);
    double mom2 = totalTrackLength / totalWeight;

    delete marlinTrk;

    return mom2;
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
