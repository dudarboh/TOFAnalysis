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
#include "DD4hep/DetElement.h"
#include "marlinutil/GeometryUtil.h"
#include "marlinutil/CalorimeterHitType.h"
#include "marlinutil/DDMarlinCED.h"
#include <UTIL/CellIDDecoder.h>
#include "UTIL/ILDConf.h"
#include "DD4hep/DetectorSelector.h"


#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "Math/Vector3D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TString.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGeoManager.h"
#include "TVirtualGeoTrack.h"
#include <TParticle.h>
#include "TPolyMarker3D.h"
#include "TGeoTrack.h"
#include "TGeoTube.h"
#include "TGeoHelix.h"
#include "TView.h"
#include "TSystem.h"
#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackImpl.h>

#include "streamlog/streamlog.h"
#include <ostream>

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
    DDMarlinCED::init(this) ;


    std::cout.precision(7);
    _nEvent = 0;
    _file.reset( new TFile(_outputFileName.c_str(), "RECREATE") );
    _tree.reset( new TTree("SETAnalysis", "SETAnalysis") );

    _tree->Branch("pdg", &_pdg);
    _tree->Branch("has_set_hit", &_hasSetHit);
    _tree->Branch("n_ecal_hits", &_nEcalHits);
    _tree->Branch("n_fit_hits", &_nFitHits);

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
    _tree->Branch("track_length_ecal", &_trackLength["ecal"]);
    _tree->Branch("phi_curl_set", &_phiCurl["set"]);
    _tree->Branch("phi_curl_ecal", &_phiCurl["ecal"]);

    _tree->Branch("pos_set_hit", &_posSetHit);
    _tree->Branch("pos_closest", &_posClosest);

    _tree->Branch("mom_hm_set", &_mom["hmSet"]);
    _tree->Branch("mom_hm_ecal", &_mom["hmEcal"]);

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
        //calculate track length
        std::vector <TrackerHit*> trackHits = track->getTrackerHits();
        //It is probably required by the fitter...
        auto sortByRho = [](TrackerHit* a, TrackerHit* b) {
            dd4hep::rec::Vector3D posA( a->getPosition() );
            dd4hep::rec::Vector3D posB( b->getPosition() );
            return posA.rho() < posB.rho();
        };
        std::sort(trackHits.begin(), trackHits.end(), sortByRho);

        auto estimateArcLength = [](double phi1, double phi2, double omega, double z1, double z2) {
            double dPhi = std::abs(phi2-phi1);
            if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
            return std::sqrt( std::pow(dPhi/omega, 2) + std::pow(z2-z1, 2) );
        };

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
        _nFitHits = hitsInFit.size();

        vector<TrackStateImpl> trackStatesPerHit;
        trackStatesPerHit.push_back( *(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter))));
        //add track states at SET hit. If it doesn't exist, then it copies the most outer tpc hit (shouldn't affect track length)
        trackStatesPerHit.push_back( *(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );

        // loop over all hits in the order of how they were fitted (fitted backwards and hits sorted by rho).
        int nHits = hitsInFit.size();
        for(int j=0; j < nHits; ++j ){
            TrackStateImpl ts;
            double chi2Tmp;
            int ndfTmp;
            marlinTrk->getTrackState(hitsInFit[j].first, ts, chi2Tmp, ndfTmp);
            trackStatesPerHit.push_back(ts);
        }
        trackStatesPerHit.push_back(*(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));

        vector<double> phi, omega, z, arcLength, arcCurl, p, d0, z0, tanL, pWeighted;
        for(size_t j=0; j < trackStatesPerHit.size(); ++j ){
            phi.push_back( trackStatesPerHit[j].getPhi() );
            d0.push_back( trackStatesPerHit[j].getD0() );
            z0.push_back( trackStatesPerHit[j].getZ0() );
            omega.push_back( trackStatesPerHit[j].getOmega() );
            tanL.push_back( trackStatesPerHit[j].getTanLambda() );
            z.push_back(trackStatesPerHit[j].getReferencePoint()[2] + trackStatesPerHit[j].getZ0() );
            HelixClass helix;
            helix.Initialize_Canonical(phi[j], d0[j], z0[j], omega[j], tanL[j], _bField);
            p.push_back( dd4hep::rec::Vector3D( helix.getMomentum() ).r() );
            //can't calculate track length from 1 element
            if (j == 0) continue;
            arcLength.push_back( estimateArcLength( phi[j-1], phi[j], omega[j-1], z[j-1], z[j] ) );

            double deltaPhi = std::abs(phi[j] - phi[j-1]);
            if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;
            arcCurl.push_back( deltaPhi );

            pWeighted.push_back( arcLength[j] /(p[j]*p[j]) );
        }

        _trackLength["set"] = std::accumulate(arcLength.begin() + 1, arcLength.end(), 0.);
        _trackLength["ecal"] = std::accumulate(arcLength.begin(), arcLength.end(), 0.);

        _phiCurl["set"] = std::accumulate(arcCurl.begin() + 1, arcCurl.end(), 0.);
        _phiCurl["ecal"] = std::accumulate(arcCurl.begin(), arcCurl.end(), 0.);

        _mom["hmSet"] = std::sqrt(_trackLength["set"] / std::accumulate(pWeighted.begin() + 1, pWeighted.end(), 0.) );
        _mom["hmEcal"] = std::sqrt(_trackLength["ecal"] / std::accumulate(pWeighted.begin(), pWeighted.end(), 0.) );

        // if (_trackLength["ecal"] - _trackLength["set"] < 0 || _trackLength["ecal"] - _trackLength["set"] > 100.){
            if (_trackLength["ecal"] > 7000.){
            DDMarlinCED::newEvent(this);
            DDMarlinCED::drawDD4hepDetector(_theDetector, 0, std::vector<std::string>{"SIT", "SET"});
            DDCEDPickingHandler& pHandler=DDCEDPickingHandler::getInstance();
            pHandler.update(evt);
            drawPfo(track, cluster);

            cout<<"_trackLength_set="<<_trackLength["set"]<<endl;
            cout<<"_trackLength_ecal="<<_trackLength["ecal"]<<endl;
            cout<<"_phiCurl_set="<<_phiCurl["set"]<<endl;
            cout<<"_phiCurl_ecal="<<_phiCurl["ecal"]<<endl;
            cout<<"ptAtIP="<<_tsMom["ip"].rho()<<endl;

            DDMarlinCED::draw(this, 1);

        }
        // if (_phiCurl["ecal"] > M_PI) drawTrack(track);

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


void SETAnalysis::drawPfo(Track* track, Cluster* cluster){
    int nSubTracks = track->getTracks().size();
    streamlog_out( DEBUG6 ) << " -- Final LCIO Track has "<<nSubTracks<< " subtracks - will use these for displaying hits "<< std::endl ;
    if ( nSubTracks == 0 ) return;

    // std::copy( track->getTrackerHits().begin() , track->getTrackerHits().end() , std::back_inserter(  hits ) ) ;
    cout<<"Track hits size:"<<track->getTrackerHits().size()<<endl;
    cout<<"N subdetectors:"<<track->getSubdetectorHitNumbers().size()<<endl;
    cout<<"VXD used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::VXD)*2-2]<<endl;
    cout<<"VXD not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::VXD)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::VXD)*2-2]<<endl;
    cout<<"SIT used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::SIT)*2-2]<<endl;
    cout<<"SIT not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::SIT)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::SIT)*2-2]<<endl;
    cout<<"FTD used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::FTD)*2-2]<<endl;
    cout<<"FTD not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::FTD)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::FTD)*2-2]<<endl;
    cout<<"TPC used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::TPC)*2-2]<<endl;
    cout<<"TPC not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::TPC)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::TPC)*2-2]<<endl;
    cout<<"SET used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::SET)*2-2]<<endl;
    cout<<"SET not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::SET)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::SET)*2-2]<<endl;

    std::vector<int> hitsColors = {0x1700ff, 0x000fff, 0x0032ff, 0x005dff, 0x0080ff};
    for(int i=0; i < nSubTracks; ++i){
        // std::copy( subTrack->getTrackerHits().begin() , subTrack->getTrackerHits().end() , std::back_inserter(  hits ) ) ;
        const Track* subTrack = track->getTracks()[i];
        int nHits = subTrack->getTrackerHits().size();
        streamlog_out( DEBUG6 ) << " -- subTrack i= "<< i << " has " <<  nHits << " hits  "<< std::endl ;
        for (int j = 0; j < nHits; ++j) {
            TrackerHit* hit = subTrack->getTrackerHits()[j];
            float x = hit->getPosition()[0];
            float y = hit->getPosition()[1];
            float z = hit->getPosition()[2];
            ced_hit(x, y, z, 0, 3, hitsColors[i]);
        }

    }

    for ( const auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        XYZVectorF hitPos;
        hitPos.SetCoordinates( hit->getPosition() );
        ced_hit(hitPos.x(), hitPos.y(), hitPos.z(), 34, 7, 0x08ff00);
    }



    const TrackState* tsIP = track->getTrackState(TrackState::AtIP);
    ced_hit(tsIP->getReferencePoint()[0], tsIP->getReferencePoint()[1], tsIP->getReferencePoint()[2], 20, 12, 0x00e8ff);
    const TrackState* tsFirst = track->getTrackState(TrackState::AtFirstHit);
    ced_hit(tsFirst->getReferencePoint()[0], tsFirst->getReferencePoint()[1], tsFirst->getReferencePoint()[2], 20, 12, 0x00e8ff);
    const TrackState* tsLast = track->getTrackState(TrackState::AtLastHit);
    ced_hit(tsLast->getReferencePoint()[0], tsLast->getReferencePoint()[1], tsLast->getReferencePoint()[2], 20, 12, 0x00e8ff);
    const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
    ced_hit(tsEcal->getReferencePoint()[0], tsEcal->getReferencePoint()[1], tsEcal->getReferencePoint()[2], 20, 12, 0x00e8ff);


    double pt = _bField * 3e-4 / std::abs( tsIP->getOmega() );
    double charge = ( tsIP->getOmega() > 0. ?  1. : -1. );
    double Px = pt * std::cos(  tsIP->getPhi() ) ;
    double Py = pt * std::sin(  tsIP->getPhi() ) ;
    double Pz = pt * tsIP->getTanLambda() ;
    double Xs = tsIP->getReferencePoint()[0] -  tsIP->getD0() * std::sin( tsIP->getPhi() ) ;
    double Ys = tsIP->getReferencePoint()[1] +  tsIP->getD0() * std::cos( tsIP->getPhi() ) ;
    double Zs = tsIP->getReferencePoint()[2] +  tsIP->getZ0() ;
    //helix at IP
    DDMarlinCED::drawHelix(_bField, charge, Xs, Ys, Zs, Px, Py, Pz, 1, 1,
                         0xff0000, 0., 2100., 3000., 0);


    pt = _bField * 3e-4 / std::abs( tsEcal->getOmega() );
    charge = ( tsEcal->getOmega() > 0. ?  1. : -1. );
    Px = pt * std::cos(  tsEcal->getPhi() ) ;
    Py = pt * std::sin(  tsEcal->getPhi() ) ;
    Pz = pt * tsEcal->getTanLambda() ;
    Xs = tsEcal->getReferencePoint()[0] -  tsEcal->getD0() * std::sin( tsEcal->getPhi() ) ;
    Ys = tsEcal->getReferencePoint()[1] +  tsEcal->getD0() * std::cos( tsEcal->getPhi() ) ;
    Zs = tsEcal->getReferencePoint()[2] +  tsEcal->getZ0() ;
    //helix at ECal
    DDMarlinCED::drawHelix(-_bField, -charge, Xs, Ys, Zs, Px, Py, Pz, 1, 1,
                         0xffb200, 0., 2100., 3000., 0);


}
