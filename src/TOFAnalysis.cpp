#include "TOFAnalysis.hpp"
#include "TOFAnaUtils.hpp"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"
#include "HelixClass.h"
#include "IMPL/TrackStateImpl.h"
#include "IMPL/TrackImpl.h"

#include "marlinutil/DDMarlinCED.h"


using namespace TOFAnaUtils;

TOFAnalysis aTOFAnalysis;

TOFAnalysis::TOFAnalysis() : Processor("TOFAnalysis"){
    registerProcessorParameter(string("outputFile"),
                               string("Name of the output root file"),
                               _outputFileName,
                               string("TOFAnalysis_RENAME.root"));
}

void TOFAnalysis::init(){
    DDMarlinCED::init(this);

    cout.precision(7);
    _nEvent = 0;
    _file.reset( new TFile(_outputFileName.c_str(), "RECREATE") );
    _tree.reset( new TTree("TOFAnalysis", "TOFAnalysis") );

    _tree->Branch("pdg", &_pdg);
    _tree->Branch("has_set_hit", &_hasSetHit);
    _tree->Branch("n_ecal_hits", &_nEcalHits);
    _tree->Branch("n_fit_hits", &_nFitHits);

    vector<string> tsNames{"ip", "first", "last", "ecal"};
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
    _tree->Branch("n_curls_set", &_nCurls["set"]);
    _tree->Branch("n_curls_ecal", &_nCurls["ecal"]);

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
    _tpcROuter = TOFAnaUtils::getTpcR().second;
    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem( "DDKalTest", nullptr, "" ) ;

    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init() ;

}


void TOFAnalysis::processEvent(LCEvent* evt){
    ++_nEvent;
    streamlog_out(MESSAGE)<<"******Event****** "<<_nEvent<<endl;

    // set the correct configuration for the tracking system for this event
    MarlinTrk::TrkSysConfig< MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson( _trkSystem, true );
    MarlinTrk::TrkSysConfig< MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson( _trkSystem, true);
    MarlinTrk::TrkSysConfig< MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon( _trkSystem, true);
    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator pfoToMc( evt->getCollection("RecoMCTruthLink") );

    LCRelationNavigator spPointToSimHit( evt->getCollection("SETSpacePointRelations") );

    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        /////////////////////////////////// Load all info from PFO
        ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

        MCParticle* mc = TOFAnaUtils::getMcMaxWeight(pfoToMc, pfo);
        if(mc == nullptr) continue;

        const vector<Cluster*>& clusters = pfo->getClusters();
        if(clusters.size() != 1) continue;
        Cluster* cluster = clusters.at(0);
        _nEcalHits = TOFAnaUtils::getNEcalHits(cluster);
        if (_nEcalHits == 0) continue;

        const vector<Track*>& tracks = pfo->getTracks();
        if(tracks.size() != 1) continue;
        Track* track = tracks.at(0);
        _nFitHits = 0;

        TrackerHit* setHit = TOFAnaUtils::getSetHit(track, _tpcROuter);
        _hasSetHit = (setHit != nullptr);
        // if (!_hasSetHit) continue;

        ///////////////////////////////////
        _pdg = mc->getPDG();

        if (_hasSetHit) _posSetHit.SetCoordinates( setHit->getPosition() );

        // if (_hasSetHit)
        const LCObjectVec& setSimHits = spPointToSimHit.getRelatedToObjects( setHit );

        /////////////////WRITE TRACK STATES/////////////////////////
        // should return at Ip, at First, at Last, at Calo
        const vector<TrackState*>& trackStates = track->getTrackStates();
        vector<string> tsNames{"ip", "first", "last", "ecal"};

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
        auto sortByRho = [](TrackerHit* a, TrackerHit* b) {
            XYZVector posA, posB;
            posA.SetCoordinates( a->getPosition() );
            posB.SetCoordinates( b->getPosition() );
            return posA.rho() < posB.rho();
        };

        auto getArcLength = [](double phi1, double phi2, double omega, double z1, double z2) {
            double dPhi = std::abs(phi2-phi1);
            if (dPhi > M_PI) dPhi = 2*M_PI - dPhi;
            return std::sqrt( std::pow(dPhi/omega, 2) + std::pow(z2-z1, 2) );
        };

        //This is how it should look like!!!!!!!!!
        // getTracksToFit();
        vector<Track*> tracksToFit;
        tracksToFit.push_back(track);
        int nSubTracks = track->getTracks().size();
        if (nSubTracks != 1){
            int nHitsTot = track->getTrackerHits().size();
            int nHits0 = track->getTracks()[0]->getTrackerHits().size();
            int nHits1 = track->getTracks()[1]->getTrackerHits().size();
            if( ((nHitsTot == nHits0 + nHits1) &&  !_hasSetHit) || ((nHitsTot == nHits0 + nHits1 + 1) && _hasSetHit )  ){
                for(int j=2; j < nSubTracks; ++j) tracksToFit.push_back( track->getTracks()[j] );
            }
            else{
                // This means VXD+SIT subTrack is not there (bug?!)! So we need to add tracks from i=1
                for(int j=1; j < nSubTracks; ++j) tracksToFit.push_back( track->getTracks()[j] );
            }
        }
        // getTrackStatesPerHit();
        vector<TrackStateImpl> trackStatesPerHit;
        for(unsigned int j=0; j < tracksToFit.size(); ++j){
            vector <TrackerHit*> trackHits = tracksToFit[j]->getTrackerHits();
            sort(trackHits.begin(), trackHits.end(), sortByRho);
            // setup initial dummy covariance matrix
            vector<float> covMatrix(15);
            // initialize variances
            covMatrix[0]  = 1e+06; //sigma_d0^2
            covMatrix[2]  = 100.; //sigma_phi0^2
            covMatrix[5]  = 0.00001; //sigma_omega^2
            covMatrix[9]  = 1e+06; //sigma_z0^2
            covMatrix[14] = 100.; //sigma_tanl^2
            double maxChi2PerHit = 100.;
            MarlinTrk::IMarlinTrack* marlinTrk = _trkSystem->createTrack();
            TrackImpl refittedTrack;

            //Need to initialize trackState at last hit
            TrackStateImpl preFit = *(tracksToFit[j]->getTrackState(TrackState::AtLastHit));
            preFit.setCovMatrix( covMatrix );
            int errorFit = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trackHits, &refittedTrack, MarlinTrk::IMarlinTrack::backward, &preFit , _bField, maxChi2PerHit);
            if (errorFit != 0) continue;

            if (j == 0) trackStatesPerHit.push_back(*(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtIP)) ));

            //fit is finished, collect the hits
            vector<pair<TrackerHit*, double> > hitsInFit;
            marlinTrk->getHitsInFit(hitsInFit);
            _nFitHits += hitsInFit.size();
            if (j == 0){
                //do reverse, so it increases in rho for the FIRST TRACK!!!!!
                for( int k = hitsInFit.size() - 1; k >= 0 ; --k ){
                    TrackStateImpl ts;
                    double chi2Tmp;
                    int ndfTmp;
                    marlinTrk->getTrackState(hitsInFit[k].first, ts, chi2Tmp, ndfTmp);
                    trackStatesPerHit.push_back(ts);
                }
            }
            else{
                // check which hit is closer to the previous hit. and iterate starting from that
                TrackerHit* minRhoHit = hitsInFit.back().first;
                TrackerHit* maxRhoHit = hitsInFit.front().first;
                XYZVector minPos, maxPos;
                XYZVectorF prevFitPos;
                minPos.SetCoordinates(minRhoHit->getPosition());
                maxPos.SetCoordinates(maxRhoHit->getPosition());
                prevFitPos.SetCoordinates( trackStatesPerHit.back().getReferencePoint() );
                if ( (minPos - prevFitPos).r() < (maxPos - prevFitPos).r() ){
                    //iterate from minRho hit
                    for( int k = hitsInFit.size() - 1; k >= 0 ; --k ){
                        TrackStateImpl ts;
                        double chi2Tmp;
                        int ndfTmp;
                        marlinTrk->getTrackState(hitsInFit[k].first, ts, chi2Tmp, ndfTmp);
                        trackStatesPerHit.push_back(ts);
                    }
                }
                else{
                    //iterate from maxRho hit
                    for(unsigned int k = 0; k < hitsInFit.size() ; ++k ){
                        TrackStateImpl ts;
                        double chi2Tmp;
                        int ndfTmp;
                        marlinTrk->getTrackState(hitsInFit[k].first, ts, chi2Tmp, ndfTmp);
                        trackStatesPerHit.push_back(ts);
                    }
                }
            }
            if (j == tracksToFit.size() - 1){
                //add track states at SET hit. If it doesn't exist, then it copies the most outer tpc hit (shouldn't affect track length as it is dublicate)
                trackStatesPerHit.push_back( *(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtLastHit)) ) );
                trackStatesPerHit.push_back( *(dynamic_cast<const TrackStateImpl*> (refittedTrack.getTrackState(TrackState::AtCalorimeter) ) ) );
            }
            delete marlinTrk;
        }
        // getArcLengths
        vector<double> phi, omega, z, arcLength, arcCurl, p, d0, z0, tanL, pWeighted;
        for(unsigned int j=0; j < trackStatesPerHit.size(); ++j ){
            phi.push_back( trackStatesPerHit[j].getPhi() );
            d0.push_back( trackStatesPerHit[j].getD0() );
            z0.push_back( trackStatesPerHit[j].getZ0() );
            omega.push_back( trackStatesPerHit[j].getOmega() );
            tanL.push_back( trackStatesPerHit[j].getTanLambda() );
            z.push_back(trackStatesPerHit[j].getReferencePoint()[2] + trackStatesPerHit[j].getZ0() );
            HelixClass helix;
            helix.Initialize_Canonical(phi[j], d0[j], z0[j], omega[j], tanL[j], _bField);
            XYZVectorF mom;
            mom.SetCoordinates(helix.getMomentum());
            p.push_back( mom.r() );
            //can't calculate track length from 1 element
            if (j == 0) continue;
            //for this we need to check if last arc between lastHit and Ecal less than pi
            if (j == trackStatesPerHit.size() - 1){
                double dLastHitToCalo = std::abs( (z[j] - z[j-1]) / tanL[j-1]);
                if( dLastHitToCalo > M_PI / std::abs(omega[j-1]) ){
                    // we cannot calculate with delta phi formula. Let's use formula with only dz assuming constant curvature
                    arcLength.push_back(std::abs(z[j] - z[j-1]) * std::sqrt( 1. + 1./(tanL[j-1]*tanL[j-1]) ) );
                    arcCurl.push_back( dLastHitToCalo * std::abs(omega[j-1]) );
                    pWeighted.push_back( arcLength[j] /(p[j]*p[j]) );
                    continue;
                }
            }

            arcLength.push_back( getArcLength( phi[j-1], phi[j], omega[j-1], z[j-1], z[j] ) );
            double deltaPhi = std::abs(phi[j] - phi[j-1]);
            if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;
            arcCurl.push_back( deltaPhi );
            pWeighted.push_back( arcLength[j] /(p[j]*p[j]) );
        }

        // just sum
        if (trackStatesPerHit.size() != 0){
            _trackLength["set"] = std::accumulate(arcLength.begin(), arcLength.end() - 1, 0.);
            _nCurls["set"] = std::accumulate(arcCurl.begin(), arcCurl.end() - 1, 0.)/(2.*M_PI);
            _mom["hmSet"] = std::sqrt(_trackLength["set"] / std::accumulate(pWeighted.begin(), pWeighted.end() - 1, 0.) );
        }
        else{
            _trackLength["set"] = 0.;
            _nCurls["set"] = 0.;
            _mom["hmSet"] = 0.;
        }

        _trackLength["ecal"] = std::accumulate(arcLength.begin(), arcLength.end(), 0.);
        _nCurls["ecal"] = std::accumulate(arcCurl.begin(), arcCurl.end(), 0.)/(2.*M_PI);
        _mom["hmEcal"] = std::sqrt(_trackLength["ecal"] / std::accumulate(pWeighted.begin(), pWeighted.end(), 0.) );

        if (false){
            DDMarlinCED::newEvent(this);
            DDMarlinCED::drawDD4hepDetector(_theDetector, 0, vector<string>{});
            DDCEDPickingHandler& pHandler=DDCEDPickingHandler::getInstance();
            pHandler.update(evt);
            TOFAnaUtils::drawPfo(track, cluster);
            DDMarlinCED::draw(this, 1);
        }

        SimTrackerHit* setSimHitFront = nullptr;
        SimTrackerHit* setSimHitBack = nullptr;
        if (_hasSetHit){
            setSimHitFront = dynamic_cast <SimTrackerHit*> (setSimHits.at(0) );
            setSimHitBack = dynamic_cast <SimTrackerHit*> (setSimHits.at(1) );
        }

        CalorimeterHit* closestHit = TOFAnaUtils::getClosestHit( cluster, _tsPos["ecal"] );
        _posClosest.SetCoordinates( closestHit->getPosition() );

        for(unsigned int j=0; j < std::size(_smearings); ++j ){
            pair<XYZVectorF, double> fastestHit = TOFAnaUtils::getFastestHit( cluster, _smearings[j] / 1000. );
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
            _tofFrankFit[j] = TOFAnaUtils::getTofFrankFit( cluster, _tsPos["ecal"], _tsMom["ecal"], _smearings[j] / 1000. );
            _tofFrankAvg[j] = TOFAnaUtils::getTofFrankAvg( cluster, _tsPos["ecal"], _tsMom["ecal"], _smearings[j] / 1000. );
        }
        _tree->Fill();
    }

}

void TOFAnalysis::end(){
    _file->Write();
}
