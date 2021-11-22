#include "TOFAnalysis.hpp"
#include "TOFAnaUtils.hpp"

#include <chrono>

#include "EVENT/LCCollection.h"
#include "UTIL/PIDHandler.h"

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/GeometryUtil.h"
#include "MarlinTrk/Factory.h"
#include "EVENT/SimTrackerHit.h"
#include "UTIL/LCRelationNavigator.h"
#include "CLHEP/Random/Randomize.h"

#include "HelixClass.h"

using namespace TOFAnaUtils;
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::ReconstructedParticle;
using EVENT::TrackerHit;
using EVENT::Track;
using EVENT::SimTrackerHit;
using EVENT::Cluster;
using EVENT::CalorimeterHit;
using EVENT::TrackState;
using EVENT::LCObject;
using UTIL::LCRelationNavigator;
using CLHEP::RandGauss;
using dd4hep::rec::Vector3D;

TOFAnalysis aTOFAnalysis ;


TOFAnalysis::TOFAnalysis() : marlin::Processor("TOFAnalysis") {}


void TOFAnalysis::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);
    _bField = MarlinUtil::getBzAtOrigin();
    _tpcOuterR = getTPCOuterR();

    _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
    _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, true);
    _trkSystem->init();

    _file.reset( new TFile("rename.root", "RECREATE") );
    _tree.reset( new TTree("TOFAnalysis", "TOFAnalysis") );

    _tree->Branch("pdg", &_pdg);
    _tree->Branch("has_set_hit", &_hasSetHit);
    _tree->Branch("n_ecal_hits", &_nEcalHits);
    _tree->Branch("ts_pos_ecal", &_posECAL);
    _tree->Branch("ts_mom_ecal", &_momECAL);
    _tree->Branch("ts_z0_ecal", &_z0ECAL);
    _tree->Branch("ts_d0_ecal", &_d0ECAL);

    _tree->Branch("track_length_set", &_trackLengthSET);
    _tree->Branch("momentum_set", &_momentumSET);

    _tree->Branch("track_length_ecal", &_trackLengthECAL);
    _tree->Branch("momentum_ecal", &_momentumECAL);

    //for every time resolution
    for(int j=0; j < 11; ++j){
        _tree->Branch(Form("tof_closest_%dps", j*10 ), &(_tof_closest[j]) );
        _tree->Branch(Form("tof_set_%dps", j*10 ), &(_tof_set[j]) );

        //for number of layers to take
        for(int l=1; l<30; ++l){
            _tree->Branch(Form("tof_avg_%dps_%dl", j*10, l), &(_tof_avg[j][l]) );
            _tree->Branch(Form("tof_fit_%dps_%dl", j*10, l), &(_tof_fit[j][l]) );
        }
    }

}


void TOFAnalysis::processEvent(EVENT::LCEvent * evt){
    RandGauss::setTheSeed( marlin::Global::EVENTSEEDER->getSeed(this) );
    ++_nEvent;
    streamlog_out(MESSAGE)<<"************Event************ "<<_nEvent<<std::endl;

    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator pfoToMc( evt->getCollection("RecoMCTruthLink") );
    LCCollection* setRelations = evt->getCollection("SETSpacePointRelations");    
    LCRelationNavigator navigatorSET = LCRelationNavigator( setRelations );

    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        //reset all branches
        for(int j=0; j < 11; ++j){
            _tof_closest[j] = 0.;
            _tof_set[j] = 0.;
            //for number of layers to take
            for(int l=1; l<30; ++l){
                _tof_avg[j][l] = 0.;
                _tof_fit[j][l] = 0.;
            }
        }

        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );

        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        if( nClusters != 1 || nTracks != 1) continue;
        Track* track = pfo->getTracks()[0];
        const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
        double tsOmega = ts->getOmega();
        double tsTanL = ts->getTanLambda();
        double tsPhi =  ts->getPhi();
        _d0ECAL = ts->getD0();
        _z0ECAL = ts->getZ0();
        _posECAL.SetCoordinates( ts->getReferencePoint() );
        HelixClass helixEcal;
        helixEcal.Initialize_Canonical(tsPhi, _d0ECAL, _z0ECAL, tsOmega, tsTanL, _bField);
        _momECAL.SetCoordinates( helixEcal.getMomentum() );

        Cluster* cluster = pfo->getClusters()[0];
        _nEcalHits = getNEcalHits(cluster);

        MCParticle* mc = TOFAnaUtils::getMcMaxWeight(pfoToMc, pfo);
        if(mc == nullptr) continue;
        _pdg = mc->getPDG();

        ///////////////////////////////////////////////////////////////
        // This part calculates track length and momentum harmonic mean
        ///////////////////////////////////////////////////////////////
        vector<Track*> subTracks = getSubTracks(track);

        {
            vector<TrackStateImpl> trackStates = getTrackStatesPerHit(subTracks, _trkSystem, false, _bField);

            double trackLength = 0.;
            double harmonicMom = 0.;
            int nTrackStates = trackStates.size();
            for( int j=1; j < nTrackStates; ++j ){
                //we check which track length formula to use
                double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
                double arcLength;
                // we cannot calculate arc length for more than pi revolution using delta phi. Use formula with only z
                if ( nTurns <= 0.5 ) arcLength = getHelixArcLength( trackStates[j-1], trackStates[j] );
                else arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );

                Vector3D mom = getHelixMomAtTrackState( trackStates[j-1], _bField );
                trackLength += arcLength;
                harmonicMom += arcLength/mom.r2();
            }
            harmonicMom = std::sqrt(trackLength/harmonicMom);

            _trackLengthSET = trackLength;
            _momentumSET = harmonicMom;
        }
        {
            vector<TrackStateImpl> trackStates = getTrackStatesPerHit(subTracks, _trkSystem, true, _bField);

            double trackLength = 0.;
            double harmonicMom = 0.;
            int nTrackStates = trackStates.size();
            for( int j=1; j < nTrackStates; ++j ){
                //we check which track length formula to use
                double nTurns = getHelixNRevolutions( trackStates[j-1], trackStates[j] );
                double arcLength;
                // we cannot calculate arc length for more than pi revolution using delta phi. Use formula with only z
                if ( nTurns <= 0.5 ) arcLength = getHelixArcLength( trackStates[j-1], trackStates[j] );
                else arcLength = getHelixLengthAlongZ( trackStates[j-1], trackStates[j] );

                Vector3D mom = getHelixMomAtTrackState( trackStates[j-1], _bField );
                trackLength += arcLength;
                harmonicMom += arcLength/mom.r2();
            }
            harmonicMom = std::sqrt(trackLength/harmonicMom);

            _trackLengthECAL = trackLength;
            _momentumECAL = harmonicMom;
        }

        TrackerHit* hitSET = getSETHit(track, _tpcOuterR);
        _hasSetHit = ( hitSET != nullptr );

        //for every time resolution
        for(int j=0; j < 11; ++j){
            double timeResolution = j*10./1000.; // in ns
            _tof_closest[j] = getTofClosest(cluster, track, timeResolution);

            _tof_set[j] = 0.;
            if ( _hasSetHit ){
                const vector<LCObject*>& simHitsSET = navigatorSET.getRelatedToObjects( hitSET );
                if ( simHitsSET.size() >= 2 ){
                    //It must be always 2, but just in case...
                    if (simHitsSET.size() > 2) streamlog_out(WARNING)<<"Found more than two SET strip hits! Writing TOF as an average of the first two elements in the array."<<std::endl;

                    SimTrackerHit* simHitSETFront = static_cast <SimTrackerHit*>( simHitsSET[0] );
                    SimTrackerHit* simHitSETBack = static_cast <SimTrackerHit*>( simHitsSET[1] );
                    double timeFront = RandGauss::shoot(simHitSETFront->getTime(), timeResolution);
                    double timeBack = RandGauss::shoot(simHitSETBack->getTime(), timeResolution);
                        _tof_set[j] = (timeFront + timeBack)/2.;
                }
                else if (simHitsSET.size() == 1){
                    streamlog_out(WARNING)<<"Found only one SET strip hit! Writing TOF from a single strip."<<std::endl;
                    SimTrackerHit* simHitSET = static_cast <SimTrackerHit*>(simHitsSET[0]);
                    _tof_set[j] = RandGauss::shoot(simHitSET->getTime(), timeResolution);
                }
                else{
                    // this happens very rarily (0.1%). When >1 simHits associated with a single strip none simHits are written by the DDSpacePointBuilder.
                    streamlog_out(WARNING)<<"Found NO simHits associated with the found SET hit! Writing TOF as 0."<<std::endl;
                }
            }

            //for number of layers to take
            for(int l=1; l<30; ++l){
                vector<CalorimeterHit*> frankHits = selectFrankEcalHits(cluster, track, l, _bField);
                _tof_avg[j][l] = getTofFrankAvg(frankHits, track, timeResolution);
                _tof_fit[j][l] = getTofFrankFit(frankHits, track, timeResolution);
            }
        }

    _tree->Fill();
    }
}

void TOFAnalysis::end(){
    _file->Write();
}
