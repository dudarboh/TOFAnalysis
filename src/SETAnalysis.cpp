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

using ROOT::Math::XYZVector;
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
    _tree->Branch("pos_closest_sim", &_posClosestSim);
    _tree->Branch("tof_closest", &_tofClosest);
    _tree->Branch("tof_closest_sim", &_tofClosestSim);

    _tree->Branch("pos_fastest", &_posFastest);
    _tree->Branch("pos_fastest_sim", &_posFastestSim);
    _tree->Branch("tof_fastest", &_tofFastest);
    _tree->Branch("tof_fastest_sim", &_tofFastestSim);

    _tree->Branch("tof_frank_fit", &_tofFrankFit);
    _tree->Branch("tof_frank_avg", &_tofFrankAvg);

    _tree->Branch("n_set_hits", &_nSetHits);
    _tree->Branch("set_hit_pos", &_setHitPos);
    _tree->Branch("set_hit_time", &_setHitTime);
    _tree->Branch("set_pos_true", &_setPosTrue);

    _tree->Branch("track_length_set", &_trackLengthSet);
    _tree->Branch("track_length_calo", &_trackLengthCalo);
    _tree->Branch("track_length_integral", &_trackLengthIntegral);

    _bField = MarlinUtil::getBzAtOrigin();
    _tpcROuter = getTpcR().second;
}


void SETAnalysis::processEvent(LCEvent* evt){

    LCCollection* pfos = evt->getCollection("PandoraPFOs");
    LCRelationNavigator pfoToMc( evt->getCollection("RecoMCTruthLink") );
    LCRelationNavigator spPointToSimHit( evt->getCollection("SETSpacePointRelations") );
    LCRelationNavigator caloHitToSimHit( evt->getCollection("EcalBarrelRelationsSimRec") );

    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        MCParticle* mc = getMcMaxWeight(pfoToMc, pfo);
        const vector<Track*>& tracks = pfo->getTracks();
        const vector<Cluster*>& clusters = pfo->getClusters();

        if(tracks.size() != 1 || clusters.size() != 1) continue;

        Track* track = tracks[0];
        Cluster* cluster = clusters[0];

        const TrackState* tsIp = track->getTrackState(TrackState::AtIP);
        const TrackState* tsFirst = track->getTrackState(TrackState::AtFirstHit);
        const TrackState* tsLast = track->getTrackState(TrackState::AtLastHit);
        const TrackState* tsCalo = track->getTrackState(TrackState::AtCalorimeter);

        ///////////////////////////////////////////////////////////////////////////////////////////////////

        _pdg = mc->getPDG();

        vector<TrackerHit*> setHits = getSetHits(track, _tpcROuter);
        _nSetHits = setHits.size();
        _nEcalHits = getNEcalHits(cluster);
        if (_nSetHits != 1 || _nEcalHits == 0) continue;
        _setHitPos = XYZVector( setHits[0]->getPosition()[0], setHits[0]->getPosition()[1], setHits[0]->getPosition()[2] );
        _setHitTime = CLHEP::RandGauss::shoot( setHits[0]->getTime(), _smearing / 1000. ) ;

        SimTrackerHit* setSimHit = dynamic_cast <SimTrackerHit*> (spPointToSimHit.getRelatedToObjects( setHits[0] )[0] );
        _setPosTrue = XYZVector( setSimHit->getPosition()[0], setSimHit->getPosition()[1], setSimHit->getPosition()[2] );

        _tsLastOmega = tsLast->getOmega();
        _tsLastTanL = tsLast->getTanLambda();
        _tsLastPhi = tsLast->getPhi();
        _tsLastD0 = tsLast->getD0();
        _tsLastZ0 = tsLast->getZ0();
        _tsLastPos = XYZVector( tsLast->getReferencePoint()[0], tsLast->getReferencePoint()[1], tsLast->getReferencePoint()[2] );
        HelixClass helixLast;
        helixLast.Initialize_Canonical(_tsLastPhi, _tsLastD0, _tsLastZ0, _tsLastOmega, _tsLastTanL, _bField);
        _tsLastMom = XYZVector( helixLast.getMomentum()[0], helixLast.getMomentum()[1], helixLast.getMomentum()[2] );

        _tsCaloOmega = tsCalo->getOmega();
        _tsCaloTanL = tsCalo->getTanLambda();
        _tsCaloPhi = tsCalo->getPhi();
        _tsCaloD0 = tsCalo->getD0();
        _tsCaloZ0 = tsCalo->getZ0();
        _tsCaloPos = XYZVector( tsCalo->getReferencePoint()[0], tsCalo->getReferencePoint()[1], tsCalo->getReferencePoint()[2] );
        HelixClass helixCalo;
        helixCalo.Initialize_Canonical(_tsCaloPhi, _tsCaloD0, _tsCaloZ0, _tsCaloOmega, _tsCaloTanL, _bField);
        _tsCaloMom = XYZVector( helixCalo.getMomentum()[0], helixCalo.getMomentum()[1], helixCalo.getMomentum()[2] );

        _trackLengthSet = getTrackLength(track, TrackState::AtIP, TrackState::AtLastHit);
        _trackLengthCalo = getTrackLength(track, TrackState::AtIP, TrackState::AtCalorimeter);
        _trackLengthIntegral = getTrackLengthIntegral(track);
        _tofFrankFit = getTofFrankFit( cluster, _tsCaloPos, _tsCaloMom, _smearing / 1000. );
        _tofFrankAvg = getTofFrankAvg( cluster, _tsCaloPos, _tsCaloMom, _smearing / 1000. );


        CalorimeterHit* closestHit =  getClosestHit( cluster, _tsCaloPos);
        _posClosest = XYZVector( closestHit->getPosition()[0], closestHit->getPosition()[1], closestHit->getPosition()[2] );
        _tofClosest = CLHEP::RandGauss::shoot( closestHit->getTime(), _smearing / 1000. ) - ( _posClosest - _tsCaloPos ).r()/CLHEP::c_light;

        //0 SimHits in barrel collection if in hit is in the ENDCAP!
        int nClosestSimHits = caloHitToSimHit.getRelatedToObjects(closestHit).size();
        if (nClosestSimHits != 0){
            SimCalorimeterHit* closestSimHit = dynamic_cast<SimCalorimeterHit*> ( caloHitToSimHit.getRelatedToObjects(closestHit)[0] );
            _posClosestSim = getFastestContPos(closestSimHit);
            _tofClosestSim = CLHEP::RandGauss::shoot( closestHit->getTime(), _smearing / 1000. ) - ( _posClosestSim - _tsCaloPos ).r()/CLHEP::c_light;
        }
        else {
            _posClosestSim = XYZVector();
            _tofClosestSim = 0.;
        }

        CalorimeterHit* fastestHit =  getFastestHit( cluster );
        _posFastest = XYZVector( fastestHit->getPosition()[0], fastestHit->getPosition()[1], fastestHit->getPosition()[2] );
        _tofFastest = CLHEP::RandGauss::shoot( fastestHit->getTime(), _smearing / 1000. ) - ( _posFastest - _tsCaloPos ).r()/CLHEP::c_light;
        int nFastestSimHits = caloHitToSimHit.getRelatedToObjects(fastestHit).size();
        if (nFastestSimHits != 0){
            SimCalorimeterHit* fastestSimHit = dynamic_cast<SimCalorimeterHit*> ( caloHitToSimHit.getRelatedToObjects(fastestHit)[0] );
            _posFastestSim = getFastestContPos(fastestSimHit);
            _tofFastestSim = CLHEP::RandGauss::shoot( fastestHit->getTime(), _smearing / 1000. ) - ( _posFastestSim - _tsCaloPos ).r()/CLHEP::c_light;
        }
        else {
            _posFastestSim = XYZVector();
            _tofFastestSim = 0.;
        }
        _tree->Fill();
    }

}

void SETAnalysis::end(){
    _file->Write();
}


MCParticle* SETAnalysis::getMcMaxWeight(LCRelationNavigator pfoToMc, ReconstructedParticle* pfo){
    const vector <LCObject*>& mcs = pfoToMc.getRelatedToObjects(pfo);
    const vector <float>& mcWeights = pfoToMc.getRelatedToWeights(pfo);
    int maxW = std::max_element(mcWeights.begin(), mcWeights.end()) - mcWeights.begin();
    MCParticle* mc = dynamic_cast <MCParticle*> ( mcs[maxW] );
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
    for(const auto hit: hits){
        XYZVector pos;
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
        XYZVector posA = XYZVector(a->getPosition()[0], a->getPosition()[1], a->getPosition()[2]);
        XYZVector posB = XYZVector(b->getPosition()[0], b->getPosition()[1], b->getPosition()[2]);
        return posA.r() < posB.r();
    };
    sort(trackHits.begin(), trackHits.end(), sortByR);

    for (size_t j=1; j < trackHits.size(); ++j){
        XYZVector pos2 = XYZVector(trackHits[j]->getPosition()[0], trackHits[j]->getPosition()[1], trackHits[j]->getPosition()[2]);
        XYZVector pos1 = XYZVector(trackHits[j-1]->getPosition()[0], trackHits[j-1]->getPosition()[1], trackHits[j-1]->getPosition()[2]);
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


CalorimeterHit* SETAnalysis::getClosestHit( Cluster* cluster, XYZVector posTrackAtCalo){
    CalorimeterHit* closestHit = nullptr;

    double closestDistance = std::numeric_limits<double>::max();
    for ( const auto hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isECALHit = ( hitType.caloID() == CHT::ecal );
        if (! isECALHit) continue;

        XYZVector hitPos( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
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


double SETAnalysis::getTofFrankFit( Cluster* cluster, XYZVector posTrackAtCalo, XYZVector momTrackAtCalo, double smearing, unsigned int nLayers ){
    vector <double> d(nLayers);
    vector <double> time(nLayers);
    vector <double> d_err(nLayers);
    vector <double> time_err(nLayers);
    vector <double> closestDistanceToLine(nLayers);

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

            XYZVector pos( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
            double dToLine = (pos - posTrackAtCalo).Cross(momTrackAtCalo.unit()).r();
            if (dToLine < closestDistanceToLine[l]){
                closestDistanceToLine[l] = dToLine;
                time[l] = CLHEP::RandGauss::shoot(hit->getTime(), smearing);
                d[l] = (pos - posTrackAtCalo).r();
            }
        }
    }

    TGraphErrors gr(nLayers, &d[0], &time[0], &d_err[0], &time_err[0]);
    gr.Fit("pol1", "Q");
    return gr.GetFunction("pol1")->GetParameter(0);
}


double SETAnalysis::getTofFrankAvg( Cluster* cluster, XYZVector posTrackAtCalo, XYZVector momTrackAtCalo, double smearing, unsigned int nLayers ){
    vector <double> d(nLayers);
    vector <double> time(nLayers);
    vector <double> closestDistanceToLine(nLayers);

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

            XYZVector pos( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
            double dToLine = (pos - posTrackAtCalo).Cross(momTrackAtCalo.unit()).r();
            if (dToLine < closestDistanceToLine[l]){
                closestDistanceToLine[l] = dToLine;
                time[l] = CLHEP::RandGauss::shoot(hit->getTime(), smearing);
                d[l] = (pos - posTrackAtCalo).r();
            }
        }
        tof += time[l] - d[l]/CLHEP::c_light;
    }
    return tof / nLayers;
}


XYZVector SETAnalysis::getFastestContPos(SimCalorimeterHit* hit){
    XYZVector pos;
    double earliestTime = std::numeric_limits<double>::max();
    for(int i=0; i < hit->getNMCParticles(); ++i) {
        double contTime = hit->getTimeCont(i);
        if(contTime < earliestTime){
            pos = XYZVector(hit->getStepPosition(i)[0], hit->getStepPosition(i)[1], hit->getStepPosition(i)[2]);
            earliestTime = contTime;
        }
    }
    return pos;
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
