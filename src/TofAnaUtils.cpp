#include "TofAnaUtils.hpp"

#include <limits>
#include <vector>
#include <algorithm>

#include "EVENT/Track.h"
#include "EVENT/TrackState.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/CalorimeterHit.h"
#include "HelixClass.h"
#include "marlinutil/CalorimeterHitType.h"

#include "DD4hep/DetectorSelector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "DDRec/Vector3D.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <marlinutil/GeometryUtil.h>

//////////////////////////////////////
// #include "EVENT/CalorimeterHit.h"
// using EVENT::TrackerHit;
//
// using dd4hep::rec::LayeredCalorimeterData;
// using dd4hep::Detector, dd4hep::DetType, dd4hep::DetElement, dd4hep::rec::FixedPadSizeTPCData;
// using std::string, std::vectorm std::cout, std::endl;

using EVENT::ReconstructedParticle;
using EVENT::Track;
using EVENT::TrackState;
using EVENT::TrackerHit;
using EVENT::CalorimeterHit;
using EVENT::Cluster;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::FixedPadSizeTPCData;
using dd4hep::rec::Vector3D;

std::pair<double, double> getTpcR(const Detector& detector){
    // const Detector& detector = Detector::getInstance();
    const DetElement tpcDet = detector.detector("TPC");
    const FixedPadSizeTPCData* tpc = tpcDet.extension <FixedPadSizeTPCData>();
    double rInner = tpc->rMinReadout/dd4hep::mm;
    double rOuter = tpc->rMaxReadout/dd4hep::mm;
    return std::make_pair(rInner, rOuter);
}


Vector3D calculateMomentum(ReconstructedParticle* pfo, int location, double bField){
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


double calculateTrackLength(ReconstructedParticle* pfo, int location){
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


double calculateTrackLengthIntegral(ReconstructedParticle* pfo){
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



void findShowerStart(ReconstructedParticle* pfo){
    std::cout<<"okok"<<std::endl;
    // auto findIter = integrations.find( _integration_method ) ;

}


double TofClosest::calculate( ReconstructedParticle* pfo ){
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
            time = CLHEP::RandGauss::shoot(hit->getTime(), _smearing);
        }
    }
    return time - closestDistance/CLHEP::c_light;
}


double TofFastest::calculate( ReconstructedParticle* pfo ){
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
            time = CLHEP::RandGauss::shoot(hit->getTime(), _smearing);
            distance = (pos - posTrackAtCalo).r();
        }
    }

    return time - distance/CLHEP::c_light;
}


double TofFrankFit::calculate( ReconstructedParticle* pfo ){
    //No track --- no correction for distance is possible time
    if ( pfo->getTracks().size() == 0 ) return 0;

    // Get Track position at ECAL
    const Track* track = pfo->getTracks()[0];
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D posTrackAtCalo( ts->getReferencePoint() );


    double bField = MarlinUtil::getBzAtOrigin();
    Vector3D momAtECAL = calculateMomentum(pfo, TrackState::AtCalorimeter, bField);

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
                time[l] = CLHEP::RandGauss::shoot(hit->getTime(), _smearing);
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


double TofFrankAvg::calculate( ReconstructedParticle* pfo ){
    //No track --- no correction for distance is possible time
    if ( pfo->getTracks().size() == 0 ) return 0;

    // Get Track position at ECAL
    const Track* track = pfo->getTracks()[0];
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    Vector3D posTrackAtCalo( ts->getReferencePoint() );


    double bField = MarlinUtil::getBzAtOrigin();
    Vector3D momAtECAL = calculateMomentum(pfo, TrackState::AtCalorimeter, bField);

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
                time[l] = CLHEP::RandGauss::shoot(hit->getTime(), _smearing);
                d[l] = (pos - posTrackAtCalo).r();
            }
        }
        tof += (time[l] - d[l]/CLHEP::c_light) / nLayers;
    }
    return tof;
}

double TofSet::calculate( ReconstructedParticle* pfo, double rTpcOuter ){
    if (pfo->getTracks().size() == 0){
        return 0.;
    }
    const Track* track = pfo->getTracks()[0];
    double tof = std::numeric_limits<double>::max();
    for ( const auto& hit : track->getTrackerHits() ){
        Vector3D pos( hit->getPosition() );
        bool isSETHit = pos.rho() > rTpcOuter;

        if (isSETHit){
            if (hit->getTime() < tof) tof = hit->getTime();
        }
    }
    return tof;
}
