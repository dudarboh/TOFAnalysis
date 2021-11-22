#ifndef TOFAnaUtils_h
#define TOFAnaUtils_h 1

#include <vector>
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "IMPL/TrackStateImpl.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "DDRec/Vector3D.h"
#include "UTIL/LCRelationNavigator.h"

namespace TOFAnaUtils{

    EVENT::MCParticle* getMcMaxWeight(UTIL::LCRelationNavigator pfoToMc, EVENT::ReconstructedParticle* pfo);

    int getNEcalHits(EVENT::Cluster* cluster);

    bool sortByRho(EVENT::TrackerHit* a, EVENT::TrackerHit* b);

    IMPL::TrackStateImpl getTrackStateAtHit(MarlinTrk::IMarlinTrack* marlinTrk, EVENT::TrackerHit* hit);

    dd4hep::rec::Vector3D getHelixMomAtTrackState(const EVENT::TrackState& ts, double bField);

    double getHelixArcLength(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    double getHelixLengthAlongZ(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    double getHelixNRevolutions(const EVENT::TrackState& ts1, const EVENT::TrackState& ts2);

    double getTPCOuterR();

    EVENT::TrackerHit* getSETHit(EVENT::Track* track, double tpcOuterR);

    std::vector<EVENT::CalorimeterHit*> selectFrankEcalHits( EVENT::Cluster* cluster, EVENT::Track* track, int maxEcalLayer, double bField );

    std::vector<EVENT::Track*> getSubTracks(EVENT::Track* track);

    std::vector<IMPL::TrackStateImpl> getTrackStatesPerHit(std::vector<EVENT::Track*> tracks, MarlinTrk::IMarlinTrkSystem* trkSystem, bool extrapolateToEcal, double bField);

    double getTofClosest( EVENT::Cluster* cluster, EVENT::Track* track, double timeResolution);

    double getTofFrankAvg( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution);

    double getTofFrankFit( std::vector<EVENT::CalorimeterHit*> selectedHits, EVENT::Track* track, double timeResolution);

}



#endif
