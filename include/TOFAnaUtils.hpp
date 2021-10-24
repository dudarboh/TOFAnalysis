#ifndef TOFAnaUtils_h
#define TOFAnaUtils_h 1

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "UTIL/LCRelationNavigator.h"
#include <UTIL/ILDConf.h>

#include "marlinutil/DDMarlinCED.h"
#include "marlinutil/CalorimeterHitType.h"
#include "marlinutil/GeometryUtil.h"

#include "IMPL/TrackStateImpl.h"


#include "DDRec/DetectorData.h"
#include "Math/Vector3D.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "CLHEP/Random/Randomize.h"

namespace TOFAnaUtils{
    using std::cout, std::endl, std::vector, std::string, std::pair, std::numeric_limits;
    using ROOT::Math::XYZVector;
    using ROOT::Math::XYZVectorF;
    using dd4hep::Detector;
    using dd4hep::DetElement;
    using dd4hep::rec::FixedPadSizeTPCData;
    using EVENT::ReconstructedParticle;
    using EVENT::MCParticle;
    using EVENT::LCObject;
    using EVENT::Track;
    using EVENT::TrackerHit;
    using EVENT::Cluster;
    using EVENT::CalorimeterHit;
    using UTIL::LCRelationNavigator;

    MCParticle* getMcMaxWeight(LCRelationNavigator pfoToMc, ReconstructedParticle* pfo);
    pair<double, double> getTpcR();
    TrackerHit* getSetHit(Track* track, double tpcROuter);

    CalorimeterHit* getClosestHit(Cluster* cluster, XYZVectorF posTrackAtCalo);
    pair<XYZVectorF, double> getFastestHit(Cluster* cluster, double smearing=0.);
    double getTofFrankFit(Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing=0., unsigned int nLayers=10);
    double getTofFrankAvg(Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing=0., unsigned int nLayers=10);

    int getNEcalHits(Cluster* cluster);

    void drawPfo(Track* track, Cluster* cluster, const TrackStateImpl tsEcal);

}

#endif
