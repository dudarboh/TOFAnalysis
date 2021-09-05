#include "TOFAnaUtils.hpp"

namespace TOFAnaUtils{

    MCParticle* getMcMaxWeight(LCRelationNavigator pfoToMc, ReconstructedParticle* pfo){
        MCParticle* mc = nullptr;
        const vector <LCObject*>& mcs = pfoToMc.getRelatedToObjects(pfo);
        const vector <float>& mcWeights = pfoToMc.getRelatedToWeights(pfo);
        if (mcs.size() == 0) return mc;
        int maxW = std::max_element(mcWeights.begin(), mcWeights.end()) - mcWeights.begin();
        mc = dynamic_cast <MCParticle*> ( mcs[maxW] );
        return mc;
    }


    pair<double, double> getTpcR(){
        const Detector& detector = Detector::getInstance();
        const DetElement tpcDet = detector.detector("TPC");
        const FixedPadSizeTPCData* tpc = tpcDet.extension <FixedPadSizeTPCData>();

        double rInner = tpc->rMinReadout/dd4hep::mm;
        double rOuter = tpc->rMaxReadout/dd4hep::mm;
        return std::make_pair(rInner, rOuter);
    }

    TrackerHit* getSetHit(Track* track, double tpcROuter){
        const vector<TrackerHit*>& hits = track->getTrackerHits();
        //Performance: loop from the end. SET hits at the end!
        for (int i=hits.size()-1; i>=0; --i){
            XYZVector pos;
            pos.SetCoordinates( hits[i]->getPosition() );
            if ( pos.rho() > tpcROuter ) return hits[i];
        }
        return nullptr;
    }


    CalorimeterHit* getClosestHit( Cluster* cluster, XYZVectorF posTrackAtCalo){
        CalorimeterHit* closestHit = nullptr;

        double closestDistance = numeric_limits<double>::max();
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


    pair<XYZVectorF, double> getFastestHit( Cluster* cluster, double smearing ){
        pair<XYZVectorF, double> earliestHit{};

        earliestHit.second = numeric_limits<double>::max();
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


    double getTofFrankFit( Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing, unsigned int nLayers ){
        vector <double> d;
        vector <double> time;
        vector <double> d_err;
        vector <double> time_err;

        for (unsigned int l=0; l < nLayers; ++l){
            double closestDistance = numeric_limits<double>::max();
            double closestTime = numeric_limits<double>::max();
            double closestDistanceToLine = numeric_limits<double>::max();

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
            if ( closestDistanceToLine == numeric_limits<double>::max() ) continue;
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


    double getTofFrankAvg( Cluster* cluster, XYZVectorF posTrackAtCalo, XYZVectorF momTrackAtCalo, double smearing, unsigned int nLayers ){
        double tof = 0.;
        int nHits = 0;

        for (unsigned int l=0; l < nLayers; ++l){
            double closestDistance = numeric_limits<double>::max();
            double closestTime = numeric_limits<double>::max();
            double closestDistanceToLine = numeric_limits<double>::max();

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
            if ( closestDistanceToLine == numeric_limits<double>::max() ) continue;

            tof += closestTime - closestDistance/CLHEP::c_light;
            ++nHits;
        }
        if (nHits == 0) return 0.;

        return tof / nHits;
    }



    int getNEcalHits(Cluster* cluster){
        int nHits = 0;
        for ( const auto hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if (isECALHit) ++nHits;
        }
        return nHits;
    }


    void drawPfo(Track* track, Cluster* cluster){
        vector<Track*> subTracks;
        subTracks.push_back(track);
        int nSubTracks = track->getTracks().size();
        if (nSubTracks != 1){
            int nHitsTot = track->getTrackerHits().size();
            int nHits0 = track->getTracks()[0]->getTrackerHits().size();
            int nHits1 = track->getTracks()[1]->getTrackerHits().size();
            if( ((nHitsTot == nHits0 + nHits1)) || ((nHitsTot == nHits0 + nHits1 + 1))  ){
                for(int j=2; j < nSubTracks; ++j) subTracks.push_back( track->getTracks()[j] );
            }
            else{
                // This means VXD+SIT subTrack is not there (bug?!)! So we need to add tracks from i=1
                for(int j=1; j < nSubTracks; ++j) subTracks.push_back( track->getTracks()[j] );
            }
        }

        // streamlog_out( DEBUG6 ) << " -- Final LCIO Track has "<<nSubTracks<< " subtracks - will use these for displaying hits "<< std::endl ;
        // if ( nSubTracks == 0 ) return;
        //
        // // std::copy( track->getTrackerHits().begin() , track->getTrackerHits().end() , std::back_inserter(  hits ) ) ;
        // cout<<"Track hits size:"<<track->getTrackerHits().size()<<endl;
        // cout<<"N subdetectors:"<<track->getSubdetectorHitNumbers().size()<<endl;
        // cout<<"VXD used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::VXD)*2-2]<<endl;
        // cout<<"VXD not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::VXD)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::VXD)*2-2]<<endl;
        // cout<<"SIT used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::SIT)*2-2]<<endl;
        // cout<<"SIT not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::SIT)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::SIT)*2-2]<<endl;
        // cout<<"FTD used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::FTD)*2-2]<<endl;
        // cout<<"FTD not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::FTD)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::FTD)*2-2]<<endl;
        // cout<<"TPC used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::TPC)*2-2]<<endl;
        // cout<<"TPC not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::TPC)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::TPC)*2-2]<<endl;
        // cout<<"SET used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::SET)*2-2]<<endl;
        // cout<<"SET not used: "<<track->getSubdetectorHitNumbers()[(ILDDetID::SET)*2-1] - track->getSubdetectorHitNumbers()[(ILDDetID::SET)*2-2]<<endl;

        for(unsigned int i=0; i < subTracks.size(); ++i){
            // std::copy( subTrack->getTrackerHits().begin() , subTrack->getTrackerHits().end() , std::back_inserter(  hits ) ) ;
            const Track* subTrack = subTracks[i];
            int nHits = subTrack->getTrackerHits().size();
            streamlog_out( DEBUG6 ) << " -- subTrack i= "<< i << " has " <<  nHits << " hits and "<< endl;
            for (int j = 0; j < nHits; ++j) {
                TrackerHit* hit = subTrack->getTrackerHits()[j];
                XYZVector hitPos;
                hitPos.SetCoordinates( hit->getPosition() );
                ced_hit(hitPos.x(), hitPos.y(), hitPos.z(), 0, 3, 0x0400ff);
            }

            const TrackState* ts = subTrack->getTrackState(TrackState::AtLastHit);
            XYZVectorF pos;
            pos.SetCoordinates( ts->getReferencePoint() );
            ced_hit(pos.x(), pos.y(), pos.z(), 0, 12, 0x00ff42);


            if (i == 0){
                const TrackState* tsIp = subTrack->getTrackState(TrackState::AtIP);
                XYZVectorF posIp;
                posIp.SetCoordinates( tsIp->getReferencePoint() );
                ced_hit(posIp.x(), posIp.y(), posIp.z(), 0, 12, 0x00ff42);

                double bField = MarlinUtil::getBzAtOrigin();
                double pt = bField * 3e-4 / std::abs( tsIp->getOmega() );
                double charge = ( tsIp->getOmega() > 0. ?  1. : -1. );
                double Px = pt * std::cos(  tsIp->getPhi() );
                double Py = pt * std::sin(  tsIp->getPhi() );
                double Pz = pt * tsIp->getTanLambda() ;
                double Xs = tsIp->getReferencePoint()[0] -  tsIp->getD0() * std::sin( tsIp->getPhi() ) ;
                double Ys = tsIp->getReferencePoint()[1] +  tsIp->getD0() * std::cos( tsIp->getPhi() ) ;
                double Zs = tsIp->getReferencePoint()[2] +  tsIp->getZ0() ;

                DDMarlinCED::drawHelix(bField, charge, Xs, Ys, Zs, Px, Py, Pz, 1, 2, 0xff0000, 0., 1855., 2450., 0);

            }
            if (i == subTracks.size() - 1){
                const TrackState* tsEcal = subTrack->getTrackState(TrackState::AtCalorimeter);
                XYZVectorF posEcal;
                posEcal.SetCoordinates( tsEcal->getReferencePoint() );
                ced_hit(posEcal.x(), posEcal.y(), posEcal.z(), 0, 12, 0x00ff42);

                double bField = MarlinUtil::getBzAtOrigin();
                double pt = bField * 3e-4 / std::abs( tsEcal->getOmega() );
                double charge = ( tsEcal->getOmega() > 0. ?  1. : -1. );
                double Px = pt * std::cos(  tsEcal->getPhi() ) ;
                double Py = pt * std::sin(  tsEcal->getPhi() ) ;
                double Pz = pt * tsEcal->getTanLambda() ;
                double Xs = tsEcal->getReferencePoint()[0] -  tsEcal->getD0() * std::sin( tsEcal->getPhi() ) ;
                double Ys = tsEcal->getReferencePoint()[1] +  tsEcal->getD0() * std::cos( tsEcal->getPhi() ) ;
                double Zs = tsEcal->getReferencePoint()[2] +  tsEcal->getZ0() ;
                //helix at ECal
                DDMarlinCED::drawHelix(-bField, -charge, Xs, Ys, Zs, Px, Py, Pz, 1, 2, 0xff9300, 0., 1855., 2450., 0);


            }

        }

        for ( const auto hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if (! isECALHit) continue;

            XYZVectorF hitPos;
            hitPos.SetCoordinates( hit->getPosition() );
            ced_hit(hitPos.x(), hitPos.y(), hitPos.z(), 0, 7, 0xff0000);
        }


    }


}
