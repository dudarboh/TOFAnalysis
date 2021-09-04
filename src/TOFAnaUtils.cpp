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


    std::pair<double, double> getTpcR(){
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

        vector<int> hitsColors = {0x1700ff, 0x000fff, 0x0032ff, 0x005dff, 0x0080ff};
        for(int i=0; i < nSubTracks; ++i){
            // std::copy( subTrack->getTrackerHits().begin() , subTrack->getTrackerHits().end() , std::back_inserter(  hits ) ) ;
            const Track* subTrack = track->getTracks()[i];
            int nHits = subTrack->getTrackerHits().size();
            streamlog_out( DEBUG6 ) << " -- subTrack i= "<< i << " has " <<  nHits << " hits and "<< endl;
            for (int j = 0; j < nHits; ++j) {
                TrackerHit* hit = subTrack->getTrackerHits()[j];
                XYZVector hitPos;
                hitPos.SetCoordinates( hit->getPosition() );
                ced_hit(hitPos.x(), hitPos.y(), hitPos.z(), 0, 3, hitsColors[i]);
            }

        }

        for ( const auto hit : cluster->getCalorimeterHits() ){
            CHT hitType( hit->getType() );
            bool isECALHit = ( hitType.caloID() == CHT::ecal );
            if (! isECALHit) continue;

            XYZVectorF hitPos;
            hitPos.SetCoordinates( hit->getPosition() );
            ced_hit(hitPos.x(), hitPos.y(), hitPos.z(), 0, 7, 0x08ff00);
        }



        const TrackState* tsIP = track->getTrackState(TrackState::AtIP);
        ced_hit(tsIP->getReferencePoint()[0], tsIP->getReferencePoint()[1], tsIP->getReferencePoint()[2], 0, 12, 0x00e8ff);
        const TrackState* tsFirst = track->getTrackState(TrackState::AtFirstHit);
        ced_hit(tsFirst->getReferencePoint()[0], tsFirst->getReferencePoint()[1], tsFirst->getReferencePoint()[2], 0, 12, 0x00e8ff);
        const TrackState* tsLast = track->getTrackState(TrackState::AtLastHit);
        ced_hit(tsLast->getReferencePoint()[0], tsLast->getReferencePoint()[1], tsLast->getReferencePoint()[2], 0, 12, 0x00e8ff);
        const TrackState* tsEcal = track->getTrackState(TrackState::AtCalorimeter);
        ced_hit(tsEcal->getReferencePoint()[0], tsEcal->getReferencePoint()[1], tsEcal->getReferencePoint()[2], 0, 12, 0x00e8ff);


        double bField = MarlinUtil::getBzAtOrigin();
        double pt = bField * 3e-4 / std::abs( tsIP->getOmega() );
        double charge = ( tsIP->getOmega() > 0. ?  1. : -1. );
        double Px = pt * std::cos(  tsIP->getPhi() ) ;
        double Py = pt * std::sin(  tsIP->getPhi() ) ;
        double Pz = pt * tsIP->getTanLambda() ;
        double Xs = tsIP->getReferencePoint()[0] -  tsIP->getD0() * std::sin( tsIP->getPhi() ) ;
        double Ys = tsIP->getReferencePoint()[1] +  tsIP->getD0() * std::cos( tsIP->getPhi() ) ;
        double Zs = tsIP->getReferencePoint()[2] +  tsIP->getZ0() ;
        //helix at IP
        DDMarlinCED::drawHelix(bField, charge, Xs, Ys, Zs, Px, Py, Pz, 1, 1,
                             0xff0000, 0., 2100., 3000., 0);


        pt = bField * 3e-4 / std::abs( tsEcal->getOmega() );
        charge = ( tsEcal->getOmega() > 0. ?  1. : -1. );
        Px = pt * std::cos(  tsEcal->getPhi() ) ;
        Py = pt * std::sin(  tsEcal->getPhi() ) ;
        Pz = pt * tsEcal->getTanLambda() ;
        Xs = tsEcal->getReferencePoint()[0] -  tsEcal->getD0() * std::sin( tsEcal->getPhi() ) ;
        Ys = tsEcal->getReferencePoint()[1] +  tsEcal->getD0() * std::cos( tsEcal->getPhi() ) ;
        Zs = tsEcal->getReferencePoint()[2] +  tsEcal->getZ0() ;
        //helix at ECal
        DDMarlinCED::drawHelix(-bField, -charge, Xs, Ys, Zs, Px, Py, Pz, 1, 1,
                             0xffb200, 0., 2100., 3000., 0);


    }


}
