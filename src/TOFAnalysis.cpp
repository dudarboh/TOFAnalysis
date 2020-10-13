#include "TOFAnalysis.h"

TOFAnalysis aTOFAnalysis;

TOFAnalysis::TOFAnalysis() : Processor("myTOFAnalysis"),
_TPCindex((ILDDetID::TPC)*2-2){
    registerProcessorParameter(string("output_filename"), string("Name of the output root file"), _outputFileName, string("TOFAnalysis.root"));
}

TOFAnalysis::~TOFAnalysis(){
    delete _tree;
    delete _file;
}

void TOFAnalysis::init(){
    _nEvt = 0;
    _start = chrono::system_clock::now();

    _file = new TFile(_outputFileName.c_str(), "RECREATE");
    _tree = new TTree("ana_tree", "Tree description");

    //Record only PFO with 1 track and 1 cluster
    //PFO parameters
    _tree->Branch("p", &_p, "p/D");
    _tree->Branch("charge", &_charge, "charge/F");
    //Track parameters
    _tree->Branch("d0", &_d0, "d0/F");
    _tree->Branch("phi", &_phi, "phi/F");
    _tree->Branch("omega", &_omega, "omega/F");
    _tree->Branch("z0", &_z0, "z0/F");
    _tree->Branch("tanL", &_tanL, "tanL/F");
    _tree->Branch("chi2", &_chi2, "chi2/F");
    _tree->Branch("ndf", &_ndf, "ndf/I");
    _tree->Branch("dEdX", &_dEdX, "dEdX/F");

    //Track states
    _tree->Branch("pFirstState", &_pFirstState, "pFirstState/F");
    _tree->Branch("d0FirstState", &_d0FirstState, "d0FirstState/F");
    _tree->Branch("phiFirstState", &_phiFirstState, "phiFirstState/F");
    _tree->Branch("omegaFirstState", &_omegaFirstState, "omegaFirstState/F");
    _tree->Branch("z0FirstState", &_z0FirstState, "z0FirstState/F");
    _tree->Branch("tanLFirstState", &_tanLFirstState, "tanLFirstState/F");
    _tree->Branch("xRefFirstState", &_xRefFirstState, "xRefFirstState/F");
    _tree->Branch("yRefFirstState", &_yRefFirstState, "yRefFirstState/F");
    _tree->Branch("zRefFirstState", &_zRefFirstState, "zRefFirstState/F");

    _tree->Branch("pLastState", &_pLastState, "pLastState/F");
    _tree->Branch("d0LastState", &_d0LastState, "d0LastState/F");
    _tree->Branch("phiLastState", &_phiLastState, "phiLastState/F");
    _tree->Branch("omegaLastState", &_omegaLastState, "omegaLastState/F");
    _tree->Branch("z0LastState", &_z0LastState, "z0LastState/F");
    _tree->Branch("tanLLastState", &_tanLLastState, "tanLLastState/F");
    _tree->Branch("xRefLastState", &_xRefLastState, "xRefLastState/F");
    _tree->Branch("yRefLastState", &_yRefLastState, "yRefLastState/F");
    _tree->Branch("zRefLastState", &_zRefLastState, "zRefLastState/F");

    _tree->Branch("pCalState", &_pCalState, "pCalState/F");
    _tree->Branch("d0CalState", &_d0CalState, "d0CalState/F");
    _tree->Branch("phiCalState", &_phiCalState, "phiCalState/F");
    _tree->Branch("omegaCalState", &_omegaCalState, "omegaCalState/F");
    _tree->Branch("z0CalState", &_z0CalState, "z0CalState/F");
    _tree->Branch("tanLCalState", &_tanLCalState, "tanLCalState/F");
    _tree->Branch("xRefCalState", &_xRefCalState, "xRefCalState/F");
    _tree->Branch("yRefCalState", &_yRefCalState, "yRefCalState/F");
    _tree->Branch("zRefCalState", &_zRefCalState, "zRefCalState/F");

    _tree->Branch("nTrHits", &_nTrHits, "nTrHits/I");
    _tree->Branch("nTPCHits", &_nTPCHits, "nTPCHits/I");
    _tree->Branch("xTrHit", &_xTrHit);
    _tree->Branch("yTrHit", &_yTrHit);
    _tree->Branch("zTrHit", &_zTrHit);
    _tree->Branch("tTrHit", &_tTrHit);

    _tree->Branch("nCalHits", &_nCalHits, "nCalHits/I");
    _tree->Branch("xCalHit", &_xCalHit);
    _tree->Branch("yCalHit", &_yCalHit);
    _tree->Branch("zCalHit", &_zCalHit);
    _tree->Branch("tCalHit", &_tCalHit);
    _tree->Branch("layerCalHit", &_layerCalHit);
    _tree->Branch("dToLineCalHit", &_dToLineCalHit);
    _tree->Branch("dToRefPointCalHit", &_dToRefPointCalHit);
}

void TOFAnalysis::processEvent(LCEvent* evt){
    //Pring status
    ++_nEvt;
    if(_nEvt%10 == 0){
        double elapsedTime = chrono::duration<double>(chrono::system_clock::now() - _start).count();
        cout<<"Event: "<<_nEvt<<"   Elapsed Time: "<<elapsedTime<<" sec     Avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    }

    //Get collection of PFOs for this event
    LCCollection* col = nullptr;
    try{
        col = evt->getCollection("PandoraPFOs");
    }
    catch (...){
        cout<<"Event "<<_nEvt<<" has no PandoarPFOs collection. Skip event"<<endl;
        return;
    }

    for (int p=0; p<col->getNumberOfElements(); ++p){
        const ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col->getElementAt(p));
        // Look only at PFOs with 1 cluster and 1 track
        if(pfo->getClusters().size() != 1 || pfo->getTracks().size() != 1) continue;
        _charge = pfo->getCharge();
        const double* mom = pfo->getMomentum();
        _p = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);

        //Track parameters
        const Track* track = pfo->getTracks()[0];
        _d0 = track->getD0();
        _phi = track->getPhi();
        _omega = track->getOmega();
        _z0 = track->getZ0();
        _tanL = track->getTanLambda();
        _chi2 = track->getChi2();
        _ndf = track->getNdf();
        _dEdX = track->getdEdx();
        _nTPCHits = track->getSubdetectorHitNumbers()[_TPCindex];

        //First hit state
        const TrackState* trackFirstState = track->getTrackState(TrackState::AtFirstHit);
        _d0FirstState = trackFirstState->getD0();
        _phiFirstState = trackFirstState->getPhi();
        _omegaFirstState = trackFirstState->getOmega();
        _z0FirstState = trackFirstState->getZ0();
        _tanLFirstState = trackFirstState->getTanLambda();

        HelixClass helixFirst;
        helixFirst.Initialize_Canonical(_phiFirstState, _d0FirstState, _z0FirstState, _omegaFirstState, _tanLFirstState, 3.5);
        float* momHelixFirst = helixFirst.getMomentum();
        _pFirstState = sqrt(momHelixFirst[0]*momHelixFirst[0] + momHelixFirst[1]*momHelixFirst[1] + momHelixFirst[2]*momHelixFirst[2]);

        const float* firstPos = trackFirstState->getReferencePoint();
        _xRefFirstState = firstPos[0];
        _yRefFirstState = firstPos[1];
        _zRefFirstState = firstPos[2];

        //Last hit state
        const TrackState* trackLastState = track->getTrackState(TrackState::AtLastHit);
        _d0LastState = trackLastState->getD0();
        _phiLastState = trackLastState->getPhi();
        _omegaLastState = trackLastState->getOmega();
        _z0LastState = trackLastState->getZ0();
        _tanLLastState = trackLastState->getTanLambda();

        HelixClass helixLast;
        helixLast.Initialize_Canonical(_phiLastState, _d0LastState, _z0LastState, _omegaLastState, _tanLLastState, 3.5);
        float* momHelixLast = helixLast.getMomentum();
        _pLastState = sqrt(momHelixLast[0]*momHelixLast[0] + momHelixLast[1]*momHelixLast[1] + momHelixLast[2]*momHelixLast[2]);

        const float* lastPos = trackLastState->getReferencePoint();
        _xRefLastState = lastPos[0];
        _yRefLastState = lastPos[1];
        _zRefLastState = lastPos[2];

        //Calo hit state
        const TrackState* trackCalState = track->getTrackState(TrackState::AtCalorimeter);
        _d0CalState = trackCalState->getD0();
        _phiCalState = trackCalState->getPhi();
        _omegaCalState = trackCalState->getOmega();
        _z0CalState = trackCalState->getZ0();
        _tanLCalState = trackCalState->getTanLambda();

        HelixClass helixCal;
        helixCal.Initialize_Canonical(_phiCalState, _d0CalState, _z0CalState, _omegaCalState, _tanLCalState, 3.5);
        float* momHelixCal = helixCal.getMomentum();
        _pCalState = sqrt(momHelixCal[0]*momHelixCal[0] + momHelixCal[1]*momHelixCal[1] + momHelixCal[2]*momHelixCal[2]);

        const float* calPos = trackCalState->getReferencePoint();
        _xRefCalState = calPos[0];
        _yRefCalState = calPos[1];
        _zRefCalState = calPos[2];

        //Track hits
        const TrackerHitVec& trHits = track->getTrackerHits();
        _nTrHits = trHits.size();
        for (int i = 0; i < _nTrHits; ++i) {
            const TrackerHit* hit = trHits[i];
            const double* trHitPos = hit->getPosition();
            _xTrHit.push_back(trHitPos[0]);
            _yTrHit.push_back(trHitPos[1]);
            _zTrHit.push_back(trHitPos[2]);
            _tTrHit.push_back(hit->getTime());
        }
        // End of Tracker variables
        // Writing calorimeter hits

        const Cluster* cluster = pfo->getClusters()[0];
        const CalorimeterHitVec& calHits = cluster->getCalorimeterHits();
        //Count only ECAl hits with time information
        _nCalHits = 0;
        for (size_t i = 0; i < calHits.size(); ++i){
            const CalorimeterHit* hit = calHits[i];
            int hitType = hit->getType();
            bool isEcal = (CHT(hitType).caloID() == CHT::ecal);
            float time = hit->getTime();

            if (!(isEcal && (time>1.e-3))) continue;
            ++_nCalHits;

            _layerCalHit.push_back(CHT(hitType).layer());

            const float* calHitPos = hit->getPosition();
            _xCalHit.push_back(calHitPos[0]);
            _yCalHit.push_back(calHitPos[1]);
            _zCalHit.push_back(calHitPos[2]);
            _tCalHit.push_back(time);

            // Compute distance to point of Impact from the tracker (ref point)
            Vector3D caloHit(calHitPos);
            Vector3D refPoint(calPos);
            Vector3D unitDir(1., _phiCalState, atan(1. / _tanLCalState), Vector3D::spherical);
            _dToRefPointCalHit.push_back((caloHit - refPoint).r());
            _dToLineCalHit.push_back((caloHit - refPoint).cross(unitDir).r());
        }
        _tree->Fill();

        //Clear all the vectors before next PFO
        _xTrHit.clear();
        _yTrHit.clear();
        _zTrHit.clear();
        _tTrHit.clear();
        _xCalHit.clear();
        _yCalHit.clear();
        _zCalHit.clear();
        _tCalHit.clear();
        _layerCalHit.clear();
        _dToLineCalHit.clear();
        _dToRefPointCalHit.clear();
    } // end of PFOs loop
}

void TOFAnalysis::end(){
    double elapsedTime = chrono::duration<double>(chrono::system_clock::now() - _start).count();
    cout<<"Finished"<<endl;
    cout<<"Total number of events: "<<_nEvt<<endl;
    cout<<"Total elapsed time: "<<elapsedTime<<" sec"<<endl;
    cout<<"Total avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    _file->Write();
    cout<<_file->GetName()<<"   file is written in the current directory"<<endl;
}
