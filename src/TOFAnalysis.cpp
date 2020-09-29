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
    gInterpreter->GenerateDictionary("vector <vector<int> >", "vector");
    gInterpreter->GenerateDictionary("vector <vector<float> >", "vector");
    gInterpreter->GenerateDictionary("vector <vector<double> >", "vector");
    _nEvt = 0;
    _start = chrono::system_clock::now();

    _file = new TFile(_outputFileName.c_str(), "RECREATE");
    _tree = new TTree("ana_tree", "Tree description");

    //Record only PFO with 1 track and 1 cluster
    //PFO parameters
    _tree->Branch("nPFOs", &_nPFOs, "nPFOs/I");
    _tree->Branch("nGoodPFOs", &_nGoodPFOs, "nGoodPFOs/I");
    _tree->Branch("p", &_p);
    _tree->Branch("charge", &_charge);
    //Track parameters
    _tree->Branch("d0", &_d0);
    _tree->Branch("phi", &_phi);
    _tree->Branch("omega", &_omega);
    _tree->Branch("z0", &_z0);
    _tree->Branch("tanL", &_tanL);
    _tree->Branch("chi2", &_chi2);
    _tree->Branch("ndf", &_ndf);
    _tree->Branch("dEdX", &_dEdX);
    _tree->Branch("length", &_length);

    //Track states
    _tree->Branch("d0Last", &_d0Last);
    _tree->Branch("phiLast", &_phiLast);
    _tree->Branch("omegaLast", &_omegaLast);
    _tree->Branch("z0Last", &_z0Last);
    _tree->Branch("tanLLast", &_tanLLast);
    _tree->Branch("xRefLast", &_xRefLast);
    _tree->Branch("yRefLast", &_yRefLast);
    _tree->Branch("zRefLast", &_zRefLast);

    _tree->Branch("d0Calo", &_d0Calo);
    _tree->Branch("phiCalo", &_phiCalo);
    _tree->Branch("omegaCalo", &_omegaCalo);
    _tree->Branch("z0Calo", &_z0Calo);
    _tree->Branch("tanLCalo", &_tanLCalo);
    _tree->Branch("xRefCalo", &_xRefCalo);
    _tree->Branch("yRefCalo", &_yRefCalo);
    _tree->Branch("zRefCalo", &_zRefCalo);

    _tree->Branch("nHitsTrack", &_nHitsTrack);
    _tree->Branch("nHitsTPC", &_nHitsTPC);
    _tree->Branch("xHit", &_xHit);
    _tree->Branch("yHit", &_yHit);
    _tree->Branch("zHit", &_zHit);
    _tree->Branch("tHit", &_tHit);

    _tree->Branch("xCluster", &_xCluster);
    _tree->Branch("yCluster", &_yCluster);
    _tree->Branch("zCluster", &_zCluster);
    _tree->Branch("phiCluster", &_phiCluster);
    _tree->Branch("thetaCluster", &_thetaCluster);

    _tree->Branch("nHitsCluster", &_nHitsCluster);
    _tree->Branch("xHitCluster", &_xHitCluster);
    _tree->Branch("yHitCluster", &_yHitCluster);
    _tree->Branch("zHitCluster", &_zHitCluster);
    _tree->Branch("tHitCluster", &_tHitCluster);
    _tree->Branch("layerHitCluster", &_layerHitCluster);
    _tree->Branch("dToLineHitCluster", &_dToLineHitCluster);
    _tree->Branch("dToRefPointHitCluster", &_dToRefPointHitCluster);
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


    // Get ID of the algorithm and its parameters
    _nPFOs = col->getNumberOfElements();
    _nGoodPFOs = 0;
    for (int p=0; p<_nPFOs; ++p){
        ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(col->getElementAt(p));
        // Look only at PFOs with 1 cluster and 1 track
        if(pfo->getClusters().size() != 1 || pfo->getTracks().size() != 1) continue;
        ++_nGoodPFOs;

        _charge.push_back(pfo->getCharge());

        double pX = pfo->getMomentum()[0];
        double pY = pfo->getMomentum()[1];
        double pZ = pfo->getMomentum()[2];
        double momentum = sqrt(pX*pX + pY*pY + pZ*pZ);
        _p.push_back(momentum);

        const Track* track = pfo->getTracks()[0];
        _d0.push_back(track->getD0());
        _phi.push_back(track->getPhi());
        _omega.push_back(track->getOmega());
        _z0.push_back(track->getZ0());
        _tanL.push_back(track->getTanLambda());
        _chi2.push_back(track->getChi2());
        _ndf.push_back(track->getNdf());
        _dEdX.push_back(track->getdEdx());
        _nHitsTPC.push_back(track->getSubdetectorHitNumbers()[_TPCindex]);

        const TrackState* trackAtLast = track->getTrackState(TrackState::AtLastHit);
        _d0Last.push_back(trackAtLast->getD0());
        _phiLast.push_back(trackAtLast->getPhi());
        _omegaLast.push_back(trackAtLast->getOmega());
        _z0Last.push_back(trackAtLast->getZ0());
        _tanLLast.push_back(trackAtLast->getTanLambda());
        _xRefLast.push_back(trackAtLast->getReferencePoint()[0]);
        _yRefLast.push_back(trackAtLast->getReferencePoint()[1]);
        _zRefLast.push_back(trackAtLast->getReferencePoint()[2]);

        const TrackState* trackAtCalo = track->getTrackState(TrackState::AtCalorimeter);
        _d0Calo.push_back(trackAtCalo->getD0());
        _phiCalo.push_back(trackAtCalo->getPhi());
        _omegaCalo.push_back(trackAtCalo->getOmega());
        _z0Calo.push_back(trackAtCalo->getZ0());
        _tanLCalo.push_back(trackAtCalo->getTanLambda());
        _xRefCalo.push_back(trackAtCalo->getReferencePoint()[0]);
        _yRefCalo.push_back(trackAtCalo->getReferencePoint()[1]);
        _zRefCalo.push_back(trackAtCalo->getReferencePoint()[2]);

        float length = abs((trackAtCalo->getPhi() - track->getPhi())/track->getOmega())*sqrt(1. + track->getTanLambda()*track->getTanLambda());
        _length.push_back(length);


        const TrackerHitVec& hits = track->getTrackerHits();
        int nHits = hits.size();
        _nHitsTrack.push_back(nHits);
        vector <double> x;
        vector <double> y;
        vector <double> z;
        vector <float> t;
        for (int i = 0; i < nHits; ++i) {
            TrackerHit* hit = hits[i];
            x.push_back(hit->getPosition()[0]);
            y.push_back(hit->getPosition()[1]);
            z.push_back(hit->getPosition()[2]);
            t.push_back(hit->getTime());
        }
        _xHit.push_back(x);
        _yHit.push_back(y);
        _zHit.push_back(z);
        _tHit.push_back(t);
        // End of Tracker variables
        // Writing calorimeter hits
        Cluster* cluster = pfo->getClusters()[0];

        _xCluster.push_back(cluster->getPosition()[0]);
        _yCluster.push_back(cluster->getPosition()[1]);
        _zCluster.push_back(cluster->getPosition()[2]);
        _phiCluster.push_back(cluster->getIPhi());
        _thetaCluster.push_back(cluster->getITheta());

        const CalorimeterHitVec& clusterHits = cluster->getCalorimeterHits();
        nHits = clusterHits.size();

        int nHitsCluster = 0;
        vector <float> xHitCluster;
        vector <float> yHitCluster;
        vector <float> zHitCluster;
        vector <float> tHitCluster;
        vector <int> layerHitCluster;
        vector <float> dToLineHitCluster;
        vector <float> dToRefPointHitCluster;

        for (int i = 0; i < nHits; ++i){
            CalorimeterHit* hit = clusterHits[i];

            int hitType = hit->getType();
            bool isEcal = (CHT(hitType).caloID() == CHT::ecal);
            float time = hit->getTime();

            if (!(isEcal && (time>1.e-3))) continue;
            ++nHitsCluster;

            int layer = CHT(hitType).layer();
            layerHitCluster.push_back(layer);

            const float* hitPos = hit->getPosition();
            xHitCluster.push_back(hitPos[0]);
            yHitCluster.push_back(hitPos[1]);
            zHitCluster.push_back(hitPos[2]);
            tHitCluster.push_back(time);

            // Compute distance to point of Impact from the tracker (ref point)
            Vector3D caloHit = hitPos;
            Vector3D refPoint = trackAtCalo->getReferencePoint();
            Vector3D unitDir = Vector3D(1., trackAtCalo->getPhi(), atan(1. / trackAtCalo->getTanLambda()), Vector3D::spherical);
            dToRefPointHitCluster.push_back((caloHit - refPoint).r());
            dToLineHitCluster.push_back((caloHit - refPoint).cross(unitDir).r());
        }
        _nHitsCluster.push_back(nHitsCluster);
        _xHitCluster.push_back(xHitCluster);
        _yHitCluster.push_back(yHitCluster);
        _zHitCluster.push_back(zHitCluster);
        _tHitCluster.push_back(tHitCluster);
        _layerHitCluster.push_back(layerHitCluster);
        _dToLineHitCluster.push_back(dToLineHitCluster);
        _dToRefPointHitCluster.push_back(dToRefPointHitCluster);

    } // end of PFOs loop
    _tree->Fill();

    //Clear all the vectors before next event
    _p.clear();
    _charge.clear();

    _d0.clear();
    _phi.clear();
    _omega.clear();
    _z0.clear();
    _tanL.clear();
    _chi2.clear();
    _ndf.clear();
    _dEdX.clear();
    _length.clear();

    _d0Last.clear();
    _phiLast.clear();
    _omegaLast.clear();
    _z0Last.clear();
    _tanLLast.clear();
    _xRefLast.clear();
    _yRefLast.clear();
    _zRefLast.clear();

    _d0Calo.clear();
    _phiCalo.clear();
    _omegaCalo.clear();
    _z0Calo.clear();
    _tanLCalo.clear();
    _xRefCalo.clear();
    _yRefCalo.clear();
    _zRefCalo.clear();

    _nHitsTrack.clear();
    _nHitsTPC.clear();
    _xHit.clear();
    _yHit.clear();
    _zHit.clear();
    _tHit.clear();

    _xCluster.clear();
    _yCluster.clear();
    _zCluster.clear();
    _phiCluster.clear();
    _thetaCluster.clear();

    _nHitsCluster.clear();
    _xHitCluster.clear();
    _yHitCluster.clear();
    _zHitCluster.clear();
    _tHitCluster.clear();
    _layerHitCluster.clear();
    _dToLineHitCluster.clear();
    _dToRefPointHitCluster.clear();
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
