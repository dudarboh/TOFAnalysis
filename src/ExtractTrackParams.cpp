#include "ExtractTrackParams.h"

// #include <iostream>
using std::cout, std::endl;

//TPC radii
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
using dd4hep::Detector, dd4hep::tesla;

#include "EVENT/LCCollection.h"
using EVENT::LCCollection;
#include "EVENT/ReconstructedParticle.h"
using EVENT::ReconstructedParticle;

#include "HelixClass.h"

#include <UTIL/LCRelationNavigator.h>

ExtractTrackParams aExtractTrackParams;

ExtractTrackParams::ExtractTrackParams() : Processor("ExtractTrackParams"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("TrackParams.root"));
    registerProcessorParameter(string("massAssumption"), string("Write track parameters of a fit that assumed this mass"), _massAssumption, string("Pion"));
}

ExtractTrackParams::~ExtractTrackParams(){
    delete _tree;
    delete _file;
}

void ExtractTrackParams::init(){
    _nEvt = 0;
    _start = system_clock::now();

    string fileName = _massAssumption + _outputFileName;
    _file = new TFile(fileName.c_str(), "RECREATE");

    string treeName = _massAssumption + "TrackParams";
    _tree = new TTree(treeName.c_str(), "Track fit parameters");

    _tree->Branch("chi2", &_chi2, "chi2/F");
    _tree->Branch("ndf", &_ndf, "ndf/I");
    _tree->Branch("dEdX", &_dEdX, "dEdX/F");

    _tree->Branch("p", &_p[0], "p/F");
    _tree->Branch("pt", &_pt[0], "pt/F");
    _tree->Branch("d0", &_d0[0], "d0/F");
    _tree->Branch("phi", &_phi[0], "phi/F");
    _tree->Branch("omega", &_omega[0], "omega/F");
    _tree->Branch("z0", &_z0[0], "z0/F");
    _tree->Branch("tanL", &_tanL[0], "tanL/F");

    //Track states
    _tree->Branch("pFirstHit", &_p[1], "pFirstHit/F");
    _tree->Branch("ptFirstHit", &_pt[1], "ptFirstHit/F");
    _tree->Branch("d0FirstHit", &_d0[1], "d0FirstHit/F");
    _tree->Branch("phiFirstHit", &_phi[1], "phiFirstHit/F");
    _tree->Branch("omegaFirstHit", &_omega[1], "omegaFirstHit/F");
    _tree->Branch("z0FirstHit", &_z0[1], "z0FirstHit/F");
    _tree->Branch("tanLFirstHit", &_tanL[1], "tanLFirstHit/F");
    _tree->Branch("xRefFirstHit", &_xRef[1], "xRefFirstHit/F");
    _tree->Branch("yRefFirstHit", &_yRef[1], "yRefFirstHit/F");
    _tree->Branch("zRefFirstHit", &_zRef[1], "zRefFirstHit/F");

    _tree->Branch("pLastHit", &_p[2], "pLastHit/F");
    _tree->Branch("ptLastHit", &_pt[2], "ptLastHit/F");
    _tree->Branch("d0LastHit", &_d0[2], "d0LastHit/F");
    _tree->Branch("phiLastHit", &_phi[2], "phiLastHit/F");
    _tree->Branch("omegaLastHit", &_omega[2], "omegaLastHit/F");
    _tree->Branch("z0LastHit", &_z0[2], "z0LastHit/F");
    _tree->Branch("tanLLastHit", &_tanL[2], "tanLLastHit/F");
    _tree->Branch("xRefLastHit", &_xRef[2], "xRefLastHit/F");
    _tree->Branch("yRefLastHit", &_yRef[2], "yRefLastHit/F");
    _tree->Branch("zRefLastHit", &_zRef[2], "zRefLastHit/F");

    _tree->Branch("pCalState", &_p[3], "pCalState/F");
    _tree->Branch("ptCalState", &_pt[3], "ptCalState/F");
    _tree->Branch("d0CalState", &_d0[3], "d0CalState/F");
    _tree->Branch("phiCalState", &_phi[3], "phiCalState/F");
    _tree->Branch("omegaCalState", &_omega[3], "omegaCalState/F");
    _tree->Branch("z0CalState", &_z0[3], "z0CalState/F");
    _tree->Branch("tanLCalState", &_tanL[3], "tanLCalState/F");
    _tree->Branch("xRefCalState", &_xRef[3], "xRefCalState/F");
    _tree->Branch("yRefCalState", &_yRef[3], "yRefCalState/F");
    _tree->Branch("zRefCalState", &_zRef[3], "zRefCalState/F");

    const Detector& detector = Detector::getInstance();
    double bTmpField[3];
    detector.field().magneticField({0., 0., 0.}, bTmpField);
    _bField = bTmpField[2]/tesla;
}

void ExtractTrackParams::processEvent(LCEvent* evt){
    ++_nEvt;
    if(_nEvt%10 == 0){
        double elapsedTime = duration<double>(system_clock::now() - _start).count();
        cout << "Event: "<<_nEvt<<"   Elapsed Time: "<<elapsedTime<<" sec     Avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    }

    LCCollection* colPFO = evt->getCollection("PandoraPFOs");

    LCCollection* colTrackPion = evt->getCollection("MarlinTrkTracks");

    LCCollection* colTrack = nullptr;
    if (_massAssumption == "Kaon") colTrack = evt->getCollection("MarlinTrkTracksKaon");
    else if (_massAssumption == "Proton") colTrack = evt->getCollection("MarlinTrkTracksProton");

    for (int i=0; i<colPFO->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        // Only simple cases of PFOs
        if( nClusters != 1 || nTracks > 1) continue;

        if (nTracks == 0){
            _chi2 = 0.;
            _ndf = 0;
            _dEdX = 0.;
            for (int j = 0; j < _nTrackStates; ++j){
                _d0[j] = 0.;
                _phi[j] = 0.;
                _omega[j] = 0.;
                _z0[j] = 0.;
                _tanL[j] = 0.;
                _xRef[j] = 0.;
                _yRef[j] = 0.;
                _zRef[j] = 0.;
                _p[j] = 0.;
                _pt[j] = 0.;
            }
            _tree->Fill();
            continue;
        }

        //Get track from PFO or collection with refitted mass assumption
        const Track* pfoTrack = pfo->getTracks()[0];

        int colTrackIdx = -1;
        if (_massAssumption == "Kaon" || _massAssumption == "Proton"){
            for (int j = 0; j < colTrackPion->getNumberOfElements(); ++j) {
                const Track* tmpTrack = dynamic_cast<Track*>(colTrackPion->getElementAt(j));
                if (pfoTrack == tmpTrack) {colTrackIdx = j; break;}
            }
        }

        const Track* track = nullptr;
        if (_massAssumption == "Pion") track = pfoTrack;
        else if (_massAssumption == "Kaon") track = dynamic_cast<Track*>(colTrack->getElementAt(colTrackIdx));
        else if (_massAssumption == "Proton") track = dynamic_cast<Track*>(colTrack->getElementAt(colTrackIdx));

        _chi2 = track->getChi2();
        _ndf = track->getNdf();
        _dEdX = track->getdEdx();

        for (int j = 0; j < _nTrackStates; ++j) {
            if(_ndf == 0){
                _d0[j] = 0.;
                _phi[j] = 0.;
                _omega[j] = 0.;
                _z0[j] = 0.;
                _tanL[j] = 0.;
                _xRef[j] = 0.;
                _yRef[j] = 0.;
                _zRef[j] = 0.;
                _p[j] = 0.;
                _pt[j] = 0.;
                continue;
            }
            const TrackState* ts = track->getTrackState(_trackStates[j]);
            _d0[j] = ts->getD0();
            _phi[j] = ts->getPhi();
            _omega[j] = ts->getOmega();
            _z0[j] = ts->getZ0();
            _tanL[j] = ts->getTanLambda();
            const float* pos = ts->getReferencePoint();
            _xRef[j] = pos[0];
            _yRef[j] = pos[1];
            _zRef[j] = pos[2];

            HelixClass helix;
            helix.Initialize_Canonical(_phi[j], _d0[j], _z0[j], _omega[j], _tanL[j], _bField);
            const float* mom = helix.getMomentum();
            _p[j] = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
            _pt[j] = sqrt(mom[0]*mom[0] + mom[1]*mom[1]);
        }
        _tree->Fill();
    } // end of PFOs loop
}

void ExtractTrackParams::end(){
    double elapsedTime = duration<double>(system_clock::now() - _start).count();
    cout<<"Finished writing Track parameters"<<endl;
    cout<<"Total number of events: "<<_nEvt<<endl;
    cout<<"Total elapsed time: "<<elapsedTime<<" sec"<<endl;
    cout<<"Total avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    _file->Write();
    cout<<_file->GetName()<<"   file is written in the current directory"<<endl;
}
