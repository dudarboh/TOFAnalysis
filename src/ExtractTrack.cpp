#include "ExtractTrack.hpp"

//TPC radii
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
using dd4hep::Detector, dd4hep::tesla;

#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "HelixClass.h"
using EVENT::LCCollection;
using EVENT::ReconstructedParticle;
using std::cout, std::endl;


ExtractTrack aExtractTrack;

ExtractTrack::ExtractTrack() : Processor("ExtractTrack"){
    registerProcessorParameter(string("outputFile"), string("Name of the output root file"), _outputFileName, string("Track.root"));
    registerProcessorParameter(string("massAssumption"), string("Write track parameters of a fit that assumed this mass"), _massAssumption, string("Pion"));
}

void ExtractTrack::init(){
    _nEvt = 0;
    _start = system_clock::now();

    string fileName = _massAssumption + _outputFileName;
    _file.reset( new TFile(fileName.c_str(), "RECREATE") );
    string treeName = _massAssumption + "Track";
    _tree.reset( new TTree(treeName.c_str(), "Track fit parameters") );

    _tree->Branch("chi2", &_chi2);
    _tree->Branch("ndf", &_ndf);
    _tree->Branch("dEdX", &_dEdX);

    //Track states
    _tree->Branch("pIP", &_p[0]);
    _tree->Branch("refIP", &_ref[0]);
    _tree->Branch("d0IP", &_d0[0]);
    _tree->Branch("phiIP", &_phi[0]);
    _tree->Branch("omegaIP", &_omega[0]);
    _tree->Branch("z0IP", &_z0[0]);
    _tree->Branch("tanLIP", &_tanL[0]);

    _tree->Branch("pFirst", &_p[1]);
    _tree->Branch("refFirst", &_ref[1]);
    _tree->Branch("d0First", &_d0[1]);
    _tree->Branch("phiFirst", &_phi[1]);
    _tree->Branch("omegaFirst", &_omega[1]);
    _tree->Branch("z0First", &_z0[1]);
    _tree->Branch("tanLFirst", &_tanL[1]);

    _tree->Branch("pLast", &_p[2]);
    _tree->Branch("refLast", &_ref[2]);
    _tree->Branch("d0Last", &_d0[2]);
    _tree->Branch("phiLast", &_phi[2]);
    _tree->Branch("omegaLast", &_omega[2]);
    _tree->Branch("z0Last", &_z0[2]);
    _tree->Branch("tanLLast", &_tanL[2]);

    _tree->Branch("pCalo", &_p[3]);
    _tree->Branch("refCalo", &_ref[3]);
    _tree->Branch("d0Calo", &_d0[3]);
    _tree->Branch("phiCalo", &_phi[3]);
    _tree->Branch("omegaCalo", &_omega[3]);
    _tree->Branch("z0Calo", &_z0[3]);
    _tree->Branch("tanLCalo", &_tanL[3]);

    const Detector& detector = Detector::getInstance();
    double bTmpField[3];
    detector.field().magneticField({0., 0., 0.}, bTmpField);
    _bField = bTmpField[2]/tesla;
}

void ExtractTrack::processEvent(LCEvent* evt){
    ++_nEvt;

    LCCollection* colPFO = evt->getCollection("PandoraPFOs");
    LCCollection* colTrackPion = evt->getCollection("MarlinTrkTracks");

    string colTrackName = "MarlinTrkTracks";
    if (_massAssumption == "Kaon" || _massAssumption == "Proton") colTrackName = Form("MarlinTrkTracks%s", _massAssumption.c_str());
    LCCollection* colTrack = evt->getCollection(colTrackName);

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
            }
            _tree->Fill();
            continue;
        }

        //Get track from PFO or collection with refitted mass assumption
        const Track* track = pfo->getTracks()[0];

        if (_massAssumption == "Kaon" || _massAssumption == "Proton"){
            for (int j = 0; j < colTrackPion->getNumberOfElements(); ++j) {
                const Track* tmpTrack = dynamic_cast<Track*>(colTrackPion->getElementAt(j));
                if (track == tmpTrack) {
                    track = dynamic_cast<Track*>(colTrack->getElementAt(j));
                    break;
                }
            }
        }

        _chi2 = track->getChi2();
        _ndf = track->getNdf();
        _dEdX = track->getdEdx();

        //God bless
        if(_ndf == 0){
            for (int j = 0; j < _nTrackStates; ++j){
                _d0[j] = 0.;
                _phi[j] = 0.;
                _omega[j] = 0.;
                _z0[j] = 0.;
                _tanL[j] = 0.;
            }
            _tree->Fill();
            continue;
        }

        for (int j = 0; j < _nTrackStates; ++j) {
            const TrackState* ts = track->getTrackState(_trackStates[j]);
            _d0[j] = ts->getD0();
            _phi[j] = ts->getPhi();
            _omega[j] = ts->getOmega();
            _z0[j] = ts->getZ0();
            _tanL[j] = ts->getTanLambda();
            const float* pos = ts->getReferencePoint();
            _ref[j] = XYZPoint(pos[0], pos[1], pos[2]);

            HelixClass helix;
            helix.Initialize_Canonical(_phi[j], _d0[j], _z0[j], _omega[j], _tanL[j], _bField);
            const float* mom = helix.getMomentum();
            _p[j] = XYZVector(mom[0], mom[1], mom[2]);
        }
        _tree->Fill();
    } // end of PFOs loop
}

void ExtractTrack::end(){
    double elapsedTime = duration<double>(system_clock::now() - _start).count();
    cout<<"Finished writing Track parameters"<<endl;
    cout<<"Total number of events: "<<_nEvt<<endl;
    cout<<"Total elapsed time: "<<elapsedTime<<" sec"<<endl;
    cout<<"Total avg speed: "<<_nEvt/elapsedTime<<" evt/sec"<<endl;
    _file->Write();
    cout<<_file->GetName()<<"   file is written in the current directory"<<endl;
}
