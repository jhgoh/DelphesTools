//#ifdef __CLING__
#ifdef __CINT__
R__LOAD_LIBRARY(libDelphes)
#endif
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TH2D.h"
#include <iostream>

void delphes2Image(const std::string finName, const std::string foutName)
{
  // Prepare output tree
  TFile* fout = TFile::Open(foutName.c_str(), "recreate");
  TTree* tree = new TTree("tree", "tree");

  unsigned short b_run = 1;
  unsigned int b_event = 0;
  float b_weight = 0;

  float b_MET_pt, b_MET_phi;

  const unsigned short Muon_N = 100;
  unsigned short b_nMuon;
  float b_Muon_pt[Muon_N], b_Muon_eta[Muon_N], b_Muon_phi[Muon_N], b_Muon_m[Muon_N];
  short b_Muon_q[Muon_N];
  float b_Muon_relIso[Muon_N];

  const unsigned short Electron_N = 100;
  unsigned short b_nElectron;
  float b_Electron_pt[Electron_N], b_Electron_eta[Electron_N], b_Electron_phi[Electron_N], b_Electron_m[Electron_N];
  short b_Electron_q[Electron_N];
  float b_Electron_relIso[Electron_N];

  const unsigned short Jet_N = 100;
  unsigned short b_nJet;
  float b_Jet_pt[Jet_N], b_Jet_eta[Jet_N], b_Jet_phi[Jet_N], b_Jet_m[Jet_N];
  short b_Jet_flav[Jet_N];
  float b_Jet_bTag[Jet_N];

  tree->Branch("run", &b_run, "run/s");
  tree->Branch("event", &b_event, "event/i");
  tree->Branch("weight", &b_weight, "weight/F");

  tree->Branch("MET_pt", &b_MET_pt, "MET_pt/F");
  tree->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  tree->Branch("nMuon", &b_nMuon, "nMuon/s");
  tree->Branch("Muon_pt", b_Muon_pt, "Muon_pt[nMuon]/F");
  tree->Branch("Muon_eta", b_Muon_eta, "Muon_eta[nMuon]/F");
  tree->Branch("Muon_phi", b_Muon_phi, "Muon_phi[nMuon]/F");
  tree->Branch("Muon_m", b_Muon_m, "Muon_m[nMuon]/F");
  tree->Branch("Muon_q", b_Muon_q, "Muon_q[nMuon]/S");
  tree->Branch("Muon_relIso", b_Muon_relIso, "Muon_relIso[nMuon]/F");

  tree->Branch("nElectron", &b_nElectron, "nElectron/s");
  tree->Branch("Electron_pt", b_Electron_pt, "Electron_pt[nElectron]/F");
  tree->Branch("Electron_eta", b_Electron_eta, "Electron_eta[nElectron]/F");
  tree->Branch("Electron_phi", b_Electron_phi, "Electron_phi[nElectron]/F");
  tree->Branch("Electron_m", b_Electron_m, "Electron_m[nElectron]/F");
  tree->Branch("Electron_q", b_Electron_q, "Electron_q[nElectron]/S");
  tree->Branch("Electron_relIso", b_Electron_relIso, "Electron_relIso[nElectron]/F");

  tree->Branch("nJet", &b_nJet, "nJet/s");
  tree->Branch("Jet_pt", b_Jet_pt, "Jet_pt[nJet]/F");
  tree->Branch("Jet_eta", b_Jet_eta, "Jet_eta[nJet]/F");
  tree->Branch("Jet_phi", b_Jet_phi, "Jet_phi[nJet]/F");
  tree->Branch("Jet_m", b_Jet_m, "Jet_m[nJet]/F");
  tree->Branch("Jet_flav", b_Jet_flav, "Jet_flav[nJet]/S");
  tree->Branch("Jet_bTag", b_Jet_bTag, "Jet_bTag[nJet]/F");

  const int nx = 64, ny = 64;
  const double minX = -3, maxX = 3, minY = -TMath::Pi(), maxY = TMath::Pi();
  TH2D hTrckPt("hTrckPt", "hTrckPt", nx, minX, maxX, ny, minY, maxY);
  TH2D hEcalPt("hEcalPt", "hEcalPt", nx, minX, maxX, ny, minY, maxY);
  TH2D hHcalPt("hHcalPt", "hHcalPt", nx, minX, maxX, ny, minY, maxY);
  TH2D hTrckN ("hTrckN" , "hTrckN" , nx, minX, maxX, ny, minY, maxY);
  TH2D hEcalN ("hEcalN" , "hEcalN" , nx, minX, maxX, ny, minY, maxY);
  TH2D hHcalN ("hHcalN" , "hHcalN" , nx, minX, maxX, ny, minY, maxY);

  const unsigned int b_NXY = nx*ny;
  unsigned int b_nXY = b_NXY;
  float b_hTrck_pt[b_NXY], b_hTrck_n[b_NXY];
  float b_hEcal_pt[b_NXY], b_hEcal_n[b_NXY];
  float b_hHcal_pt[b_NXY], b_hHcal_n[b_NXY];
  tree->Branch("nXY", &b_nXY, "nXY/i");
  tree->Branch("hTrck_pt", b_hTrck_pt, "hTrck_pt[nXY]/F");
  tree->Branch("hEcal_pt", b_hEcal_pt, "hEcal_pt[nXY]/F");
  tree->Branch("hHcal_pt", b_hHcal_pt, "hHcal_pt[nXY]/F");
  tree->Branch("hTrck_n" , b_hTrck_n , "hTrck_n[nXY]/F" );
  tree->Branch("hEcal_n" , b_hEcal_n , "hEcal_n[nXY]/F" );
  tree->Branch("hHcal_n" , b_hHcal_n , "hHcal_n[nXY]/F" );

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(finName.c_str());

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchGen = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    std::cout << entry << '/' << numberOfEntries << '\r';

    hTrckPt.Reset();
    hEcalPt.Reset();
    hHcalPt.Reset();
    hTrckN.Reset();
    hEcalN.Reset();
    hHcalN.Reset();

    const HepMCEvent* event = (const HepMCEvent*)branchEvent->At(0);
    b_event = event->Number;
    b_weight = event->Weight;

    if ( branchMET->GetEntries() > 0 ) {
      const MissingET* met = (const MissingET*)branchMET->At(0);
      b_MET_pt = met->MET;
      b_MET_phi = met->Phi;
    }
    else {
      b_MET_pt = b_MET_phi = 0;
    }

    b_nMuon = 0;
    for ( int i=0; i<branchMuon->GetEntries(); ++i ) {
      const Muon* muon = (const Muon*) branchMuon->At(i);
      const TLorentzVector p4 = muon->P4();

      b_Muon_pt[b_nMuon] = muon->PT;
      b_Muon_eta[b_nMuon] = muon->Eta;
      b_Muon_phi[b_nMuon] = muon->Phi;
      b_Muon_m[b_nMuon] = p4.M();
      b_Muon_q[b_nMuon] = muon->Charge;

      b_Muon_relIso[b_nMuon] = muon->IsolationVar;///muon->PT;

      ++b_nMuon;
      if ( b_nMuon >= Muon_N ) break;
    }

    b_nElectron = 0;
    for ( int i=0; i<branchElectron->GetEntries(); ++i ) {
      const Electron* electron = (const Electron*) branchElectron->At(i);
      const TLorentzVector p4 = electron->P4();

      b_Electron_pt[b_nElectron] = electron->PT;
      b_Electron_eta[b_nElectron] = electron->Eta;
      b_Electron_phi[b_nElectron] = electron->Phi;
      b_Electron_m[b_nElectron] = p4.M();
      b_Electron_q[b_nElectron] = electron->Charge;

      b_Electron_relIso[b_nElectron] = electron->IsolationVarRhoCorr;///electron->PT;

      ++b_nElectron;
      if ( b_nElectron >= Electron_N ) break;
    }

    b_nJet = 0;
    for ( int i=0; i<branchJet->GetEntries(); ++i ) {
      const Jet* jet = (const Jet*) branchJet->At(i);
      const TLorentzVector p4 = jet->P4();

      b_Jet_pt[b_nJet] = jet->PT;
      b_Jet_eta[b_nJet] = jet->Eta;
      b_Jet_phi[b_nJet] = jet->Phi;
      b_Jet_m[b_nJet] = p4.M();

      b_Jet_flav[b_nJet] = jet->Flavor;
      b_Jet_bTag[b_nJet] = jet->BTag;

      ++b_nJet;
      if ( b_nJet >= Jet_N ) break;
    }

    // Read tracks and towers to fill 2D images
    for ( int i=0; i<branchTrack->GetEntries(); ++i ) {
      const Track* track = (const Track*) branchTrack->At(i);

      hTrckPt.Fill(track->Eta, track->Phi, track->PT);
      hTrckN.Fill(track->Eta, track->Phi);
    }

    for ( int i=0; i<branchTower->GetEntries(); ++i ) {
      const Tower* tower = (const Tower*) branchTower->At(i);
      hEcalPt.Fill(tower->Eta, tower->Phi, tower->Eem);
      hHcalPt.Fill(tower->Eta, tower->Phi, tower->Ehad);
      if ( tower->Eem  > 1e-6 ) hEcalN.Fill(tower->Eta, tower->Phi);
      if ( tower->Ehad > 1e-6 ) hHcalN.Fill(tower->Eta, tower->Phi);
    }

    // Save 2D image into arrays
    for ( int j=0; j<ny; ++j ) {
      for ( int i=0; i<nx; ++i ) {
        b_hTrck_pt[i+j*nx] = hTrckPt.GetBinContent(i+1,j+1);
        b_hEcal_pt[i+j*nx] = hEcalPt.GetBinContent(i+1,j+1);
        b_hHcal_pt[i+j*nx] = hHcalPt.GetBinContent(i+1,j+1);
        b_hTrck_n[i+j*nx]  = hTrckN.GetBinContent(i+1,j+1);
        b_hEcal_n[i+j*nx]  = hEcalN.GetBinContent(i+1,j+1);
        b_hHcal_n[i+j*nx]  = hHcalN.GetBinContent(i+1,j+1);
      }
    }

    tree->Fill();
  }

  tree->Write();
  fout->Close();

}

