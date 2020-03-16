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
#include <iostream>

//------------------------------------------------------------------------------
// Global options to switch on/off output branches
const bool dojetDau = true;
const bool doPrintDebug = false;

//------------------------------------------------------------------------------
int getLast(TClonesArray* branch, const int iGen)
{
  const GenParticle* p = (const GenParticle*)branch->At(iGen);
  if ( p->D1 == -1 or p->D2 == -1 ) return iGen;

  for ( int i=p->D1; i<=p->D2; ++i ) {
    const GenParticle* dau = (const GenParticle*)branch->At(i);
    if ( p->PID == dau->PID ) return getLast(branch, i);
  }
  return iGen;
}

void findPartonAncestors(const GenParticle* p, const TClonesArray* branchGen,
                         const std::map<const GenParticle*, int>& dauPtrToIdx,
                         std::set<const GenParticle*>& matches)
{
  if ( p == nullptr or p->M1 == -1 ) return;

  for ( int i=p->M1, n=std::max(p->M1, p->M2); i<=n; ++i ) {
    const GenParticle* mo = (const GenParticle*)branchGen->At(i);
    auto match = dauPtrToIdx.find(mo);

    if ( match != dauPtrToIdx.end() ) matches.insert(match->first);
    else findPartonAncestors(mo, branchGen, dauPtrToIdx, matches);
  }
}

void delphes2FlatTuple(const std::string finName, const std::string foutName)
{
  // Prepare output tree
  TFile* fout = TFile::Open(foutName.c_str(), "recreate");
  TTree* tree = new TTree("Events", "Events");

  unsigned int b_run = 1;
  unsigned long b_event = 0;
  float b_weight = 0;

  float b_MET_pt, b_MET_phi;

  const unsigned int Muon_N = 100;
  unsigned int b_nMuon;
  float b_Muon_pt[Muon_N], b_Muon_eta[Muon_N], b_Muon_phi[Muon_N], b_Muon_m[Muon_N];
  int b_Muon_q[Muon_N];
  float b_Muon_relIso[Muon_N];

  const unsigned int Electron_N = 100;
  unsigned int b_nElectron;
  float b_Electron_pt[Electron_N], b_Electron_eta[Electron_N], b_Electron_phi[Electron_N], b_Electron_m[Electron_N];
  int b_Electron_q[Electron_N];
  float b_Electron_relIso[Electron_N];

  const unsigned int Jet_N = 100;
  unsigned int b_nJet;
  float b_Jet_pt[Jet_N], b_Jet_eta[Jet_N], b_Jet_phi[Jet_N], b_Jet_m[Jet_N];
  int b_Jet_flav[Jet_N];
  float b_Jet_bTag[Jet_N];
  int b_Jet_partonIdx[Jet_N];

  const unsigned int jetDau_N = 10000;
  unsigned int b_njetDau;
  float b_jetDau_pt[jetDau_N], b_jetDau_eta[jetDau_N], b_jetDau_phi[jetDau_N];
  int b_jetDau_q[jetDau_N], b_jetDau_pdgId[jetDau_N];
  unsigned int b_jetDau_jetIdx[jetDau_N];

  const unsigned int GenParton_N = 1000;
  unsigned int b_nGenParton;
  float b_GenParton_pt[GenParton_N], b_GenParton_eta[GenParton_N], b_GenParton_phi[GenParton_N], b_GenParton_m[GenParton_N];
  int b_GenParton_pdgId[GenParton_N], b_GenParton_q3[GenParton_N];
  int b_GenParton_mother[GenParton_N], b_GenParton_dau1[GenParton_N], b_GenParton_dau2[GenParton_N];

  const unsigned int GenJet_N = 100;
  unsigned int b_nGenJet;
  float b_GenJet_pt[GenJet_N], b_GenJet_eta[GenJet_N], b_GenJet_phi[GenJet_N], b_GenJet_m[GenJet_N];
  int b_GenJet_partonIdx[GenJet_N];

  tree->Branch("run", &b_run, "run/i");
  tree->Branch("event", &b_event, "event/l");
  tree->Branch("weight", &b_weight, "weight/F");

  tree->Branch("MET_pt", &b_MET_pt, "MET_pt/F");
  tree->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  tree->Branch("nMuon", &b_nMuon, "nMuon/i");
  tree->Branch("Muon_pt", b_Muon_pt, "Muon_pt[nMuon]/F");
  tree->Branch("Muon_eta", b_Muon_eta, "Muon_eta[nMuon]/F");
  tree->Branch("Muon_phi", b_Muon_phi, "Muon_phi[nMuon]/F");
  tree->Branch("Muon_m", b_Muon_m, "Muon_m[nMuon]/F");
  tree->Branch("Muon_q", b_Muon_q, "Muon_q[nMuon]/I");
  tree->Branch("Muon_relIso", b_Muon_relIso, "Muon_relIso[nMuon]/F");

  tree->Branch("nElectron", &b_nElectron, "nElectron/i");
  tree->Branch("Electron_pt", b_Electron_pt, "Electron_pt[nElectron]/F");
  tree->Branch("Electron_eta", b_Electron_eta, "Electron_eta[nElectron]/F");
  tree->Branch("Electron_phi", b_Electron_phi, "Electron_phi[nElectron]/F");
  tree->Branch("Electron_m", b_Electron_m, "Electron_m[nElectron]/F");
  tree->Branch("Electron_q", b_Electron_q, "Electron_q[nElectron]/I");
  tree->Branch("Electron_relIso", b_Electron_relIso, "Electron_relIso[nElectron]/F");

  tree->Branch("nJet", &b_nJet, "nJet/i");
  tree->Branch("Jet_pt", b_Jet_pt, "Jet_pt[nJet]/F");
  tree->Branch("Jet_eta", b_Jet_eta, "Jet_eta[nJet]/F");
  tree->Branch("Jet_phi", b_Jet_phi, "Jet_phi[nJet]/F");
  tree->Branch("Jet_m", b_Jet_m, "Jet_m[nJet]/F");
  tree->Branch("Jet_flav", b_Jet_flav, "Jet_flav[nJet]/I");
  tree->Branch("Jet_bTag", b_Jet_bTag, "Jet_bTag[nJet]/F");
  tree->Branch("Jet_partonIdx", b_Jet_partonIdx, "Jet_partonIdx[nJet]/I");

  if ( dojetDau ) {
    tree->Branch("njetDau", &b_njetDau, "njetDau/i");
    tree->Branch("jetDau_pt", b_jetDau_pt, "jetDau_pt[njetDau]/F");
    tree->Branch("jetDau_eta", b_jetDau_eta, "jetDau_eta[njetDau]/F");
    tree->Branch("jetDau_phi", b_jetDau_phi, "jetDau_phi[njetDau]/F");
    tree->Branch("jetDau_q", b_jetDau_q, "jetDau_q[njetDau]/I");
    tree->Branch("jetDau_pdgId", b_jetDau_pdgId, "jetDau_pdgId[njetDau]/I");
    tree->Branch("jetDau_jetIdx", b_jetDau_jetIdx, "jetDau_jetIdx[njetDau]/I");
  }

  tree->Branch("nGenParton", &b_nGenParton, "nGenParton/i");
  tree->Branch("GenParton_pt", b_GenParton_pt, "GenParton_pt[nGenParton]/F");
  tree->Branch("GenParton_eta", b_GenParton_eta, "GenParton_eta[nGenParton]/F");
  tree->Branch("GenParton_phi", b_GenParton_phi, "GenParton_phi[nGenParton]/F");
  tree->Branch("GenParton_m", b_GenParton_m, "GenParton_m[nGenParton]/F");
  tree->Branch("GenParton_pdgId", b_GenParton_pdgId, "GenParton_pdgId[nGenParton]/I");
  tree->Branch("GenParton_q3", b_GenParton_q3, "GenParton_q3[nGenParton]/I");
  tree->Branch("GenParton_mother", b_GenParton_mother, "GenParton_mother[nGenParton]/I");
  tree->Branch("GenParton_dau1", b_GenParton_dau1, "GenParton_dau1[nGenParton]/I");
  tree->Branch("GenParton_dau2", b_GenParton_dau2, "GenParton_dau2[nGenParton]/I");

  tree->Branch("nGenJet", &b_nGenJet, "nGenJet/i");
  tree->Branch("GenJet_pt", b_GenJet_pt, "GenJet_pt[nGenJet]/F");
  tree->Branch("GenJet_eta", b_GenJet_eta, "GenJet_eta[nGenJet]/F");
  tree->Branch("GenJet_phi", b_GenJet_phi, "GenJet_phi[nGenJet]/F");
  tree->Branch("GenJet_m", b_GenJet_m, "GenJet_m[nGenJet]/F");
  tree->Branch("GenJet_partonIdx", b_GenJet_partonIdx, "GenJet_partonIdx[nGenJet]/I");

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
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    std::cout << entry << '/' << numberOfEntries << '\r';

    const HepMCEvent* event = (const HepMCEvent*)branchEvent->At(0);
    b_event = event->Number;
    b_weight = event->Weight;

    // Build gen particle collection, for the ttbar decays
    std::vector<std::vector<int> > GenParton_topDaus;
    b_nGenParton = 0;
    if ( branchGen ) {
      for ( int i=0; i<branchGen->GetEntries(); ++i ) {
        const GenParticle* p = (const GenParticle*)branchGen->At(i);
        const int pid = p->PID, absId = abs(p->PID);
        if ( absId != 6 ) continue; // top only.
        if ( p->D1 == -1 or p->D2 == -1 ) continue; // should have valid daughters

        bool isDupl = false;
        for ( int j=p->D1; j<=p->D2; ++j ) {
          const GenParticle* dau  = (const GenParticle*)branchGen->At(j);
          if ( pid == dau->PID ) { isDupl = true; break; }
        }
        if ( isDupl ) continue;

        GenParton_topDaus.emplace_back();
        for ( int j=p->D1; j<=p->D2; ++j ) GenParton_topDaus.back().push_back(j);

        // Fill top quarks
        b_GenParton_pt[b_nGenParton] = p->PT;
        b_GenParton_eta[b_nGenParton] = p->Eta;
        b_GenParton_phi[b_nGenParton] = p->Phi;
        b_GenParton_m[b_nGenParton] = p->Mass;
        b_GenParton_pdgId[b_nGenParton] = p->PID;
        b_GenParton_q3[b_nGenParton] = p->Charge*3;
        b_GenParton_dau1[b_nGenParton] = b_GenParton_dau2[b_nGenParton] = -1;
        b_GenParton_mother[b_nGenParton] = -1;

        ++b_nGenParton;
        if ( b_nGenParton >= GenParton_N ) break;
      }
    }

    std::map<const GenParticle*, int> dauPtrToIdx;
    for ( int i=0, n=GenParton_topDaus.size(); i<n; ++i ) {
      const auto& dauIdxs = GenParton_topDaus.at(i);
      if ( dauIdxs.empty() ) continue; // no daughters

      b_GenParton_dau1[i] = b_nGenParton;
      b_GenParton_dau2[i] = b_nGenParton+dauIdxs.size()-1;

      for ( auto j : dauIdxs ) {
        const GenParticle* dau = (const GenParticle*)branchGen->At(j);

        // Fill top quark daughters
        b_GenParton_pt[b_nGenParton] = dau->PT;
        b_GenParton_eta[b_nGenParton] = dau->Eta;
        b_GenParton_phi[b_nGenParton] = dau->Phi;
        b_GenParton_m[b_nGenParton] = dau->Mass;
        b_GenParton_pdgId[b_nGenParton] = dau->PID;
        b_GenParton_q3[b_nGenParton] = dau->Charge*3;
        b_GenParton_dau1[b_nGenParton] = b_GenParton_dau2[b_nGenParton] = -1;
        b_GenParton_mother[b_nGenParton] = i;

        const int iDau = b_nGenParton;
        dauPtrToIdx[dau] = b_nGenParton;
        ++b_nGenParton;
        if ( b_nGenParton >= GenParton_N ) break;

/*
        if ( abs(dau->PID) < 23 or abs(dau->PID) > 25 ) continue;
        if ( dau->D1 == -1 or dau->D2 == -1 ) continue;
        dau = (const GenParticle*)branchGen->At(getLast(branchGen, j));

        int ngdau = 0;
        for ( int k=dau->D1; k<=dau->D2; ++k ) {
          const GenParticle* gdau = (const GenParticle*)branchGen->At(k);

          // Fill W/Z/H daughters
          b_GenParton_pt[b_nGenParton] = gdau->PT;
          b_GenParton_eta[b_nGenParton] = gdau->Eta;
          b_GenParton_phi[b_nGenParton] = gdau->Phi;
          b_GenParton_m[b_nGenParton] = gdau->Mass;
          b_GenParton_pdgId[b_nGenParton] = gdau->PID;
          b_GenParton_q3[b_nGenParton] = gdau->Charge*3;
          b_GenParton_dau1[b_nGenParton] = b_GenParton_dau2[b_nGenParton] = -1;
          b_GenParton_mother[b_nGenParton] = iDau;

          ++ngdau;
          dauPtrToIdx[gdau] = b_nGenParton;
          ++b_nGenParton;
          if ( b_nGenParton >= GenParton_N ) break;

          const int iGdau = b_nGenParton;
          // For the tau decays
          if ( abs(gdau->PID) != 15 ) continue;
          if ( gdau->D1 == -1 or gdau->D2 == -1 ) continue;
          gdau = (const GenParticle*)branchGen->At(getLast(branchGen, k));

          int nggdau = 0;
          for ( int l=gdau->D1; l<=gdau->D2; ++l ) {
            const GenParticle* ggdau = (const GenParticle*)branchGen->At(l);
            const int absId = abs(ggdau->PID);
            if ( absId < 11 or absId >= 15 ) continue;

            // Fill W/Z/H daughters
            b_GenParton_pt[b_nGenParton] = ggdau->PT;
            b_GenParton_eta[b_nGenParton] = ggdau->Eta;
            b_GenParton_phi[b_nGenParton] = ggdau->Phi;
            b_GenParton_m[b_nGenParton] = ggdau->Mass;
            b_GenParton_pdgId[b_nGenParton] = ggdau->PID;
            b_GenParton_q3[b_nGenParton] = ggdau->Charge*3;
            b_GenParton_dau1[b_nGenParton] = b_GenParton_dau2[b_nGenParton] = -1;
            b_GenParton_mother[b_nGenParton] = iGdau;

            ++nggdau;
            ++b_nGenParton;
            dauPtrToIdx[ggdau] = b_nGenParton;
            if ( b_nGenParton >= GenParton_N ) break;
          }
          b_GenParton_dau1[iGdau] = iGdau+1;
          b_GenParton_dau2[iGdau] = iGdau+nggdau;
          if ( b_nGenParton >= GenParton_N ) break;
        }
        b_GenParton_dau1[iDau] = iDau+1;
        b_GenParton_dau2[iDau] = iDau+ngdau;

        if ( b_nGenParton >= GenParton_N ) break;
*/
      }
      if ( b_nGenParton >= GenParton_N ) break;
    }

    if ( branchMET and branchMET->GetEntries() > 0 ) {
      const MissingET* met = (const MissingET*)branchMET->At(0);
      b_MET_pt = met->MET;
      b_MET_phi = met->Phi;
    }
    else {
      b_MET_pt = b_MET_phi = 0;
    }

    b_nMuon = 0;
    if ( branchMuon ) {
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
    }

    b_nElectron = 0;
    if ( branchElectron ) {
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
    }

    b_nJet = b_njetDau = 0;
    if ( branchJet ) {
      for ( int i=0; i<branchJet->GetEntries(); ++i ) {
        const Jet* jet = (const Jet*) branchJet->At(i);
        //const TLorentzVector p4 = jet->P4();

        b_Jet_pt[b_nJet] = jet->PT;
        b_Jet_eta[b_nJet] = jet->Eta;
        b_Jet_phi[b_nJet] = jet->Phi;
        b_Jet_m[b_nJet] = jet->Mass;
        b_Jet_flav[b_nJet] = jet->Flavor;
        b_Jet_bTag[b_nJet] = jet->BTag;
        b_Jet_partonIdx[b_nJet] = -1;

        // Keep the subjet particles
        TRefArray cons = jet->Particles; //Constituents;
        std::set<const GenParticle*> matches;
        for ( int j=0; j<cons.GetEntriesFast(); ++j ) {
          if ( b_njetDau > jetDau_N ) break;

          const TObject* obj = cons.At(j);
          if ( !obj ) continue;

          if ( dojetDau ) {
            const Track* track = dynamic_cast<const Track*>(obj);
            const Tower* tower = dynamic_cast<const Tower*>(obj);
            const Muon* muon = dynamic_cast<const Muon*>(obj);
            const Electron* elec = dynamic_cast<const Electron*>(obj);
            const Photon* phot = dynamic_cast<const Photon*>(obj);
            const GenParticle* gen = dynamic_cast<const GenParticle*>(obj);

            if ( track ) {
              const GenParticle* p = dynamic_cast<const GenParticle*>(track->Particle.GetObject());
              findPartonAncestors(p, branchGen, dauPtrToIdx, matches);
              b_jetDau_pt[b_njetDau] = track->PT;
              b_jetDau_eta[b_njetDau] = track->Eta;
              b_jetDau_phi[b_njetDau] = track->Phi;
              b_jetDau_q[b_njetDau] = track->Charge;
              //b_jetDau_pdgId[b_njetDau] = track->Charge*211;
              b_jetDau_pdgId[b_njetDau] = p->PID;
            }
            else if ( tower ) {
              b_jetDau_pt[b_njetDau] = tower->ET;
              b_jetDau_eta[b_njetDau] = tower->Eta;
              b_jetDau_phi[b_njetDau] = tower->Phi;
              b_jetDau_q[b_njetDau] = 0;
              //const bool isPhoton = ( tower->Eem > tower->Ehad ); // Crude estimation
              //if ( isPhoton ) b_jetDau_pdgId[b_njetDau] = 22; // photons
              //else b_jetDau_pdgId[b_njetDau] = 2112; // set as neutron
              int nPhoton = 0;
              TRefArray ps = tower->Particles;
              for ( int k=0; k<ps.GetEntries(); ++k ) {
                const GenParticle* p = dynamic_cast<const GenParticle*>(ps.At(k));
                findPartonAncestors(p, branchGen, dauPtrToIdx, matches);

                if ( p->PID == 22 ) ++nPhoton;
              }
              b_jetDau_pdgId[b_njetDau] = nPhoton > 0 ? 22 : 2112; // set as neutron if no photon found in this tower
            }
            else if ( muon ) {
              b_jetDau_pt[b_njetDau] = muon->PT;
              b_jetDau_eta[b_njetDau] = muon->Eta;
              b_jetDau_phi[b_njetDau] = muon->Phi;
              b_jetDau_q[b_njetDau] = muon->Charge;
              b_jetDau_pdgId[b_njetDau] = muon->Charge*-13;
            }
            else if ( elec ) {
              b_jetDau_pt[b_njetDau] = elec->PT;
              b_jetDau_eta[b_njetDau] = elec->Eta;
              b_jetDau_phi[b_njetDau] = elec->Phi;
              b_jetDau_q[b_njetDau] = elec->Charge;
              b_jetDau_pdgId[b_njetDau] = track->Charge*-11;
            }
            else if ( phot ) {
              b_jetDau_pt[b_njetDau] = phot->PT;
              b_jetDau_eta[b_njetDau] = phot->Eta;
              b_jetDau_phi[b_njetDau] = phot->Phi;
              b_jetDau_q[b_njetDau] = 0;
              b_jetDau_pdgId[b_njetDau] = 22;
            }
            else if ( gen ) {
              b_jetDau_pt[b_njetDau] = gen->PT;
              b_jetDau_eta[b_njetDau] = gen->Eta;
              b_jetDau_phi[b_njetDau] = gen->Phi;
              b_jetDau_q[b_njetDau] = gen->Charge;
              //b_jetDau_pdgId[b_njetDau] = track->Charge*211;
              b_jetDau_pdgId[b_njetDau] = gen->PID;
            }
            else {
              if ( doPrintDebug ) std::cout << obj->IsA()->GetName() << endl;
              continue;
            }

            b_jetDau_jetIdx[b_njetDau] = b_nJet;
            ++b_njetDau;
          }
        }

        double dR0 = 1e9;
        TLorentzVector jetLVec;
        jetLVec.SetPtEtaPhiM(b_Jet_pt[i], b_Jet_eta[i], b_Jet_phi[i], b_Jet_m[i]);
        for ( auto p : matches ) {
          const int idx1 = dauPtrToIdx[p];
          const int dau1 = b_GenParton_dau1[idx1];
          if ( dau1 != -1 ) continue;
          const int pdgId = b_GenParton_pdgId[idx1];
          if ( abs(pdgId) > 5 ) continue; // match quarks only

          const int idx0 = b_Jet_partonIdx[b_nJet];
          if ( idx0 == -1 ) {
            TLorentzVector lvec0;
            lvec0.SetPtEtaPhiM(b_GenParton_pt[idx0], b_GenParton_eta[idx0], b_GenParton_phi[idx0], b_GenParton_m[idx0]);
            dR0 = jetLVec.DeltaR(lvec0);
            b_Jet_partonIdx[b_nJet] = idx1;
            continue;
          }

          TLorentzVector lvec1;
          lvec1.SetPtEtaPhiM(b_GenParton_pt[idx1], b_GenParton_eta[idx1], b_GenParton_phi[idx1], b_GenParton_m[idx1]);
          const double dR1 = jetLVec.DeltaR(lvec1);
          if ( dR0 > dR1 ) {
            dR0 = dR1;
            b_Jet_partonIdx[b_nJet] = idx1;
            if ( doPrintDebug ) {
              std::cout << "Multiple matches detected on Jet" << b_nJet << " (" << b_Jet_pt[i] << "," << b_Jet_eta[i] << "," << b_Jet_phi[i] << "). Replacing best match:"
                << "\n  <- PdgId=" << b_GenParton_pdgId[idx0] << " (" << b_GenParton_pt[idx0] << "," << b_GenParton_eta[idx0] << "," << b_GenParton_phi[idx0] << ")"
                << "\n  -> PdgId=" << b_GenParton_pdgId[idx1] << " (" << b_GenParton_pt[idx1] << "," << b_GenParton_eta[idx1] << "," << b_GenParton_phi[idx1] << ")"
                << endl;
            }
          }
        }

        ++b_nJet;
        if ( b_nJet >= Jet_N ) break;
      }
    }

    b_nGenJet = 0; //b_nSubGenJet = 0;
    if ( branchGenJet ) {
      for ( int i=0; i<branchGenJet->GetEntries(); ++i ) {
        const Jet* jet = (const Jet*) branchGenJet->At(i);
        //const TLorentzVector p4 = jet->P4();

        b_GenJet_pt[b_nGenJet] = jet->PT;
        b_GenJet_eta[b_nGenJet] = jet->Eta;
        b_GenJet_phi[b_nGenJet] = jet->Phi;
        b_GenJet_m[b_nGenJet] = jet->Mass;
        b_GenJet_partonIdx[b_nGenJet] = -1;

        // Loop over the constituents to find its origin
        TRefArray cons = jet->Constituents;
        std::set<const GenParticle*> matches;
        for ( int j=0; j<cons.GetEntriesFast(); ++j ) {
          const TObject* obj = cons.At(j);
          if ( !obj ) continue;

          const GenParticle* p = dynamic_cast<const GenParticle*>(obj);
          if ( !p ) continue;

          findPartonAncestors(p, branchGen, dauPtrToIdx, matches);
        }

        double dR0 = 1e9;
        TLorentzVector genJetLVec;
        genJetLVec.SetPtEtaPhiM(b_GenJet_pt[i], b_GenJet_eta[i], b_GenJet_phi[i], b_GenJet_m[i]);
        for ( auto p : matches ) {
          const int idx1 = dauPtrToIdx[p];
          const int dau1 = b_GenParton_dau1[idx1];
          if ( dau1 != -1 ) continue;
          const int pdgId = b_GenParton_pdgId[idx1];
          if ( abs(pdgId) > 5 ) continue; // match quarks only

          const int idx0 = b_GenJet_partonIdx[b_nGenJet];
          if ( idx0 == -1 ) {
            TLorentzVector lvec0;
            lvec0.SetPtEtaPhiM(b_GenParton_pt[idx0], b_GenParton_eta[idx0], b_GenParton_phi[idx0], b_GenParton_m[idx0]);
            dR0 = genJetLVec.DeltaR(lvec0);
            b_GenJet_partonIdx[b_nGenJet] = idx1;
            continue;
          }

          TLorentzVector lvec1;
          lvec1.SetPtEtaPhiM(b_GenParton_pt[idx1], b_GenParton_eta[idx1], b_GenParton_phi[idx1], b_GenParton_m[idx1]);
          const double dR1 = genJetLVec.DeltaR(lvec1);
          if ( dR0 > dR1 ) {
            dR0 = dR1;
            b_GenJet_partonIdx[b_nGenJet] = idx1;
            if ( doPrintDebug ) {
              std::cout << "Multiple matches detected on GenJet" << b_nGenJet << " (" << b_GenJet_pt[i] << "," << b_GenJet_eta[i] << "," << b_GenJet_phi[i] << "). Replacing best match:"
                << "\n  <- PdgId=" << b_GenParton_pdgId[idx0] << " (" << b_GenParton_pt[idx0] << "," << b_GenParton_eta[idx0] << "," << b_GenParton_phi[idx0] << ")"
                << "\n  -> PdgId=" << b_GenParton_pdgId[idx1] << " (" << b_GenParton_pt[idx1] << "," << b_GenParton_eta[idx1] << "," << b_GenParton_phi[idx1] << ")"
                << endl;
            }
          }
          b_GenJet_partonIdx[b_nGenJet] = idx1;
        }

        ++b_nGenJet;
        if ( b_nGenJet >= GenJet_N ) break;
      }
    }

    tree->Fill();
  }

  tree->Write();
  fout->Close();

}

