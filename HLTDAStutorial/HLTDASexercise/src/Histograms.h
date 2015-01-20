#ifndef Histograms_H
#define Histograms_H

#ifndef ROOTANALYSIS
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#endif

#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"

#include <iostream>
using namespace std; 

class HistoKin {
public:

#ifndef ROOTANALYSIS
  HistoKin(std::string name, TFileService& fs) : theName(name) {
   

    float bins[43];
    float previous = 0.;
    float step = 0.5;
    for(int id = 0; id != 41; ++id) {
//       cout << "id: " << id << " prev: " << previous << endl;
      bins[id] = previous;
      previous += step;

    }
    bins[41] = 30.;
    bins[42] = 50.;
    hPt        = fs.make<TH1F>(theName+"_hPt",       "Transv. Momentum; P_{T} (GeV/c^{2}); # events", 42, bins);
    hPtForEff  = fs.make<TH1F>(theName+"_hPtForEff", "Transv. Momentum; P_{T} (GeV/c^{2}); # events", 42, bins);

    hEta       = fs.make<TH1F>(theName+"_hEta",      "Pseudo rapidity; #eta; # events",               50,  -5.0, 5.0);
    hPhi       = fs.make<TH1F>(theName+"_hPhi",      "Phi; #phi (rad); # events",                    100, -3.14, 3.14);
    hMass      = fs.make<TH1F>(theName+"_hMass",     "Mass; Mass (GeV/c^{2}); # events",             100,     0, 10.0);
    hNObj      = fs.make<TH1F>(theName+"_hNObj",     "# objects",                                     50,     0, 50);
    hPtEta     = fs.make<TH2F>(theName+"_hPtEta",    "Pt vs Eta; #eta; p_{T} (GeV/c^{2})",            50,  -5.0, 5.0, 100, 0, 50);  

    hPt->Sumw2();
    hEta->Sumw2();
    hPhi->Sumw2();
    hMass->Sumw2();
  }

#endif
  
  
  HistoKin(std::string name) : theName(name) {
    float bins[43];
    float previous = 0.;
    float step = 0.5;
    for(int id = 0; id != 41; ++id) {
      bins[id] = previous;
      previous += step;
    }
    bins[41] = 30.;
    bins[42] = 50.;

    hPt        = new TH1F(theName+"_hPt",       "Transv. Momentum; P_{T} (GeV/c^{2}); # events", 42, bins);
    hPtForEff  = new TH1F(theName+"_hPtForEff", "Transv. Momentum; P_{T} (GeV/c^{2}); # events", 42, bins);

    hEta       = new TH1F(theName+"_hEta",      "Pseudo rapidity; #eta; # events",               50,  -5.0, 5.0);
    hPhi       = new TH1F(theName+"_hPhi",      "Phi; #phi (rad); # events",                    100, -3.14, 3.14);
    hMass      = new TH1F(theName+"_hMass",     "Mass; Mass (GeV/c^{2}); # events",             100,     0, 10.0);
    hNObj      = new TH1F(theName+"_hNObj",     "# objects",                                     50,     0, 50);
    hPtEta     = new TH2F(theName+"_hPtEta",    "Pt vs Eta; #eta; p_{T} (GeV/c^{2})",            50,  -5.0, 5.0, 100, 0, 50);  

    hPt->Sumw2();
    hEta->Sumw2();
    hPhi->Sumw2();
    hMass->Sumw2();
    hNObj->Sumw2();
  }





  HistoKin() : theName("") {
    hPt = 0;
    hPtForEff = 0;
    hEta = 0;
    hPhi = 0;
    hMass = 0;
    hPtEta = 0;
    hNObj = 0;
  }



  HistoKin(std::string name, std::string dir, TFile *file) : theName(name) {
    if(dir != "") {
      dir = dir + "/";
    }

    hPt = (TH1F *) file->Get(TString(dir) + theName+"_hPt");
    hPtForEff = (TH1F *) file->Get(TString(dir) + theName+"_hPtForEff");
    hEta = (TH1F *) file->Get(TString(dir) + theName+"_hEta");
    hPhi = (TH1F *) file->Get(TString(dir) + theName+"_hPhi");
    hPtEta = (TH2F *) file->Get(TString(dir) + theName+"_hPtEta");
    hMass = (TH1F *) file->Get(TString(dir) + theName+"_hMass");
    hNObj = (TH1F *) file->Get(TString(dir) + theName+"_hNObj");

  }


  HistoKin * Clone(std::string name) {
    HistoKin *ret = new HistoKin();
    ret->theName = name;

    if(hPt != 0) hPt->Clone((ret->theName+"_hPt").Data());
    if(hPtForEff != 0) hPt->Clone((ret->theName+"_hPt").Data());
    if(hEta != 0) hEta->Clone((ret->theName+"_hEta").Data());
    if(hPhi != 0) hPhi->Clone((ret->theName+"_hPhi").Data());
    if(hPtEta != 0) hPtEta->Clone((ret->theName+"_hPtEta").Data());
    if(hMass != 0) hPhi->Clone((ret->theName+"_hMass").Data());
    if(hNObj != 0) hNObj->Clone((ret->theName+"_hNObj").Data());

    return ret;
  }



  void Add(const HistoKin* histSet) {
    if(hPt != 0) hPt->Add(histSet->hPt);
    if(hPtForEff != 0) hPt->Add(histSet->hPt);

    if(hEta != 0) hEta->Add(histSet->hEta);
    if(hPhi != 0) hPhi->Add(histSet->hPhi);
    if(hPtEta != 0) hPtEta->Add(histSet->hPtEta);
    if(hMass != 0) hPhi->Add(histSet->hMass);
    if(hNObj != 0) hPhi->Add(histSet->hNObj);

  }



  void Scale(double scaleFact) {
    if(hPt != 0) hPt->Scale(scaleFact);
    if(hPtForEff != 0) hPt->Scale(scaleFact);
    if(hEta != 0) hEta->Scale(scaleFact);
    if(hPhi != 0) hPhi->Scale(scaleFact);
    if(hPtEta != 0) hPtEta->Scale(scaleFact);
    if(hMass != 0) hMass->Scale(scaleFact);
    if(hNObj != 0) hNObj->Scale(scaleFact);

  }


  void Write() {
    if(hPt != 0) hPt->Write();
    if(hPtForEff != 0) hPt->Write();
    if(hEta != 0) hEta->Write();
    if(hPhi != 0) hPhi->Write();
    if(hPtEta != 0) hPtEta->Write();
    if(hMass != 0) hMass->Write();
    if(hNObj != 0) hNObj->Write();

  }
  
  void Fill(double pt, double eta, double phi, double mass, double weight) {
    hPt->Fill(pt, weight);
    hPtForEff->Fill(pt, weight);
    hEta->Fill(eta, weight);
    hPhi->Fill(phi, weight);
    hPtEta->Fill(eta,pt,weight);
    hMass->Fill(mass, weight);
  }

  void FillNObj(int nObj, double weight) {
    hNObj->Fill(nObj, weight);
  }
  /// Destructor
  virtual ~HistoKin() {}

  // Operations
  TString theName;

  TH1F *hPt;
  TH1F *hPtForEff;

  TH1F *hEta;
  TH1F *hPhi;
  TH2F *hPtEta;
  TH1F *hMass;
  TH1F *hNObj;
};









class HistoKinPair {
public:

#ifndef ROOTANALYSIS
  HistoKinPair(std::string name, float massMin, float massMax, TFileService& fs) : theName(name) {
    hPt1VsPt2     = fs.make<TH2F>(theName+"_hPt1VsPt2","Transv. Momentum; P^{1}_{T} (GeV); P^{2}_{T} ", 300, 0, 150, 300, 0, 150);
    hDEta         = fs.make<TH1F>(theName+"_hDEta","Delta #eta; #Delta #eta; #events",50,0,5);
    hDPhi         = fs.make<TH1F>(theName+"_hDPhi","Delta #phi; #Delta #phi (rad); # events", 100, 0, 3.14);
    hDPhiVsDEta   = fs.make<TH2F>(theName+"_hDPhiVsDEta","Delta #phi vs Delta #eta; #Delta #eta; #Delta #phi (rad)", 50,0,5, 100,0,3.14);
    hDR           = fs.make<TH1F>(theName+"_hDR","Delta R; #Delta R; # events", 100, 0, 5);
    hMass         = fs.make<TH1F>(theName+"_hMass","Mass (GeV)", 300, massMin, massMax);
    hPt1VsPt2->Sumw2();
    hDEta->Sumw2();
    hDPhi->Sumw2();
    hDR->Sumw2();
    hMass->Sumw2();
  }

#endif
  
  


  HistoKinPair() : theName("") {
    hPt1VsPt2=0;
    hDEta=0;
    hDPhi=0;
    hDPhiVsDEta=0;
    hDR=0;
    hMass=0;

  }



  HistoKinPair(std::string name, std::string dir, TFile *file) : theName(name) {
    if(dir != "") {
      dir = dir + "/";
    }

    hPt1VsPt2   = (TH2F *) file->Get(dir + theName+"_hPt1VsPt2");
    hDEta       = (TH1F *) file->Get(dir + theName+"_hDEta");
    hDPhi       = (TH1F *) file->Get(dir + theName+"_hDPhi");
    hDPhiVsDEta = (TH2F *) file->Get(dir + theName+"_hDPhiVsDEta");
    hDR         = (TH1F *) file->Get(dir + theName+"_hDR");
    hMass       = (TH1F *) file->Get(dir + theName+"_hMass");


  }


  void Fill(double pt1, double pt2, double deltaPhi, double deltaEta, double deltaR, double mass, double weight) {
    hPt1VsPt2->Fill(pt2, pt1, weight);
    hDEta->Fill(deltaEta, weight);
    hDPhi->Fill(deltaPhi, weight);
    hDPhiVsDEta->Fill(deltaEta, deltaPhi, weight);
    hDR->Fill(deltaR, weight);
    hMass->Fill(mass, weight);
  }

  /// Destructor
  virtual ~HistoKinPair() {}

  // Operations
  TString theName;

  TH2F *hPt1VsPt2;
  TH1F *hDEta;
  TH1F *hDPhi;
  TH2F *hDPhiVsDEta;
  TH1F *hDR;
  TH1F *hMass;

};






class HistoImpactParam {
public:

#ifndef ROOTANALYSIS
  HistoImpactParam(std::string name, TFileService& fs) : theName(name) {
    hD0 = fs.make<TH1F>(theName+"_hD0","D0;d0 (cm);#events", 1000, 0, 1);
    hD0Err = fs.make<TH1F>(theName+"_hD0Err","D0;d0 (cm);#events", 100, 0, 0.1);
    hD0Sig = fs.make<TH1F>(theName+"_hD0Sig","D0;d0 (cm);#events", 1000, 0, 100);
  }

#endif
  
  


  HistoImpactParam() : theName("") {
    hD0 = 0;
    hD0Err = 0;
    hD0Sig = 0;
  }



  HistoImpactParam(std::string name, std::string dir, TFile *file) : theName(name) {
    if(dir != "") {
      dir = dir + "/";
    }

    hD0          = (TH1F *) file->Get(dir + theName+"_hD0");
    hD0Err       = (TH1F *) file->Get(dir + theName+"_hD0Err");
    hD0Sig       = (TH1F *) file->Get(dir + theName+"_hD0Sig");
  }

  
  void Fill(double d0, double d0Err, double weight) {
    hD0->Fill(d0, weight);
    hD0Err->Fill(d0Err, weight);
    hD0Sig->Fill(d0/d0Err, weight);
  }

  /// Destructor
  virtual ~HistoImpactParam() {}

  // Operations
  TString theName;

  TH1F *hD0;
  TH1F *hD0Err;
  TH1F *hD0Sig;

};

class HistoVertex {
public:

#ifndef ROOTANALYSIS
  HistoVertex(std::string name, TFileService& fs) : theName(name) {
    hNormChi2 = fs.make<TH1F>(theName+"_hNormChi2","Norm #Chi^{2};Norm #Chi^{2};#events", 40, 0, 40);
    hVtxProb = fs.make<TH1F>(theName+"_hVtxProb", "Vertex prob.;Vtx prob; # events", 100,0,1);
    hLxy = fs.make<TH1F>(theName+"_hLxy", "Transv. displ w.r.t BS;Lxy (cm); # events", 50,0,0.2);
    hLxySig = fs.make<TH1F>(theName+"_hLxysig", "Sig. Transv. displ w.r.t BS;LxySig; # events", 100,-5,30);
    hLxySigVsPtDs = fs.make<TH2F>(theName+"_hLxysigVsPtDs", "Sig. Transv. displ w.r.t BS;P_{T}()D_{s}; LxySig;", 100, 0, 200, 100,-5,30);
    hCosPointAngle = fs.make<TH1F>(theName+"_hCosPointAngle", "Cosine pointing angle;Cos pointing angle; # events", 300,-1,1);
    hMass = fs.make<TH1F>(theName+"_hMass", "Mass;Mass (GeV); # events", 200,0,5);
    hNormChi2->Sumw2();
    hVtxProb->Sumw2();
    hLxy->Sumw2();
    hLxySig->Sumw2();
    hLxySigVsPtDs->Sumw2();
    hCosPointAngle->Sumw2();
    hMass->Sumw2();
  }

#endif

  


  HistoVertex() : theName("") {
    hNormChi2 = 0;
    hVtxProb = 0;
    hLxy = 0;
    hLxySig = 0;
    hLxySigVsPtDs = 0;
    hCosPointAngle = 0;
    hMass = 0;
  }



  HistoVertex(std::string name, std::string dir, TFile *file) : theName(name) {
    if(dir != "") {
      dir = dir + "/";
    }

    hNormChi2 = (TH1F *) file->Get(dir + theName+"_hNormChi2");
    hVtxProb = (TH1F *) file->Get(dir + theName+"_hVtxProb");
    hLxy = (TH1F *) file->Get(dir + theName+"_hLxy");
    hLxySig = (TH1F *) file->Get(dir + theName+"_hLxySig");
    hLxySigVsPtDs = (TH2F *) file->Get(dir + theName+"_hLxySigVsPtDs");
    hCosPointAngle = (TH1F *) file->Get(dir + theName+"_hCosPointAngle");
    hMass =  (TH1F *) file->Get(dir + theName+"_hMass");
  }

  
  void Fill(double normChi2, double vtxProb, double lxy, double lxySig, double ptDs, double cosPA, double mass, double weight) {
    hNormChi2->Fill(normChi2,weight);
    hVtxProb->Fill(vtxProb,weight);
    hLxy->Fill(lxy,weight);
    hLxySig->Fill(lxySig,weight);
    hLxySigVsPtDs->Fill(ptDs, lxySig,weight);
    hCosPointAngle->Fill(cosPA, weight);
    hMass->Fill(mass, weight);
}

  /// Destructor
  virtual ~HistoVertex() {}

  // Operations
  TString theName;

  TH1F *hNormChi2;
  TH1F *hVtxProb;
  TH1F *hLxy;
  TH1F *hLxySig;
  TH2F *hLxySigVsPtDs;
  TH1F *hCosPointAngle;
  TH1F *hMass;
};




#endif
