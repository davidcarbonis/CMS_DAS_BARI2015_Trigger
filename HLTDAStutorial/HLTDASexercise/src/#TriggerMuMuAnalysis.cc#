#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "HLTDAStutorial/HLTDASexercise/src/Histograms.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"

#include "Math/SMatrix.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TMatrixD.h>
#include <iostream>
#include <map>
#include <set>
#include <TROOT.h>

using namespace std;
using namespace reco;
using namespace edm;
using namespace trigger;


//-----------------------------------------------------------------


// generically maximum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}


//-----------------------------------------------------------------
class TriggerMuMuAnalysis : public edm::EDAnalyzer {
public:
  explicit TriggerMuMuAnalysis(const edm::ParameterSet&);
  ~TriggerMuMuAnalysis();

private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void endRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool muoncuts(MuonCollection::const_iterator themu, const reco::BeamSpot& beamspot);
  bool triggerfired(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname);
  bool triggerfound(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname);
  unsigned int triggerIndex(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname);


  HLTConfigProvider hltConfig_; // to get configuration for L1s/Pre
  edm::InputTag hlTriggerResults_;    // Input tag for TriggerResults
  edm::InputTag hlTriggerEvent_;      // Input tag for TriggerEvent

  std::string triggerName_;
  edm::InputTag collectionName_;

  edm::InputTag beamSpotTag_;
  edm::InputTag offlineMuonsTag_;      // Input tag for Offline (reco) muons
  edm::InputTag L3muCandLabel_;
  edm::InputTag L3muDisplVtxCandLabel_;

  double maxDRmatch_;
  bool debug;

  HistoVertex *hRecoDiMuonVertex;
  HistoVertex *hRecoDiMuonVertex_den;

  HistoKin *hL3MuonPt1;
  HistoKin *hL3MuonPt2;
  HistoKin *hL3MuonAll;
  HistoKinPair *hL3DiMu;

  HistoKin *hRecoMuonPt1;
  HistoKin *hRecoMuonPt2;
  HistoKin *hRecoMuonAll;
  HistoKinPair *hRecoDiMu;

  HistoKin *hRecoMuonPt1_noV;
  HistoKin *hRecoMuonPt2_noV;
  HistoKin *hRecoMuonAll_noV;
  HistoKinPair *hRecoDiMu_noV;

  HistoKin *hRecoMuonPt1_den;
  HistoKin *hRecoMuonPt2_den;
  HistoKin *hRecoMuonAll_den;
  HistoKinPair *hRecoDiMu_den;

  TH1D* hNprimVtx_den;
  TH1D* hNprimVtx_num;

  int counter;
  int countInAccepted;
  int countInTriggered;
  int countRecoMatchTrig;
  int countHasVertex;
  int counterMoreThan3Muons;
  int counterMoreThanOneds;
  int counterTaus;
  int countDiMuValidVtx;
  int hltL1s_, hltPre_, hlt2mu_, hltVtx_, hlAccept_;

};

//
//

TriggerMuMuAnalysis::TriggerMuMuAnalysis(const edm::ParameterSet& iConfig) {

  counter = 0;
  countInAccepted = 0;
  countInTriggered = 0;
  countHasVertex=0;
  counterMoreThan3Muons = 0;
  counterMoreThanOneds = 0;
  counterTaus = 0;
  countRecoMatchTrig=0;
  countDiMuValidVtx=0;
  hltL1s_ = 0; hltPre_ = 0; hlt2mu_ = 0; hltVtx_ = 0; hlAccept_ = 0;

  hlTriggerResults_  = iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultsTag", edm::InputTag("TriggerResults", "", "HLT"));
  hlTriggerEvent_    = iConfig.getUntrackedParameter<edm::InputTag>("TriggerEventTag", edm::InputTag("hltTriggerSummaryAOD", "", "HLT"));
  beamSpotTag_       = iConfig.getUntrackedParameter<edm::InputTag>("BeamSpotTag",edm::InputTag("offlineBeamSpot","","RECO"));
  offlineMuonsTag_    = iConfig.getUntrackedParameter<edm::InputTag>("OfflineMuonsTag",edm::InputTag("muons","","RECO"));
  
  triggerName_ = iConfig.getUntrackedParameter<std::string>("PathName","HLT_Dimuon0_Jpsi_v");
  collectionName_ = iConfig.getUntrackedParameter<edm::InputTag>("CollectionName",edm::InputTag("hltL3MuonCandidates","","HLT"));
  L3muCandLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("L3MuonFilter",edm::InputTag("hltJpsiL3Filtered","","HLT"));
  L3muDisplVtxCandLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("DiMuonVertexTag",edm::InputTag("hltVertexmumuFilterJpsi","","HLT"));

  maxDRmatch_ = iConfig.getUntrackedParameter<double>("maxDRforRecoMatching",0.1);
  
  debug = iConfig.getUntrackedParameter<bool>("debug",false);
}


TriggerMuMuAnalysis::~TriggerMuMuAnalysis() {}


bool TriggerMuMuAnalysis::muoncuts(MuonCollection::const_iterator themu, const reco::BeamSpot& beamspot){
  bool glb = themu->isGlobalMuon();
  //bool trk = themu->isTrackerMuon();
  if (!(glb)) return false;
  const reco::TrackRef mutrackref = themu->outerTrack();
  int nstations = themu->numberOfMatches();
  //int nsahits =   mutrackref->numberOfValidHits();
  double sachi2 =    mutrackref->chi2()/mutrackref->ndof();
  //if (nstations<2 || nsahits<1 || sachi2>10) return false;
  if (nstations<2 || sachi2>10) return false;
  const reco::TrackRef glbmutrackref = themu->innerTrack();
  //int pixhits = glbmutrackref->hitPattern().numberOfValidPixelHits();
  //int trkhits = glbmutrackref->hitPattern().numberOfValidTrackerHits();
  //if (pixhits<1 || trkhits<6) return false;
  double glbchi2 = glbmutrackref->chi2()/glbmutrackref->ndof();
  double dxy= glbmutrackref->dxy(beamspot.position());
  double pt=  sqrt(themu->px()*themu->px()+themu->py()*themu->py());
  if (pt<3.5 || fabs(dxy)>0.2 || glbchi2>10) return false;
  return true;
}

bool TriggerMuMuAnalysis::triggerfired(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname){
  const edm::TriggerNames TrigNames_ = ev.triggerNames(*triggerResultsHandle_);
  const unsigned int ntrigs = triggerResultsHandle_->size();
  for (unsigned int itr=0; itr<ntrigs; itr++){
    TString trigName=TrigNames_.triggerName(itr);
    if (!triggerResultsHandle_->accept(itr)) continue;
    if(trigName.Contains(trigname))      return true;
  }
  return false;
}

bool TriggerMuMuAnalysis::triggerfound(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname){
  const edm::TriggerNames TrigNames_ = ev.triggerNames(*triggerResultsHandle_);
  const unsigned int ntrigs = triggerResultsHandle_->size();
  for (unsigned int itr=0; itr<ntrigs; itr++){
    TString trigName=TrigNames_.triggerName(itr);
    if(trigName.Contains(trigname))      return true;
  }
  return false;
}

unsigned int TriggerMuMuAnalysis::triggerIndex(const edm::Event& ev, edm::Handle<edm::TriggerResults> triggerResultsHandle_, TString trigname){
  const edm::TriggerNames TrigNames_ = ev.triggerNames(*triggerResultsHandle_);
  const unsigned int ntrigs = triggerResultsHandle_->size();
  unsigned int itr;
  for (itr=0; itr<ntrigs; itr++){
    TString trigName=TrigNames_.triggerName(itr);
    if(trigName.Contains(trigname))      return itr;
  }
  return itr;
}



//
// member functions
//

// ------------ method called to for each event  ------------

void TriggerMuMuAnalysis::analyze(const edm::Event& ev, const edm::EventSetup& iSetup) {

  float mumass = 0.1057;
  int pdgIdMu = 13;

  counter++;
  if (debug)  cout<<"==== Event "<<counter<<" ===="<<endl;

  edm::Handle<reco::VertexCollection> privtxs;
  ev.getByLabel(edm::InputTag("offlinePrimaryVertices","","RECO"), privtxs);
  int nvtx=(*privtxs).size();

  // *** get Handle to the TriggerEvent
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  ev.getByLabel(hlTriggerEvent_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    cout << "Error in getting TriggerEvent product from Event!" << endl;
    return;
  }

  // *** get Handle to the TriggerResults
  edm::Handle<TriggerResults> HLTR; 
  ev.getByLabel(hlTriggerResults_, HLTR);
  if (!HLTR.isValid()) {
    if (debug) cout << "HLT TriggerResults with label " << hlTriggerResults_ << " not found!";
    return;
  }

  // get offline reco muons
  edm::Handle<MuonCollection> muons;
  ev.getByLabel(offlineMuonsTag_,muons);

  // get beamspot
  edm::Handle<BeamSpot> beamSpoth;
  ev.getByLabel(beamSpotTag_, beamSpoth);
  const reco::BeamSpot& vertexBeamSpot = *beamSpoth;


  // Only events in which the path actually fired had stored the filter results and products:	  

  bool triggerFound = triggerfound(ev,HLTR,triggerName_);
  if (triggerFound) countInAccepted++;
  bool triggerFired = triggerfired(ev,HLTR,triggerName_);

  const unsigned int numberOfHltPaths = HLTR->size();
  const unsigned int numberOfHltFilters = triggerEventHandle_->sizeFilters();

  unsigned int pathIndex = triggerIndex(ev,HLTR,triggerName_);
  if (pathIndex>=numberOfHltPaths) {
    cout << " WARNING: path " << triggerName_ << " out of range in TriggerResults" << endl;
    return;
  }

  if (HLTR->wasrun(pathIndex)) {
    if (!triggerFound) cout << " WARNING: path exists in HLTR but not found in TriggerResults" << endl;
  }
  else {
    if (triggerFound) cout << " WARNING: path found in TriggerResults but it does not exist in HLTR" << endl;
  }


   int posL1s_ = -1;
   int posPre_ = -1;
   int pos2mu_ = -1;
   int posVtx_ = -1;
   const std::vector<std::string> & moduleLabels(hltConfig_.moduleLabels(pathIndex));
   for (unsigned int j = 0; j < moduleLabels.size(); ++j) {
     const std::string & label = hltConfig_.moduleType(moduleLabels[j]);
     if (label == "HLTLevel1GTSeed") {
       posL1s_ = j;
       if (debug) cout << " *** " << moduleLabels[j] << " is a " << label << endl;
     }
     else if (label == "HLTPrescaler") {
       posPre_ = j;
       if (debug) cout << " *** " << moduleLabels[j] << " is a " << label << endl;
     }
     else if (label == "HLTMuonDimuonL3Filter") {
       pos2mu_ = j;
       if (debug) cout << " *** " << moduleLabels[j] << " is a " << label << endl;
     }
     else if (label == "HLTDisplacedmumuFilter") {
       posVtx_ = j;
       if (debug) cout << " *** " << moduleLabels[j] << " is a " << label << endl;
     }
   }
   
  const int index(static_cast<int>(HLTR->index(pathIndex)));
  if (HLTR->accept(pathIndex)) {
    hlAccept_++;
    hltL1s_++; hltPre_++; hlt2mu_++; hltVtx_++;
    if (debug) cout << " Path accepted, with state " << HLTR->state(pathIndex) 
		    << " , and index = " << index << " vs " << posL1s_ << " / " << posPre_ << " / " << pos2mu_ << " / " << posVtx_ << endl;
  } else {
    if (index >  posL1s_) hltL1s_++;
    if (index >  posPre_) hltPre_++;
    if (index >  pos2mu_) hlt2mu_++;
    if (index >  posVtx_) hltVtx_++;
    if (debug) cout << " Path not accepted, with state " << HLTR->state(pathIndex) 
		    << " , but index = " << index << " vs " << posL1s_ << " / " << posPre_ << " / " << pos2mu_ << " / " << posVtx_ << endl;
  }


  trigger::size_type collIndex(0);
  collIndex = triggerEventHandle_->collectionIndex(collectionName_);
  if (collIndex < triggerEventHandle_->sizeCollections()) {
    const trigger::Keys& Keys(triggerEventHandle_->collectionKeys());
    const trigger::size_type n0 (collIndex == 0? 0 : Keys.at(collIndex-1));
    const trigger::size_type n1 (Keys.at(collIndex));
    for (trigger::size_type i = n0; i != n1; ++i) {
      const trigger::TriggerObject& obj( triggerEventHandle_->getObjects().at(i) );
      if (abs(obj.id()) == 13) {
	TLorentzVector L3muon_4mom(0., 0., 0., 0.);
	L3muon_4mom.SetPtEtaPhiM(obj.pt(), obj.eta(), obj.phi(), obj.mass());
	hL3MuonAll->Fill(L3muon_4mom.Pt(),L3muon_4mom.Eta(),L3muon_4mom.Phi(),mumass,1);
      }
    }
  }

  if (triggerFired) {
    countInTriggered++;
    const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());

    // hltL3MuonCandidates

    const unsigned int L3filterIndex(triggerEventHandle_->filterIndex(L3muCandLabel_));
    vector<TLorentzVector> trigmus;
    if (debug) cout << "- Looking for L3 muon candidates; L3filterIndex = " << L3filterIndex<< " out of " << numberOfHltFilters << " filters" << endl;
    if (L3filterIndex < numberOfHltFilters) {
      const trigger::Vids& VIDS(triggerEventHandle_->filterIds(L3filterIndex));
      const trigger::Keys& KEYS(triggerEventHandle_->filterKeys(L3filterIndex));
      const trigger::size_type nI(VIDS.size()); const trigger::size_type nK(KEYS.size());
      trigger::size_type n;      if (nI>nK) n=nI; else n=nK;
      if (triggerFired && debug) cout << "Objects that fired the Trigger:" << endl;
      for (size_type i=0; i!=n; ++i) {
	const TriggerObject& TO(TOC[KEYS[i]]);
	if (abs(TO.id())==pdgIdMu && TO.pt()!=0){
	  TLorentzVector tmu=TLorentzVector(0.,0.,0.,0.); 
	  tmu.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),mumass);
	  trigmus.push_back(tmu);
	  if (triggerFired && debug) cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
					  << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass() << " "  << " " << TO.energy() 
					  << endl;
	  hL3MuonAll->Fill(tmu.Pt(),tmu.Eta(),tmu.Phi(),mumass,1);
	}
      }

      // fill the L3 monitoring plots
      if (trigmus.size()==2){
	TLorentzVector mu1; TLorentzVector mu2;
	if (trigmus[0].Pt()<trigmus[1].Pt()){mu1=trigmus[1]; mu2=trigmus[0];} else {mu1=trigmus[0]; mu2=trigmus[1];}
	
	hL3MuonPt1->Fill(mu1.Pt(),mu1.Eta(),mu1.Phi(),mumass,1);
	hL3MuonPt2->Fill(mu2.Pt(),mu2.Eta(),mu2.Phi(),mumass,1);
	TLorentzVector DiMu = mu1+mu2;
	hL3DiMu->Fill(mu1.Pt(),mu2.Pt(), mu1.DeltaPhi(mu2),fabs(mu1.Eta()-mu2.Eta()), mu1.DeltaR(mu2),DiMu.M(),1);
	
	// now match offline muons to the trigger muons
	int recomatch[2]={-1,-1}; 
	int it=-1;
	vector<TLorentzVector> matchedRecoMu;
	vector<MuonCollection::const_iterator> matched_mu_ptrs;
	for (vector<TLorentzVector>::const_iterator trigMu = trigmus.begin(); trigMu!=trigmus.end(); ++trigMu){
	  it++;
	  double bestDR=maxDRmatch_; int ir=-1;
	  TLorentzVector offmu(0,0,0,0);
	  MuonCollection::const_iterator mitemp;
	  if (triggerFired && debug) cout<< "- Matching one RECO muons to one L3 muon"<<endl;
	  for (MuonCollection::const_iterator recoMu = muons->begin(); recoMu!=muons->end(); ++recoMu){
	    if (!muoncuts(recoMu,vertexBeamSpot)) continue;
	    ir++;
	    TLorentzVector offmu2(0,0,0,0);
	    offmu2.SetPtEtaPhiM(recoMu->pt(),recoMu->eta(),recoMu->phi(),mumass);
	    double dr=trigMu->DeltaR(offmu2);
	    if (dr<bestDR && (ir!=recomatch[it])){
	      bestDR=dr;
	      recomatch[it]=ir;
	      mitemp=recoMu;
	      offmu=offmu2;
	    }
	  }
	  if (recomatch[it]!=-1){
	    matchedRecoMu.push_back(offmu);
	    matched_mu_ptrs.push_back(mitemp);
	  }
	}
	
	if (matched_mu_ptrs.size()==2) {
	  if (debug) {
	    for (int itr=0; itr<2; itr++){
	      cout<<" ---> Matched L3 muon:";
	      trigmus[itr].Print();
	      cout<<"      with RECO muon:";
	      matchedRecoMu[itr].Print();
	    }
	  }
	  countRecoMatchTrig++;
	  TLorentzVector rmu1; TLorentzVector rmu2;
	  if (matchedRecoMu[0].Pt()<matchedRecoMu[1].Pt()){rmu1=matchedRecoMu[1]; rmu2=matchedRecoMu[0];} else {rmu1=matchedRecoMu[0]; rmu2=matchedRecoMu[1];}
	  hRecoMuonPt1_noV->Fill(rmu1.Pt(),rmu1.Eta(),rmu1.Phi(),mumass,1);
	  hRecoMuonPt2_noV->Fill(rmu2.Pt(),rmu2.Eta(),rmu2.Phi(),mumass,1);
	  hRecoMuonAll_noV->Fill(rmu1.Pt(),rmu1.Eta(),rmu1.Phi(),mumass,1);
	  hRecoMuonAll_noV->Fill(rmu2.Pt(),rmu2.Eta(),rmu2.Phi(),mumass,1);
	  TLorentzVector rDiMu = rmu1+rmu2;
	  hRecoDiMu_noV->Fill(rmu1.Pt(),rmu2.Pt(), rmu1.DeltaPhi(rmu2),fabs(rmu1.Eta()-rmu2.Eta()), rmu1.DeltaR(rmu2),rDiMu.M(),1);
	  
	  // now prepare tracks to build vertexes
	  edm::ESHandle<TransientTrackBuilder> Builder;
	  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", Builder); 
	  vector<TransientTrack> tts;
	  KalmanVertexFitter vfitter;
	  reco::TrackRef intrack1 = matched_mu_ptrs[0]->innerTrack();
	  reco::TrackRef intrack2 = matched_mu_ptrs[1]->innerTrack();
	  if ( (intrack1.isNonnull() && intrack1.isAvailable()) && (intrack2.isNonnull() && intrack2.isAvailable()) ) {
	    TransientTrack t1=Builder->build(intrack1);
	    TransientTrack t2=Builder->build(intrack2);
	    tts.push_back(t1); tts.push_back(t2);
	    // build vertex and check
	    TransientVertex recomuvtx;
	    recomuvtx=vfitter.vertex(tts);
	    if (debug) cout<<"- Checking if it triggered also the HLT dimuon vertex flavor: ";
	    if (recomuvtx.isValid()){
	      countDiMuValidVtx++; bool hltvtx=false;
	      if (triggerFired){
		countHasVertex++; hltvtx=true;
		if (debug) cout << "YES" << endl;
	      }
	      else	      {
		if (debug) cout << "NO" << endl;
	      }
	      // build dimu vertex from offline muons
	      GlobalPoint secondaryVertex = recomuvtx.position();
	      GlobalError err = recomuvtx.positionError();
	      //calculate decay length  significance w.r.t. the beamspot
	      GlobalPoint displacementFromBeamspot( -1*((vertexBeamSpot.x0()-secondaryVertex.x()) + (secondaryVertex.z()-vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()),
						    -1*((vertexBeamSpot.y0()-secondaryVertex.y()) + (secondaryVertex.z()-vertexBeamSpot.z0()) * vertexBeamSpot.dydz()),
						    0);
	      double lxy = displacementFromBeamspot.perp();
	      double lxyerr = sqrt(err.rerr(displacementFromBeamspot));
	      double normChi2 = recomuvtx.normalisedChiSquared();
	      
	      // fill reco mu plot, vertexed flavor (denominator for HLT vtx eff)
	      hRecoMuonPt1_den->Fill(rmu1.Pt(),rmu1.Eta(),rmu1.Phi(),mumass,1);
	      hRecoMuonPt2_den->Fill(rmu2.Pt(),rmu2.Eta(),rmu2.Phi(),mumass,1);
	      hRecoMuonAll_den->Fill(rmu1.Pt(),rmu1.Eta(),rmu1.Phi(),mumass,1);
	      hRecoMuonAll_den->Fill(rmu2.Pt(),rmu2.Eta(),rmu2.Phi(),mumass,1);
	      hRecoDiMu_den->Fill(rmu1.Pt(),rmu2.Pt(), rmu1.DeltaPhi(rmu2), fabs(rmu1.Eta()-rmu2.Eta()), rmu1.DeltaR(rmu2), rDiMu.M(),1);
	      hRecoDiMuonVertex_den->Fill(normChi2, TMath::Prob(normChi2*2,2), lxy, lxy/lxyerr, rDiMu.Pt(),0,rDiMu.M(),1);
	      hNprimVtx_den->Fill(nvtx);
	      if (hltvtx){
		hRecoMuonPt1->Fill(rmu1.Pt(),rmu1.Eta(),rmu1.Phi(),mumass,1);
		hRecoMuonPt2->Fill(rmu2.Pt(),rmu2.Eta(),rmu2.Phi(),mumass,1);
		hRecoMuonAll->Fill(rmu1.Pt(),rmu1.Eta(),rmu1.Phi(),mumass,1);
		hRecoMuonAll->Fill(rmu2.Pt(),rmu2.Eta(),rmu2.Phi(),mumass,1);
		hRecoDiMu->Fill(rmu1.Pt(),rmu2.Pt(), rmu1.DeltaPhi(rmu2),fabs(rmu1.Eta()-rmu2.Eta()), rmu1.DeltaR(rmu2),rDiMu.M(),1);
		hRecoDiMuonVertex->Fill(normChi2, TMath::Prob(normChi2*2,2), lxy, lxy/lxyerr, rDiMu.Pt(),0,rDiMu.M(),1);
		hNprimVtx_num->Fill(nvtx);
	      }
	    } else {if (debug) cout << "Vertex built with reco muons not valid!" << endl;}
	  } else {if (debug) cout << "There weren't two inner tracks to build the muon vertex" << endl;}
	} else {if (debug) cout << "Not enough offline muons matched to trigger muons!" << endl;}
      } else {if (debug) cout << "Less than 2 L3 muons were found in the trigger object!" << endl;} 
    } else {if (debug) cout << "Module index out of HLT modules index range ..." << endl;}    
  }  
  
}

// ------------ method called once each job just before starting event loop  ------------
void TriggerMuMuAnalysis::beginJob() {
  cout << "begin job" << endl;
  edm::Service<TFileService> fs;
  hL3MuonPt1=new HistoKin("L3MuonPt1",*fs);
  hL3MuonPt2=new HistoKin("L3MuonPt2",*fs);
  hL3MuonAll=new HistoKin("L3MuonAll",*fs);
  hL3DiMu=new HistoKinPair("L3DiMuon", 0, 150,*fs);

  hRecoMuonPt1=new HistoKin("RecoMuonPt1",*fs);
  hRecoMuonPt2=new HistoKin("RecoMuonPt2",*fs);
  hRecoMuonAll=new HistoKin("RecoMuonAll",*fs);
  hRecoDiMu=new HistoKinPair("RecoDiMuon", 0, 150,*fs);

  hRecoDiMuonVertex = new HistoVertex("RecoDiMuKVFVertex",*fs);

  hRecoMuonPt1_den=new HistoKin("RecoMuonPt1_den",*fs);
  hRecoMuonPt2_den=new HistoKin("RecoMuonPt2_den",*fs);
  hRecoMuonAll_den=new HistoKin("RecoMuonAll_den",*fs);
  hRecoDiMu_den=new HistoKinPair("RecoDiMuon_den", 0, 150,*fs);

  hRecoDiMuonVertex_den = new HistoVertex("RecoDiMuKVFVertex_den",*fs);

  hRecoMuonPt1_noV=new HistoKin("RecoMuonPt1_noV",*fs);
  hRecoMuonPt2_noV=new HistoKin("RecoMuonPt2_noV",*fs);
  hRecoMuonAll_noV=new HistoKin("RecoMuonAll_noV",*fs);
  hRecoDiMu_noV=new HistoKinPair("RecoDiMuon_noV", 0, 150,*fs);

  hNprimVtx_num = fs->make<TH1D>("hNprimVtx_num","# hNprimVtx_num; # vertices; # events",25,0,50);
  hNprimVtx_den = fs->make<TH1D>("hNprimVtx_den","# hNprimVtx_den; # vertices; # events",25,0,50);
  hNprimVtx_num->Sumw2();
  hNprimVtx_den->Sumw2();
}

void TriggerMuMuAnalysis::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
   bool changed = true;
   bool configured = hltConfig_.init(iRun, iSetup, hlTriggerResults_.process(), changed);
   if (!configured) {
     cout << "WARNING: HLTConfigProvider cannot be configured" << endl;
   }
   if (changed) {
     if (debug) cout << "--> Starting a new HLT table" << endl;
   }
}

void TriggerMuMuAnalysis::endRun(edm::Run const & run, edm::EventSetup const & setup) { }



// ------------ method called once per job, just after ending the event loop  ------------
void 
TriggerMuMuAnalysis::endJob() {
  cout << "Total # events: " << counter << endl;
  cout << "Total # events triggered by " << triggerName_ << " : " << countInTriggered << endl;
  cout << "Total # events triggered by " << triggerName_ << " with offline/l3 muon matchings: " << countRecoMatchTrig << endl;
  cout << "Total # events where offline muon make a vertex: " << countDiMuValidVtx << endl;
  cout << "Total # events triggered by " << triggerName_ << " with vertex: " << countHasVertex << endl;
  cout << "Fraction of " << triggerName_ << " : " << double(countInTriggered)/double(countInAccepted) << endl;

  cout << "Filtered events: " << counter << " ==> " 
       << hltL1s_ << " -> " << hltPre_ << " -> " << hlt2mu_ << " -> " << hltVtx_ << " ==> " << hlAccept_;

}


#include "FWCore/Framework/interface/MakerMacros.h"  
DEFINE_FWK_MODULE( TriggerMuMuAnalysis );
