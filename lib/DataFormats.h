// event info
typedef struct {
  int run;
  int lumi;
  long long int event;
  int bx;
  int orbit;
} _EventInfo;


// vertex info
typedef struct{
  int nVertices;
  int isValid[20];
  float x[20];	
  float y[20];	
  float z[20];	
  float xErr[20];	
  float yErr[20];	
  float zErr[20];	
  float chi2[20];
  int ndf[20];
  float normChi2[20];
} _VertexInfo;


// basic track info
typedef struct {
  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;
} _TrackInfo;

// muon info
typedef struct {

  int isTracker;
  int isStandAlone;
  int isGlobal;

  int isTightMuon;
  int isMediumMuon;
  int isLooseMuon;

  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;

  float trkPt;
  float trkPtErr;
  float trketa;
  float trkPhi;
  
  float normChiSquare;
  float d0_BS;
  float dz_BS;

  float d0_PV;
  float dz_PV;
  
  // +++ should get rid of these now that we have the isTight, isMedium, and isLoose ID's built into CMSSW
  int numPixelLayers;   //number of pixel layers with valid hits
  int numTrackerLayers; //number of tracker layers with valid hits 
  int numStripLayers;   //number of strip layers with valid hits

  float validFracTracker; //valid fraction of tracker hits

  int numValidMuonHits;
  int numValidPixelHits;
  int numValidTrackerHits;
  int numValidStripHits;
  int numSegmentMatches;
  int numOfMatchedStations;
  // +++

  float trackIsoSumPt;
  float trackIsoSumPtCorr;
  float hcalIso;
  float ecalIso;
  float relCombIso;

  // PF information
  int isPFMuon;

  float pfPt;
  float pfEta;
  float pfPhi;
  
  float sumChargedHadronPtR03; // sum-pt of charged Hadron 
  float sumChargedParticlePtR03; // sum-pt of charged Particles(inludes e/mu) 
  float sumNeutralHadronEtR03;  // sum pt of neutral hadrons
  float sumPhotonEtR03;  // sum pt of PF photons
  float sumPUPtR03;  // sum pt of charged Particles not from PV  (for Pu corrections)

  float sumChargedHadronPtR04; 
  float sumChargedParticlePtR04;
  float sumNeutralHadronEtR04;  
  float sumPhotonEtR04;
  float sumPUPtR04;

  int isHltMatched[3];
  float hltPt[3];
  float hltEta[3];
  float hltPhi[3];

  // +++ should get rid of these now that we have the isTight, isMedium, and isLoose ID's built into CMSSW
  float segmentCompatibility;
  float combinedQualityChi2LocalPosition;
  float combinedQualityTrkKink;
  // +++
  
} _MuonInfo;

//MET
typedef struct {
  float px;
  float py;
  float pt;
  float phi;
  float sumEt;
} _MetInfo;


// pf Jets
typedef struct {
  int   nJets;
  float px[10];
  float py[10];
  float pz[10];
  float pt[10];
  float eta[10];
  float phi[10];
  float mass[10];
  int   charge[10];
  int   partonFlavour[10];
  /////// Energy Fractions //////
  //Charged Hadron
  float chf[10];
  //NeutralHadron
  float nhf[10];
  //Charged EM
  float cef[10];
  //Neutral EM
  float nef[10];
  //Mu
  float muf[10];
  // HF Hadron Fraction
  float hfhf[10];
  // HF EM Fraction
  float hfef[10];
  /////// Multiplicities //////
  // Total Charged Mult
  int cm[10];
  //Charged Hadron Mult
  int chm[10];
  //NeutralHadron Mult
  int nhm[10];
  //Charged EM Mult
  int cem[10];
  //Neutral EM Mult
  int nem[10];
  //Mu Mult
  int mum[10];
  // HF Hadron Mult
  int hfhm[10];
  // HF EM Mult
  int hfem[10];
  // Jet Correction Factor--Above momentum is already corrected!!
  // This factor will return momentum to uncorrected value!!
  float jecFactor[10];
  // Jet Energy Correction Uncertainty
  float jecUnc[10];
  // b-Tag
  float csv[10];
  // Gen Jet Values
  bool genMatched[10];
  float genPx[10];
  float genPy[10];
  float genPz[10];
  float genPt[10];
  float genEta[10];
  float genPhi[10];
  float genMass[10];
  ///// Gen Jet Energy Fractions ///////
  // EM Fraction
  float genEMF[10];
  // Had Fraction
  float genHadF[10];
  // Invisible Fraction
  float genInvF[10];
  // Auxiliary Fraction (Undecayed Sigmas, etc.)
  float genAuxF[10];
  // PUID
  float puid[10];
} _PFJetInfo;

// generator level jets
typedef struct {
  int nJets;
  float px[10];
  float py[10];
  float pz[10];
  float pt[10];
  float eta[10];
  float phi[10];
  float mass[10];
  int   charge[10];
} _GenJetInfo;

// generator level composite Candidate
typedef struct {
  float mass;
  float pt;
  float eta;  // pseudo rapidity
  float y;    // rapidity
  float phi;  // phi
} _genPartInfo;
