// The size of the different arrays
const unsigned int N_TRIGGER_INFO = 6;
const unsigned int N_VERTEX_INFO = 20;

const unsigned int N_MU_INFO = 10;
const unsigned int N_ELECTRON_INFO = 10;
const unsigned int N_TAU_INFO = 10;
const unsigned int N_JET_INFO = 10;


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
  int isValid[N_VERTEX_INFO];
  float x[N_VERTEX_INFO];	
  float y[N_VERTEX_INFO];	
  float z[N_VERTEX_INFO];	
  float xErr[N_VERTEX_INFO];	
  float yErr[N_VERTEX_INFO];	
  float zErr[N_VERTEX_INFO];	
  float chi2[N_VERTEX_INFO];
  int ndf[N_VERTEX_INFO];
  float normChi2[N_VERTEX_INFO];
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

  int nMuons;
  int nSelectedMuons;
  int nMuonPairs;

  int isTracker[N_MU_INFO];
  int isStandAlone[N_MU_INFO];
  int isGlobal[N_MU_INFO];

  int isTightMuon[N_MU_INFO];
  int isMediumMuon[N_MU_INFO];
  int isLooseMuon[N_MU_INFO];

  int charge[N_MU_INFO];
  float pt[N_MU_INFO];
  float ptErr[N_MU_INFO];
  float eta[N_MU_INFO];
  float phi[N_MU_INFO];

  float trkPt[N_MU_INFO];
  float trkPtErr[N_MU_INFO];
  float trketa[N_MU_INFO];
  float trkPhi[N_MU_INFO];

  float d0_BS[N_MU_INFO];
  float dz_BS[N_MU_INFO];

  float d0_PV[N_MU_INFO];
  float dz_PV[N_MU_INFO];

  float trackIsoSumPt[N_MU_INFO];
  float trackIsoSumPtCorr[N_MU_INFO];
  float hcalIso[N_MU_INFO];
  float ecalIso[N_MU_INFO];
  float relCombIso[N_MU_INFO];

  // PF information
  int isPFMuon[N_MU_INFO];

  float pfPt[N_MU_INFO];
  float pfEta[N_MU_INFO];
  float pfPhi[N_MU_INFO];
  
  float sumChargedHadronPtR03[N_MU_INFO];   // sum-pt of charged Hadron 
  float sumChargedParticlePtR03[N_MU_INFO]; // sum-pt of charged Particles(inludes e/mu) 
  float sumNeutralHadronEtR03[N_MU_INFO];   // sum pt of neutral hadrons
  float sumPhotonEtR03[N_MU_INFO];          // sum pt of PF photons
  float sumPUPtR03[N_MU_INFO];              // sum pt of charged Particles not from PV  (for Pu corrections)

  float sumChargedHadronPtR04[N_MU_INFO]; 
  float sumChargedParticlePtR04[N_MU_INFO];
  float sumNeutralHadronEtR04[N_MU_INFO];  
  float sumPhotonEtR04[N_MU_INFO];
  float sumPUPtR04[N_MU_INFO];

  int isHltMatched[N_MU_INFO][N_TRIGGER_INFO];
  float hltPt[N_MU_INFO][N_TRIGGER_INFO];
  float hltEta[N_MU_INFO][N_TRIGGER_INFO];
  float hltPhi[N_MU_INFO][N_TRIGGER_INFO];

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
  float px[N_JET_INFO];
  float py[N_JET_INFO];
  float pz[N_JET_INFO];
  float pt[N_JET_INFO];
  float eta[N_JET_INFO];
  float phi[N_JET_INFO];
  float mass[N_JET_INFO];
  int   charge[N_JET_INFO];
  int   partonFlavour[N_JET_INFO];
  /////// Energy Fractions //////
  //Charged Hadron
  float chf[N_JET_INFO];
  //NeutralHadron
  float nhf[N_JET_INFO];
  //Charged EM
  float cef[N_JET_INFO];
  //Neutral EM
  float nef[N_JET_INFO];
  //Mu
  float muf[N_JET_INFO];
  // HF Hadron Fraction
  float hfhf[N_JET_INFO];
  // HF EM Fraction
  float hfef[N_JET_INFO];
  /////// Multiplicities //////
  // Total Charged Mult
  int cm[N_JET_INFO];
  //Charged Hadron Mult
  int chm[N_JET_INFO];
  //NeutralHadron Mult
  int nhm[N_JET_INFO];
  //Charged EM Mult
  int cem[N_JET_INFO];
  //Neutral EM Mult
  int nem[N_JET_INFO];
  //Mu Mult
  int mum[N_JET_INFO];
  // HF Hadron Mult
  int hfhm[N_JET_INFO];
  // HF EM Mult
  int hfem[N_JET_INFO];
  // Jet Correction Factor--Above momentum is already corrected!!
  // This factor will return momentum to uncorrected value!!
  float jecFactor[N_JET_INFO];
  // Jet Energy Correction Uncertainty
  float jecUnc[N_JET_INFO];
  // b-Tag
  float csv[N_JET_INFO];
  // Gen Jet Values
  bool genMatched[N_JET_INFO];
  float genPx[N_JET_INFO];
  float genPy[N_JET_INFO];
  float genPz[N_JET_INFO];
  float genPt[N_JET_INFO];
  float genEta[N_JET_INFO];
  float genPhi[N_JET_INFO];
  float genMass[N_JET_INFO];
  ///// Gen Jet Energy Fractions ///////
  // EM Fraction
  float genEMF[N_JET_INFO];
  // Had Fraction
  float genHadF[N_JET_INFO];
  // Invisible Fraction
  float genInvF[N_JET_INFO];
  // Auxiliary Fraction (Undecayed Sigmas, etc.)
  float genAuxF[N_JET_INFO];
  // PUID
  float puid[N_JET_INFO];
} _PFJetInfo;

// generator level jets
typedef struct {
  int nJets;
  float px[N_JET_INFO];
  float py[N_JET_INFO];
  float pz[N_JET_INFO];
  float pt[N_JET_INFO];
  float eta[N_JET_INFO];
  float phi[N_JET_INFO];
  float mass[N_JET_INFO];
  int   charge[N_JET_INFO];
} _GenJetInfo;

// generator level composite Candidate
typedef struct {
  int charge;
  float mass;
  float pt;
  float eta;  // pseudo rapidity
  float y;    // rapidity
  float phi;  // phi
} _genPartInfo;
