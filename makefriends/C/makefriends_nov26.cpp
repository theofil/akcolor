// nov26 fork of makefriends.cpp: inclusive friends (every input event is written, the dijet
// gate is a flag, not a cut), <= 2 jets with constituents, CMS-PF-like particle type per
// constituent, per-jet beam connection (bc). `diff makefriends.cpp makefriends_nov26.cpp`
// lists the nov26 changes; the spec is in nov26.md.
// make makefriends_nov26
// ./makefriends_nov26 VBFZ_herwig.root --genWeight evweight --xs 1.099 --boson Z --totEve 1000

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2Poly.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "ROOT/RDataFrame.hxx"

#include "fastjet/ClusterSequence.hh"

static const int    MAXNPART          = 10000;
static const int    NJETMAX           = 2;     // nov26: the two leading jets only (no 3rd jet)
static const double JET_R             = 0.4;
static const int    NC_MAX            = 80;   // max constituents for pull vector (-1 = all); if >-1, implies pt-ordering
static const bool   USE_KT_JETS       = false; // if true, use kt algorithm with inclusive_jets() instead of anti-kT
static const double JC_PT_MIN         = 0.5;   // pt cut on particles before clustering and on jcs for pull vector
static const double JC_ETA_MAX        = 4.7;   // |eta| cut on particles before clustering
static const double JET_PT_MIN        = 30.0;  // inclusive jet pt threshold
static const double JET_ETA_MAX       = 3.0;   // |eta| cut on reconstructed jets
static const double MJJ_MIN           = 0.0;   // dijet gate: invariant mass threshold (no explicit cut)
static const double BC_DY             = 0.4;   // beam connection: region starts BC_DY in y beyond the jet, toward its own beam
static const double BC_DPHI           = 0.4;   // beam connection: half-width in phi of the "in" strip around the jet
static const double BC_YMAX           = 4.7;   // beam connection: |y| acceptance of the region
static const double PV_A              = 1.0;   // exponent on z = pt_jc/pt_jet in pull-vector weight
static const double PV_B              = 1.0;   // exponent on r in pull-vector weight
static const double PV_C              = 0.0;   // exponent on JET_R in pull-vector weight denominator
static const bool   SAVE_JCS          = true; // if true, save per-jet constituent 4-vectors in tree

// ── helpers ──────────────────────────────────────────────────────────────────

inline double deltaPhi(double p1, double p2) {
    double r = p1 - p2;
    while (r >  M_PI) r -= 2*M_PI;
    while (r < -M_PI) r += 2*M_PI;
    return r;
}

// CMS reco::PFCandidate::ParticleType codes: 1 h (charged hadron), 2 e, 3 mu, 4 gamma,
// 5 h0 (neutral hadron); 0 is kept for empty constituent slots. charge = objects[5] (-1/0/+1).
inline int pfType(int pid, int charge) {
    switch (std::abs(pid)) {
        case 11: return 2;
        case 13: return 3;
        case 22: return 4;
    }
    return (charge != 0) ? 1 : 5;
}

// Beam connection (bc) of a jet: activity between the jet and its own beam.
// Region R = particles i with sign(y_j)*(y_i - y_j) > BC_DY && |y_i| < BC_YMAX (sign(0) = +1)
//   in   = sum pT_i over R with |dphi(i,j)| <  BC_DPHI
//   out  = sum pT_i over R with |dphi(i,j)| >= BC_DPHI   (in + out = scalar sum over R)
//   met  = pT of the vector sum over R
//   dphi = Phi_mpi_pi(phi_j - phi(vector sum)) in [-pi, pi]; -99 if R is empty
struct JetBC { double in = 0., out = 0., met = 0., dphi = -99.; };

JetBC beamConnection(const fastjet::PseudoJet& jet, const std::vector<fastjet::PseudoJet>& parts) {
    JetBC bc;
    const double yj = jet.rapidity(), phij = jet.phi();
    const double s  = (yj < 0) ? -1. : 1.;
    double sx = 0., sy = 0.;
    bool any = false;
    for (const auto& p : parts) {
        const double y = p.rapidity();
        if (!(s * (y - yj) > BC_DY && std::abs(y) < BC_YMAX)) continue;
        any = true;
        if (std::abs(TVector2::Phi_mpi_pi(p.phi() - phij)) < BC_DPHI) bc.in  += p.pt();
        else                                                          bc.out += p.pt();
        sx += p.px();
        sy += p.py();
    }
    bc.met = std::hypot(sx, sy);
    if (any) bc.dphi = TVector2::Phi_mpi_pi(phij - std::atan2(sy, sx));
    return bc;
}

struct JetPV { double pvm = -99., pva = -99., spva = -99.; };

JetPV fillPV(const fastjet::PseudoJet& jet,
             const std::vector<fastjet::PseudoJet>& jcs,
             double a = 1.0, double b = 1.0, double c = 0.0) {
    double pv0 = 0., pv1 = 0.;
    for (const auto& jc : jcs) {
        if (jc.pt() < 1e-3) continue;
        double dY   = jc.rapidity() - jet.rapidity();
        double dPhi = deltaPhi(jc.phi(), jet.phi());
        double r    = std::sqrt(dY*dY + dPhi*dPhi);
        double w    = std::pow(jc.pt() / jet.pt(), a) * std::pow(r / std::pow(JET_R, c), b);
        pv0 += w * dY;
        pv1 += w * dPhi;
    }
    double r = std::sqrt(pv0*pv0 + pv1*pv1);
    JetPV pv; pv.pvm = r;
    if (r > 0) {
        double th = std::atan2(pv1/r, pv0/r);
        pv.pva  = th;
        pv.spva = (jet.rapidity() < 0) ? std::atan2(pv1/r, -pv0/r) : th;
    }
    return pv;
}

TLorentzVector sumP4(const std::vector<fastjet::PseudoJet>& jets, int n) {
    TLorentzVector v(0,0,0,0);
    for (int i = 0; i < n && i < (int)jets.size(); ++i)
        v += TLorentzVector(jets[i].px(), jets[i].py(), jets[i].pz(), jets[i].e());
    return v;
}

// ── argument parsing ─────────────────────────────────────────────────────────

struct Args {
    std::string inputfile, genWeight = "genWeight", output = "", boson = "";
    double xs = 1.0;
    long   totEve = -1, skip = 0, jobId = -1;
    int    debug  = 0;
};

Args parseArgs(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " file.root [--genWeight W] [--xs X] [--totEve N]"
                     " [--skip N] [--jobId N] [--output out.root] [--boson H|W|Z] [--debug N]\n";
        exit(1);
    }
    Args a; a.inputfile = argv[1];
    for (int i = 2; i < argc-1; ++i) {
        std::string f = argv[i];
        if      (f == "--genWeight") a.genWeight = argv[++i];
        else if (f == "--xs"       ) a.xs        = std::stod(argv[++i]);
        else if (f == "--totEve"   ) a.totEve    = std::stol(argv[++i]);
        else if (f == "--skip"     ) a.skip      = std::stol(argv[++i]);
        else if (f == "--jobId"    ) a.jobId     = std::stol(argv[++i]);
        else if (f == "--output"   ) a.output    = argv[++i];
        else if (f == "--boson"    ) a.boson     = argv[++i];
        else if (f == "--debug"    ) a.debug     = std::stoi(argv[++i]);
    }
    if (!a.boson.empty() && a.boson != "H" && a.boson != "W" && a.boson != "Z") {
        std::cerr << "--boson must be H, W or Z (got " << a.boson << ")\n";
        exit(1);
    }
    return a;
}

// ── main ──────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    Args args = parseArgs(argc, argv);

    // ── open input ───────────────────────────────────────────────────────────
    TFile* ifile = TFile::Open(args.inputfile.c_str());
    if (!ifile || ifile->IsZombie()) { std::cerr << "Cannot open " << args.inputfile << "\n"; return 1; }

    const std::vector<std::string> knownTrees = {"Events","events","ntuple/tree","tree","Data"};
    TTree* ttree = nullptr;
    for (const auto& tn : knownTrees) {
        ttree = dynamic_cast<TTree*>(ifile->Get(tn.c_str()));
        if (ttree) { std::cout << "Found tree: " << tn << "\n"; break; }
    }
    if (!ttree) { std::cerr << "No known TTree found\n"; return 1; }

    Long64_t nentries = ttree->GetEntries();
    std::cout << args.inputfile << " with " << nentries << " entries\n";

    // ── sum of weights ────────────────────────────────────────────────────────
    ROOT::RDataFrame df(ttree->GetName(), args.inputfile.c_str());
    double sumWtot = df.Sum<double>(args.genWeight.c_str()).GetValue();
    std::cout << "sumWtot = " << sumWtot << "\n";

    // ── output file + tree ────────────────────────────────────────────────────
    std::string outname = args.inputfile;
    size_t sl = outname.rfind('/'); if (sl != std::string::npos) outname = outname.substr(sl+1);
    if (outname.size() >= 5) outname = outname.substr(0, outname.size()-5);
    outname += "slimfriend";
    if (args.jobId != -1) outname += "_" + std::to_string(args.jobId);
    outname += ".root";
    if (!args.output.empty()) outname = args.output;

    TFile* ofile = TFile::Open(outname.c_str(), "RECREATE");
    TTree* otree = new TTree("events", "events");
    TH1D*  hSumW = new TH1D("hSumW","hSumW",1,0,1); hSumW->Sumw2();

    // leading-jet |SPVA| (20 bins, 0–π), filled for dijet-gate events (gate == 1)
    TH1F* hLeadJetSPVA = new TH1F("hLeadJetSPVA","hLeadJetSPVA", 20, 0., M_PI);
    hLeadJetSPVA->Sumw2();

    // leading-jet PVM |t⃗| (20 bins, 0–0.06), matching pvm.py linspace(0,0.06,21)
    TH1F* hLeadJetSPVM = new TH1F("hLeadJetSPVM","hLeadJetSPVM", 20, 0., 0.06);
    hLeadJetSPVM->Sumw2();
    hLeadJetSPVM->GetXaxis()->SetTitle("|#vec{t}|");

    TH1D* hNJets = new TH1D("hNJets","hNJets", 10, 0, 10);
    hNJets->Sumw2();
    hNJets->GetXaxis()->SetTitle("N_{jets}");

    TH1F* hJetPt1 = new TH1F("hJetPt1","hJetPt1", 100, 0., 500.);
    hJetPt1->Sumw2();
    hJetPt1->GetXaxis()->SetTitle("leading jet p_{T} [GeV]");

    // TH2Poly jetShape — 64 phi × 30 rho annular sectors
    TH2Poly* jetShape = new TH2Poly("jetShape","jetShape",-0.5,0.5,-0.5,0.5);
    jetShape->Sumw2();
    {
        const int nArc = 12, nPhi = 128, nRho = 60;
        for (int iPhi = 0; iPhi < nPhi; ++iPhi) {
            double p1 = 2*M_PI * iPhi / nPhi;
            double p2 = 2*M_PI * (iPhi+1) / nPhi;
            for (int iRho = 0; iRho < nRho; ++iRho) {
                double r1 = 0.5 * iRho / nRho;
                double r2 = 0.5 * (iRho+1) / nRho;
                std::vector<double> xv, yv;
                for (int k = 0; k < nArc; ++k) {
                    double phi = p1 + (p2-p1)*k/(nArc-1);
                    xv.push_back(r2*std::cos(phi)); yv.push_back(r2*std::sin(phi));
                }
                for (int k = nArc-1; k >= 0; --k) {
                    double phi = p1 + (p2-p1)*k/(nArc-1);
                    xv.push_back(r1*std::cos(phi)); yv.push_back(r1*std::sin(phi));
                }
                jetShape->AddBin((int)xv.size(), xv.data(), yv.data());
            }
        }
    }

    // TH2Poly jetShapeRaw — same binning as jetShape, weight = pt_i/pt_jet (no radial factor)
    TH2Poly* jetShapeRaw = new TH2Poly("jetShapeRaw","jetShapeRaw",-0.5,0.5,-0.5,0.5);
    jetShapeRaw->Sumw2();
    {
        const int nArc = 12, nPhi = 128, nRho = 60;
        for (int iPhi = 0; iPhi < nPhi; ++iPhi) {
            double p1 = 2*M_PI * iPhi / nPhi;
            double p2 = 2*M_PI * (iPhi+1) / nPhi;
            for (int iRho = 0; iRho < nRho; ++iRho) {
                double r1 = 0.5 * iRho / nRho;
                double r2 = 0.5 * (iRho+1) / nRho;
                std::vector<double> xv, yv;
                for (int k = 0; k < nArc; ++k) {
                    double phi = p1 + (p2-p1)*k/(nArc-1);
                    xv.push_back(r2*std::cos(phi)); yv.push_back(r2*std::sin(phi));
                }
                for (int k = nArc-1; k >= 0; --k) {
                    double phi = p1 + (p2-p1)*k/(nArc-1);
                    xv.push_back(r1*std::cos(phi)); yv.push_back(r1*std::sin(phi));
                }
                jetShapeRaw->AddBin((int)xv.size(), xv.data(), yv.data());
            }
        }
    }

    // output branches
    Long64_t o_evIdx;
    Int_t   o_nJets, o_nJetsAll;
    Bool_t  o_gate;
    Float_t o_jetPt[NJETMAX], o_jetEta[NJETMAX], o_jetPhi[NJETMAX], o_jetM[NJETMAX];
    Float_t o_jetSPVA[NJETMAX], o_jetPVA[NJETMAX], o_jetPVM[NJETMAX];
    Float_t o_jetBCin[NJETMAX], o_jetBCout[NJETMAX], o_jetBCmet[NJETMAX], o_jetBCdphi[NJETMAX];
    Int_t   o_jetNC[NJETMAX];
    Float_t o_kWeight;
    Int_t   o_leadJetIndex;
    Float_t o_mjj, o_dYjj, o_dPhijj, o_ptjj;
    Float_t o_bosonM, o_bosonY, o_bosonEta, o_bosonPt, o_bosonPhi;

    otree->Branch("evIdx",    &o_evIdx,     "evIdx/L");
    otree->Branch("nJets",    &o_nJets,     "nJets/I");
    otree->Branch("nJetsAll", &o_nJetsAll,  "nJetsAll/I");
    otree->Branch("gate",     &o_gate,      "gate/O");
    otree->Branch("jetPt0",   &o_jetPt[0],  "jetPt0/F");
    otree->Branch("jetPt1",   &o_jetPt[1],  "jetPt1/F");
    otree->Branch("jetEta0",  &o_jetEta[0], "jetEta0/F");
    otree->Branch("jetEta1",  &o_jetEta[1], "jetEta1/F");
    otree->Branch("jetPhi0",  &o_jetPhi[0], "jetPhi0/F");
    otree->Branch("jetPhi1",  &o_jetPhi[1], "jetPhi1/F");
    otree->Branch("jetM0",    &o_jetM[0],   "jetM0/F");
    otree->Branch("jetM1",    &o_jetM[1],   "jetM1/F");
    otree->Branch("jetSPVA0", &o_jetSPVA[0],"jetSPVA0/F");
    otree->Branch("jetSPVA1", &o_jetSPVA[1],"jetSPVA1/F");
    otree->Branch("jetPVA0",  &o_jetPVA[0], "jetPVA0/F");
    otree->Branch("jetPVA1",  &o_jetPVA[1], "jetPVA1/F");
    otree->Branch("jetPVM0",  &o_jetPVM[0], "jetPVM0/F");
    otree->Branch("jetPVM1",  &o_jetPVM[1], "jetPVM1/F");
    otree->Branch("jetNC0",   &o_jetNC[0],  "jetNC0/I");
    otree->Branch("jetNC1",   &o_jetNC[1],  "jetNC1/I");
    otree->Branch("jetBCin0",   &o_jetBCin[0],   "jetBCin0/F");
    otree->Branch("jetBCin1",   &o_jetBCin[1],   "jetBCin1/F");
    otree->Branch("jetBCout0",  &o_jetBCout[0],  "jetBCout0/F");
    otree->Branch("jetBCout1",  &o_jetBCout[1],  "jetBCout1/F");
    otree->Branch("jetBCmet0",  &o_jetBCmet[0],  "jetBCmet0/F");
    otree->Branch("jetBCmet1",  &o_jetBCmet[1],  "jetBCmet1/F");
    otree->Branch("jetBCdphi0", &o_jetBCdphi[0], "jetBCdphi0/F");
    otree->Branch("jetBCdphi1", &o_jetBCdphi[1], "jetBCdphi1/F");
    otree->Branch("kWeight",      &o_kWeight,      "kWeight/F");
    otree->Branch("leadJetIndex", &o_leadJetIndex, "leadJetIndex/I");
    otree->Branch("mjj",          &o_mjj,          "mjj/F");
    otree->Branch("dYjj",         &o_dYjj,         "dYjj/F");
    otree->Branch("dPhijj",       &o_dPhijj,       "dPhijj/F");
    otree->Branch("ptjj",         &o_ptjj,         "ptjj/F");
    otree->Branch("bosonM",       &o_bosonM,       "bosonM/F");
    otree->Branch("bosonY",       &o_bosonY,       "bosonY/F");
    otree->Branch("bosonEta",     &o_bosonEta,     "bosonEta/F");
    otree->Branch("bosonPt",      &o_bosonPt,      "bosonPt/F");
    otree->Branch("bosonPhi",     &o_bosonPhi,     "bosonPhi/F");

    Float_t o_jcsPt     [NJETMAX][NC_MAX];
    Float_t o_jcsDEta   [NJETMAX][NC_MAX];
    Float_t o_jcsDEtaRaw[NJETMAX][NC_MAX];
    Float_t o_jcsDPhi   [NJETMAX][NC_MAX];
    Float_t o_jcsM   [NJETMAX][NC_MAX];
    Float_t o_jcsW   [NJETMAX][NC_MAX];
    Char_t  o_jcsPFType[NJETMAX][NC_MAX];  // pfType() code, 0 = empty slot
    Char_t  o_jcsCharge[NJETMAX][NC_MAX];  // -1/0/+1
    if (SAVE_JCS) {
        std::string nc = std::to_string(NC_MAX);
        otree->Branch("jcsPt0",   o_jcsPt[0],   ("jcsPt0["   + nc + "]/F").c_str());
        otree->Branch("jcsPt1",   o_jcsPt[1],   ("jcsPt1["   + nc + "]/F").c_str());
        otree->Branch("jcsDEta0",    o_jcsDEta[0],    ("jcsDEta0["    + nc + "]/F").c_str());
        otree->Branch("jcsDEta1",    o_jcsDEta[1],    ("jcsDEta1["    + nc + "]/F").c_str());
        otree->Branch("jcsDEtaRaw0", o_jcsDEtaRaw[0], ("jcsDEtaRaw0[" + nc + "]/F").c_str());
        otree->Branch("jcsDEtaRaw1", o_jcsDEtaRaw[1], ("jcsDEtaRaw1[" + nc + "]/F").c_str());
        otree->Branch("jcsDPhi0",    o_jcsDPhi[0],    ("jcsDPhi0["    + nc + "]/F").c_str());
        otree->Branch("jcsDPhi1", o_jcsDPhi[1], ("jcsDPhi1[" + nc + "]/F").c_str());
        otree->Branch("jcsM0",    o_jcsM[0],    ("jcsM0["    + nc + "]/F").c_str());
        otree->Branch("jcsM1",    o_jcsM[1],    ("jcsM1["    + nc + "]/F").c_str());
        otree->Branch("jcsW0",    o_jcsW[0],    ("jcsW0["    + nc + "]/F").c_str());
        otree->Branch("jcsW1",    o_jcsW[1],    ("jcsW1["    + nc + "]/F").c_str());
        otree->Branch("jcsPFType0", o_jcsPFType[0], ("jcsPFType0[" + nc + "]/B").c_str());
        otree->Branch("jcsPFType1", o_jcsPFType[1], ("jcsPFType1[" + nc + "]/B").c_str());
        otree->Branch("jcsCharge0", o_jcsCharge[0], ("jcsCharge0[" + nc + "]/B").c_str());
        otree->Branch("jcsCharge1", o_jcsCharge[1], ("jcsCharge1[" + nc + "]/B").c_str());
    }

    // input branches
    Int_t    numparticles;
    Double_t objects_arr[8][MAXNPART];
    Double_t evweight;

    ttree->SetBranchAddress("numparticles",       &numparticles);
    ttree->SetBranchAddress("objects",             objects_arr);
    ttree->SetBranchAddress(args.genWeight.c_str(), &evweight);

    // fastjet setup
    fastjet::JetDefinition jetdef(USE_KT_JETS ? fastjet::kt_algorithm
                                                     : fastjet::antikt_algorithm, JET_R);

    long alleve = 0, ngate = 0;
    double sumW = 0., sumW2 = 0.;
    long start = args.skip;
    long end   = (args.totEve < 0) ? nentries : std::min((long)nentries, start + args.totEve);

    // ── event loop ───────────────────────────────────────────────────────────
    for (long iev = start; iev < end; ++iev) {
        ttree->GetEntry(iev);
        if (iev % 10000 == 0) std::cout << "event " << iev << "\n";

        // reset output (zero-pad so unused jet/constituent slots are 0)
        o_evIdx = iev;
        o_nJets = 0;
        o_nJetsAll = 0;
        o_gate = false;
        o_leadJetIndex = -1;
        o_mjj = o_dYjj = o_dPhijj = o_ptjj = -99.f;
        memset(o_jetPt,   0, sizeof(o_jetPt));
        memset(o_jetEta,  0, sizeof(o_jetEta));
        memset(o_jetPhi,  0, sizeof(o_jetPhi));
        memset(o_jetM,    0, sizeof(o_jetM));
        memset(o_jetSPVA, 0, sizeof(o_jetSPVA));
        memset(o_jetPVA,  0, sizeof(o_jetPVA));
        memset(o_jetPVM,  0, sizeof(o_jetPVM));
        memset(o_jetBCin,   0, sizeof(o_jetBCin));
        memset(o_jetBCout,  0, sizeof(o_jetBCout));
        memset(o_jetBCmet,  0, sizeof(o_jetBCmet));
        memset(o_jetBCdphi, 0, sizeof(o_jetBCdphi));
        memset(o_jetNC,   0, sizeof(o_jetNC));
        if (SAVE_JCS) {
            memset(o_jcsPt,   0, sizeof(o_jcsPt));
            memset(o_jcsDEta,    0, sizeof(o_jcsDEta));
            memset(o_jcsDEtaRaw, 0, sizeof(o_jcsDEtaRaw));
            memset(o_jcsDPhi,    0, sizeof(o_jcsDPhi));
            memset(o_jcsM,    0, sizeof(o_jcsM));
            memset(o_jcsW,    0, sizeof(o_jcsW));
            memset(o_jcsPFType,  0, sizeof(o_jcsPFType));
            memset(o_jcsCharge,  0, sizeof(o_jcsCharge));
        }
        o_kWeight = (Float_t)(evweight * args.xs / sumWtot);
        sumW  += o_kWeight;
        sumW2 += o_kWeight * o_kWeight;
        hSumW->Fill(0., o_kWeight);

        // particle loop → build momtocluster (+ collect boson decay candidates)
        std::vector<fastjet::PseudoJet> momfj;
        momfj.reserve(numparticles);
        std::vector<TLorentzVector> b_higgs, b_tauM, b_tauP, b_nuTau, b_nuTauBar, b_nuE, b_nuEBar;
        for (int mm = 0; mm < numparticles; ++mm) {
            double E  = objects_arr[0][mm];
            double px = objects_arr[1][mm];
            double py = objects_arr[2][mm];
            double pz = objects_arr[3][mm];
            int    pid = (int)objects_arr[4][mm];
            if (!args.boson.empty()) {
                if      (pid ==  25) b_higgs.emplace_back(px, py, pz, E);
                else if (pid ==  15) b_tauM.emplace_back(px, py, pz, E);
                else if (pid == -15) b_tauP.emplace_back(px, py, pz, E);
                else if (pid ==  16) b_nuTau.emplace_back(px, py, pz, E);
                else if (pid == -16) b_nuTauBar.emplace_back(px, py, pz, E);
                else if (pid ==  12) b_nuE.emplace_back(px, py, pz, E);
                else if (pid == -12) b_nuEBar.emplace_back(px, py, pz, E);
            }
            if (std::abs(pid) == 25 || std::abs(pid) == 12 || std::abs(pid) == 14 ||
                             std::abs(pid) == 16 || std::abs(pid) == 15) continue;
            fastjet::PseudoJet pj(px, py, pz, E);
            pj.set_user_index(mm);  // constituent -> objects[4][mm] (pid), objects[5][mm] (charge)
            if (pj.perp() > JC_PT_MIN && std::abs(pj.eta()) < JC_ETA_MAX)
                momfj.push_back(pj);
        }

        // boson 4-vector: H = stable pid 25; W = flavor-matched (τ,ν_τ) pair,
        // Z = (ν_e,ν̄_e) pair — each with mass closest to the pole
        o_bosonM = o_bosonY = o_bosonEta = o_bosonPt = o_bosonPhi = -99.f;
        if (!args.boson.empty()) {
            TLorentzVector bp4;
            bool  found = false;
            double bestDiff = 1e18;
            auto tryPairs = [&](const std::vector<TLorentzVector>& A,
                                const std::vector<TLorentzVector>& B, double pole) {
                for (const auto& a : A) for (const auto& b : B) {
                    TLorentzVector s = a + b;
                    double d = std::abs(s.M() - pole);
                    if (d < bestDiff) { bestDiff = d; bp4 = s; found = true; }
                }
            };
            if (args.boson == "H") {
                if (!b_higgs.empty()) { bp4 = b_higgs[0]; found = true; }
            } else if (args.boson == "W") {
                tryPairs(b_tauM, b_nuTauBar, 80.379);   // W⁻ → τ⁻ ν̄_τ
                tryPairs(b_tauP, b_nuTau,    80.379);   // W⁺ → τ⁺ ν_τ
            } else if (args.boson == "Z") {
                tryPairs(b_nuE, b_nuEBar, 91.1876);     // Z → ν_e ν̄_e
            }
            if (found) {
                o_bosonM   = (Float_t)std::round(bp4.M()        * 100) / 100.f;
                o_bosonY   = (Float_t)std::round(bp4.Rapidity() * 1000) / 1000.f;
                o_bosonEta = (Float_t)std::round(bp4.Eta()      * 1000) / 1000.f;
                o_bosonPt  = (Float_t)std::round(bp4.Pt()       * 100) / 100.f;
                o_bosonPhi = (Float_t)std::round(bp4.Phi()      * 1000) / 1000.f;
            }
        }

        // cluster
        fastjet::ClusterSequence cs(momfj, jetdef);
        std::vector<fastjet::PseudoJet> jets;
        if (USE_KT_JETS) {
            for (const auto& j : fastjet::sorted_by_pt(cs.inclusive_jets()))
                if (std::abs(j.eta()) < JET_ETA_MAX) jets.push_back(j);
        } else {
            for (const auto& j : fastjet::sorted_by_pt(cs.inclusive_jets(JET_PT_MIN)))
                if (std::abs(j.eta()) < JET_ETA_MAX) jets.push_back(j);
        }

        // pull vectors for all jets
        std::vector<JetPV> pvs(jets.size());
        std::vector<JetBC> bcs(jets.size());
        std::vector<std::vector<fastjet::PseudoJet>> all_jcs(jets.size());
        for (int ij = 0; ij < (int)jets.size(); ++ij) {
            std::vector<fastjet::PseudoJet> jcs;
            if (USE_KT_JETS) {
                fastjet::PseudoJet sub1, sub2;
                if (jets[ij].has_parents(sub1, sub2))
                    jcs = {sub1, sub2};
                else
                    jcs = jets[ij].constituents(); // fallback: single-particle jet has no parents
            } else {
                for (const auto& jc : jets[ij].constituents())
                    if (jc.pt() > JC_PT_MIN) jcs.push_back(jc);
                if (NC_MAX > -1) {
                    std::sort(jcs.begin(), jcs.end(),
                        [](const fastjet::PseudoJet& a, const fastjet::PseudoJet& b){ return a.pt() > b.pt(); });
                    if ((int)jcs.size() > NC_MAX) jcs.resize(NC_MAX);
                }
            }
            if (SAVE_JCS) all_jcs[ij] = jcs;
            pvs[ij] = fillPV(jets[ij], jcs, PV_A, PV_B, PV_C);
            if (ij < NJETMAX) bcs[ij] = beamConnection(jets[ij], momfj);  // clustering particles
        }

        // per-event inclusive histograms (no dijet selection)
        hNJets->Fill((float)jets.size(), o_kWeight);
        if (!jets.empty()) hJetPt1->Fill(jets[0].pt(), o_kWeight);

        // nov26: every event is written. Dijet variables need >= 2 jets (else -99); the dijet
        // gate (>= 2 jets, mjj > MJJ_MIN, opposite-eta leading pair) is a flag, not a cut.
        o_nJetsAll = (Int_t)jets.size();
        bool passNJets = (int)jets.size() >= 2;
        if (passNJets) {
            TLorentzVector jj = sumP4(jets, 2);
            o_mjj    = (Float_t)jj.M();
            o_dYjj   = (Float_t)(jets[0].rapidity() - jets[1].rapidity());
            o_dPhijj = (Float_t)deltaPhi(jets[0].phi(), jets[1].phi());
            o_ptjj   = (Float_t)jj.Pt();
            o_gate   = jj.M() > MJJ_MIN && jets[0].eta() * jets[1].eta() < 0.;
        }

        if (o_gate) {
            ++ngate;

            // |SPVA| and PVM for leading and subleading jet — two fills per event
            for (int ij = 0; ij < 2; ++ij) {
                if (pvs[ij].spva > -90.)
                    hLeadJetSPVA->Fill(std::abs(pvs[ij].spva), o_kWeight);
                if (pvs[ij].pvm > -90.)
                    hLeadJetSPVM->Fill(pvs[ij].pvm, o_kWeight);
            }

            // fill jetShape from both leading jets (sign-flip dY for backward jet)
            for (int ij = 0; ij < 2; ++ij) {
                for (const auto& jc : jets[ij].constituents()) {
                    double dY   = jc.rapidity() - jets[ij].rapidity();
                    double dPhi = deltaPhi(jc.phi(), jets[ij].phi());
                    double r    = std::sqrt(dY*dY + dPhi*dPhi);
                    double w    = std::pow(jc.pt() / jets[ij].pt(), PV_A) * std::pow(r / std::pow(JET_R, PV_C), PV_B);
                    double wRaw = jc.pt() / jets[ij].pt();
                    double dYs  = (jets[ij].rapidity() < 0) ? -dY : dY;
                    jetShape->Fill(dYs, dPhi, w);
                    jetShapeRaw->Fill(dYs, dPhi, wRaw);
                }
            }

        }
        {
            // fill output tree once per event (every event; the two leading jets only)
            int nmax = std::min((int)jets.size(), NJETMAX);
            o_nJets = nmax;
            // slot randomization when both slots are filled: odd iev → leading at [0], even →
            // leading at [1]; a single jet stays in slot 0; leadJetIndex = -1 without jets
            bool swap = (nmax == 2) && (iev % 2 == 0);
            if (nmax > 0) o_leadJetIndex = swap ? 1 : 0;
            for (int i = 0; i < nmax; ++i) {
                int si = swap ? 1 - i : i;
                o_jetPt[si]   = (Float_t)std::round(jets[i].pt()  * 100) / 100.f;
                o_jetEta[si]  = (Float_t)std::round(jets[i].eta() * 1000) / 1000.f;
                o_jetPhi[si]  = (Float_t)std::round(jets[i].phi() * 1000) / 1000.f;
                o_jetM[si]    = (Float_t)std::round(jets[i].m()   * 100) / 100.f;
                o_jetSPVA[si] = (Float_t)pvs[i].spva;
                o_jetPVA[si]  = (Float_t)pvs[i].pva;
                o_jetPVM[si]  = (Float_t)pvs[i].pvm;
                o_jetNC[si]   = (Int_t)jets[i].constituents().size();
                o_jetBCin[si]   = (Float_t)bcs[i].in;
                o_jetBCout[si]  = (Float_t)bcs[i].out;
                o_jetBCmet[si]  = (Float_t)bcs[i].met;
                o_jetBCdphi[si] = (Float_t)bcs[i].dphi;
                if (SAVE_JCS) {
                    const auto& jcs_i = all_jcs[i];
                    int nk = (int)jcs_i.size();
                    for (int k = 0; k < nk; ++k) {
                        double dEta_raw = jcs_i[k].eta() - jets[i].eta();
                        double dEtas    = (jets[i].eta() < 0) ? -dEta_raw : dEta_raw;
                        double dPhi     = deltaPhi(jcs_i[k].phi(), jets[i].phi());
                        double dY_raw   = jcs_i[k].rapidity() - jets[i].rapidity();
                        double r        = std::sqrt(dY_raw*dY_raw + dPhi*dPhi);
                        double w        = std::pow(jcs_i[k].pt() / jets[i].pt(), PV_A)
                                        * std::pow(r / std::pow(JET_R, PV_C), PV_B);
                        o_jcsPt  [si][k] = (Float_t)jcs_i[k].pt();
                        o_jcsDEta   [si][k] = (Float_t)dEtas;
                        o_jcsDEtaRaw[si][k] = (Float_t)dEta_raw;
                        o_jcsDPhi[si][k] = (Float_t)dPhi;
                        o_jcsM   [si][k] = (Float_t)jcs_i[k].m();
                        o_jcsW   [si][k] = (Float_t)w;
                        int idx = jcs_i[k].user_index();  // input particle index (particle loop)
                        if (idx >= 0) {
                            int q = (int)objects_arr[5][idx];
                            o_jcsPFType[si][k] = (Char_t)pfType((int)objects_arr[4][idx], q);
                            o_jcsCharge[si][k] = (Char_t)q;
                        }
                    }
                }
            }
            otree->Fill();
            ++alleve;
        }
    }

    std::cout << "hSumW integral   " << hSumW->Integral() << "\n";
    std::cout << "events in tree   " << alleve            << "\n";
    std::cout << "gate events      " << ngate             << "\n";
    std::cout << "sumW             " << sumW              << "\n";
    std::cout << "sumW2            " << sumW2             << "\n";

    ofile->cd();
    otree->Write();
    hSumW->Write();
    jetShape->Write();
    jetShapeRaw->Write();
    hLeadJetSPVA->Write();
    hLeadJetSPVM->Write();
    ofile->Write();

    ofile->Close();
    ifile->Close();
    return 0;
}
