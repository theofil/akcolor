// g++ $(root-config --cflags) $(fastjet-config --cxxflags) -O2 -o makefriends makefriends.cpp $(root-config --libs) $(fastjet-config --libs)
// ./makefriends VBFHtoInv.root --genWeight evweight --xs 3.901 --totEve 1000

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
#include "ROOT/RDataFrame.hxx"

#include "fastjet/ClusterSequence.hh"

static const int    MAXNPART          = 10000;
static const int    NJETMAX           = 3;
static const double JET_R             = 0.4;
static const int    NC_MAX            = 80;   // max constituents for pull vector (-1 = all); if >-1, implies pt-ordering
static const bool   USE_KT_JETS       = false; // if true, use kt algorithm with inclusive_jets() instead of anti-kT
static const double JC_PT_MIN         = 0.0;   // pt cut on particles before clustering and on jcs for pull vector
static const double JC_ETA_MAX        = 20.0;  // |eta| cut on particles before clustering
static const double JET_PT_MIN        = 30.0;  // inclusive jet pt threshold
static const double JET_ETA_MAX       = 3.0;   // |eta| cut on reconstructed jets
static const double MJJ_MIN           = 400.0; // dijet invariant mass threshold
static const double DYJJ_MIN          = 2.5;   // |dYjj| threshold (leading pair rapidity separation)
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

    // leading-jet |SPVA| (20 bins, 0–π), filled after mjj>MJJ_MIN + opposite-eta selection
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
    Int_t   o_nJets;
    Float_t o_jetPt[NJETMAX], o_jetEta[NJETMAX], o_jetPhi[NJETMAX], o_jetM[NJETMAX];
    Float_t o_jetSPVA[NJETMAX], o_jetPVA[NJETMAX], o_jetPVM[NJETMAX];
    Int_t   o_jetNC[NJETMAX];
    Float_t o_kWeight;
    Int_t   o_leadJetIndex;
    Float_t o_mjj, o_dYjj, o_dPhijj, o_ptjj;
    Float_t o_bosonM, o_bosonY, o_bosonEta, o_bosonPt, o_bosonPhi;

    otree->Branch("nJets",    &o_nJets,     "nJets/I");
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
    otree->Branch("jetPt2",   &o_jetPt[2],  "jetPt2/F");
    otree->Branch("jetEta2",  &o_jetEta[2], "jetEta2/F");
    otree->Branch("jetPhi2",  &o_jetPhi[2], "jetPhi2/F");
    otree->Branch("jetM2",    &o_jetM[2],   "jetM2/F");
    otree->Branch("jetSPVA2", &o_jetSPVA[2],"jetSPVA2/F");
    otree->Branch("jetPVA2",  &o_jetPVA[2], "jetPVA2/F");
    otree->Branch("jetPVM2",  &o_jetPVM[2], "jetPVM2/F");
    otree->Branch("jetNC2",   &o_jetNC[2],  "jetNC2/I");
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
        otree->Branch("jcsPt2",      o_jcsPt[2],      ("jcsPt2["      + nc + "]/F").c_str());
        otree->Branch("jcsDEta2",    o_jcsDEta[2],    ("jcsDEta2["    + nc + "]/F").c_str());
        otree->Branch("jcsDEtaRaw2", o_jcsDEtaRaw[2], ("jcsDEtaRaw2[" + nc + "]/F").c_str());
        otree->Branch("jcsDPhi2",    o_jcsDPhi[2],    ("jcsDPhi2["    + nc + "]/F").c_str());
        otree->Branch("jcsM2",       o_jcsM[2],       ("jcsM2["       + nc + "]/F").c_str());
        otree->Branch("jcsW2",       o_jcsW[2],       ("jcsW2["       + nc + "]/F").c_str());
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

    long alleve = 0;
    double sumW = 0., sumW2 = 0.;
    long start = args.skip;
    long end   = (args.totEve < 0) ? nentries : std::min((long)nentries, start + args.totEve);

    // ── event loop ───────────────────────────────────────────────────────────
    for (long iev = start; iev < end; ++iev) {
        ttree->GetEntry(iev);
        if (iev % 10000 == 0) std::cout << "event " << iev << "\n";

        // reset output (zero-pad so unused jet/constituent slots are 0)
        o_nJets = 0;
        memset(o_jetPt,   0, sizeof(o_jetPt));
        memset(o_jetEta,  0, sizeof(o_jetEta));
        memset(o_jetPhi,  0, sizeof(o_jetPhi));
        memset(o_jetM,    0, sizeof(o_jetM));
        memset(o_jetSPVA, 0, sizeof(o_jetSPVA));
        memset(o_jetPVA,  0, sizeof(o_jetPVA));
        memset(o_jetPVM,  0, sizeof(o_jetPVM));
        memset(o_jetNC,   0, sizeof(o_jetNC));
        if (SAVE_JCS) {
            memset(o_jcsPt,   0, sizeof(o_jcsPt));
            memset(o_jcsDEta,    0, sizeof(o_jcsDEta));
            memset(o_jcsDEtaRaw, 0, sizeof(o_jcsDEtaRaw));
            memset(o_jcsDPhi,    0, sizeof(o_jcsDPhi));
            memset(o_jcsM,    0, sizeof(o_jcsM));
            memset(o_jcsW,    0, sizeof(o_jcsW));
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
            if (pid == 25 || std::abs(pid) == 12 || std::abs(pid) == 14 ||
                             std::abs(pid) == 16 || std::abs(pid) == 15) continue;
            fastjet::PseudoJet pj(px, py, pz, E);
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
        }

        // per-event inclusive histograms (no dijet selection)
        hNJets->Fill((float)jets.size(), o_kWeight);
        if (!jets.empty()) hJetPt1->Fill(jets[0].pt(), o_kWeight);

        // dijet selection: mjj>MJJ_MIN, opposite-eta leading pair, |dYjj|>DYJJ_MIN
        bool passNJets = (int)jets.size() >= 2;
        double dYjj    = passNJets ? (jets[0].rapidity() - jets[1].rapidity()) : -99.;
        if (passNJets &&
            sumP4(jets, 2).M() > MJJ_MIN &&
            jets[0].eta() * jets[1].eta() < 0. &&
            std::abs(dYjj) > DYJJ_MIN) {

            TLorentzVector jj = sumP4(jets, 2);
            o_mjj    = (Float_t)jj.M();
            o_dYjj   = (Float_t)dYjj;
            o_dPhijj = (Float_t)deltaPhi(jets[0].phi(), jets[1].phi());
            o_ptjj   = (Float_t)jj.Pt();
            bool swap = (iev % 2 == 0);  // odd iev → leading at [0], even → leading at [1]
            o_leadJetIndex = swap ? 1 : 0;

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

            // fill output tree once per event
            int nmax = std::min((int)jets.size(), NJETMAX);
            o_nJets = nmax;
            for (int i = 0; i < nmax; ++i) {
                // randomization confined to jets 0/1; jet2 stays pt-ordered 3rd jet
                int si = (i < 2) ? (swap ? 1 - i : i) : i;
                o_jetPt[si]   = (Float_t)std::round(jets[i].pt()  * 100) / 100.f;
                o_jetEta[si]  = (Float_t)std::round(jets[i].eta() * 1000) / 1000.f;
                o_jetPhi[si]  = (Float_t)std::round(jets[i].phi() * 1000) / 1000.f;
                o_jetM[si]    = (Float_t)std::round(jets[i].m()   * 100) / 100.f;
                o_jetSPVA[si] = (Float_t)pvs[i].spva;
                o_jetPVA[si]  = (Float_t)pvs[i].pva;
                o_jetPVM[si]  = (Float_t)pvs[i].pvm;
                o_jetNC[si]   = (Int_t)jets[i].constituents().size();
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
                    }
                }
            }
            otree->Fill();
            ++alleve;
        }
    }

    std::cout << "hSumW integral   " << hSumW->Integral() << "\n";
    std::cout << "events in tree   " << alleve            << "\n";
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
