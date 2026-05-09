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
#include "TH2Poly.h"
#include "TLorentzVector.h"
#include "ROOT/RDataFrame.hxx"

#include "fastjet/ClusterSequence.hh"

static const int MAXNPART = 10000;
static const int NJETMAX  = 2;

// ── helpers ──────────────────────────────────────────────────────────────────

inline double deltaPhi(double p1, double p2) {
    double r = p1 - p2;
    while (r >  M_PI) r -= 2*M_PI;
    while (r < -M_PI) r += 2*M_PI;
    return r;
}

struct JetPV { double pvm = -99., pva = -99., spva = -99.; };

JetPV fillPV(const fastjet::PseudoJet& jet,
             const std::vector<fastjet::PseudoJet>& jcs) {
    double pv0 = 0., pv1 = 0.;
    for (const auto& jc : jcs) {
        if (jc.pt() < 1e-3) continue;
        double dY   = jc.rapidity() - jet.rapidity();
        double dPhi = deltaPhi(jc.phi(), jet.phi());
        double r    = std::sqrt(dY*dY + dPhi*dPhi);
        double w    = (jc.pt() / jet.pt()) * r;
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
    std::string inputfile, genWeight = "genWeight", output = "";
    double xs = 1.0;
    long   totEve = -1, skip = 0, jobId = -1;
    int    debug  = 0;
};

Args parseArgs(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " file.root [--genWeight W] [--xs X] [--totEve N]"
                     " [--skip N] [--jobId N] [--output out.root] [--debug N]\n";
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
        else if (f == "--debug"    ) a.debug     = std::stoi(argv[++i]);
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

    // leading-jet |SPVA| (20 bins, 0–π), filled after mjj>400 + opposite-eta selection
    TH1F* hLeadJetSPVA = new TH1F("hLeadJetSPVA","hLeadJetSPVA", 20, 0., M_PI);
    hLeadJetSPVA->Sumw2();

    // TH2Poly jetShape — 64 phi × 30 rho annular sectors
    TH2Poly* jetShape = new TH2Poly("jetShape","jetShape",-0.5,0.5,-0.5,0.5);
    jetShape->Sumw2();
    {
        const int nArc = 12, nPhi = 64, nRho = 30;
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

    // output branches
    Int_t   o_nJets;
    Float_t o_jetPt[NJETMAX], o_jetEta[NJETMAX], o_jetPhi[NJETMAX], o_jetM[NJETMAX];
    Float_t o_jetSPVA[NJETMAX], o_jetPVA[NJETMAX], o_jetPVM[NJETMAX];
    Int_t   o_jetNC[NJETMAX];
    Float_t o_kWeight;

    otree->Branch("nJets",   &o_nJets,  "nJets/I");
    otree->Branch("jetPt",   o_jetPt,   "jetPt[nJets]/F");
    otree->Branch("jetEta",  o_jetEta,  "jetEta[nJets]/F");
    otree->Branch("jetPhi",  o_jetPhi,  "jetPhi[nJets]/F");
    otree->Branch("jetM",    o_jetM,    "jetM[nJets]/F");
    otree->Branch("jetSPVA", o_jetSPVA, "jetSPVA[nJets]/F");
    otree->Branch("jetPVA",  o_jetPVA,  "jetPVA[nJets]/F");
    otree->Branch("jetPVM",  o_jetPVM,  "jetPVM[nJets]/F");
    otree->Branch("jetNC",   o_jetNC,   "jetNC[nJets]/I");
    otree->Branch("kWeight", &o_kWeight,"kWeight/F");

    // input branches
    Int_t    numparticles;
    Double_t objects_arr[8][MAXNPART];
    Double_t evweight;

    ttree->SetBranchAddress("numparticles",       &numparticles);
    ttree->SetBranchAddress("objects",             objects_arr);
    ttree->SetBranchAddress(args.genWeight.c_str(), &evweight);

    // fastjet setup
    fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, 0.4);
    const double jcPtMin   = 0.0, jcEtaMax  = 10.0;
    const double jetPtMin  = 30.0, jetEtaMax =  3.0;

    long alleve = 0;
    double sumW = 0., sumW2 = 0.;
    long start = args.skip;
    long end   = (args.totEve < 0) ? nentries : std::min((long)nentries, start + args.totEve);

    // ── event loop ───────────────────────────────────────────────────────────
    for (long iev = start; iev < end; ++iev) {
        ttree->GetEntry(iev);
        if (iev % 10000 == 0) std::cout << "event " << iev << "\n";

        // reset output
        o_nJets = 0;
        for (int i = 0; i < NJETMAX; ++i) {
            o_jetPt[i] = o_jetEta[i] = o_jetPhi[i] = o_jetM[i] = -99.9f;
            o_jetSPVA[i] = o_jetPVA[i] = o_jetPVM[i] = -99.9f;
            o_jetNC[i] = -99;
        }
        o_kWeight = (Float_t)(evweight * args.xs / sumWtot);
        sumW  += o_kWeight;
        sumW2 += o_kWeight * o_kWeight;
        hSumW->Fill(0., o_kWeight);

        // particle loop → build momtocluster
        std::vector<fastjet::PseudoJet> momfj;
        momfj.reserve(numparticles);
        for (int mm = 0; mm < numparticles; ++mm) {
            double E  = objects_arr[0][mm];
            double px = objects_arr[1][mm];
            double py = objects_arr[2][mm];
            double pz = objects_arr[3][mm];
            int    pid = (int)objects_arr[4][mm];
            if (pid == 25 || std::abs(pid) == 12 || std::abs(pid) == 14 ||
                             std::abs(pid) == 16 || std::abs(pid) == 15) continue;
            fastjet::PseudoJet pj(px, py, pz, E);
            if (pj.perp() > jcPtMin && std::abs(pj.eta()) < jcEtaMax)
                momfj.push_back(pj);
        }

        // cluster
        fastjet::ClusterSequence cs(momfj, jetdef);
        auto all_jets = fastjet::sorted_by_pt(cs.inclusive_jets(jetPtMin));
        std::vector<fastjet::PseudoJet> jets;
        for (const auto& j : all_jets)
            if (std::abs(j.eta()) < jetEtaMax) jets.push_back(j);

        // pull vectors for all jets
        std::vector<JetPV> pvs(jets.size());
        for (int ij = 0; ij < (int)jets.size(); ++ij) {
            auto jcs = jets[ij].constituents();
            pvs[ij] = fillPV(jets[ij], jcs);
        }

        // dijet selection: mjj>400, opposite-eta leading pair
        if ((int)jets.size() >= 2 &&
            sumP4(jets, 2).M() > 400. &&
            jets[0].eta() * jets[1].eta() < 0.) {

            // |SPVA| for leading and subleading jet — two fills per event
            for (int ij = 0; ij < 2; ++ij)
                if (pvs[ij].spva > -90.)
                    hLeadJetSPVA->Fill(std::abs(pvs[ij].spva), o_kWeight);

            // fill jetShape from both leading jets (sign-flip dY for backward jet)
            for (int ij = 0; ij < 2; ++ij) {
                for (const auto& jc : jets[ij].constituents()) {
                    double dY   = jc.rapidity() - jets[ij].rapidity();
                    double dPhi = deltaPhi(jc.phi(), jets[ij].phi());
                    double z    = jc.pt() / jets[ij].pt();
                    double w    = z * std::sqrt(dY*dY + dPhi*dPhi);
                    double dYs  = (jets[ij].rapidity() < 0) ? -dY : dY;
                    jetShape->Fill(dYs, dPhi, w);
                }

                // fill output tree (once per jet, matching Python behaviour)
                int nmax = std::min((int)jets.size(), NJETMAX);
                o_nJets = nmax;
                for (int i = 0; i < nmax; ++i) {
                    o_jetPt[i]   = (Float_t)std::round(jets[i].pt()  * 100) / 100.f;
                    o_jetEta[i]  = (Float_t)std::round(jets[i].eta() * 1000) / 1000.f;
                    o_jetPhi[i]  = (Float_t)std::round(jets[i].phi() * 1000) / 1000.f;
                    o_jetM[i]    = (Float_t)std::round(jets[i].m()   * 100) / 100.f;
                    o_jetSPVA[i] = (Float_t)pvs[i].spva;
                    o_jetPVA[i]  = (Float_t)pvs[i].pva;
                    o_jetPVM[i]  = (Float_t)pvs[i].pvm;
                    o_jetNC[i]   = (Int_t)jets[i].constituents().size();
                }
                otree->Fill();
                ++alleve;
            }
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
    hLeadJetSPVA->Write();
    ofile->Write();
    ofile->Close();
    ifile->Close();
    return 0;
}
