from __future__ import print_function
import os, re, sys
from array import array
import itertools
import numpy as np

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


# akcolor
import fastjet

# parse external parameters
import argparse

# overide default argparse behavior of the error method
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = MyParser()
parser = argparse.ArgumentParser(description='making friends')
parser.add_argument('files', nargs='+', help='needs at minimum 1 file')
parser.add_argument('--output', default = '', help = 'name of the plot directory to be created')
parser.add_argument('--ttree', default = '', help = 'name of the tree, if is inside a TDir, Dirname/TreeName, otherwise will attemp to fetch from ListOfKnownTTrees')
parser.add_argument('--xs', default = 1.0, help = 'user provided cross section for normalization', type=float)
parser.add_argument('--genWeight', default = "genWeight", help = 'MC weight from event generator', type=str)
parser.add_argument('--goFast', default = 1.0, help = 'if set, will only process fraction of the entries', type=float)
parser.add_argument('--debug', default = 0, help = 'user provided cross section for normalization', type=int)
parser.add_argument('--skip', default = 0, help = 'skip first events', type=float)
parser.add_argument('--totEve', default = -1, help = 'totEve to process', type=float)
parser.add_argument('--jobId', default = -1, help = 'jobId to append in the output', type=float)
args = parser.parse_args()

# make a simple counter
class counter:
    """counting events"""
    pass

# list of known TTrees
listOfknownTrees = ['Events','events', 'ntuple/tree', 'tree', 'Data']

# define maximum numparticles, as a trivial funciton
maxNpartcl = lambda : 10000

ptMIN = 0.0  # pt threshold for jet constituents used in pull vector calculation

# get momentum of ith particle in format  E, px, py, pz, id
def ith_momentum(obj, ith):
    arraySize = int(maxNpartcl())
    return (obj[0*arraySize+ith],obj[1*arraySize+ith],obj[2*arraySize+ith],obj[3*arraySize+ith],int(obj[4*arraySize+ith]))


# get momenta
def get_momenta(obj, numprtcl):
    mom = []
    for i in range(0,numprtcl):
        mom.append(ith_momentum(obj, i))

    return np.array(mom, dtype=[('E', 'f8'), ('px', 'f8'), ('py', 'f8'), ('pz', 'f8'), ('id', 'int')])

# convert to fastjet, but only if the pt is > minptc
def convert_tofj(momin, minptc, maxrapc, excludeNeutrinos = True):
    arrayout = []
    for mm in range(len(momin)):
        #print(momin[mm][1], momin[mm][2], momin[mm][3], momin[mm][0])
        fj = fastjet.PseudoJet(momin[mm][1], momin[mm][2], momin[mm][3], momin[mm][0])
        if excludeNeutrinos:
            if fj.perp() > minptc and abs(fj.eta()) < maxrapc:
                arrayout.append(fj)
            fj.set_user_index(int(momin[mm][4]))
    return arrayout

######## This is not used anywhere but is left for documentation
###def get_clusters(events, filename, jetalgo, jetR, maxevents=100000):
###    if maxevents > len(events):
###        maxevents = len(events)
###    print('Analyzing', maxevents, 'events from', filename)
###    # put the Higgs momenta into an array:
###    higgs = []
###    # return the cluster:
###    clusters_jets = []
###    # return the unclustered objects:
###    unclustered = []
###    # jet algorithm
###    jetdef = fastjet.JetDefinition(jetalgo, jetR)
###    # loop over events and analyze:
###    for yy in tqdm(range(0,maxevents)):
###        # put the momenta for clustering into array:
###        momtocluster = []
###        # and the rest into another array;
###        momNOcluster = []
###        # all the momenta from this event
###        momenta = events[yy]
###        #print(momenta)
###        for mm in range(0,len(momenta)):
###            if momenta[mm][4] == 25: # find a Higgs boson
###             higgs.append(momenta[mm])
###            if momenta[mm][4] != 25 and abs(momenta[mm][4]) != 12 and abs(momenta[mm][4]) != 14 and abs(momenta[mm][4]) != 16 and abs(momenta[mm][4]) != 11 and abs(momenta[mm][4]) != 13:
###                momtocluster.append(momenta[mm])
###            else:
###                momNOcluster.append(momenta[mm])
###        momfj = convert_tofj(momtocluster, jcPtMin,jcEtaMax)
###        momfj_unclustered = convert_tofj(momNOcluster, jcPtMin,jcEtaMax)
###
###        cluster = fastjet.ClusterSequence(momfj, jetdef)
###        clusters_jets.append(cluster)
###        unclustered.append(momfj_unclustered)
###    return clusters_jets, unclustered

# from https://gitlab.cern.ch/cms-sw/cmssw/blob/e303d9f2c3d4f25397db5feb7ad59d2f20c842f2/PhysicsTools/HeppyCore/python/utils/deltar.py
def deltaPhi( p1, p2):
    '''Computes delta phi, handling periodic limit conditions.'''
    res = p1 - p2
    while res > np.pi:
        res -= 2*np.pi
    while res < -np.pi:
        res += 2*np.pi
    return res

# equipe with the pull vector the selected pseudojets
def fillPV(j, jetconst):
    theta = -99
    r = -99
    spva = -99 # signed pull vector angle
    pullV = np.array([0., 0.])
    for jc in jetconst:
        if jc.pt() < 1.e-3: continue
        dY    = jc.rapidity() - j.rapidity()
        dPhi  =  deltaPhi(jc.phi(), j.phi())
        r     = (dY**2 + dPhi**2)**0.5
        pullV +=  (jc.pt()/j.pt())*r*np.array([dY, dPhi])
    r = (pullV[0]**2 + pullV[1]**2)**0.5
    if r>0:
        theta  = np.arctan2(pullV[1]/r, pullV[0]/r)
        spva   = theta
        if j.rapidity() < 0: spva = np.arctan2(pullV[1]/r, -pullV[0]/r)
    j.pvm  = r
    j.pva  = theta
    j.spva = spva


# translate a LorentzVector to TLorentzVector
toTLV = lambda x: ROOT.TLorentzVector(x.px(), x.py(), x.pz(), x.e())  


# vector sum of four vectors
def sumP4(fvecs):
    #megaJ = ROOT.Math.LorentzVector(ROOT.Math.PxPyPzE4D('double'))(0,0,0,0)
    megaJ = ROOT.TLorentzVector(0,0,0,0)
    for jet in fvecs:
        megaJ += jet.p4
    return megaJ




#  main starts here
if __name__ == "__main__":
    tfiles = [ROOT.TFile.Open(f) for f in args.files]

    # read the TTrees from files check if ttree is from list of known, or force to use external
    if args.ttree!='':listOfknownTrees = [args.ttree]

    ttrees = [tfile.Get(ttree) for tfile in tfiles  for ttree in listOfknownTrees if tfile != None and tfile.Get(ttree)!=None]
    for ttree in ttrees: ttree.tfile = ttree.GetCurrentFile()

    if len(ttrees) != len(tfiles):
        print('not all tfiles have been found with a valid ttree, exiting')
        os._exit(0)

    minimalPrint = True


    # check if all ttrees have the same name and set it to equal to args.ttree, print warning if not
    ttreeNames = [ttree.GetName() for ttree in ttrees]
    if len(set(ttreeNames))!=1 and len(ttreeNames)>0: print('Warning: not all TTrees have the same name %s'%ttreeNames)
    else: args.ttree = ttreeNames[0]


    # RDataFrame of all files
    print("Getting sum of weights for all input files")
    sumWtot = 0
    Ntot = 0
    for ii, ttree in enumerate(ttrees):
        df = ROOT.RDataFrame(ttree)
        sumW = 0
        if args.genWeight != 'genWeight': ttree.SetAlias('genWeight', args.genWeight)
        try:
            sumW = df.Sum(args.genWeight).GetValue()
        except TypeError:
            print(args.genWeight, " is not part of the given ttrees, use --genWeight branchName to fix this or assume no such branch and count each entry as 1 event")
            sumW = ttree.GetEntries()
        print('%s with %d entries and sumW = %2.1f'%(tfiles[ii].GetName(), ttree.GetEntries(), sumW))
        Ntot += ttree.GetEntries()
        sumWtot += sumW

    print("sumWtot %d   Ntot %d"%(sumWtot, Ntot))

    # one output tree for all inputs, using the name of the first by default
    outname = tfiles[0].GetName()[0:-4].split('/')[-1]+'.slimfriend'
    if args.jobId!=-1: outname += '_'+str(int(args.jobId))
    outname += '.root'

    # change the outname on user's request
    if args.output != '': outname = args.output
    ofile = ROOT.TFile(outname,"RECREATE")

    #otree = ROOT.TTree(ttree.GetName(), ttree.GetName()) # same name as original
    otree = ROOT.TTree("events", "events") # use events for the name of the output tree

    hSumW  = ROOT.TH1D("hSumW","hSumW",1, 0, 1)
    hSumW.Sumw2()

    _rhoEdges = np.linspace(0, 0.5, 31)
    _phiEdges = np.linspace(0, 2*np.pi, 65)
    _nArc     = 12
    jetShape  = ROOT.TH2Poly("jetShape",  "jetShape",  -0.5, 0.5, -0.5, 0.5)
    jetShape.Sumw2()
    for _iPhi in range(64):
        _p1, _p2 = _phiEdges[_iPhi], _phiEdges[_iPhi+1]
        for _iRho in range(30):
            _r1, _r2 = _rhoEdges[_iRho], _rhoEdges[_iRho+1]
            _arc = np.linspace(_p1, _p2, _nArc)
            _xv  = np.concatenate([_r2*np.cos(_arc), _r1*np.cos(_arc[::-1])])
            _yv  = np.concatenate([_r2*np.sin(_arc), _r1*np.sin(_arc[::-1])])
            _xv_arr = array('d', _xv.tolist())
            _yv_arr = array('d', _yv.tolist())
            jetShape.AddBin(len(_xv), _xv_arr, _yv_arr)

    # ttree variables
    nJetsMax = 2
    tvars = []
    t_nJets       = array('i', [0]); tvars += [t_nJets]
    t_jetPt       = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPt]
    t_jetEta      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetEta]
    t_jetPhi      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPhi]
    t_jetM        = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetM] 
    t_jetSPVA     = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetSPVA]     
    t_jetPVA      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPVA]     
    t_jetPVM      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPVM]
    t_jetNC       = array('i', [0   for i in range(nJetsMax)]); tvars += [t_jetNC]
    t_kWeight     = array('f', [1.]) ; tvars += [t_kWeight]
    
   #kWeight
 
    def reset():
        global tvars
        for var in tvars:
            typecode = var.typecode
            for i in range(len(var)):
                if typecode == 'f': var[i] = -99999.9
                if typecode == 'i': var[i] = int(-99999)
                if typecode == 'B': var[i] = False
    

    otree.Branch("nJets",      t_nJets,      "nJets/I")
    otree.Branch("jetPt",      t_jetPt,      "jetPt[nJets]/F")
    otree.Branch("jetEta",     t_jetEta,     "jetEta[nJets]/F")
    otree.Branch("jetPhi",     t_jetPhi,     "jetPhi[nJets]/F")
    otree.Branch("jetM",       t_jetM,       "jetM[nJets]/F")
    otree.Branch("jetSPVA",    t_jetSPVA,    "jetSPVA[nJets]/F")
    otree.Branch("jetPVA",     t_jetPVA,     "jetPVA[nJets]/F")
    otree.Branch("jetPVM",     t_jetPVM,     "jetPVM[nJets]/F")
    otree.Branch("jetNC",      t_jetNC,      "jetNC[nJets]/I")
    otree.Branch("kWeight",    t_kWeight,    "kWeight/F")
    
    def fillJets(collection):
        nmax = min(len(collection), nJetsMax)
        t_nJets[0] = nmax 
        for iobj, obj in enumerate(collection[0:nmax]):
            t_jetPt[iobj]      = round(obj.pt() ,  2)
            t_jetEta[iobj]     = round(obj.eta(),  3)
            t_jetPhi[iobj]     = round(obj.phi(),  3)
            t_jetM[iobj]       = round(obj.m(), 2)
            t_jetSPVA[iobj]    = round(obj.spva,    7)
            t_jetPVA[iobj]     = round(obj.pva,     5)
            t_jetPVM[iobj]     = round(obj.pvm,     7)
            t_jetNC[iobj]      = len(obj.constituents())

    def reset():
        global tvars
        for var in tvars:
            typecode = var.typecode
            for i in range(len(var)):
                if typecode == 'f': var[i] = -99.9
                if typecode == 'i': var[i] = int(-99)
                if typecode == 'B': var[i] = False

    count = counter()
    count.alleve   = 0
    count.sumW     = 0
    count.sumW2    = 0


    # jet algorithm
    jetR    = 0.4
    jetalgo = fastjet.antikt_algorithm
    jetdef  = fastjet.JetDefinition(jetalgo, jetR)

    # loop over the trees
    for ii, ttree in enumerate(ttrees):
        print("opening %s" % tfiles[ii].GetName())

        # Get the total number of entries in the tree
        nentries = ttree.GetEntries()

        # Set the starting point for the event loop based on args.skip
        start_entry = int(args.skip)
        if start_entry >= nentries:
            print("Skipping entire tree %s, as skip value is too large." % tfiles[ii].GetName())
            continue

        # Event loop starting directly from 'start_entry'
        for iev in range(start_entry, nentries):
            # Retrieve the specific entry
            ttree.GetEntry(iev)
            event = ttree  
            if args.totEve != -1 and (iev >= args.skip + args.totEve): break

            if args.goFast and iev >=  args.goFast*ttree.GetEntries()  : break
            if minimalPrint and iev%10000 == 0: print('event %d of %s'%(iev, tfiles[ii].GetName()))

            # magic starts here
            reset()

            # put the Higgs momenta into an array:
            higgs = []

            # return the cluster:
            clusters_jets = []

            # return the unclustered objects:
            unclustered = []

            # put the momenta for clustering into array:
            momtocluster = []

            # and the rest into another array;
            momNOcluster = []

            t_kWeight[0]      = event.genWeight*args.xs/sumWtot
            count.sumW        += t_kWeight[0]
            count.sumW2       += t_kWeight[0]*t_kWeight[0]
            hSumW.Fill(0, t_kWeight[0])


            numparticles  = getattr(event, "numparticles")
            objects = getattr(event,"objects")
            momenta = get_momenta(objects, numparticles)

            if args.debug >= 3: print('numparticles %d'%numparticles)
            if args.debug >= 3: print(momenta.size)

            ### particle loop: mm is counting the particles of each event
            for mm in range(0,len(momenta)):

                (E, px, py, pz, pid)  = momenta[mm]
                if args.debug == 4:
                    if mm == 0: print('########## event %d ###########'%iev)
                    print('(E %2.3f, px %2.3f, py %2.3f, pz %2.3f, pid %d)'%(E, px, py, pz, pid))

                if pid == 25: # find a Higgs boson
                    HiggsP4 = ROOT.TLorentzVector(px, py, pz, E)
                    higgs.append(HiggsP4)
                    if args.debug == 4:print('Higgs boson found (E %2.3f, px %2.3f, py %2.3f, pz %2.3f, pid %d)'%(E, px, py, pz, pid) , higgs)

                if pid != 25 and abs(pid) != 12 and abs(pid) != 14 and abs(pid) != 16 and abs(pid) != 15: # exclude stable Higgs, neutrinos and taus
                    momtocluster.append(momenta[mm])
                else:
                    momNOcluster.append(momenta[mm])

            # Phase-II offline cuts
            # jcPtMin, jcEtaMax      = 2.0, 4.0
            # jetPtMin, jetEtaMax    = 30.0, 3.0

            # Run2/3 offline cuts
            # jcPtMin, jcEtaMax      = 0.5, 2.4 
            # jetPtMin, jetEtaMax    = 30.0, 2.0   

            # Phase-II no cuts
            jcPtMin, jcEtaMax      = 0.0, 10.0
            jetPtMin, jetEtaMax    = 30.0, 3.0

            # Phase-II L1 cuts
            # jcPtMin, jcEtaMax      = 2.0, 2.4
            # jetPtMin, jetEtaMax    = 30.0, 2.0

            momfj = convert_tofj(momtocluster, jcPtMin,jcEtaMax)
            momfj_unclustered = convert_tofj(momNOcluster, jcPtMin,jcEtaMax)

            cluster = fastjet.ClusterSequence(momfj, jetdef)
            events_jets = fastjet.sorted_by_pt(cluster.inclusive_jets(jetPtMin))
            jets = [jet for jet in events_jets if abs(jet.eta()) < jetEtaMax]
            
            if args.debug == 2:print('len events_jets = ', len(events_jets))
            if args.debug == 2:print('len jets = ', len(jets))

            # devise the jets with pull vector and TLorentzVector accessed through p4 member
            for jet in jets:
                jcs = [jc for jc in jet.constituents() if jc.pt()>ptMIN and abs(jc.eta())<jcEtaMax]
                fillPV(jet, jcs)
                jet.p4 = toTLV(jet)

            # alias to jets for code combatibility
            genjs = jets

            # fill jetShape for mjj>400 events with forward leading jet
            if len(jets) >= 2 and sumP4(jets[0:2]).M() > 400 and jets[0].eta()*jets[1].eta() < 0:
                for jet in jets[0:2]:
                    for jc in jet.constituents():
                        dY   = jc.rapidity() - jet.rapidity()
                        dPhi = deltaPhi(jc.phi(), jet.phi())
                        r2   = dY**2 + dPhi**2
                        z    = jc.pt()/jet.pt()
                        w    = z * np.sqrt(r2)
                        dYsigned = dY
                        if jet.rapidity()<0: dYsigned = -dY
                        jetShape.Fill(dYsigned, dPhi, w)

                    # fill tree
                    fillJets  (jets)
                    otree.Fill()
                    count.alleve += 1

    print('hSumW.Integral() %2.1f'%hSumW.Integral())
    print('number of events processed %d'%count.alleve)
    print('number of entries in the TTree %d'%ttree.GetEntries())
    print('number of sumW %3.10f'%count.sumW)
    print('number of sumW2 %3.10f'%count.sumW2)

    # creating output
    ofile.cd()
    otree.Write()
    hSumW.Write()
    jetShape.Write()
    ofile.Write()
    ofile.Close()
