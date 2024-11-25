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

# https://root.cern.ch/root/html534/guides/users-guide/InputOutput.html#the-logical-root-file-tfile-and-tkey
def printListOfkeys(tfile):tfile.GetListOfKeys().Print()

# print list of possible variables of the ROOT file
def printListOfLeaves(myttree, filename = ''):
    if filename != '':
        fp = open(filename, 'w')
        print('opening %s to write the list of leaves for %s'%(filename, myttree.GetName()))
        for leave in myttree.GetListOfLeaves():
            print(leave)
            fp.write(str(leave)+'\n')
        print('closing %s'%filename)
        fp.close()
    else:
        for leave in myttree.GetListOfLeaves():
            print(leave)

# make a simple counter
class counter:
    """counting events"""
    pass

# list of known TTrees
listOfknownTrees = ['Events','events', 'ntuple/tree', 'tree', 'Data']

# define maximum numparticles, as a trivial funciton
maxNpartcl = lambda : 10000

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
    pullV = np.array([0., 0.])
    for jc in jetconst:
        if jc.pt() < 1.e-3: continue
        dY    = jc.rapidity() - j.rapidity()
        dPhi  =  deltaPhi(jc.phi(), j.phi())
        r     = (dY**2 + dPhi**2)**0.5
        jvec  = np.array([j.px(), j.py(), j.pz()])
        jcvec = np.array([jc.px(), jc.py(), jc.pz()])
        jcL   = ((jvec@jcvec)/(jvec@jvec))*jvec
        jcT   = jcvec - jcL
        ptRel = (jcT@jcT)**0.5
        pullV +=  ptRel*np.array([dY, dPhi])
        # print('PtRel = %2.2f r = %2.2f |jcT| = %2.2f'%( PtRel(jc, j), r, (jcT@jcT)**0.5)) # works!
    j.pv1 = pullV[0]
    j.pv2 = pullV[1]
    r = (j.pv1**2 + j.pv2**2)**0.5
    if r>0:
        theta  = np.arctan2(j.pv2/r, j.pv1/r)
    j.pvm = r
    j.pva = theta

def RPA(j1, j2):
    '''cos21 needs j1 and j2 in that order, do not invert'''
    theta21 = -99
    cos21 = -99
    p  = np.array([[j1.pv1], [j1.pv2]])
    v2 = np.array([[j2.rapidity()] , [j2.phi()] ])
    v1 = np.array([[j1.rapidity()] , [j1.phi()] ])
    r = v2 - v1
    mag_p = (p.T).dot(p)[0][0]**0.5
    mag_r = (r.T).dot(r)[0][0]**0.5
    if mag_p*mag_r>0:
        cos21 = (r.T).dot(p)[0][0]/(mag_p*mag_r)
        theta21 = np.arccos(cos21)*(180/np.pi)
    if args.debug == 3:
        print('p = ', p)
        print('v2 = ', v2)
        print('v1 = ', v1)
        print('r = ', r)
        print('mag_r = ', mag_r)
        print('mag_rp = ', mag_p)
        print('(r.T).dot(p)[0][0] = ', (r.T).dot(p)[0][0])
    return theta21

# equipe with the pull vector the selected pseudojets
def fillPV3(j, jetconst):
    pullV3 = 0
    for jc in jetconst:
        if jc.pt() < 1.e-3: continue
        dY    = jc.rapidity() - j.rapidity()
        pullV3 +=  (jc.pt()/j.pt())*dY
    j.pv3 = pullV3


# translate a LorentzVector to TLorentzVector
toTLV = lambda x: ROOT.TLorentzVector(x.px(), x.py(), x.pz(), x.e())  


# vector sum of four vectors
def sumP4(fvecs):
    #megaJ = ROOT.Math.LorentzVector(ROOT.Math.PxPyPzE4D('double'))(0,0,0,0)
    megaJ = ROOT.TLorentzVector(0,0,0,0)
    for jet in fvecs:
        megaJ += jet.p4
    return megaJ

# calculate sum of pt in a cone of DR = 0.3 arround a particle
def Pt03(p, particles):
    myP4 = toTLV(p)
    sumPt = 0
    for aParticle in particles:
        aParticleP4 = toTLV(aParticle)
        DR = myP4.DeltaR(aParticleP4)
        if 0 < DR < 0.3:
            sumPt +=  aParticleP4.Pt()

    return sumPt

# calculate the pt of a particle relative to its closest jet
def PtRel(p, jet):
    myP4 = toTLV(p)
    myJet = toTLV(jet)
    ptrel = myP4.Perp(myJet.Vect())
    return ptrel    
            


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
    outname = tfiles[0].GetName()[0:-4].split('/')[-1]+'friend'
    if args.jobId!=-1: outname += '_'+str(int(args.jobId))
    outname += '.root'

    # change the outname on user's request
    if args.output != '': outname = args.output
    ofile = ROOT.TFile(outname,"RECREATE")

    #otree = ROOT.TTree(ttree.GetName(), ttree.GetName()) # same name as original
    otree = ROOT.TTree("events", "events") # use events for the name of the output tree

    hSumW  = ROOT.TH1D("hSumW","hSumW",1, 0, 1)
    hSumW.Sumw2()

    # ttree variables
    nJetsMax = 5
    tvars = []
    t_nJets       = array('i', [0]); tvars += [t_nJets]
    t_jetPt       = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPt]
    t_jetEta      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetEta]
    t_jetPhi      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPhi]
    t_jetM        = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetM] 
    t_jetPV1      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPV1]
    t_jetPV2      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPV2]     
    t_jetPV3      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPV3]     
    t_jetPVA      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPVA]     
    t_jetPVM      = array('f', [0.  for i in range(nJetsMax)]); tvars += [t_jetPVM]     
    t_jetFlag     = array('i', [0   for i in range(nJetsMax)]); tvars += [t_jetFlag]    
    t_jetMuIndex  = array('i', [-1 for i in range(nJetsMax)]) ; tvars += [t_jetMuIndex] 
    t_jetBtag     = array('B', [False  for i in range(nJetsMax)]); tvars += [t_jetBtag]
    t_mjj         = array('f', [0.]) ; tvars += [t_mjj]      
    t_ptjj        = array('f', [0.]) ; tvars += [t_ptjj]     
    t_dYjj        = array('f', [0.]) ; tvars += [t_dYjj]    
    t_dPhijj      = array('f', [0.]) ; tvars += [t_dPhijj]   
    t_mbb         = array('f', [0.]) ; tvars += [t_mbb]      
    t_ptbb        = array('f', [0.]) ; tvars += [t_ptbb]     
    t_dYbb        = array('f', [0.]) ; tvars += [t_dYbb]    
    t_dPhibb      = array('f', [0.]) ; tvars += [t_dPhibb]   
    t_c21         = array('f', [-3.]); tvars += [t_c21]      
    t_c12         = array('f', [-3.]); tvars += [t_c12]      
    t_b21         = array('f', [-3.]); tvars += [t_b21]      
    t_b12         = array('f', [-3.]); tvars += [t_b12]      
    t_met         = array('f', [0.]) ; tvars += [t_met]      
    t_metphi      = array('f', [0.]) ; tvars += [t_metphi]   
    t_higgsPt     = array('f', [0.]) ; tvars += [t_higgsPt]      
    t_higgsM      = array('f', [0.]) ; tvars += [t_higgsM]      
    t_higgsY      = array('f', [0.]) ; tvars += [t_higgsY]      
    t_higgsPhi    = array('f', [0.]) ; tvars += [t_higgsPhi]      
    t_higgsEta    = array('f', [0.]) ; tvars += [t_higgsEta]      
    t_weight      = array('f', [1.]) ; tvars += [t_weight]   
    t_kWeight     = array('f', [1.]) ; tvars += [t_kWeight]
    t_kFile       = array('i', [0]); tvars += [t_kFile]
    t_kBadFlag    = array('i', [0]); tvars += [t_kBadFlag]
    
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
    otree.Branch("jetPV1",     t_jetPV1,     "jetPV1[nJets]/F")
    otree.Branch("jetPV2",     t_jetPV2,     "jetPV2[nJets]/F")
    otree.Branch("jetPV3",     t_jetPV3,     "jetPV3[nJets]/F")
    otree.Branch("jetPVA",     t_jetPVA,     "jetPVA[nJets]/F")
    otree.Branch("jetPVM",     t_jetPVM,     "jetPVM[nJets]/F")
    otree.Branch("jetFlag",    t_jetFlag,    "jetFlag[nJets]/I")
    otree.Branch("jetBtag",    t_jetBtag,    "jetBtag[nJets]/O")
    otree.Branch("jetMuIndex", t_jetMuIndex, "jetMuIndex[nJets]/I")
    otree.Branch("c21",        t_c21,        "c21/F")
    otree.Branch("c12",        t_c12,        "c12/F")
    otree.Branch("b21",        t_b21,        "b21/F")
    otree.Branch("b12",        t_b12,        "b12/F")
    otree.Branch("mjj",        t_mjj,        "mjj/F")
    otree.Branch("ptjj",       t_ptjj,       "ptjj/F")
    otree.Branch("dYjj",       t_dYjj,       "dYjj/F")
    otree.Branch("dPhijj",     t_dPhijj,     "dPhijj/F")
    otree.Branch("mbb",        t_mbb,        "mbb/F")
    otree.Branch("ptbb",       t_ptbb,       "ptbb/F")
    otree.Branch("dYbb",       t_dYbb,       "dYbb/F")
    otree.Branch("dPhibb",     t_dPhibb,     "dPhibb/F")
    otree.Branch("met",        t_met,        "met/F")
    otree.Branch("metphi",     t_metphi,     "metphi/F")
    otree.Branch("higgsPt",    t_higgsPt,    "higgsPt/F")
    otree.Branch("higgsM",     t_higgsM,     "higgsM/F")
    otree.Branch("higgsY",     t_higgsY,     "higgsY/F")
    otree.Branch("higgsPhi",   t_higgsPhi,   "higgsPhi/F")
    otree.Branch("higgsEta",   t_higgsEta,   "higgsEta/F")
    otree.Branch("weight",     t_weight,     "weight/F")
    otree.Branch("kWeight",    t_kWeight,    "kWeight/F")
    otree.Branch("kFile",      t_kFile,      "kFile/I")
    otree.Branch("kBadFlag",   t_kBadFlag,   "kBadFlag/I")
    
    def fillJets(collection):
        nmax = min(len(collection), nJetsMax)
        t_nJets[0] = nmax 
        for iobj, obj in enumerate(collection[0:nmax]):
            t_jetPt[iobj]      = round(obj.pt() ,  1)
            t_jetEta[iobj]     = round(obj.eta(),  2)
            t_jetPhi[iobj]     = round(obj.phi(),  2)
            t_jetM[iobj]       = round(obj.m(), 1)
            t_jetPV1[iobj]     = round(obj.pv1,     7)
            t_jetPV2[iobj]     = round(obj.pv2,     7)
            t_jetPV3[iobj]     = round(obj.pv3,     7)
            t_jetPVA[iobj]     = round(obj.pva,     5)
            t_jetPVM[iobj]     = round(obj.pvm,     7)
            #t_jetFlag[iobj]    = obj.flag
            #t_jetBtag[iobj]    = obj.btag
            #t_jetMuIndex[iobj] = obj.muIndex
    
    nMuonsMax = 10
    t_nMuons        = array('i', [0])
    t_muonPt        = array('f', [0 for i in range(nMuonsMax)])
    t_muonEta       = array('f', [0 for i in range(nMuonsMax)])
    t_muonPhi       = array('f', [0 for i in range(nMuonsMax)])
    t_muonM         = array('f', [0 for i in range(nMuonsMax)])
    t_muonPtRel     = array('f', [0 for i in range(nMuonsMax)])
    t_muonPt03      = array('f', [0 for i in range(nMuonsMax)])
    otree.Branch("nMuons",     t_nMuons,   "nMuons/I")
    otree.Branch("muonPt",     t_muonPt,     "muonPt[nMuons]/F")
    otree.Branch("muonEta",    t_muonEta,    "muonEta[nMuons]/F")
    otree.Branch("muonPhi",    t_muonPhi,    "muonPhi[nMuons]/F")
    otree.Branch("muonM",      t_muonM,      "muonM[nMuons]/F")
    otree.Branch("muonPtRel",  t_muonPtRel,   "muonPtRel[nMuons]/F")
    otree.Branch("muonPt03",   t_muonPt03,   "muonPt03[nMuons]/F")

    def fillMuons(collection):
        nmax = min(nMuonsMax, len(collection))
        t_nMuons[0] = nmax 
        for iobj, obj in enumerate(collection[0:nmax]):
            t_muonPt[iobj]    = round(obj.pt() ,  1)
            t_muonEta[iobj]   = round(obj.eta(),  2)
            t_muonPhi[iobj]   = round(obj.phi(),  2)
            t_muonM[iobj]     = round(obj.mass(), 2)
            t_muonPtRel[iobj] = obj.PtRel
            t_muonPt03[iobj]  = obj.Pt03
        





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
            t_kFile[0]        = ii
            count.sumW        += t_kWeight[0]
            count.sumW2       += t_kWeight[0]*t_kWeight[0]
            hSumW.Fill(0, t_kWeight[0])


            numparticles  = getattr(event, "numparticles")
            objects = getattr(event,"objects")
            momenta = get_momenta(objects, numparticles)

            # by default all events are good, i.e., have kBadFlag = 0
            t_kBadFlag[0] = 0
            if numparticles >= maxNpartcl(): t_kBadFlag[0] = 1 

            if args.debug >= 3: print('numparticles %d'%numparticles)
            if args.debug >= 3: print(momenta.size)

            METpx = 0
            METpy = 0
            invisible = [1000022, 1000012, 1000014, 1000016, 2000012, 2000014, 2000016, 1000039, 5100039, 4000012, 4000014, 4000016, 9900012, 9900014, 9900016, 39, 12, 14, 16, 25]
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

                if pid != 25 and abs(pid) != 12 and abs(pid) != 14 and abs(pid) != 16: 
                    momtocluster.append(momenta[mm])
                    METpx -= px
                    METpy -= py
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
                jcs = [jc for jc in jet.constituents() if jc.pt()>jcPtMin and abs(jc.eta())<jcEtaMax]
                fillPV(jet, jcs)
                fillPV3(jet, jcs)
                jet.p4 = toTLV(jet)

            # alias to jets for code combatibility
            genjs = jets 

            # for dijet events
            if len(genjs) >= 2:
                t_mjj[0]  = round(sumP4(genjs[0:2]).M() , 1)
                t_ptjj[0] = round(sumP4(genjs[0:2]).Pt() , 1)
                t_dYjj[0]   = 0. if len(genjs) < 2 else round(genjs[0].rapidity() - genjs[1].rapidity()  , 1)
                t_dPhijj[0] = 0. if len(genjs) < 2 else round(deltaPhi(genjs[0].phi(), genjs[1].phi())   , 1)
                j1 = genjs[0]
                j2 = genjs[1]
                theta21 = RPA(j1, j2)
                theta12 = RPA(j2, j1)
                t_c21[0]  = round(theta21 , 2)
                t_c12[0]  = round(theta12 , 2)
                if args.debug == 2: print('dijet found, calculating the pull vector for the first 2 leading jets (magnitute, angle) [%2.4f, %2.2f], [%2.4f, %2.2f] :'%(j1.pvm, j1.pva, j2.pvm, j2.pva))


            # compute gen met
            t_met[0] = np.sqrt(METpx*METpx + METpy*METpy)
            t_metphi[0] = np.arctan2(METpy, METpx)

            # store Higgs if it's there
            if (len(higgs)>0): 
                t_higgsPt[0] = higgs[0].Pt()
                t_higgsM[0] = higgs[0].M()
                t_higgsY[0] = higgs[0].Rapidity()
                t_higgsPhi[0] = higgs[0].Phi()
                t_higgsEta[0] = higgs[0].Eta()

            #clusters_jets.append(cluster)
            #unclustered.append(momfj_unclustered)

# rapidity () 

#            for ii in range(numparticles):
#                E  = event.objects[0][ii] 
#                px = event.objects[1][ii] 
#                py = event.objects[2][ii] 
#                pz = event.objects[3][ii] 
#                Id = event.objects[4][ii]
#                Ch = event.objects[5][ii] 


            # fill tree
            fillJets  (jets)
            #fillMuons (genmus)
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
    ofile.Write()
    ofile.Close()
