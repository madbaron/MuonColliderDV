from array import array
import os
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TVector3, TTree, TMath
from math import *
from optparse import OptionParser

#########################
# parameters

Bfield = 3.56  # T

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-t', '--trackCollection', help='--trackCollection SiTracks_Refitted',
                  type=str, default='SiTracks_Refitted')
parser.add_option('-o', '--outDir', help='--outDir ./',
                  type=str, default='./')
(options, args) = parser.parse_args()


def getOriginPID(mcp):
    # Look for sbottom mothers
    origin_PDGid = 0
    momVec = mcp.getParents()
    while (len(momVec) > 0 and fabs(origin_PDGid) != 1000005):
        mc_mother = momVec[0]
        origin_PDGid = mc_mother.getPDG()
        momVec = mc_mother.getParents()

    return origin_PDGid

#########################
# declare histograms


arrBins_R = array('d', (0., 10., 20., 31., 51., 74., 102.,
                        127., 150., 200., 250., 340., 450., 554.))
arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
                         3.5, 4., 5., 6., 7., 8., 10.,20.,30.,50.))
# arrBins_pT = array('d', (0., 0.5, 1., 1.5, 2., 2.5, 3.,
#                         3.5, 4., 5., 6., 7., 8., 10., 20., 30., 50., 75., 100., 250., 500., 1000., 1500.))
arrBins_theta = array('d', (30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180.))
arrBins_dz0 = array('d', (-2., -1., -0.8, -0.6, -0.5, -0.4, -
                          0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1., 2.))

h_truth_Rprod = TH1D('truth_Rprod', 'truth_Rprod',
                     len(arrBins_R)-1, arrBins_R)  # mm
h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_pT)-1, arrBins_pT)
h_truth_theta = TH1D('truth_theta', 'truth_theta',
                     len(arrBins_theta)-1, arrBins_theta)
h_truth_phi = TH1D('truth_phi', 'truth_phi', 20, -TMath.Pi(), TMath.Pi())

h_track_d0 = TH1D('track_d0', 'track_d0', 100, -5., 5.)
h_track_z0 = TH1D('track_z0', 'track_z0', 100, -20., 20.)
h_track_pT = TH1D('track_pT', 'track_pT', len(arrBins_pT)-1, arrBins_pT)
h_track_phi = TH1D('track_phi', 'track_phi', 20, -TMath.Pi(), TMath.Pi())
h_track_theta = TH1D('track_theta', 'track_theta',
                     len(arrBins_theta)-1, arrBins_theta)
h_track_nholes = TH1D('track_nholes', 'track_nholes', 20, 0., 20.)
h_track_nhits = TH1D('track_nhits', 'track_nhits', 20, 0., 20.)
h_track_chi2ndf = TH1D(
    'track_chi2ndf', 'track_chi2ndf', 100, 0., 100.)
h_track_Rprod = TH1D('track_Rprod',
                     'track_Rprod', len(arrBins_R)-1, arrBins_R)  # mm

histos_list = [h_truth_Rprod, h_truth_pT, h_truth_theta, h_truth_phi, h_track_chi2ndf, h_track_d0,
               h_track_z0, h_track_pT, h_track_phi, h_track_theta, h_track_nholes, h_track_nhits, h_track_Rprod]

for histo in histos_list:
    histo.SetDirectory(0)

#########################

tree = TTree("tracks_tree", "tracks_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
pt = array('d', [0])
pt_truth = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
d0 = array('d', [0])
z0 = array('d', [0])
z_truth = array('d', [0])
sigma_d0 = array('d', [0])
sigma_z0 = array('d', [0])
omega = array('d', [0])
chi2 = array('d', [0])
ndf = array('i', [0])
nhits = array('i', [0])
nholes = array('i', [0])
pdgID = array('i', [0])
r_truth = array('d', [0])

# create the branches and assign the fill-variables to them as doubles (D)
tree.Branch("pT",  pt,  'var/D')
tree.Branch("pTtruth",  pt_truth,  'var/D')
tree.Branch("phi", phi, 'var/D')
tree.Branch("theta", theta, 'var/D')
tree.Branch("d0", d0, 'var/D')
tree.Branch("z0", z0, 'var/D')
tree.Branch("z_truth", z_truth, 'var/D')
tree.Branch("sigma_d0", sigma_d0, 'var/D')
tree.Branch("sigma_z0", sigma_z0, 'var/D')
tree.Branch("omega", omega, 'var/D')
tree.Branch("chi2", chi2, 'var/D')
tree.Branch("ndf", ndf, 'var/I')
tree.Branch("nhits", nhits, 'var/I')
tree.Branch("nholes", nholes, 'var/I')
tree.Branch("pdgID", pdgID, 'var/I')
tree.Branch("r_truth", r_truth, 'var/D')

#########################
# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievent, event in enumerate(reader):

    if ievent % 100 == 0:
        print("Processing event " + str(ievent))

    # Look at particles from Rhad decays and their matched tracks
    mcpCollection = event.getCollection('MCParticle')
    relationCollection = event.getCollection(
        'MCParticle_'+options.trackCollection)
    relation = UTIL.LCRelationNavigator(relationCollection)

    '''
    prod = -1
    for mcp in mcpCollection:

        charge = mcp.getCharge()
        status = mcp.getGeneratorStatus()

        if fabs(charge) > 0:
            if fabs(mcp.getPDG()) == 13:
                vx = mcp.getVertex()
                vxv = TVector3(vx[0], vx[1], vx[2])
                dp3 = mcp.getMomentum()
                tlv = TVector3(dp3[0], dp3[1], dp3[2])
                prod = tlv.Dot(vxv)

    if prod < 0:
        print("Gun shot backwards. Skipping")
        continue
    '''

    # getting tracks and relative hits
    tracks = event.getCollection(options.trackCollection)

    for itrack, track in enumerate(tracks):
        pt[0] = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
        phi[0] = track.getPhi()
        theta[0] = TMath.Pi()/2-atan(track.getTanLambda())
        d0[0] = track.getD0()
        z0[0] = track.getZ0()
        sigma_d0[0] = track.getCovMatrix()[0]
        sigma_z0[0] = track.getCovMatrix()[9]
        omega[0] = track.getOmega()
        chi2[0] = track.getChi2()
        ndf[0] = track.getNdf()

        hits = track.getTrackerHits()
        numhits = len(hits)
        holes = int(track.getdEdxError())  # BADHACK
        nhits[0] = numhits
        nholes[0] = holes

        mcpvec = relation.getRelatedFromObjects(track)
        if len(mcpvec) > 0:
            mcp = mcpvec[0]
            pdgID[0] = abs(mcp.getPDG())
            vx = mcp.getVertex()
            z_truth[0] = vx[2]
            dp3 = mcp.getMomentum()
            tlv = TLorentzVector()
            tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())
            pt_truth[0] = tlv.Perp()
            rprod = sqrt(vx[0]*vx[0]+vx[1]*vx[1])
            r_truth[0] = rprod
        else:
            pt_truth[0] = -1.
            z_truth[0] = -99999.
            pdgID[0] = 0
            r_truth[0] = -1.

        tree.Fill()

    # Look at particles from Rhad decays and their matched tracks
    mcpCollection = event.getCollection('MCParticle')
    relationCollection = event.getCollection(
        'MCParticle_'+options.trackCollection)
    #relationCollection = event.getCollection('MCParticle_SelectedTracks')
    relation = UTIL.LCRelationNavigator(relationCollection)

    for mcp in mcpCollection:

        charge = mcp.getCharge()
        status = mcp.getGeneratorStatus()

        if fabs(charge) > 0:
            if fabs(mcp.getPDG()) == 13:
                vx = mcp.getVertex()
                rprod = sqrt(vx[0]*vx[0]+vx[1]*vx[1])
                dp3 = mcp.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())

                if tlv.Perp() > 1 and not mcp.isDecayedInTracker():
                    h_truth_Rprod.Fill(rprod)
                    h_truth_pT.Fill(tlv.Perp())
                    h_truth_phi.Fill(tlv.Phi())
                    h_truth_theta.Fill(tlv.Theta())

                    tracks = relation.getRelatedToObjects(mcp)
                    for track in tracks:
                        h_track_Rprod.Fill(rprod)
                        h_track_pT.Fill(tlv.Perp())
                        h_track_phi.Fill(tlv.Phi())
                        h_track_theta.Fill(tlv.Theta())
                        h_track_d0.Fill(track.getD0())
                        h_track_z0.Fill(track.getZ0())
                        h_track_nholes.Fill(int(track.getdEdxError()))
                        h_track_nhits.Fill(len(track.getTrackerHits()))
                        h_track_chi2ndf.Fill(track.getChi2()/track.getNdf())

reader.close()

# write histograms
output_file = TFile(options.outDir + "ntup_" +
                    options.trackCollection + ".root", 'RECREATE')
#tree.Write()
for histo in histos_list:
    histo.Write()
output_file.Close()
