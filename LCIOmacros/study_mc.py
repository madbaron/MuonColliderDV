from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch

#########################
# parameters

Bfield = 3.56  # T

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_truth.root',
                  type=str, default='ntup_truth.root')
(options, args) = parser.parse_args()
#########################


def find_C1_in_decay_chain(list_particles):
    result = []

    for mcp in mcpCollection:
        pdg = mcp.getPDG()
        if fabs(pdg) == 1000005:
            daughters = mcp.getDaughters()
            if len(daughters) == 0:
                print('WARNING: Found chargino with no daughters')
                continue
            else:
                for d in daughters:
                    if d.getPDG() == 1000022:
                        result.append(mcp)
    return result


#########################
# declare histograms
photon_pt = TH1D('photon_pt', 'photon_pt', 100, 0, 5000)  # GeV
stop_beta = TH1D('stop_beta', 'stop_beta', 100, 0, 1)  # GeV
stop_theta = TH1D('stop_theta', 'stop_theta', 100, 0., 180.)  # GeV
stop_gamma = TH1D('stop_gamma', 'stop_gamma', 60, 0, 6)  # GeV
stop_gammabeta = TH1D('stop_gammabeta', 'stop_gammabeta', 120, 0, 6)  # GeV
bquark_pt = TH1D('bquark_pt', 'bquark_pt', 100, 0, 5)  # GeV
mc_pt = TH1D('mc_pt', 'mc_pt', 100, 0, 5000)  # GeV
mc_pz = TH1D('mc_pz', 'mc_pz', 200, -10000, 10000)  # GeV
mc_E = TH1D('mc_E', 'mc_E', 200, -10000, 10000)  # GeV
mc_m = TH1D('mc_m', 'mc_m', 200, -10000, 10000)  # GeV
mc_charge = TH1D('mc_charge', 'mc_charge', 4, -2, 2)  # GeV
mc_genStatus = TH1D('mc_genStatus', 'mc_genStatus', 40, 0,
                    40)  # it should be in the 21-23 range
mc_endPoint = TH1D('mc_endPoint', 'mc_endPoint', 100, 0., 1500.)  # mm
B_endPoint = TH1D('B_endPoint', 'B_endPoint', 40, 0., 20.)  # mm
mc_radialEndPoint = TH1D(
    'mc_radialEndPoint', 'mc_radialEndPoint', 100, 0., 1500.)  # mm
mc_pdg = TH1D('mc_pdg', 'mc_pdg', 2, -2., 2.)
nC1perEvent = TH1D('nC1perEvent', 'nC1perEvent', 4, 0, 4)

histos_list = [photon_pt, stop_beta, stop_theta, stop_gamma, stop_gammabeta, bquark_pt, mc_pt, 
                mc_pz, mc_E, mc_m, mc_charge, mc_genStatus, mc_endPoint, B_endPoint, 
                mc_radialEndPoint, mc_pdg, nC1perEvent]

for histo in histos_list:
    histo.SetDirectory(0)

#########################

tree = TTree("truth_tree", "truth_tree")

# create 1 dimensional float arrays as fill variables, in this way the float
# array serves as a pointer which can be passed to the branch
pt_vec = array('d', [0])
phi_vec = array('d', [0])
theta_vec = array('d', [0])
z_truth_vec = array('d', [0])
r_truth_vec = array('d', [0])
pdgID_vec = array('i', [0])

# create the branches and assign the fill-variables to them as doubles (D)
tree.Branch("pT",  pt_vec,  'var/D')
tree.Branch("phi", phi_vec, 'var/D')
tree.Branch("theta", theta_vec, 'var/D')
tree.Branch("z_truth", z_truth_vec, 'var/D')
tree.Branch("r_truth", r_truth_vec, 'var/D')
tree.Branch("pdgID", pdgID_vec, 'var/I')

# create a reader and open an LCIO file
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(options.inFile)

# loop over all events in the file
for ievt, event in enumerate(reader):

    if ievt % 100 == 0:
        print("Processing event " + str(ievt))

    # find the last stops
    mcpCollection = event.getCollection('MCParticle')
    for part in mcpCollection:
        if fabs(part.getPDG()) == 22:
            p = part.getMomentum()  # GeV
            pt = sqrt(p[0]*p[0] +
                      p[1]*p[1])
            photon_pt.Fill(pt)
        if fabs(part.getPDG()) == 5:
            p = part.getMomentum()  # GeV
            pt = sqrt(p[0]*p[0] +
                      p[1]*p[1])
            bquark_pt.Fill(pt)

    result = find_C1_in_decay_chain(mcpCollection)
    nC1 = len(result)

    good_pid = [1000005, 1005321, 1000522, 1000512]

    for c in mcpCollection:
        #        for c in result:
        if fabs(c.getPDG()) in good_pid:
            end = c.getEndpoint()
            endpoint = sqrt(
                end[0]*end[0] + end[1]*end[1] + end[2]*end[2])
            pos = c.getVertex()
            print(str(c.getPDG()) + " at " + str(pos[0]) + "  " +
                  str(pos[1]) + "  " + str(pos[2]))
            print("decay at " + str(end[0]) + "  " +
                  str(end[1]) + "  " + str(end[2]))

            daughters = c.getDaughters()
            dau_list = []

            for d in daughters:
                dau_list.append(d.getPDG())
                if fabs(d.getPDG()) == 5:
                    # found the b, now we have to figure out where it's going
                    b_daughter_pgdid = 5
                    the_part = d
                    n_step = 0
                    while b_daughter_pgdid == 5:
                        #print("Step "+str(n_step))
                        b_daughters = the_part.getDaughters()
                        for k in b_daughters:
                            b_daughter_pgdid = fabs(k.getPDG())
                            end = k.getEndpoint()
                            endpoint = sqrt(
                                end[0]*end[0] + end[1]*end[1] + end[2]*end[2])
                            momentum = k.getMomentum()  # GeV
                            pt = sqrt(
                                momentum[0]*momentum[0] + momentum[1]*momentum[1])
                            mass = k.getMass()
                            # print(" " + str(k.getPDG()) +
                            #      " " + str(endpoint) + " " + str(pt) + " " + str(mass))
                            if fabs(b_daughter_pgdid) > 500 and fabs(b_daughter_pgdid) < 600:
                                B_endPoint.Fill(endpoint)
                        the_part = k
                        n_step = n_step + 1

            print("daughters " + str(dau_list))

            endpoint = c.getEndpoint()
            endpoint_l = sqrt(
                endpoint[0]*endpoint[0] + endpoint[1]*endpoint[1] + endpoint[2]*endpoint[2])
            if endpoint_l > 0:
                endpoint_r = sqrt(
                    endpoint[0]*endpoint[0] + endpoint[1]*endpoint[1])
                momentum = c.getMomentum()  # GeV
                pt = sqrt(momentum[0]*momentum[0] +
                          momentum[1]*momentum[1])
                pz = momentum[2]
                momentum = c.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(momentum[0], momentum[1],
                               momentum[2], c.getEnergy())
                stop_beta.Fill(tlv.Beta())
                stop_gamma.Fill(tlv.Gamma())
                stop_gammabeta.Fill(tlv.Gamma()*tlv.Beta())
                theta = tlv.Theta()*180./TMath.Pi()
                stop_theta.Fill(theta)
                mc_pt.Fill(pt)
                mc_pz.Fill(pz)
                mc_E.Fill(c.getEnergy())
                mc_m.Fill(c.getMass())
                mc_charge.Fill(c.getCharge())
                mc_genStatus.Fill(c.getGeneratorStatus())
                mc_radialEndPoint.Fill(endpoint_r)  # mm
                mc_endPoint.Fill(endpoint_l)  # mm
                c_pdg = c.getPDG()
                if c_pdg > 0:
                    mc_pdg.Fill(1)
                else:
                    mc_pdg.Fill(-1)
        
        vx = c.getVertex()
        dp3 = c.getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], c.getEnergy())

        pt_vec[0] = tlv.Perp()
        phi_vec[0] = tlv.Phi()
        theta_vec[0] = tlv.Theta()
        z_truth_vec[0] = vx[2] 
        r_truth_vec[0] = sqrt(vx[0]*vx[0]+vx[1]*vx[1])
        pdgID_vec[0] = c.getPDG()

        tree.Fill()

    nC1perEvent.Fill(nC1)

reader.close()

# write outputs
output_file = TFile(options.outFile, 'RECREATE')
tree.Write()
for histo in histos_list:
    histo.Write()
output_file.Close()
