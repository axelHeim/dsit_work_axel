""" Attempt to reconstruct decay mode 3 (of the 6 modes in use on May 20th '21) """

import basf2
from basf2 import *

import sys
import basf2 as b2
import modularAnalysis as ma
import stdV0s
import variables.collections as vc
import variables.utils as vu




main = b2.create_path()

ma.inputMdstList('default', [], main)


ma.fillParticleList(
    "pi-:FSP",
    "pionID > 0.1 and dr < 0.5 and abs(dz) < 2 and thetaInCDCAcceptance",
    path=main,
)
ma.fillParticleList(
    "K-:FSP",
    "pionID > 0.1 and dr < 0.5 and abs(dz) < 2 and thetaInCDCAcceptance",
    path=main,
)

ma.reconstructDecay(
    "D-:K_2pi -> K+:FSP pi-:FSP pi-:FSP",
    cut="",
    path=main,
)

ma.reconstructDecay(
    "B+ -> D-:K_2pi pi+:FSP pi+:FSP",
    cut="",
    path=main,
)


# Create list of variables to save into the output file
b_vars = []

standard_vars = vc.kinematics + vc.mc_kinematics + vc.mc_truth
b_vars += vc.deltae_mbc
b_vars += standard_vars

# Variables for final states (electrons, positrons, pions)
fs_vars = vc.pid + vc.track + vc.track_hits + standard_vars
b_vars += vu.create_aliases_for_selected(
    fs_vars,
    "B+ -> [D- -> ^K+ ^pi- ^pi-] ^pi+ ^pi+",
    prefix=["K_p", "pim_1", "pim_2", "pip_1", "pip_2"],
)

# match reconstructed with MC particles
ma.matchMCTruth("B+", path=main)










# Save variables to an output file (ntuple)
ma.variablesToNtuple(
    "B+",
    variables=b_vars,
    filename="B+2_D-_2pi+.root",
    treename="tree",
    path=main,
)




process(main, max_event=1000)


# print out the summary
print(b2.statistics)