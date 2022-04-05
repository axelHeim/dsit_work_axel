# +
import basf2
from basf2 import *

import sys
import basf2 as b2
import modularAnalysis as ma
import stdV0s
import variables.collections as vc
import variables.utils as vu



main = create_path()

inputMdstList('default', [], main)

process(path, max_event=10000)
