"""This script shows how to calculate the number of reaches upstream of a target"""

import os
import sys
sys.path.append('attribution_and_accumulation')

from attribution_and_accumulation import Navigator

REGION_ID = '07'
TOPOLOGY_FILE = "upstream_{}.npz".format(REGION_ID)
TOPOLOGY_FILE = os.path.join("WatershedTopology", TOPOLOGY_FILE)

REGION = Navigator(TOPOLOGY_FILE)

TEST_REACH = 4867727

N = len(REGION.all_upstream(TEST_REACH))

print(N)
