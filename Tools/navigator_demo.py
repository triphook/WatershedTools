from Tools.attribution_and_accumulation import Navigator

region_id = '07'
topology_file = r"..\WatershedTopology\upstream_{}.npz".format(region_id)

region = Navigator(topology_file)

test_reach = 4867727

n = len(region.all_upstream(test_reach))

print(n)