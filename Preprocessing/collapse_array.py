import numpy as np
import os
import pickle
from array_tools import get_nhd, read_dbf

def collapse_array(paths, a2g):
    out_paths = []
    for i, row in enumerate(paths):
        path_ids = row[np.where(row > 0)[0]]
        path_gridcodes = list(map(lambda x: a2g.get(x, 0), path_ids))
        out_paths.append(path_gridcodes)
    return out_paths


for region, region_dir in get_nhd().items():

    # Initialize paths
    g2f_table = os.path.join(region_dir, "NHDPlusCatchment", "Catchment.dbf")
    nhd_obj = os.path.join(r"T:\CA3T\bin", "upstream_{}.npz".format(region))
    out_file = nhd_obj.rstrip(".npz") + "_gc.p"

    # Generate gridcode to featureid array
    gridcodes, comids = np.array(read_dbf(g2f_table, ["GRIDCODE", "FEATUREID"])).T
    gridcode_array = np.zeros(gridcodes.max() + 1, dtype=[('comid', 'i4'), ('alias', 'i4')])
    gridcode_array['comid'][gridcodes] = comids

    # Generate gridcode to alias array
    alias_to_comid = np.load(nhd_obj)['conversion_array']
    comid_to_alias = dict(zip(alias_to_comid, np.arange(alias_to_comid.size)))
    gridcode_array['alias'] = np.vectorize(lambda x: comid_to_alias.get(x, 0))(gridcode_array['comid'])
    alias_to_gridcode = dict(zip(gridcode_array['alias'], np.arange(gridcode_array['alias'].size)))

    arrays = np.load(nhd_obj)
    paths, times, path_map, alias_to_reach = \
        map(lambda x: arrays[x], ['paths', 'times', 'path_map', 'conversion_array'])


    paths = collapse_array(paths, alias_to_gridcode)

    with open(out_file, 'wb') as g:
        pickle.dump([paths, np.int32(path_map)], g)
        print(out_file)
