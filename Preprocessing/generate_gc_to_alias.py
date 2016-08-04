import os
import numpy as np
import collections
import pickle

# Reads the contents of a dbf table
def read_dbf(dbf_file, out_fields=None, include_fields=False):
    import ogr

    # Initialize file
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(dbf_file)
    layer = data_source.GetLayer(0)

    # Set fields
    ld = layer.GetLayerDefn()
    fields = {ld.GetFieldDefn(i).GetName() for i in range(ld.GetFieldCount())}
    if out_fields:
        missing_fields = set(out_fields) - fields
        remedies = [lambda x: x.upper(), lambda x: x.capitalize()]
        while missing_fields and remedies:
            new_fields = map(remedies.pop(0), out_fields)
            out_fields = [nf if nf in fields else out_fields[i] for i, nf in enumerate(new_fields)]
            missing_fields = set(out_fields) - fields
        if missing_fields:
            print("Fields {} not found in {}".format(out_fields, dbf_file))
        fields = [field for field in out_fields if not field in missing_fields]

    # Read data
    if len(fields) > 1:
        if include_fields:
            table = [{f: row.GetField(f) for f in fields} for row in layer]
        else:
            table = [[row.GetField(f) for f in fields] for row in layer]
    else:
        table = [row.GetField(list(fields)[0]) for row in layer]
    # noinspection PyUnusedLocal
    data_source = None
    return table


# Assembles a dictionary of NHD Plus directory structure indexed by region
def get_nhd(nhd_dir=r"T:\NationalData\NHDPlusV2", region_filter='all'):
    from collections import OrderedDict
    # Get catchment grid, gridcode to comid translation files, and flow tables
    regions = {'01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09',
               '10U', '10L', '11', '12', '13', '14', '15', '16', '17', '18'}
    all_paths = collections.defaultdict()
    region_dirs = {"NHDPlus{}".format(region) for region in regions}
    for root_dir, sub_dirs, _ in os.walk(nhd_dir):
        if set(sub_dirs) & region_dirs:
            for sub_dir in sub_dirs:
                region = sub_dir.lstrip("NHDPlus")
                if region in regions:
                    all_paths[sub_dir.lstrip("NHDPlus")] = os.path.join(root_dir, sub_dir)
    return OrderedDict(sorted(all_paths.items()))


for region, region_dir in get_nhd().items():
    print(region)

    # Initialize paths
    g2f_table = os.path.join(region_dir, "NHDPlusCatchment", "Catchment.dbf")
    nhd_obj = os.path.join(r"T:\CA3T\bin", "upstream_{}.npz".format(region))
    outfile = r"T:\CA3T\bin\id_lookup_{}.npz".format(region)

    # Generate gridcode to featureid array
    gridcodes, comids = np.array(read_dbf(g2f_table, ["GRIDCODE", "FEATUREID"])).T
    gridcode_array = np.zeros(gridcodes.max() + 1, dtype=[('comid', 'i4'), ('alias', 'i4')])
    gridcode_array['comid'][gridcodes] = comids

    # Generate gridcode to alias array
    alias_to_comid = np.load(nhd_obj)['conversion_array']
    comid_to_alias = dict(zip(alias_to_comid, np.arange(alias_to_comid.size)))
    gridcode_array['alias'] = np.vectorize(lambda x: comid_to_alias.get(x, 0))(gridcode_array['comid'])

    np.savez_compressed(outfile, gc_lookup=gridcode_array)
