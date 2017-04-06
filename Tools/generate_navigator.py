import os
import numpy as np
import ogr
from collections import defaultdict, OrderedDict


def read_dbf(dbf_file, out_fields=None):
    """ Reads the contents of a dbf table """

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
        table = [[row.GetField(f) for f in fields] for row in layer]
    else:
        table = [row.GetField(list(fields)[0]) for row in layer]
    # noinspection PyUnusedLocal
    data_source = None
    return table


# Assembles a dictionary of NHD Plus directory structure indexed by region
def get_nhd(nhd_dir):
    nhd_states = OrderedDict((('01', {"ME", "NH", "VT", "MA", "CT", "RI", "NY"}),
                              ('02', {"VT", "NY", "PA", "NJ", "MD", "DE", "WV", "DC", "VA"}),
                              ('03N', {"VA", "NC", "SC", "GA"}),
                              ('03S', {"FL", "GA"}),
                              ('03W', {"FL", "GA", "TN", "AL", "MS"}),
                              ('04', {"WI", "MN", "MI", "IL", "IN", "OH", "PA", "NY"}),
                              ('05', {"IL", "IN", "OH", "PA", "WV", "VA", "KY", "TN"}),
                              ('06', {"VA", "KY", "TN", "NC", "GA", "AL", "MS"}),
                              ('07', {"MN", "WI", "SD", "IA", "IL", "MO"}),
                              ('08', {"MO", "KY", "TN", "AR", "MS", "LA"}),
                              ('09', {"ND", "MN", "SD"}),
                              ('10U', {"MT", "ND", "WY", "SD", "MN", "NE", "IA"}),
                              ('10L', {"CO", "WY", "MN", "NE", "IA", "KS", "MO"}),
                              ('11', {"CO", "KS", "MO", "NM", "TX", "OK", "AR", "LA"}),
                              ('12', {"NM", "TX", "LA"}),
                              ('13', {"CO", "NM", "TX"}),
                              ('14', {"WY", "UT", "CO", "AZ", "NM"}),
                              ('15', {"NV", "UT", "AZ", "NM", "CA"}),
                              ('16', {"CA", "OR", "ID", "WY", "NV", "UT"}),
                              ('17', {"WA", "ID", "MT", "OR", "WY", "UT", "NV"}),
                              ('18', {"OR", "NV", "CA"})))

    all_paths = defaultdict()
    regions = list(nhd_states.keys())
    region_dirs = {"NHDPlus{}".format(region) for region in regions}
    for root_dir, sub_dirs, _ in os.walk(nhd_dir):
        if set(sub_dirs) & region_dirs:
            for sub_dir in sub_dirs:
                region = sub_dir.lstrip("NHDPlus")
                if region in regions:
                    all_paths[sub_dir.lstrip("NHDPlus")] = os.path.join(root_dir, sub_dir)
    return OrderedDict(sorted(all_paths.items()))


def generate_navigator(region, region_dir, output_directory, include_time=False, include_gridcode=False,
                       filename="upstream", overwrite=True):
    def convert_attributions(cd, attrs):
        return [(name, {cd.get(k, -1): v for k, v in a.items()}) for name, a in attrs.items()]

    def dict_to_array(cd):
        max_alias = max(cd.values())
        out_array = np.zeros(max_alias + 1)
        for comid, alias in cd.items():
            out_array[alias] = comid
        return out_array

    def map_paths(paths):

        # Get starting row and column for each value
        column_numbers = np.tile(np.arange(paths.shape[1]) + 1, (paths.shape[0], 1)) * (paths > 0)
        path_begins = np.argmax(column_numbers > 0, axis=1)
        max_reach = np.max(paths)
        path_map = np.zeros((max_reach + 1, 3))
        n_paths = paths.shape[0]
        for i, path in enumerate(paths):
            for j, val in enumerate(path):
                if val:
                    if i == n_paths:
                        end_row = 0
                    else:
                        next_row = (path_begins[i + 1:] <= j)
                        if next_row.any():
                            end_row = np.argmax(next_row)
                        else:
                            end_row = n_paths - i - 1
                    values = np.array([i, i + end_row + 1, j])
                    path_map[val] = values

        return path_map

    def get_gridcode(gridcode_table):
        return dict(read_dbf(gridcode_table, ["GRIDCODE", "FEATUREID"]))

    def get_nodes(flow_table, vaa_table, type_table):

        # Get a set of nodes
        flows = read_dbf(flow_table, ["TOCOMID", "FROMCOMID"])  # Get all nodes
        streamcalc = dict(read_dbf(vaa_table, ["ComID", "StreamCalc"]))
        divergence = dict(read_dbf(vaa_table, ["ComID", "Divergence"]))
        ftype = dict(read_dbf(type_table, ["ComID", "FCode"]))
        nodes = list(filter(all, flows))  # Filter out nodes where from or to value is zero
        nodes = list(filter(lambda x: ftype.get(x[1]) != 56600, nodes))

        # Convert to indices
        unique_comids = np.unique(nodes)
        conversion_dict = {comid: i + 1 for i, comid in enumerate(unique_comids)}
        nodes = np.vectorize(lambda x: conversion_dict.get(x))(nodes)
        streamcalc = dict(zip(map(conversion_dict.get, streamcalc.keys()), streamcalc.values()))
        divergence = dict(zip(map(conversion_dict.get, divergence.keys()), divergence.values()))

        # Find connections to sever
        sc_nodes = np.array([bool(streamcalc.get(n)) for n in nodes.flat]).reshape(nodes.shape)
        div_nodes = np.array([bool(divergence.get(n) == 2) for n in nodes.flat]).reshape(nodes.shape)
        sever = (((sc_nodes[:, 1] == True) & (sc_nodes[:, 0] == False)) |  # Streamcalc 0 -> Streamcalc 1
                 ((div_nodes[:, 1] == False) & (div_nodes[:, 0] == True)))  # Divergence 0 -> Divergence 2
        active_nodes = nodes[~sever]

        return active_nodes, conversion_dict

    def get_outlets(flow_table, flow_lines, vaa_table, conversion_dict):

        # Identify reaches that empty into another NHD region
        in_region = set(read_dbf(flow_lines, ["ComID"]))
        nodes = read_dbf(flow_table, ["FROMCOMID", "TOCOMID"])
        _, ds_reaches = zip(*nodes)
        basin_outlets = set(ds_reaches) - in_region
        basin_outlets = {x[0] for x in nodes if x[1] in basin_outlets}

        # Identify all the reaches that are designated as a terminal path (has to be converted from HydroSeq)
        c, h, t = zip(*read_dbf(vaa_table, ["ComID", "HydroSeq", "TerminalPa"]))
        terminal_paths = set(map(dict(zip(h, c)).get, set(t) & set(h)))
        outlets = (terminal_paths | basin_outlets) - {0}

        return set(filter(None, map(conversion_dict.get, outlets)))

    def upstream_trace(nodes, outlets, attributes=(), max_length=3000, max_paths=500000):

        n_attributes = len(attributes)

        # Output arrays
        all_paths = np.zeros((max_paths, max_length), dtype=np.int32)
        all_attributes = np.zeros((n_attributes, max_paths, max_length), dtype=np.float)

        # Bounds
        path_cursor = 0
        longest_path = 0

        for index, start_node in enumerate(outlets):

            # Reset everything except the master. Trace is done separately for each outlet
            # Holders
            queue = np.zeros((nodes.shape[0], 2))
            active_reach = np.zeros(max_length, dtype=np.int32)
            active_attributes = np.zeros((n_attributes, max_length), dtype=np.float)

            # Cursors
            start_cursor = 0
            queue_cursor = 0
            active_reach_cursor = 0
            active_node = start_node

            while True:
                upstream = nodes[nodes[:, 0] == active_node]
                active_reach[active_reach_cursor] = active_node
                for i, attribute_dict in enumerate(attributes):
                    active_attributes[i, active_reach_cursor] = attribute_dict.get(active_node, 0.0)
                active_reach_cursor += 1
                if active_reach_cursor > longest_path:
                    longest_path = active_reach_cursor
                if upstream.size:
                    active_node = upstream[0][1]
                    for j in range(1, upstream.shape[0]):
                        queue[queue_cursor] = upstream[j]
                        queue_cursor += 1
                else:
                    all_paths[path_cursor, start_cursor:] = active_reach[start_cursor:]
                    for i in range(n_attributes):
                        all_attributes[i, path_cursor] = np.cumsum(active_attributes[i]) * (all_paths[path_cursor] > 0)
                    queue_cursor -= 1
                    path_cursor += 1
                    last_node, active_node = queue[queue_cursor]
                    if not any((last_node, active_node)):
                        break
                    active_reach_cursor = np.where(active_reach == last_node)[0][0] + 1
                    start_cursor = active_reach_cursor
                    active_reach[active_reach_cursor:] = 0.
                    for i in range(n_attributes):
                        active_attributes[i, active_reach_cursor:] = 0.

        return all_paths[:path_cursor, :longest_path + 1], all_attributes[:, :path_cursor, :longest_path + 1]

    def write_outfile(outfile, paths, path_map, conversion_array, attributes, names):
        attribute_dict = {name: attributes[i] for i, name in enumerate(names)}
        np.savez_compressed(outfile, paths=paths, path_map=path_map, alias_index=conversion_array, **attribute_dict)

    def collapse_array(paths, attribute_matrix):
        out_paths = []
        out_attributes = [[] for _ in range(attribute_matrix.shape[0])]
        path_starts = []
        for i, row in enumerate(paths):
            active_path = (row > 0)
            path_starts.append(np.argmax(active_path))
            out_paths.append(row[active_path])
            for j in range(attribute_matrix.shape[0]):
                out_attributes[j].append(attribute_matrix[j][i][active_path])
        out_attributes = np.array(out_attributes)
        return np.array(out_paths), np.array(out_attributes), np.array(path_starts)

    def fetch_times():
        erom_table = os.path.join(region_dir, "EROMExtension", "EROM_MA0001.dbf")
        velocity_dict = table_to_dict(erom_table, "V0001E")
        length_dict = table_to_dict(vaa_table, "LengthKM")
        time_dict = {}
        for comid in set(velocity_dict.keys()) & set(length_dict.keys()):
            try:
                time_dict[comid] = velocity_dict.get(comid, 0.0) / length_dict.get(comid, 0.0)
            except ZeroDivisionError:
                time_dict[comid] = 0.0
        return time_dict

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Designate paths
    flow_table = os.path.join(region_dir, "NHDPlusAttributes", "PlusFlow.dbf")
    flow_lines = os.path.join(region_dir, "NHDSnapshot", "Hydrography", "NHDFlowline.dbf")
    vaa_table = os.path.join(region_dir, "NHDPlusAttributes", "PlusFlowlineVAA.dbf")
    gridcode_table = os.path.join(region_dir, "NHDPlusCatchment", "featureidgridcode.dbf")

    outfile = os.path.join(output_directory, filename + "_{}.npz".format(region))

    if not os.path.exists(outfile) or overwrite:
        # Do the work
        all_nodes, conversion_dict = get_nodes(flow_table, vaa_table, flow_lines)

        attributions = {}
        if include_time:
            travel_times = fetch_times()
            attributions['time'] = travel_times
        if include_gridcode:
            attributions['gridcode'] = get_gridcode(gridcode_table)

        attributions = convert_attributions(conversion_dict, attributions)

        if attributions:
            names, attribution_arrays = zip(*attributions)
            attribution_array = np.array(attribution_arrays)

        outlets = get_outlets(flow_table, flow_lines, vaa_table, conversion_dict)

        paths, attribute_matrix = upstream_trace(all_nodes, outlets, attribution_array)
        path_map = map_paths(paths)

        collapse = True
        if collapse:
            paths, attribute_matrix, start_cols = collapse_array(paths, attribute_matrix)

        conversion_array = dict_to_array(conversion_dict)

        write_outfile(outfile, paths, path_map, conversion_array, attribute_matrix, names)


def table_to_dict(table, field, index_field="ComID"):
    out_dict = dict(read_dbf(table, [index_field, field]))
    out_dict.pop(None, None)
    return out_dict


def main():
    # Output path
    output_folder = "../WatershedTopology/Gridcode_and_Time"
    include_time = True
    include_gridcode = True

    for region, region_dir in get_nhd(nhd_dir=r"T:\NationalData\NHDPlusV2").items():
        print(region)
        # Attribution tables (if any: format is {"Name": attribution dictionary...}
        generate_navigator(region, region_dir, output_folder, include_time, include_gridcode, overwrite=False)


main()
