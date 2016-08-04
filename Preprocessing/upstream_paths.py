import numpy as np
import os


from read_tools import get_nhd, read_dbf

def snapshot(nodes, outlets, time_dict, max_length=3000, max_paths=500000):

    # Output arrays
    all_paths = np.zeros((max_paths, max_length), dtype=np.int32)
    all_times = np.zeros((max_paths, max_length), dtype=np.float)

    # Bounds
    path_cursor = 0
    longest_path = 0

    for index, start_node in enumerate(outlets):

        # Reset everything except the master. Trace is done separately for each outlet
        # Holders
        queue = np.zeros((nodes.shape[0], 2))
        active_reach = np.zeros(max_length, dtype=np.int32)
        active_times = np.zeros(max_length, dtype=np.float)

        # Cursors
        start_cursor = 0
        queue_cursor = 0
        active_reach_cursor = 0

        active_node = start_node

        while True:
            upstream = nodes[nodes[:,0]==active_node]
            active_reach[active_reach_cursor] = active_node
            active_times[active_reach_cursor] = time_dict.get(active_node, 0.0)
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
                all_times[path_cursor] = np.cumsum(active_times) * (all_paths[path_cursor] > 0)
                queue_cursor -= 1
                path_cursor += 1
                last_node, active_node = queue[queue_cursor]
                if not any((last_node, active_node)):
                    break
                active_reach_cursor = np.where(active_reach == last_node)[0][0] + 1
                start_cursor = active_reach_cursor
                active_reach[active_reach_cursor:] = 0.
                active_times[active_reach_cursor:] = 0.

    return all_paths[:path_cursor, :longest_path + 1], all_times[:path_cursor, :longest_path + 1]


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
    sever = (((sc_nodes[:,1] == True) & (sc_nodes[:, 0] == False)) |   # Streamcalc 0 -> Streamcalc 1
             ((div_nodes[:,1] == False) & (div_nodes[:, 0] == True)))  # Divergence 0 -> Divergence 2
    active_nodes = nodes[~sever]

    return active_nodes, conversion_dict


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


def get_times(all_nodes, q_table, vaa_table, conversion_dict):

    # Read in velocity and length data into lookup tables
    velocities = tuple(zip(*read_dbf(q_table, ["COMID", "V0001E"])))
    lengths = tuple(zip(*read_dbf(vaa_table, ["ComID", "LengthKM"])))

    velocity_dict = dict(zip(map(conversion_dict.get, velocities[0]), velocities[1]))
    length_dict = dict(zip(map(conversion_dict.get, lengths[0]), lengths[1]))

    velocity_dict.pop(None, None)
    length_dict.pop(None, None)

    # Create lookup dictionary with times for each reach
    time_dict = {}
    all_reaches = {reach for node in all_nodes for reach in node if reach}
    for reach in all_reaches:  # loop through all reaches
        length = length_dict.get(reach, 0) * 1000.  # km -> m
        velocity = velocity_dict.get(reach, 0) * 0.3048  # ft/sec -> m/sec
        if length and velocity:
            time_dict[reach] = (length / velocity) / 86400.  # seconds -> days

    return time_dict


def write_outfile(outfile, paths, times, path_map, conversion_array):
    np.savez_compressed(outfile, paths=paths, times=times, path_map=path_map, conversion_array=conversion_array)


def dict_to_array(cd):
    max_alias = max(cd.values())
    out_array = np.zeros(max_alias + 1)
    for comid, alias in cd.items():
        out_array[alias] = comid
    return out_array


def create_upstream_paths(region, region_dir, output_directory, region_filter='all'):

    print("Processing region {}...".format(region))

    # Designate paths
    flow_table = os.path.join(region_dir, "NHDPlusAttributes", "PlusFlow.dbf")
    flow_lines = os.path.join(region_dir, "NHDSnapshot", "Hydrography", "NHDFlowline.dbf")
    vaa_table = os.path.join(region_dir, "NHDPlusAttributes", "PlusFlowlineVAA.dbf")
    q_table = os.path.join(region_dir, "EROMExtension", "EROM_MA0001.dbf")

    outfile = os.path.join(output_directory, "upstream_{}.npz".format(region))

    if not os.path.exists(outfile):
        # Do the work
        all_nodes, conversion_dict = get_nodes(flow_table, vaa_table, flow_lines)

        time_dict = get_times(all_nodes, q_table, vaa_table, conversion_dict)

        outlets = get_outlets(flow_table, flow_lines, vaa_table, conversion_dict)

        paths, times = snapshot(all_nodes, outlets, time_dict, cd=conversion_dict)

        path_map = map_paths(paths)

        conversion_array = dict_to_array(conversion_dict)

        write_outfile(outfile, paths, times, path_map, conversion_array)


def main():
    output_directory = r"T:\CA3T\bin\UpstreamPaths"
    if not os.path.isdir(output_directory): os.mkdir(output_directory)
    nhd_dir = r"T:\NationalData\NHDPlusV2"
    region_filter = 'all'
    for region, region_dir in get_nhd(nhd_dir).items():
        if region == region_filter or region_filter == 'all':
            create_upstream_paths(region, region_dir, output_directory)


def gis_it(comids, field="COMID", lookup_dict=None):

    if lookup_dict:
        comids = filter(None, map(lookup_dict.get, comids))

    print("\"{}\" = ".format(field) + " OR \"{}\" = ".format(field).join(map(str, comids)))

if __name__ == "__main__":
    time_it = False
    if time_it:
        import cProfile
        cProfile.run('main()')
    else:
        main()
