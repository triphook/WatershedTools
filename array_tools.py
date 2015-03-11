__author__ = 'Trip Hook'

import numpy as np
import gdal, ogr
import sys
import collections
import os
import csv

# Object representing a simple bounding rectangle
class Envelope(object):
    def __init__(self, left, right, bottom, top):
        self.left = left
        self.right = right
        self.bottom = bottom
        self.top = top

    def overlap(self, r2):
        def range_overlap(a_min, a_max, b_min, b_max):
            return (a_min <= b_max) and (b_min <= a_max)
        if not all((range_overlap(self.left, self.right, r2.left, r2.right),
                    range_overlap(self.bottom, self.top, r2.bottom, r2.top))):
            return None
        else:
            left, right = sorted([self.left, self.right, r2.left, r2.right])[1:3]
            bottom, top = sorted([self.bottom, self.top, r2.bottom, r2.top])[1:3]
        return Envelope(left, right, bottom, top)

    @property
    def area(self):
        return abs(self.top - self.bottom) * abs(self.right - self.left)

    def __repr__(self):
        return "Rectangle(left: {}, right: {}, top: {}, bottom: {}".format(self.left, self.right, self.top, self.bottom)

    def __eq__(self, other):
        if (self.left, self.right, self.bottom, self.top) == (other.left, other.right, other.bottom, other.top):
            return True
        return False


# Wrapper for an ESRI raster grid and function for reading into an array
class Raster(object):
    def __init__(self, path, precision=None):
        self.path = path
        self.obj = gdal.Open(path)
        gt = self.obj.GetGeoTransform()
        self.cell_size = gt[1]
        self.tl = (gt[0], gt[3])  # top left
        self.size = (self.obj.RasterXSize * self.cell_size, self.obj.RasterYSize * self.cell_size)
        self.shape = Envelope(self.tl[0], self.tl[0] + self.size[0], self.tl[1] - self.size[1], self.tl[1])
        self.max_val = int(self.obj.GetRasterBand(1).GetMaximum())
        if not precision:
            precision = len(str(int(self.max_val)))
        self.precision = 10 ** precision
        self.envelope = None
        self._array = np.array([])

    def array(self, envelope, zero_min=True):
        if not self._array.size or envelope != self.envelope:
            offset_x = (envelope.left - self.shape.left)
            offset_y = (self.shape.top - envelope.top)
            x_max = (envelope.right - envelope.left)
            y_max = (envelope.top - envelope.bottom)
            bounds = map(lambda x: int(x / self.cell_size), (offset_x, offset_y, x_max, y_max))
            self._array = self.obj.ReadAsArray(*bounds)
            if zero_min:
                self._array[self._array < 0] = 0
            if self.precision > 10000:
                self._array = np.int64(self._array)
            self.envelope = envelope
        return self._array


# Object representing a stream reach.  Used to construct a hierarchy for accumulation function
class Reach(object):
    def __init__(self):
        self.name = None
        self.parent = None
        self.children = set()
        self._upstream = set()
        self._paths = []

        self.length = 0  # meters
        self.velocity = 0  # meters
        self._time = -1  # seconds

    def __repr__(self):
        return "Reach({})".format(self.name)

    @property
    def time(self):
        if self._time != -1:
            return self._time
        else:
            if self.length and self.velocity:
                self._time = self.length / self.velocity
            else:
                self._time = 0
            return self._time

    # Recursively search upstream in the drainage network, returns a set of all upstream reaches
    @property
    def upstream(self):
        if self._upstream:
            return self._upstream
        else:
            for child in self.children:
                self._upstream.add(child)
                self._upstream |= child.upstream
            return self._upstream

    # Recursively search upstream in the drainage network, returns a set of all paths upstream
    def upstream_paths(self, h):
        def look_upstream(n):
            if not h[n].children:
                self._paths.append(current_path[:])
            for child in h[n].children:
                current_path.append(child)
                look_upstream(child.name)
            if current_path:
                current_path.pop()
        if self._paths:
            return self._paths
        else:
            current_path = [self]
            look_upstream(self.name)
            return self._paths

    def time_maps(self, h):
        paths = self.upstream_paths(h)
        time_maps = [zip(path, np.cumsum([reach.time for reach in path])) for path in paths]
        items = {s for time_map in time_maps for s in time_map}
        return items


# Accumulates class areas for all reaches upstream of each reach
def accumulate(allocated_data, flow_table, translate=None, excl_params=None, excl_field=None, excl_vals=None):
    sys.setrecursionlimit(50000)
    reaches = build_hierarchy(flow_table, excl_params, excl_field, excl_vals, translate)
    accumulated = collections.defaultdict(dict)
    for reach in iter(reaches):
        us_reaches = set([reach.name] + [r.name for r in reach.upstream])
        accumulated[reach.name] = sum([allocated_data.get(r, 0.0) for r in us_reaches])
    return accumulated


# Allocates raster classes to a set of overlapping zones
def allocate(allocation_raster, zone_rasters, tile='max'):

    # Accept single raster as input
    if not type(zone_rasters) in (list, set):
        zone_rasters = [zone_rasters]

    # Overlap rasters and create envelope covering common areas
    overlap_area = allocation_raster.shape
    for raster in zone_rasters:
        overlap_area = raster.shape.overlap(overlap_area)
        if not overlap_area:
            sys.exit("Zone and allocation rasters do not overlap")

    finished = np.array([])
    tiles = make_tiles(overlap_area, tile)
    for tile in tiles:
        # Generate combined array and set precision adjustments
        all_rasters = sorted(zone_rasters + [allocation_raster], key=lambda x: x.precision, reverse=True)
        combined_array = np.int64(allocation_raster.array(tile))
        for i, raster in enumerate(zone_rasters):
            raster.adjust = np.prod([r.precision for r in all_rasters[i + 1:]])
            combined_array += np.int64(raster.array(tile) * raster.adjust)
        allocation_raster.adjust = 1

        # Break down combined array to get zones and classes
        zones, counts = np.unique(combined_array.flat, return_counts=True)  # Cells in each zone and class
        zones = np.array([(zones / r.adjust) - ((zones / r.adjust) / r.precision) * r.precision for r in all_rasters])
        counts *= (allocation_raster.cell_size ** 2)  # Convert cell count to area
        final = np.vstack((zones, counts))[:, (zones[0] > 0) & (zones[1] > 0)]
        if finished.size:
            finished = np.hstack((finished, final))
        else:
            finished = final

    return finished


def build_hierarchy(flow_table, from_field="FROMCOMID", to_field="TOCOMID", id_field="COMID",
                    excl_params=None, excl_field=None, excl_vals=None, translate=None):

    # Read from/to relationships in flow table
    flows = read_dbf(flow_table, [from_field, to_field])

    # Exclude flows that meet certain criteria (e.g., exclude streams where FCODE is 56600 (coastline)
    if all((excl_params, excl_field, excl_vals)):
        params = read_dbf(excl_params, [id_field, excl_field])
        # This bit is necessary because some NHD plus regions have different capitalization patterns
        if not params:
            params = dict(read_dbf(excl_params, ['ComID', "FCode"]))
        params = dict(params)
        flows = [flow for flow in flows if not set(excl_vals) & set(map(lambda x: params.get(x), flow))]

    # Translate COMID into another field for processing (usually GRIDCODE for raster catchments)
    if translate:
        flows = [map(lambda x: translate.get(x), row) for row in flows]

    # Build hierarchy
    hierarchy = collections.defaultdict(lambda: Reach())
    for reach_id, parent in flows:
        if reach_id:
            hierarchy[reach_id].name = reach_id
            hierarchy[parent].name = parent
            hierarchy[reach_id].parent = hierarchy[parent]
            hierarchy[parent].children.add(hierarchy[reach_id])

    return filter(lambda reach: bool(reach.name), hierarchy.values())


# Divides a bounding envelope into smaller tiles for processing on less powerful computers
def make_tiles(envelope, tile_size):
    if tile_size == 'max':
        return [envelope]
    else:
        h = range(int(envelope.left), int(envelope.right), tile_size) + [envelope.right]
        v = range(int(envelope.bottom), int(envelope.top), tile_size) + [envelope.top]
        return [Envelope(h[i], h[i + 1], v[j], v[j + 1]) for i in range(len(h) - 1) for j in range(len(v) - 1)]


# Returns the average of a sequence
def mean(iterable):
    if iterable:
        return float(sum(iterable) / len(iterable))
    else:
        return 0.0


# Allows the initialization of a nested dictionary
def nested_dict():
    return collections.defaultdict(nested_dict)


# Reads the contents of a dbf table
def read_dbf(dbf_file, out_fields=None):
    # Initialize file
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(dbf_file)
    layer = data_source.GetLayer(0)

    # Set fields
    ld = layer.GetLayerDefn()
    fields = {ld.GetFieldDefn(i).GetName() for i in range(ld.GetFieldCount())}
    if out_fields:
        missing_fields = set(out_fields) - fields
        if missing_fields:
            print "Specified field(s) {} not found in dbf file {}".format(",".join(missing_fields), dbf_file)
            return None
        fields = [field for field in out_fields if not field in missing_fields]
    # Read data
    if len(fields) > 1:
        table = [[row.GetField(f) for f in fields] for row in layer]
    else:
        table = [row.GetField(list(fields)[0]) for row in layer]
    # noinspection PyUnusedLocal
    data_source = None
    return table


# Reads a dbf file and returns dictionary with one field as a key and another as a value
def read_translation(tf, from_field, to_field, flip=False):
    tt = read_dbf(tf, [from_field, to_field])
    comid_to_gridcode = dict(tt)
    if flip:
        gridcode_to_comid = dict([(row[1], row[0]) for row in tt])
        return comid_to_gridcode, gridcode_to_comid
    else:
        return comid_to_gridcode


# Breaks down an allocation array into a nested dictionary
def make_histogram(output_array):
    output_dict = nested_dict()
    for row in output_array.T:
        row = list(row)
        top = row.pop(0)
        active_d = row.pop()
        for i in range(len(row)):
            active_d = {row.pop(): active_d}
        output_dict[top].update(active_d)
    return output_dict


# Sums all values for each group
def sum_by_group(groups, values, oper='sum'):
    order = np.argsort(groups)
    groups = groups[order]
    values = values[order]
    if oper == 'sum':
        values.cumsum(out=values)
    elif oper == 'prod':
        values.cumprod(out=values)
    index = np.ones(len(groups), 'bool')
    index[:-1] = groups[1:] != groups[:-1]
    values = values[index]
    groups = groups[index]
    values[1:] = values[1:] - values[:-1]
    return np.vstack((groups, values))


# Write dictionary or NumPy array to csv file
def write_to_file(output_dict, outfile_name, translate=None, id_field='COMID'):

    def is_num(x):
        try:
            float(x)
            return True
        except:
            return False

    def format_keys(r):
        return map(lambda x: "class_{}".format(x) if is_num(x) else x, r)

    if not os.path.exists(os.path.dirname(outfile_name)):
        os.mkdir(os.path.dirname(outfile_name))

    with open(outfile_name, 'wb') as f:
        classes = set([_class for row in output_dict.values() for _class in row.keys()])  # all classes in allocation
        header = format_keys([id_field] + sorted(classes))
        writer = csv.DictWriter(f, header)
        writer.writeheader()
        for comid, row in output_dict.iteritems():
            row = dict(zip(format_keys(row.keys()), row.values()))
            comid = translate.get(comid, -comid) if translate else comid
            row[id_field] = comid
            if row:
                writer.writerow(row)