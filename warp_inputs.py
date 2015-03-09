import numpy as np
import gdal, ogr
import os
import sys
import collections
import csv
import urllib2
import cPickle


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
        self._array = "null"

    def array(self, envelope, zero_min=True):
        if self._array == "null" or envelope != self.envelope:
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

    def __repr__(self):
        return "Reach({})".format(self.name)

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


# Accumulates class areas for all reaches upstream of each reach
def accumulate(allocated_data, flow_table, translate=None, excl_params=None, excl_field=None, excl_vals=None):

    def build_hierarchy():
        # Read from/to relationships in flow table
        flows = read_dbf(flow_table, ['FROMCOMID', 'TOCOMID'])

        # Exclude flows that meet certain criteria (e.g., exclude streams where FCODE is 56600 (coastline)
        if all((excl_params, excl_field, excl_vals)):
            print excl_params
            params = read_dbf(excl_params, ['COMID', excl_field])
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

    def accumulate_data():
        accumulated = collections.defaultdict(dict)
        for reach in iter(reaches):
            us_reaches = set([reach.name] + [r.name for r in reach.upstream])
            accumulated[reach.name] = sum([allocated_data.get(r, 0.0) for r in us_reaches])
        return accumulated

    sys.setrecursionlimit(50000)
    reaches = build_hierarchy()
    out_array = accumulate_data()
    return out_array


# Allocates raster classes to a set of overlapping zones 
def allocate(allocation_raster, zone_rasters, tile='max'):

    # Overlap rasters and create envelope covering common areas
    overlap_area = allocation_raster.shape
    for raster in zone_rasters:
        overlap_area = raster.shape.overlap(overlap_area)
        if not overlap_area:
            sys.exit("Zone and allocation rasters do not overlap")

    for envelope in make_tiles(overlap_area, tile):
        # Generate combined array and set precision adjustments
        all_rasters = sorted(zone_rasters + [allocation_raster], key=lambda r: r.precision, reverse=True)
        combined_array = np.int64(allocation_raster.array(envelope))
        for i, raster in enumerate(zone_rasters):
            raster.adjust = np.prod([r.precision for r in all_rasters[i + 1:]])
            combined_array += np.int64(raster.array(envelope) * raster.adjust)
        allocation_raster.adjust = 1

        # Break down combined array to get zones and classes
        zones, counts = np.unique(combined_array.flat, return_counts=True)  # Cells in each zone and class
        zones = np.array([(zones / r.adjust) - ((zones / r.adjust) / r.precision) * r.precision for r in all_rasters])
        counts *= (allocation_raster.cell_size ** 2)  # Convert cell count to area
        finished = np.vstack((zones, counts))[:, (zones[0] > 0) & (zones[1] > 0)]
        return finished


# Divides a bounding envelope into smaller tiles for processing on less powerful computers
def make_tiles(envelope, tile_size):
    if tile_size == 'max':
        return [envelope]
    else:
        h = range(int(envelope.left), int(envelope.right), tile_size) + [envelope.right]
        v = range(int(envelope.bottom), int(envelope.top), tile_size) + [envelope.top]
        return [Envelope(h[i], h[i + 1], v[j + 1], v[j]) for i in range(len(h) - 1) for j in range(len(v) - 1)]


# Returns the average of a sequence
def mean(iterable):
    if iterable:
        return float(sum(iterable) / len(iterable))
    else:
        return 0.0


# Allows the initialization of a nested dictionary
def nested_dict():
    return collections.defaultdict(dict)


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


# Sums all values for each group
def sum_by_group(groups, values):
    order = np.argsort(groups)
    groups = groups[order]
    values = values[order]
    values.cumsum(out=values)
    index = np.ones(len(groups), 'bool')
    index[:-1] = groups[1:] != groups[:-1]
    values = values[index]
    groups = groups[index]
    values[1:] = values[1:] - values[:-1]
    return np.vstack((groups, values))


# Gets pesticide use data from usgs.gov
def pull_use_data(use_pickle, year='all'):  # Pull use data from usgs.gov and store as a serialized dictionary
    if not os.path.exists(use_pickle):
        low_range, high_range = range(8, 15), range(1, 8)
        pairs = [('high', i) for i in high_range] + [('low', i) for i in low_range]
        use_data = collections.defaultdict(nested_dict)
        for level, table_num in pairs:
            url = r"http://pubs.usgs.gov/ds/752/EPest.{}.county.estimates.table{}.txt".format(level, table_num)
            print "Reading {}...".format(url)
            lines = urllib2.urlopen(url).read().split("\n")
            header = lines.next()
            for line in lines:
                if year == 'all' or str(year) in line:
                    vals = line.split()
                    if len(vals) != len(header):
                        vals = [-1 for _ in header]
                    l = dict(zip(header, vals))
                    fips = str(l['STATE_FIPS_CODE']).zfill(2) + str(l['COUNTY_FIPS_CODE']).zfill(3)
                    use_data[l['COMPOUND']][fips][level] = l['KG']
        with open(use_pickle, 'wb') as f:
            cPickle.dump(use_data, f)
    else:
        with open(use_pickle) as f:
            use_data = cPickle.load(f)
    return use_data


# Writes accumulated data to a csv file
def write_to_file(output_dict, outfile_name, translate=None, id_field='COMID'):

    if not os.path.exists(os.path.dirname(outfile_name)):
        os.mkdir(os.path.dirname(outfile_name))
    with open(outfile_name, 'wb') as f:
        header = [id_field, "LOAD"]
        writer = csv.DictWriter(f, header)
        writer.writeheader()
        for comid, load in output_dict.iteritems():
            if translate:
                comid = translate.get(comid)
            row = {id_field: comid, "LOAD": load}
            writer.writerow(row)


overlaps = {'01': {"ME", "NH", "VT", "MA", "CT", "RI", "NY"},
            '02': {"VT", "NY", "PA", "NJ", "MD", "DE", "WV", "DC", "VA"},
            '03N': {"VA", "NC", "SC", "GA"},
            '03S': {"FL", "GA"},
            '03W': {"FL", "GA", "TN", "AL", "MS"},
            '04': {"WI", "MN", "MI", "IL", "IN", "OH", "PA", "NY"},
            '05': {"IL", "IN", "OH", "PA", "WV", "VA", "KY", "TN"},
            '06': {"VA", "KY", "TN", "NC", "GA", "AL", "MS"},
            '07': {"MN", "WI", "SD", "IA", "IL", "MO"},
            '08': {"MO", "KY", "TN", "AR", "MS", "LA"},
            '09': {"ND", "MN", "SD"},
            '10U': {"MT", "ND", "WY", "SD", "MN", "NE", "IA"},
            '10L': {"CO", "WY", "MN", "NE", "IA", "KS", "MO"},
            '11': {"CO", "KS", "MO", "NM", "TX", "OK", "AR", "LA"},
            '12': {"NM", "TX", "LA"},
            '13': {"CO", "NM", "TX"},
            '14': {"WY", "UT", "CO", "AZ", "NM"},
            '15': {"NV", "UT", "AZ", "NM", "CA"},
            '16': {"CA", "OR", "ID", "WY", "NV", "UT"},
            '17': {"WA", "ID", "MT", "OR", "WY", "UT", "NV"},
            '18': {"OR", "NV", "CA"}}

if __name__ == "__main__":
    print "This is just a library. Run another script"