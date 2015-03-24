import numpy as np
import gdal, ogr
import sys
import collections
import os
import csv

# Object representing a simple bounding rectangle, used primarily to measure raster overlap
class Envelope(object):
    def __init__(self, left, right, bottom, top):
        self.left = left
        self.right = right
        self.bottom = bottom
        self.top = top

    #  Returns the rectangle corresponding to the overlap with another Envelope object
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


# A modified dictionary which holds a set of "reach" objects representing NHD Plus stream reaches
# indexed by COMID or another user designated "id_field"
class Hierarchy(dict):

    def __init__(self, region_dir, region,
                 init_table='PlusFlow', from_field="FROMCOMID", to_field="TOCOMID", id_field="ComID",
                 exclusions=[("PlusFlowlineVAA", "Fcode", 56600)],
                 velocity_table='EROM_Annual', velocity_field='V0001E'):  # (table, field, values)

        self.region = region
        self.dir = region_dir
        self.from_field = from_field
        self.to_field = to_field
        self.id_field = id_field
        self.velocity_table = velocity_table
        self.velocity_field = velocity_field
        self.init_table = init_table
        self.exclusions = set()

        # Initialize paths to NHD Plus files
        nhd_files = [('CatchmentRaster', r"NHDPlusCatchment\cat"),
                     ('CatchmentTable', r"NHDPlusCatchment\featureidgridcode.dbf"),
                     ('PlusFlowlineVAA', r"NHDPlusAttributes\PlusFlowlineVAA.dbf"),
                     ('Flowline', r"NHDSnapshot\Hydrography\NHDFlowline.shp"),
                     ('PlusFlow', r"NHDPlusAttributes\PlusFlow.dbf"),
                     ('EROM_Annual', r"EROMExtension\EROM_MA0001.dbf")]
        nhd_files += [('EROM_{}'.format(m), r"EROMExtension\EROM_MA000{}.dbf".format(m)) for m in range(1, 13)]
        self.paths = collections.defaultdict(lambda: "No path specified")
        not_found = []
        for name, path in nhd_files:
            full_path = os.path.join(self.dir, path)
            if os.path.exists(full_path):
                self.paths[name] = full_path
            else:
                not_found.append(name)
        if not_found:
            print "Files not found for: " + ", ".join(not_found)

        self.initialize()
        self.exclude(exclusions)


    def __missing__(self, key):
        self[key] = Reach()
        return self[key]

    def initialize(self):
        # Read length and velocity attributes into dictionaries
        length_dict = dict(read_dbf(self.paths['PlusFlowlineVAA'], [self.id_field, "LengthKM"]))

        velocity_dict = {}
        if self.velocity_table:
            velocity_dict = dict(read_dbf(self.paths[self.velocity_table], [self.id_field, self.velocity_field]))
        for reach_id, parent in read_dbf(self.paths[self.init_table], [self.from_field, self.to_field]):
            if all((reach_id, parent, not reach_id in self.exclusions, not parent in self.exclusions)):
                self[reach_id].name = reach_id
                self[parent].name = parent
                self[reach_id].parent = self[parent]
                self[parent].children.add(self[reach_id])
                self[reach_id].length = length_dict.get(reach_id, 0) * 1000  # meters
                self[reach_id].velocity = velocity_dict.get(reach_id, -1) * 0.3048  # meters per second

    def exclude(self, exclusions):
        remove_lines = []
        for table, field, value in sorted(exclusions):
            params = dict(read_dbf(self.paths[table], [self.id_field, field]))
            remove_lines += filter(lambda x: params.get(x) == value, self.keys())

        for line in remove_lines:
            del self[line]


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

        # Recursively search upstream in the drainage network, returns a set of all upstream reaches
        @property
        def all_upstream(self):
            if not self._upstream:
                map(lambda x: self._upstream.update(x.all_upstream | {x}), self.children)
                for child in self.children:
                    self._upstream.add(child)
                    self._upstream |= child.all_upstream
            return self._upstream

        @property
        def time(self):
            if self._time == -1 and all((self.length, self.velocity)):
                self._time = self.length / self.velocity
            return self._time

        # Recursively search upstream in the drainage network, returns a set of all paths upstream
        @property
        def upstream_paths(self):
            if not self._paths:
                for child in self.children:
                    if child.upstream_paths:
                        self._paths.extend([child] + path for path in child.upstream_paths)
                    else:
                        self._paths.append([child])
            return self._paths


# Accumulates class areas for all reaches upstream of each reach
def accumulate(allocated_data, hierarchy, translate=True):
    sys.setrecursionlimit(50000)

    # Translate from GRIDCODE to COMID for NHD Plus catchment allocations
    if translate:
        translation_dict = dict(read_dbf(hierarchy.paths['CatchmentTable'], ["GRIDCODE", "FEATUREID"]))
        translated_comids = map(lambda x: translation_dict.get(x, -1), allocated_data[0])
        allocated_data = np.vstack([translated_comids, allocated_data[1:]])

    # Accumulate data
    accumulated_data = np.empty([1, 3])  # Initialize blank array
    all_classes = np.unique(allocated_data[1])  # Get a set of all unique subclasses
    for reach_id, reach in hierarchy.iteritems():  # Iterate through all reaches in hierarchy
        us_reaches = np.array(list(set([reach_id] + [r.name for r in reach.all_upstream])))
        for _class in all_classes:
            in_class = allocated_data[:,(allocated_data[1] == _class)]
            upstream = in_class[2:, np.in1d(in_class[0], us_reaches)]
            accumulated_data = np.vstack((accumulated_data, np.array([reach_id, _class, upstream.sum()])))
    return accumulated_data.T


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


# Assembles a dictionary of NHD Plus directory structure indexed by region
def get_nhd(nhd_dir, region_filter='all'):
    # Get catchment grid, gridcode to comid translation files, and flow tables
    regions = {'01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09',
               '10U', '10L', '11', '12', '13', '14', '15', '16', '17', '18'}
    all_paths = collections.defaultdict()
    region_dirs = {"NHDPlus{}".format(region) for region in regions}
    for root_dir, sub_dirs, _ in os.walk(nhd_dir):
        if set(sub_dirs) & region_dirs:
            for sub_dir in sub_dirs:
                all_paths[sub_dir.lstrip("NHDPlus")] = os.path.join(root_dir, sub_dir)
    return all_paths


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


# Returns the average of a sequence
def mean(iterable):
    if iterable:
        return float(sum(iterable) / len(iterable))
    else:
        return 0.0


# Allows the initialization of a nested dictionary
def nested_dict():
    return collections.defaultdict(nested_dict)


# Divides a bounding envelope into smaller tiles for processing on less powerful computers
def make_tiles(envelope, tile_size):
    if tile_size == 'max':
        return [envelope]
    else:
        h = range(int(envelope.left), int(envelope.right), tile_size) + [envelope.right]
        v = range(int(envelope.bottom), int(envelope.top), tile_size) + [envelope.top]
        return [Envelope(h[i], h[i + 1], v[j], v[j + 1]) for i in range(len(h) - 1) for j in range(len(v) - 1)]


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
        remedies = [lambda x: x.upper(), lambda x: x.capitalize()]
        while missing_fields and remedies:
            new_fields = map(remedies.pop(0), out_fields)
            out_fields = [nf if nf in fields else out_fields[i] for i, nf in enumerate(new_fields)]
            missing_fields = set(out_fields) - fields
        if missing_fields:
            print "Fields {} not found in {}".format(out_fields, dbf_file)
        fields = [field for field in out_fields if not field in missing_fields]

    # Read data
    if len(fields) > 1:
        table = [[row.GetField(f) for f in fields] for row in layer]
    else:
        table = [row.GetField(list(fields)[0]) for row in layer]
    # noinspection PyUnusedLocal
    data_source = None
    return table


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
def write_to_file(output_dict, outfile_name, id_field='COMID'):

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
            row[id_field] = comid
            if row:
                writer.writerow(row)