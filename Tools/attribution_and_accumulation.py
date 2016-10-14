import gdal
import numpy as np
import os


class Envelope(object):
    # Object representing a simple bounding rectangle, used primarily to measure raster overlap

    def __init__(self, left, right, bottom, top):
        self.left, self.right, self.bottom, self.top = map(float, (left, right, bottom, top))

    # Returns the rectangle corresponding to the overlap with another Envelope object
    def overlap(self, r2):
        if all(((self.left <= r2.right), (r2.left <= self.right), (self.bottom <= r2.top), (r2.bottom <= self.top))):
            left, right = sorted([self.left, self.right, r2.left, r2.right])[1:3]
            bottom, top = sorted([self.bottom, self.top, r2.bottom, r2.top])[1:3]
            return Envelope(left, right, bottom, top)
        else:
            print("Rasters do not overlap: {}\nAllocation: {}".format(self.shape, r2))

    def tiles(self, tile_size='max', progress=True):
        if tile_size == 'max':
            yield self
        else:
            x = list(range(int(self.left), int(self.right), tile_size)) + [self.right]
            y = list(range(int(self.bottom), int(self.top), tile_size)) + [self.top]
            total_tiles = (len(x) - 1) * (len(y) - 1)
            counter = iter(range(total_tiles))
            for i in range(len(x) - 1):
                for j in range(len(y) - 1):
                    if progress and total_tiles > 1:
                        print("Processing tile {} of {}".format(next(counter) + 1, total_tiles))
                    yield Envelope(*map(float, (x[i], x[i + 1], y[j], y[j + 1])))

    @property
    def area(self):
        return abs(self.top - self.bottom) * abs(self.right - self.left)

    def __eq__(self, other):
        return (self.left, self.right, self.bottom, self.top) == (other.left, other.right, other.bottom, other.top)


class Navigator(object):
    def __init__(self, paths_file):
        data = np.load(paths_file)
        self.paths, self.path_map, self.start_cols, self.alias_to_comid = \
            data['paths'], data['path_map'], data['start_cols'], data['conversion_array']
        self.n_reaches = self.alias_to_comid.size
        self.comid_to_alias = dict(zip(self.alias_to_comid, np.arange(self.n_reaches)))

        a = (self.start_cols == 0)
        self.last_outlet = np.where(a)[0][a.cumsum() - 1]

    def all_upstream(self, reach, mode='comid'):
        reach = self._format_input(reach, mode)
        start_row, end_row, col = map(int, self.path_map[reach])
        start_col = list(self.paths[start_row]).index(reach)
        upstream_reaches = list(self.paths[start_row:end_row])
        upstream_reaches.append(self.paths[start_row][start_col:])
        output = np.concatenate(upstream_reaches)
        return self._format_output(output, mode)

    def all_downstream(self, reach, mode='comid'):
        reach = self._format_input(reach, mode)
        start_row, _, _ = map(int, self.path_map[reach])
        last_outlet = self.last_outlet[start_row]
        start_cols = self.start_cols[last_outlet:start_row + 1]
        active_paths = np.where(start_cols == np.minimum.accumulate(start_cols[::-1])[::-1])[0] + last_outlet
        output = np.zeros((start_cols[-1] + len(self.paths[start_row])) * 1.5, dtype=np.int32)
        for i, start, end in zip(active_paths, self.start_cols[active_paths], self.start_cols[active_paths][1:]):
            output[start:end] = self.paths[i][:end - start]
        return self._format_output(output[output > 0], mode)

    def upstream_paths(self, reach, mode='comid'):
        reach = self._format_input(reach, mode)
        start_row, end_row, _ = map(int, self.path_map[reach])
        path = np.zeros(max(map(lambda x: self.start_cols[x] + self.paths[x].size, range(start_row, end_row))),
                        dtype=np.int32)
        baseline = self.start_cols[start_row]
        for i in range(start_row, end_row - 1):
            stub = self.paths[i]
            start_col = self.start_cols[i] - baseline
            path[start_col:] = 0
            path[start_col:start_col + stub.size] = stub
            output = path[path != 0]
            yield self._format_output(output, mode)

    def _format_input(self, reach_id, mode):
        return reach_id if mode == 'alias' else self.comid_to_alias[reach_id]

    def _format_output(self, output, mode):
        return output if mode == 'alias' else self.alias_to_comid[output]


class Raster(object):
    def __init__(self, path, no_data=255):
        self.no_data = no_data
        self.obj = gdal.Open(path)
        left, self.cell_size, _, top, *_ = self.obj.GetGeoTransform()
        x_size, y_size = np.array([self.obj.RasterXSize, self.obj.RasterYSize]) * self.cell_size
        self.shape = Envelope(left, left + x_size, top - y_size, top)
        self.max_val = int(self.obj.GetRasterBand(1).ComputeRasterMinMax(1)[1])
        self.precision = 10 ** (int(np.log10(self.max_val)) + 1)
        self.values = np.where(np.array(self.obj.GetRasterBand(1).GetHistogram()) > 0)[0]
        self._array = np.array([])
        self.n_classes = len(self.values)
        self.array_envelope = self.shape

    def array(self, envelope=None, zero_min=True):
        if not self._array.size or envelope != self.array_envelope:
            offset_x, offset_y = (envelope.left - self.shape.left), (self.shape.top - envelope.top)
            x_max, y_max = (envelope.right - envelope.left), (envelope.top - envelope.bottom)
            bounds = map(lambda x: int(x / self.cell_size), (offset_x, offset_y, x_max, y_max))
            self._array = self.obj.ReadAsArray(*bounds)
            self._array[self._array == self.no_data] = 0
            if zero_min:
                self._array[self._array < 0] = 0
            if self.precision > 10000:
                self._array = np.int64(self._array)
            self.array_envelope = envelope
        return self._array


def allocate(allocation_raster, zone_raster, region, tile_size='max'):
    overlap_area = allocation_raster.shape.overlap(zone_raster.shape)
    assert overlap_area

    # Perform allocation
    allocation = np.zeros(allocation_raster.precision * zone_raster.precision, dtype=np.int64)
    for tile in overlap_area.tiles(tile_size):
        combined = np.int64(allocation_raster.array(tile)) + np.int64(
            zone_raster.array(tile) * np.int64(allocation_raster.precision))
        allocation += np.int32(np.bincount(combined.flat, minlength=allocation.size))
    allocation *= (int(allocation_raster.cell_size) ** 2)  # Convert cell count to area

    # Rearrange allocation to be indexed by alias instead of gridcode
    allocation = np.reshape(allocation, (zone_raster.precision, allocation_raster.precision))
    gridcodes = np.where(allocation.sum(axis=1) > 0)[0]
    aliases = region.gridcode_to_alias[gridcodes]
    new_allocation = np.zeros((region.n_reaches + 1, allocation.shape[1]))
    new_allocation[aliases] = allocation[gridcodes]

    return new_allocation


def histogram(allocation, watershed, progress=True):
    upstream = np.zeros(allocation.shape)
    for reach in range(1, watershed.path_map.shape[0]):
        if progress and not reach % 10000:
            print("Accumulating reach {} of {}".format(reach, watershed.path_map.shape[0]))
        upstream_reaches = watershed.all_upstream(reach)
        if upstream_reaches.any():
            upstream[reach] = np.take(allocation, upstream_reaches, axis=0).sum(axis=0)
    return upstream


def write_to_file(results, out_file, zone_ids, out_format="csv"):
    print("Writing to {}...".format(out_file))
    classes = np.arange(results.shape[1])
    if out_format == "npz":
        np.savez_compressed(out_file, results=results, header=zone_ids)
    elif out_format == "csv":
        results = np.vstack([zone_ids, results[1:].T]).T
        np.savetxt(out_file, results, delimiter=',', header="Zone/Class," + ",".join(map(str, classes)))
    else:
        print("Invalid output format {}".format(out_format))


def main():
    # Specify paths here
    allocation_raster = r"T:\NationalData\NLCD_2011\nlcd_2011_landcover_2011_edition_2014_10_10.img"
    nhd_dir = r"T:\NationalData\NHDPlusV2"

    for region_id in ['01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09',
                      '10U', '10L', '11', '12', '13', '14', '15', '16', '17', '18']:
        region_dir = os.path.join(nhd_dir, r"NHDPlus{}".format(region_id))
        output_file_root = r"..\Output\region_{}".format(region_id)

        catchment_raster = os.path.join(region_dir, "NHDPlusCatchment", "cat")
        topology_file = r"..\Preprocessed\WatershedTopology\upstream_{}.npz".format(region_id)

        region = Navigator(topology_file)
        allocation = allocate(Raster(allocation_raster), Raster(catchment_raster), region, tile_size=300000)
        accumulation = histogram(allocation, region)
        write_to_file(allocation, output_file_root + "_alloc.csv", region.alias_to_reach)
        write_to_file(accumulation, output_file_root + "_accum.csv", region.alias_to_reach)


if __name__ == "__main__":
    time_it = True
    if time_it:
        import cProfile

        cProfile.run("main()")
    else:
        main()
