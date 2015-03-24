from array_tools import *
import os
import time

def ca3t(nhd_dir, allocation_file, output_dir, overwrite):
    nhd_files = get_nhd(nhd_dir)
    basename = os.path.basename(allocation_file)
    for region, path in sorted(nhd_files.iteritems()):
        region_start = time.time()
        print "\tProcessing region {}...".format(region)
        output_basename = os.path.join(output_dir, "{}_{}.csv".format(basename, region))
        if not os.path.exists(output_basename) or overwrite:
            catchment_raster = Raster(os.path.join(path, "NHDPlusCatchment", "cat"))
            allocation_array = allocate(Raster(allocation_file), catchment_raster, tile=1000000)
            if accumulate:
                region_hierarchy = Hierarchy(path, region)
                allocation_array = accumulate(allocation_array, region_hierarchy)
            allocation_dict = make_histogram(allocation_array)
            write_to_file(allocation_dict, output_basename)
            print "\tRegion complete in {} seconds".format(int(time.time() - region_start))
        else:
            print "\t...already processed"

def main():
    nhd_dir = r"G:\NationalData\NHDPlusV2"
    attribution_file = r"G:\NationalData\CDL\CDL_finals\cdlx_10"
    output_dir = r"G:\CA3T\TestOutput"
    overwrite = True
    ca3t(nhd_dir, attribution_file, output_dir, overwrite)

main()