from array_tools import *
import os
import time

def nhd_batch(allocation_dir, nhd_dir, output_dir, accumulate=True, overwrite=True, region_filter='all', allocation_filter='all'):
    # Fetches necessary NHD Plus files for running allocation and accumulation
    def assemble_batch():
        # Assemble rasters for allocation
        all_files = []
        for f in os.listdir(allocation_dir):
            if f in allocation_filter or (allocation_filter == 'all' and not "." in f):
                all_files.append(os.path.join(allocation_dir, f))

        # Get catchment grid, gridcode to comid translation files, and flow tables
        h_files = {}
        regions = [d.lstrip("NHDPlus") for d in os.listdir(nhd_dir)
                   if os.path.isdir(os.path.join(nhd_dir, d)) and d.startswith("NHDPlus")]

        for r in regions:
            if r in region_filter or region_filter == 'all':
                catchment_path = r"NHDPlusCatchment\cat"
                translation_path = r"NHDPlusCatchment\featureidgridcode.dbf"
                flow_path = r"NHDPlusAttributes\PlusFlow.dbf"
                nhd_path = r"NHDSnapshot\Hydrography\NHDFlowline.dbf"
                h_files[r] = [os.path.join(nhd_dir, r"NHDPlus{}".format(r), f)
                              for f in catchment_path, flow_path, translation_path, nhd_path]
        return all_files, h_files

    allocation_files, hydro_files = assemble_batch()
    for allocation_file in allocation_files:
        if not allocation_file.endswith("info"):
            file_start = time.time()
            basename = os.path.basename(allocation_file)
            print "Processing file {}...".format(allocation_file)
            for region in sorted(hydro_files.keys()):
                if not region == "08":
                    region_start = time.time()
                    print "\tProcessing region {}...".format(region)
                    output_basename = os.path.join(output_dir, "{}_{}.csv".format(basename, region))
                    if not os.path.exists(output_basename) or overwrite:
                        catchment_raster, flow_table, translation_file, nhd_file = hydro_files[region]
                        translate_f, translate_b = read_translation(translation_file, "FEATUREID", "GRIDCODE", flip=True)
                        allocation_dict = allocate(Raster(allocation_file), Raster(catchment_raster), tile=1000000)
                        if accumulate:
                            excl_f, excl_v = 'FCODE', [56600]
                            allocation_dict = accumulate(allocation_dict, flow_table, translate_f, nhd_file, excl_f, excl_v)
                        output = make_histogram(allocation_dict)
                        write_to_file(output, output_basename, translate_b)
                        print "\tRegion complete in {} seconds".format(int(time.time() - region_start))
                    else:
                        print "\t...already processed"
            print "File complete in {} seconds".format(int(time.time() - file_start))

def main():
    allocation_directory = r"G:\WARP\Data\Projected"
    nhd_directory = r"G:\NationalData\NHDPlusV2"
    output_directory = r"G:\WARP\Outputs"
    accumulate = False
    overwrite_output = False
    region_filter = 'all'
    allocation_filter = 'all'
    nhd_batch(allocation_directory, nhd_directory, output_directory, accumulate, overwrite_output, region_filter, allocation_filter)

main()