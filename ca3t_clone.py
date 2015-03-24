from array_tools import *
import os
import time

# Sample script designed to replicate the functionality of the
# NHD Plus Catchment Attribution Allocation and Accumulation Tool (CA3T2)

# For this to function correctly, the attribution raster MUST have the same projection and cell size as the NHD
# plus catchment rasters

# Loops through the NHD Plus regions and performs attribution, allocation and accumulation
def ca3t(nhd_dir, allocation_file, output_dir, overwrite):

    # Assemble NHD file path structure
    nhd_files = get_nhd(nhd_dir)

    # Loop through NHD regions
    for region, path in sorted(nhd_files.iteritems()):

        region_start = time.time()  # starts function timer
        print "\tProcessing region {}...".format(region)

        # Set output file name
        output_filename = os.path.join(output_dir, "{}_{}.csv".format(os.path.basename(allocation_file), region))

        # Skips regions that have already been processed unless overwrite is on
        if not os.path.exists(output_filename) or overwrite:

            # Initialize the catchment (zone) raster
            catchment_raster = Raster(os.path.join(path, "NHDPlusCatchment", "cat"))

            # Perform allocation
            allocation_array = allocate(Raster(allocation_file), catchment_raster, tile=1000000)

            # Perform accumulation
            if accumulate:
                region_hierarchy = Hierarchy(path, region)
                allocation_array = accumulate(allocation_array, region_hierarchy)

            # Convert allocation array into histogram in dictionary format {zone1: zone1count, zone2: zone2count...}
            allocation_dict = make_histogram(allocation_array)

            # Write dictionary to output file
            write_to_file(allocation_dict, output_filename)

            print "\tRegion complete in {} seconds".format(int(time.time() - region_start))
        else:
            print "\t...already processed"

#  Set parameters here
def main():
    nhd_dir = r"G:\NationalData\NHDPlusV2"  # The parent directory containing all NHD files
    attribution_file = r"G:\NationalData\CDL\CDL_finals\cdlx_10"  # The grid file to be allocated
    output_dir = r"G:\CA3T\TestOutput"  # Output directory
    overwrite = True  # Overwrite output
    ca3t(nhd_dir, attribution_file, output_dir, overwrite)

main()