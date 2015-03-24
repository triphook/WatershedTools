from array_tools import *
import collections


def timesheds(nhd_dir, threshold, outfile):
    threshold *= 24 * 60 * 60 # days to seconds
    nhd_paths = get_nhd(nhd_dir)
    for region, path in sorted(nhd_paths.items()):
        print region
        h = Hierarchy(path, region)
        for path in h[5205118].upstream_paths:
            print sum(reach.time for reach in path)


def main():
    nhd_dir = r"G:\NationalData\NHDPlusV2"
    threshold = 5  #days
    outfile = r"G:\SAM\Outputs\test.txt"
    timesheds(nhd_dir, threshold, outfile)


main()