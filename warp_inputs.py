import time
import urllib2
import cPickle
from array_tools import *

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

def catchment_loads(nhd_dir, state_dir, allocation_grid, use_table, compound, out_dir, overwrite=True):

    # Get all NHD regions and their corresponding catchment ID rasters and flow tables
    nhd_data = ((region,
                 Raster(os.path.join(nhd_dir, "NHDPlus{}".format(region), r"NHDPlusCatchment", "cat")),
                 os.path.join(nhd_dir, "NHDPlus{}".format(region), r"NHDPlusAttributes", "PlusFlow.dbf"),
                 os.path.join(nhd_dir, "NHDPlus{}".format(region), r"NHDSnapshot", "Hydrography", "NHDFlowline.dbf"),
                 os.path.join(nhd_dir, "NHDPlus{}".format(region), r"NHDPlusCatchment", "featureidgridcode.dbf"))
                 for region in sorted(overlaps.keys()))

    # Initialize land use raster for allocation
    land_use_raster = Raster(allocation_grid)

    # Load use data dictionary
    use_data = pull_use_data(use_table)[compound] # {'16073': {'high': '5', 'low': '5'}, '16077': {'high': '0.4'}, ...}

    # Loop through regions and run allocation/accumulation
    for region, nhd_raster, flow_table, nhd_atts, lookup_table in nhd_data:

        out_file = os.path.join(out_dir, "r{}_loads.csv".format(region))

        if overwrite or not os.path.exists(out_file):

            region_start = time.time()
            print "Processing region {}...".format(region)

            # Build dictionaries to translate from gridcode to comid and back
            comid_to_gridcode, gridcode_to_comid = read_translation(lookup_table, "FEATUREID", "GRIDCODE", True)

            # Get rasters for each state that the region overlaps
            state_rasters = [Raster(os.path.join(state_dir, "{}_cty_30m".format(state)), 6) for state in overlaps[region]]

            # Run through all states and count cells by county by watershed
            raw = np.hstack([allocate(land_use_raster, [nhd_raster, state_r]) for state_r in state_rasters])

            # Filter out pixels that will not be assigned any load
            filtered = raw[:, raw[2] > 0]

            # Get county use rates
            counties, area_by_county = sum_by_group(filtered[1], filtered[-1])
            use_by_county = map(lambda x: mean(map(float, use_data.get(str(x), {}).values())), counties)
            county_use_rates = dict(zip(counties, (use_by_county / area_by_county)))

            # Weight allocation by county use table weights
            weighted = np.vstack((filtered, filtered[-1] * map(lambda x: county_use_rates.get(x, 0.0), filtered[1])))

            # Dissolve watersheds with counts in multiple counties
            dissolved = sum_by_group(weighted[0], weighted[-1])
            dissolved = dict(zip(map(int, dissolved[0]), dissolved[1])) # {1835008: 4.5999999999767169 ...}

            # Accumulate all upstream values
            #accumulated = accumulate(dissolved, flow_table, comid_to_gridcode, nhd_atts, "FCODE", [56600])
            accumulated = dissolved

            # Write output to a csv file
            write_to_file(accumulated, out_file, gridcode_to_comid)

            region_time = int(time.time() - region_start)
            print "\tFinished region in {} minutes, {} seconds".format(int(region_time / 60.0), int(region_time % 60))


def main():
    year = 2009
    nhd_dir = r"G:\NationalData\NHDPlusV2"
    state_dir = r"G:\NationalData\County_Grids"
    allocation_grid = r"G:\NationalData\nlcd_2011_landcover_2011_edition_2014_10_10\nlcd_2011_ag"
    use_table = r"G:\WARP\UseTables\raw_{}.p".format(year)
    compound = "DIAZINON"
    out_dir = r"G:\WARP\Outputs\non_accum"

    catchment_loads(nhd_dir, state_dir, allocation_grid, use_table, compound, out_dir, False)

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

main()
