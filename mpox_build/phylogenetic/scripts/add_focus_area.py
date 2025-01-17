import pandas as pd
import argparse


def load_metadata(file):
    '''
    Loads metadata tsv as df.
    '''
    with open(file) as tfile:
        metadata = pd.read_csv(tfile, sep = '\t', parse_dates= ['date'])
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description= "adds focus area to mpox sequences" ,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--metadata', type=str, required=True, help= "downloaded and unzipped metadata file ")
    parser.add_argument('--output', type=str, required=True, help= "new metadata file" )

    args = parser.parse_args()


metadata = load_metadata(args.metadata)


metadata_for_build = metadata.copy()
metadata_for_build["focus_areas"] = metadata_for_build.region
metadata_for_build.focus_areas[metadata_for_build.focus_areas != "North America"] = "Global"
metadata_for_build.focus_areas[(metadata_for_build.institution.str.contains("Los Angeles")) & (metadata_for_build.country == "USA") ] = "Los Angeles County"
metadata_for_build.focus_areas[metadata_for_build.division == "New York" ] = "New York City"
metadata_for_build.focus_areas[(metadata_for_build.institution.str.contains("UW")) & (metadata_for_build.division == "Washington") ] = "Washington"
metadata_for_build.focus_areas[(metadata_for_build.institution.str.contains("CDPH")) & (metadata_for_build.country == "USA") ] = "Other California"
metadata_for_build.focus_areas[metadata_for_build.location == "Cook County IL"] = "Cook County"


with open(args.output, 'w') as f:
    metadata_for_build.to_csv(f, sep = '\t')
