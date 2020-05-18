#!/usr/bin/env python3

"""
Pipeline to add pulsar parameters
- Input: DL1/2 data file.
- Output: DL1/2 data file with barycenter corrected times and phases into the parameters dataframe

Usage:

$> python lstchain_add_pulsar_parameters.py 
--input-file dl2_LST-1.Run02033.0137.h5
--eph-file

"""

import argparse
import numpy as np
import pandas as pd
import astropy.units as u
import pint.toa as toa

from astropy.time import Time
import pint.models as models
from tables import open_file
import os

np.set_printoptions(precision=18)

from lstchain.io import (
    get_dataset_keys,
    write_dl2_dataframe
)
from lstchain.reco.utils import filter_events
from lstchain.io.io import (
    dl2_params_lstcam_key,
)


parser = argparse.ArgumentParser(description="DL1 to DL2")

# Required arguments
parser.add_argument('--input-file', '-f', type=str,
                    dest='input_file',
                    help='path to a DL2 HDF5 file',
                    default=None, required=True)

parser.add_argument('--path-models', '-p', action='store', type=str,
                     dest='path_models',
                     help='Path where to find the trained RF',
                     default='./trained_models')

# Optional arguments
parser.add_argument('--output-dir', '-o', action='store', type=str,
                     dest='output_dir',
                     help='Path where to store the reco dl2 events',
                     default='./dl2_data')


parser.add_argument('--config', '-c', action='store', type=str,
                    dest='config_file',
                    help='Path to a configuration file. If none is given, a standard configuration is applied',
                    default=None, required=False)



args = parser.parse_args()

def main():

    custom_config = {}
    if args.config_file is not None:
        try:
            custom_config = read_configuration_file(os.path.abspath(args.config_file))
        except("Custom configuration could not be loaded !!!"):
            pass

    config = replace_config(standard_config, custom_config)

    data = pd.read_hdf(args.input_file, key=dl2_params_lstcam_key)
    data = filter_events(data, filters=config["events_filters"])

    os.makedirs(args.output_dir, exist_ok=True)
    output_file = os.path.join(args.output_dir, 
                               os.path.splitext(os.path.basename(args.input_file))[0] + '_corr.h5')
    if os.path.exists(output_file):
        raise IOError(output_file + ' exists, exiting.')

    string_array = []
    
    data = pd.read_hdf(args.input_file, key=dl2_params_lstcam_key)
    timestamp = data['dragon_time']
    utc = Time(timestamp, format='unix', scale='utc')
    mjd = utc.mjd
    
    for i in range(0,len(mjd)):
        string_array.append(f'{mjd[i]}.000004.1.000.000.8y.x.ff 1432.00000000 {mjd[i]} 3.72100 lst')
    
    with open(f'{output_file}.tim', 'w') as f:
        f.write('FORMAT 1\n')
        f.write('MODE 1\n')
        for item in string_array:
            f.write("%s\n" % item)
            
    t = toa.get_TOAs(f'{output_file}.tim', usepickle=False)

    data['tdb'] = m.get_barycentric_toas(t)
    data['phase'] = m.phase(t)[1]

    dl2_keys = get_dataset_keys(args.input_file)
    if dl2_params_lstcam_key in dl2_keys:
        dl2_keys.remove(dl2_params_lstcam_key)

    with open_file(args.input_file, 'r') as h5in:
        with open_file(output_file, 'a') as h5out:

            # Write the selected DL2 info
            for k in dl2_keys:
                if not k.startswith('/'):
                    k = '/' + k

                path, name = k.rsplit('/', 1)
                if path not in h5out:
                    grouppath, groupname = path.rsplit('/', 1)
                    g = h5out.create_group(
                        grouppath, groupname, createparents=True
                        )
                else:
                    g = h5out.get_node(path)

                h5in.copy_node(k, g, overwrite=True)

    write_dl2_dataframe(data, output_file)
    os.remove(f'{output_file}.tim')

if __name__ == '__main__':
    main()
