#!/usr/bin/env python3

"""
Pipeline for the reconstruction of Energy, disp and gamma/hadron
separation of events stored in a simtelarray file.
- Input: DL2 file.
- Output: Gamma candidate list.

Usage:

$> python lstchain_dl2_to_dl3.py 
--input-file dl2_LST-1.Run02033.0137.h5

"""

import argparse
import astropy.units as u
import numpy as np
import os
import pandas as pd
from tables import open_file
import joblib
from lstchain.reco.utils import filter_events
from lstchain.reco import dl1_to_dl2
from lstchain.io import (
    read_configuration_file,
    standard_config,
    replace_config,
    write_dl2_dataframe,
    get_dataset_keys,
)
from lstchain.io.io import (
    dl1_params_lstcam_key,
    dl1_params_src_dep_lstcam_key,
    dl1_images_lstcam_key,
    dl2_params_lstcam_key,
)


parser = argparse.ArgumentParser(description="DL2 to DL3")

# Required arguments
parser.add_argument('--input-file', '-f', type=str,
                    dest='input_file',
                    help='path to a DL2 HDF5 file',
                    default=None, required=True)

# Optional arguments
parser.add_argument('--output-dir', '-o', action='store', type=str,
                     dest='output_dir',
                     help='Path where to store the reco dl3 events',
                     default='./dl3_data')

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
    dl2_filtered = filter_events(data, filters=config["events_filters"])

    os.makedirs(args.output_dir, exist_ok=True)
    output_file = os.path.join(args.output_dir, os.path.basename(args.input_file).replace('dl2','dl3'))

    if os.path.exists(output_file):
        raise IOError(output_file + ' exists, exiting.')

    dl2_keys = get_dataset_keys(args.input_file)

    if dl2_params_lstcam_key in dl2_keys:
        dl2_keys.remove(dl2_params_lstcam_key)

    with open_file(args.input_file, 'r') as h5in:
        with open_file(output_file, 'a') as h5out:

            # Write the selected DL1 info
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

    write_dl2_dataframe(dl2_filtered, output_file)

if __name__ == '__main__':
    main()
