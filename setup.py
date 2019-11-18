#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# import sys
import setuptools
import lstchain
import os

def list_scripts():
    script_dir = 'scripts'
    script_list = [f'{os.path.join(script_dir, f)}' for f in os.listdir(script_dir) if f.startswith('lstchain_')]
    return script_list

entry_points = {}
entry_points['console_scripts'] = [
    'calc_camera_calibration = lstchain.tools.calc_camera_calibration:main',
    'onsite_create_drs4_pedestal_file = scripts.onsite.onsite_create_drs4_pedestal_file:main',
    'onsite_create_calibration_file = scripts.onsite.onsite_create_calibration_file:main',
    'lstchain_data_create_pedestal_file = scripts.lstchain_data_create_pedestal_file:main'
]


setuptools.setup(name='lstchain',
                 version=lstchain.__version__,
                 description="DESCRIPTION",  # these should be minimum list of what is needed to run
                 packages=setuptools.find_packages(),
                 install_requires=['h5py',
                                   'seaborn'
                                   ],
                 package_data={'lstchain': ['data/lstchain_standard_config.json']},
                 tests_require=['pytest', 'pytest-ordering'],
                 author='LST collaboration',
                 author_email='',
                 license='',
                 url='https://github.com/cta-observatory/cta-lstchain',
                 long_description='',
                 classifiers=[],
                 entry_points=entry_points,
                 scripts=list_scripts()
                 )
