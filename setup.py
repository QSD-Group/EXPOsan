#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from setuptools import setup

setup(
    name='exposan',
    packages=['exposan'],
    version='1.4.1',
    license='UIUC',
    author='Quantitative Sustainable Design Group',
    author_email='quantitative.sustainable.design@gmail.com',
    description='Exposition of sanitation and resource recovery systems',
    long_description=open('README.rst', encoding='utf-8').read(),
    url='https://github.com/QSD-Group/EXPOsan',
    project_urls={
        'Homepage': 'https://qsdsan.com',
        'Documentation': 'https://qsdsan.readthedocs.io',
        'Repository': 'https://github.com/QSD-Group/EXPOsan',
    },
    install_requires=['qsdsan',],
    package_data=
        {'exposan': [
            'adm/*',
            'adm/data/*',
            'asm/*',
            'asm/data/*',
            'biogenic_refinery/*',
            'biogenic_refinery/data/*',
            'bsm1/*',
            'bsm1/data/*',
            'bsm2/*',
            'bsm2/data/*',
            'bwaise/*',
            'bwaise/data/*',
            'cas/*',
            'eco_san/*',
            'eco_san/data/*',
            'htl/*',
            'htl/data/*',
            'metab/*',
            'metab/data/*',
            'metab/utils/*',
            'metro/*',
            'new_generator/*',
            'new_generator/data/*',
            'pm2_batch/*',
            'pm2_batch/data/*',
            'pm2_ecorecover/*',
            'pm2_ecorecover/data/*',
            'pou_disinfection/*',
            'pou_disinfection/data/*',
            'reclaimer/*',
            'reclaimer/data/*',
            'saf/*',
            'saf/data/*',
        ]},
    classifiers=[
        'License :: OSI Approved :: University of Illinois/NCSA Open Source License',
        'Environment :: Console',
        'Topic :: Education',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Manufacturing',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        ],
    keywords=[
        'quantitative sustainable design',
        'sanitation',
        'resource recovery',
        'techno-economic analysis',
        'life cycle assessment',
    ],
)