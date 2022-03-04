#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from setuptools import setup

setup(
    name='exposan',
    packages=['exposan'],
    version='1.1.1',
    license='UIUC',
    author='Quantitative Sustainable Design Group',
    author_email='quantitative.sustainable.design@gmail.com',
    description='Exposition of sanitation and resource recovery systems',
    long_description=open('README.rst', encoding='utf-8').read(),
    url='https://github.com/QSD-Group/EXPOsan',
    project_urls={
        'QSDsan': 'https://github.com/QSD-Group/QSDsan',
        'QSDsan documentation': 'https://qsdsan.readthedocs.io/',
    },
    install_requires=['qsdsan>=1.1.3',],
    package_data=
        {'exposan': [
                    'bsm1/*',
                    'bsm1/data/*',
                    'bwaise/*',
                    'bwaise/data/*',
                    'cas/*',
                    ]},
    platforms=['Windows', 'Mac', 'Linux'],
    classifiers=['License :: OSI Approved :: University of Illinois/NCSA Open Source License',
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
                 'Operating System :: MacOS',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: POSIX :: BSD',
                 'Operating System :: POSIX :: Linux',
                 'Operating System :: Unix',
                 'Programming Language :: Python :: 3.8',
                 'Programming Language :: Python :: 3.9',
                 ],
    keywords=['quantitative sustainable design', 'sanitation', 'resource recovery', 'techno-economic analysis', 'life cycle assessment'],
)