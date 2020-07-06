"""
BAP - Bond Angle Printer 
"""

from os.path import abspath, dirname
from setuptools import setup, find_packages

project_dir = abspath(dirname(__file__))

setup(
    name='bap',
    version='1.0.0',
    description='quick way of getting bond angles out of a VASP POSCAR',
    url="https://github.com/badw/bap",
    author="Benjamin A. D. Williamson",
    author_email="benjamin.williamson@ntnu.no",
    license='MIT',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics'
        ],
    keywords='chemistry pymatgen physics dft vasp',
    packages=find_packages(),
    install_requires=['tabulate','numpy', 'pandas', 'pymatgen','argparse'],
    entry_points={'console_scripts': [
                      'bap-bap = bap.bond_angles:main',
                      'bap-tilt = bap.tilt_angles:main']},
    )
