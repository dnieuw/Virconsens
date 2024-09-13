from setuptools import setup, find_packages
from virconsens import __version__

setup(name='virconsens',
    version=__version__,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'virconsens=virconsens.virconsens:main',
        ]
    },
    install_requires=[
        'biopython>=1.70',
        'pysam>=0.20'
    ],
    description='Tool to create a consensus sequence from mapped virus Nanopore data',
    url='https://github.com/dnieuw/virconsens',
    author='David F. Nieuwenhuijse',
    author_email='d.nieuwenhuijse@erasmusmc.nl',
    license='BSD 3-Clause',
    zip_safe=False)
