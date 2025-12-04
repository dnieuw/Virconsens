import os
import re
from setuptools import setup, find_packages

# Read the version from __init__.py without importing the package
with open(os.path.join("virconsens", "__init__.py"), "r") as f:
    version_file = f.read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        __version__ = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string.")

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
