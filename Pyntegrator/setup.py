import os
from setuptools import setup, find_packages
import sys

try:
    import numpy
except ImportError:
    raise ImportError(
        "Please install NumPy first, or use the Anaconda Python Distribution "
        "(https://store.continuum.io/cshop/anaconda/) which comes with NumPy "
        "installed."
    )



version_py = os.path.join(os.path.dirname(__file__), 'Pyntegrator', 'version.py')
version = open(version_py).read().split('=')[-1].strip().replace('"','')

requirements = open(os.path.join(os.path.dirname(__file__), 'requirements.txt')).readlines()


long_description = open('README.rst').read()
setup(
        name='Pyntegrator',
        version=version,
        description="Integrative analysis of high-thoughput sequencing data",
        #long_description=long_description,
        license="MIT",
        install_requires=requirements,
        packages=find_packages(),
        package_data={
            'Pyntegrator':[
                'test/data/gdc*',
                'test/data/make_examples_from_pybedtools.py',
                'test/data/x.*',
            ]
        },
        scripts=[
            'Pyntegrator/scripts/download_metaseq_example_data.py',
            'Pyntegrator/scripts/metaseq-cli',
            'Pyntegrator/scripts/speedtest.py',
        ],
        author='Ryan Dale',
        author_email='dalerr@niddk.nih.gov',
        url='./Pyntegrator',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Intended Audience :: Developers',
            'Intended Audience :: System Administrators',
            'Operating System :: POSIX',
            'Operating System :: MacOS :: MacOS X',
            'Environment :: Console',
            'License :: OSI Approved :: MIT License',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Medical Science Apps.', 
        ]
    )
