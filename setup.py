import sys
from setuptools import setup, find_packages

import ez_setup
ez_setup.use_setuptools()

long_description = """
A traits based Melodic (ICA) results viewer for neuroscientists
"""
requires = [
    'setuptools',
    'matplotlib >= 1.3',
    'nibabel >= 2.1',
    'numpy >= 1.8',
    'scipy >= 0.13'
]
if sys.version < (3, 2, 0):
    requires.append('backports.functools_lru_cache >= 1.3')

setup(
    name='melview',
    version='1.0.1',
    description='Melodic ICA data viewer',
    author='Dave Flitney',
    author_email='flitney@fmrib.ox.ac.uk',
    packages=find_packages(),
    long_description=long_description,
    classifiers=[
        "Licence :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
    ],
    entry_points={
        'console_scripts': [
            'melview=melview.melodic_traits:main',
            'mv2fix=melview.mv2fix:main',
        ]
    },
    keywords='melodic viewer',
    license='GPL',
    install_requires=requires
)
