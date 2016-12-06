from setuptools import setup, find_packages

import ez_setup
ez_setup.use_setuptools()

setup(
    name='melview',
    version='1.0.1',
    description='Melodic ICA data viewer',
    author='Dave Flitney',
    author_email='flitney@fmrib.ox.ac.uk',
    packages=find_packages(),
    long_description="""\
melview is a traits based Melodic (ICA) results viewer for neuroscientists ...
    """,
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
    install_requires=[
        'setuptools',
        'matplotlib',
        'nibabel',
        'numpy',
        'scipy',
    ],
)
