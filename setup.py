try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(
    name = "vase",
    packages = ["vase"],
    version = "0.2.5",
    description = "Variant Annotation, Segregation and Exclusion",
    author = "David A. Parry",
    author_email = "david.parry@igmm.ed.ac.uk",
    url = "https://github.com/david-a-parry/vase",
    download_url = 'https://github.com/david-a-parry/vase/archive/0.2.5.tar.gz',
    license='MIT',
    install_requires=[
          'pysam',
          'parse_vcf>=0.2.6',
          'natsort',
    ],
    extras_require={
        'BGZIP': ['biopython'],
        'REPORTER': ['xlsxwriter', 'requests'],
        'MYGENEINFO': ['mygene'],
    },
    scripts = ["bin/vase", "bin/burden_test_vase", "bin/vase_reporter",
               "bin/coordinates_from_genes", "bin/filter_gts",
               "bin/phase_by_transmission"],
    include_package_data=True,
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)
