try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

test_requirements = ['nose', 'xlrd']
setup(
    name="vase",
    packages=["vase"],
    version="0.4.2",
    description="Variant Annotation, Segregation and Exclusion",
    author="David A. Parry",
    author_email="david.parry@igmm.ed.ac.uk",
    url="https://github.com/david-a-parry/vase",
    download_url='https://github.com/david-a-parry/vase/archive/0.4.2.tar.gz',
    license='MIT',
    install_requires=[
          'pysam>=0.14',
          'natsort',
          'numpy',
    ],
    tests_require=test_requirements,
    extras_require={
        'BGZIP': ['biopython'],
        'REPORTER': ['xlsxwriter', 'requests'],
        'MYGENEINFO': ['mygene'],
        'tests': test_requirements,
    },
    scripts=["bin/vase", "bin/burden_test_vase", "bin/vase_reporter",
             "bin/coordinates_from_genes", "bin/filter_gts",
             "bin/phase_by_transmission", "bin/remove_info_fields",
            ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)
