try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(
    name = "vase",
    packages = ["vase"],
    version = "0.1.1a",
    description = "Variant Annotation, Segregation and Exclusion",
    author = "David A. Parry",
    author_email = "david.parry@igmm.ed.ac.uk",
    url = "https://github.com/gantzgraf/vase",
    download_url = 'https://github.com/gantzgraf/vase/archive/v0.1.1a.tar.gz',
    license='MIT',
    install_requires=[
          'pysam',
          'parse_vcf>=0.2.1',
          'natsort',
      ],
    scripts = ["bin/vase", "bin/burden_test_vase", "bin/vase_reporter"],
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
