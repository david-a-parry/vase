try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(
    name = "vase",
    packages = ["vase"],
    version = "0.1.0",
    description = "Variant Annotation, Segregation and Exclusion",
    author = "David A. Parry",
    author_email = "david.parry@igmm.ed.ac.uk",
    url = "https://github.com/gantzgraf/vase",
    license='GPLv3',
    install_requires=[
          'pysam',
          'parse_vcf>=0.2.1',
      ],
    scripts = ["bin/vase", "bin/burden_test_vase," "bin/vase_reporter"],
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)
