from distutils.core import setup
setup(
    name = "vape",
    packages = ["vape"],
    version = "0.0.1",
    description = "TO DO",
    author = "David A. Parry",
    author_email = "david.parry@igmm.ed.ac.uk",
    url = "https://github.com/gantzgraf/vape",
    license='GPLv3',
    install_requires=[
          'pysam',
          'parse_vcf',
      ],
    scripts = ["bin/vape"],
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        ],
)
