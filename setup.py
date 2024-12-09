"""
@file setup.py
@brief the PyPI required setup.py
"""

import os
from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

version = "0.1.0"

with open('requirements.txt') as f:
    install_requires = f.readlines()

setup(name='julia2-tool',
      version=version,
      license='MIT',
      author="Benjamin Glick",
      author_email='glick@glick.cloud',
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: Apache Software License",
          "Natural Language :: English", "Operating System :: OS Independent",
          "Programming Language :: Python :: 3",
          "Topic :: Scientific/Engineering"
      ],
      python_requires=">=3.6.0",
      packages=find_packages(),
      url='https://github.com/benhg/julia2-tool',
      keywords=["RNA", "Sequencing", "Index Hopping", "Bioinformatics"],
      entry_points={'console_scripts': ['julia2=julia2.julia2:main']},
      install_requires=install_requires,
      long_description=long_description,
      long_description_content_type='text/markdown'
)
