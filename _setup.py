import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.1.0'
PACKAGE_NAME = 'pyND'
AUTHOR = 'J. Christopher Howk'
AUTHOR_EMAIL = 'jhowk@nd.edu'
URL = 'https://github.com/jchowk/pyND'
github_project = 'jchowk/pyND'

LICENSE = 'GNU GENERAL PUBLIC LICENSE'
DESCRIPTION = 'Codes used by the Notre Dame CGM group, focused on absorption line studies.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

PYTHON_REQUIRES=">=3.6"

INSTALL_REQUIRES = [
      'numpy',
      'astropy',
      'scipy',
      'matplotlib',
      'linetools'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      packages=find_packages()
      )
