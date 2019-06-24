import os
import sys
from setuptools import setup, find_packages

directory = os.path.abspath(os.path.dirname(__file__))
if sys.version_info >= (3, 0):
    with open(os.path.join(directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
else:
    with open(os.path.join(directory, 'README.md')) as f:
        long_description = f.read()

setup(name='S1_ARD',
      packages=find_packages(),
      include_package_data=True,
      version='0.9',
      description='A common workspace for our publication in the MDPI Data special issue on data cubes',
      classifiers=[
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python',
      ],
      install_requires=['numpy',
                        'matplotlib',
                        'pandas',
                        'scipy',
                        'jupyter',
                        'IPython',
                        'ipywidgets',
                        'scikit-learn',
                        'pyroSAR==0.9',
                        'astropy',
                        'spatialist==0.2.9',
                        'pathos'],
      extras_require={
          'docs': ['sphinx'],
      },
      url='https://github.com/johntruckenbrodt/S1_ARD',
      author='John Truckenbrodt',
      author_email='john.truckenbrodt@uni-jena.de',
      zip_safe=False,
      long_description=long_description,
      long_description_content_type='text/markdown')
