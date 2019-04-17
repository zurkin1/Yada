#Run using:
#python setup.py sdist bdist_wheel
#python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/yada-0.1*
#(pass is Feb..9x2!)
#Manage at: https://test.pypi.org/manage/projects/
#Install with:
#echo 'y'| pip uninstall yada
#pip install -i https://test.pypi.org/simple/ --no-cache-dir yada
from setuptools import find_packages, setup

with open("README.md", "r") as fh:
    long_description = fh.read()

package_data={'yada': ['data/*.*']}

setup(name='yada',
      version='0.4',
      description='An ML package for gene deconvolution',
      long_description=long_description,
      url='http://github.com/zurkin1/yada',
      author='Dani Livne',
      author_email='dani.livne@yahoo.com',
      license='MIT',
      #packages=['yada'],
      packages=find_packages(),
      zip_safe=False,
      include_package_data=True,
      classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            ]
      )