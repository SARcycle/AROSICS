#!/bin/bash
INSTALLDIR=/home/localadmin/.venv/arosicspy_test
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo snap refresh
sudo apt update -y
sudo apt upgrade -y
sudo apt autoremove -y
sudo apt install -y build-essential 
sudo apt install -y python3.10-dev python3.10-venv gdal-bin libgdal-dev python3-gdal
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal
python3.10 -m venv $INSTALLDIR
$INSTALLDIR/bin/python3 -m pip install --upgrade pip
$INSTALLDIR/bin/pip3 install --no-cache-dir "numpy>1.0.0"
$INSTALLDIR/bin/pip3 install --no-cache-dir wheel
$INSTALLDIR/bin/pip3 install --no-cache-dir "setuptools>=67"
$INSTALLDIR/bin/pip3 install --no-cache-dir --no-cache --force-reinstall gdal[numpy]=="$(gdal-config --version).*"
#$INSTALLDIR/bin/pip3 install --no-cache-dir GDAL==$(gdal-config --version)
$INSTALLDIR/bin/pip3 install --no-cache-dir cartopy
$INSTALLDIR/bin/pip3 install --no-cache-dir geopandas
$INSTALLDIR/bin/pip3 install --no-cache-dir "joblib>=1.3.0"
$INSTALLDIR/bin/pip3 install --no-cache-dir matplotlib
$INSTALLDIR/bin/pip3 install --no-cache-dir numpy
$INSTALLDIR/bin/pip3 install --no-cache-dir pandas
#$INSTALLDIR/bin/pip3 install --no-cache-dir pyfftw
$INSTALLDIR/bin/pip3 install --no-cache-dir pykrige
$INSTALLDIR/bin/pip3 install --no-cache-dir "pyproj>2.2.0"
$INSTALLDIR/bin/pip3 install --no-cache-dir "scikit-image>=0.21.0"
$INSTALLDIR/bin/pip3 install --no-cache-dir shapely
$INSTALLDIR/bin/pip3 install --no-cache-dir tqdm
$INSTALLDIR/bin/pip3 install --no-cache-dir matplotlib
$INSTALLDIR/bin/pip3 install --no-cache-dir wget
$INSTALLDIR/bin/pip3 install --no-cache-dir arosics
$INSTALLDIR/bin/pip3 install --no-cache-dir rasterio
$INSTALLDIR/bin/pip3 install --no-cache-dir xarray
$INSTALLDIR/bin/pip3 install --no-cache-dir rioxarray
