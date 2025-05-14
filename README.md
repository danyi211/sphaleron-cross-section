# Sphaleron theory cross section

## Install LHAPDF
The installation steps are from the [official documents](https://www.lhapdf.org/install.html).
```
wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.X.Y.tar.gz -O LHAPDF-6.X.Y.tar.gz
# ^ or use a web browser to download, which will get the filename correct
tar xf LHAPDF-6.X.Y.tar.gz
cd LHAPDF-6.X.Y
./configure --prefix=/path/for/installation
make
make install
```
If you aren't installing LHAPDF into a standard, system-wide location, you will probably need to set some environment variables for it to work at runtime (and perhaps even during the installation).
```
export PATH=$PATH:/foo/lhapdf/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/foo/lhapdf/lib
export PYTHONPATH=$PYTHONPATH:/foo/lhapdf/lib/python3.9/site-packages
```
Here the `/foo` is usually the `/path/for/installation` in the configure step above.

## Integrate partonic cross-section with PDF
The partonic cross-section is Eq. (2.5) from [arXiv:1910.04761](https://arxiv.org/abs/1910.04761).
