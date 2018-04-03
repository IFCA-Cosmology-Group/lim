# lim

lim is a python application designed to analytically compute various statistics of line intensity maps using a wide variety of models.  This code is a work in progress, so it may change significantly and there may be undetected bugs.

### Prerequisites

lim requires several packages which should be familiar to astronomers working with python, including numpy, scipy, and astropy.  It also makes substantial use of Steven Murray's [hmf](https://www.github.com/steven-murray/hmf) package, which can be installed along with its dependencies with the command

```
pip install hmf
```

Finally, lim uses the hmf matter power spectrum, which uses the python camb wrapper if available, and the Eisenstein-Hu transfer function if not.

### Quickstart

In the folder containing the lim functions, you can quickly get the default CO power spectrum by running in an interpreter

```
from lim import LineModel
m = LineModel()
m.Pk
```

All parameters have default values, which can be changed either when creating the model or using the built-in update() method.  For example, to change the observing frequency from the default you could either run

```
m = LineModel(nuObs=15*u.GHz)
m.z
```

or

```
m = LineModel()
m.update(nuObs=15*u.GHz)
m.z
```

## DocTests

To quickly check several of the parameters, lim includes doctests.  In a terminal, simply run

```
python lim.py
```

Note that the expected power spectra for the doctests were computed assuming the camb module is installed.  If you do not have camb installed, i.e. if

```
import camb
```
gives an error, hmf will use the EH transfer function and the doctest on Pk will give incorrect numbers

## Authors

* **Patrick C. Breysse**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Code based on matlab routines originally developed with Ely Kovetz
* Thanks to Dongwoo Chung and George Stein for debugging help


