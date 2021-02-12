# Functions for Kernel Densit Estimates

A small collection of functions written in pure R for calculating
kernel density estimates and working with them.  The functions use
simple approximate methods, but are adequate for most analysis or
visualisation purposes.

* kde() computes a simple kernel density estimate with a gaussian kernel similar to stats::density()
* kde_2d() computes a 2-dimensional KDE similar to MASS::kde2d
* p_fun(), q_fun(), d_fun() and r_fun() can be used to produce p-, q-, d- and r-prefix functions for working with distributions defined by KDEs, as per the convention used for standard statistical distributions in R.

The purpose of the package is mainly to illustrate some simple programming techniques, but may be uesful as a lightweight suite of tools for visualisation.
