# mc2err

A small C library with a Python interface to calculate averages and error bars for Markov chain Monte Carlo data streams.
It ingests ordered data in arbitrarily sized segments and updates its estimate of averages and error bars when they are requested.
The analysis backend estimates an equilibration point to remove unequilibrated data
 and truncates the autocorrection function to ignore small correlations in which the signal is overwhelmed by noise.

This library is still in a pre-alpha state,
and the main purpose of the repository at the moment is to test/demonstrate the Python interface for the library.

`setup.py` to compile and install the C library and its Python bindings (requires the C Foreign Function Interface package).
 
`example.py` is a basic demonstration of how the library is intended to be used.

Preliminary functionality should be operational in early 2021.
This library will be presented in a talk at the 2021 APS March Meeting,
 and a first release should occur well in advance of the conference.
