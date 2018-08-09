# This *should* allow us to run if DISPLAY is not set...
# ... however, Linux builds on Travis CI still fail.
#
# Instead of setting up a headless X server as in
# https://stackoverflow.com/questions/35403127/testing-matplotlib-based-plots-in-travis-ci
# I will set MPLBACKEND=agg
# (https://matplotlib.org/faq/environment_variables_faq.html)
#
# The following is left because (1) harmless and (2) documentation
import matplotlib
matplotlib.use('agg')
