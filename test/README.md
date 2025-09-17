# Testing

A very basic test is provided in the Python script:

     python test_dp.py

which performs an integration of the DP-Langevin equation 
on a small rectangular 2D grid. The grid-averaged mean density 
is written out as a time series into a text file named "solution.txt" for 
all simulation time steps ("epochs"). 
Time slices of the density grid are printed for selected epochs.

A similar test is provided in the Jupyter notebook 

    TestDP.ipynb

but the final density grid is also plotted as an image, and the mean-density time series is graphed. Both plots are exported to PNG files.

If you build from source, you will need to point Python to your local build of `dplvn`. 
Comment out the lines:

    import sys, os
    sys.path.insert(0, os.path.join(os.path.pardir, "build"))

and check that Python is using your local build:

    print(dplvn.__file__)
