# Testing

A very basic test is provided in the Python script:

     python test_demo.py

which performs an integration of the DP-Langevin equation 
on a 2D grid of size 64x64. The grid-averaged mean density 
is written out as a time series into a text file named "solution.txt" for all 10001 epochs.

A similar test is provided in the Jupyter notebook 

    test_demo.ipynb

but the resulting time series is instead unpacked from the returned `numpy` 
array and graphed; the results plot is exported to a PNG file.