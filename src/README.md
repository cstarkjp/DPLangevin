# Build notes

Full deployment to PyPI is in progress. For now, if you want to use this
package, first set up the Python environment (see `environment.yml` for 
a `conda` install), then do a build using `meson-python`. 

Either build and deploy to your local Python environment (note: do this from the cloned repo root directory, not in `src/`)

    rm -rf build; pip install .

or just build in-place

    rm -rf build; meson setup build; meson compile -C build   

Then run either the test Python script or the test Jupyter/Python notebook.