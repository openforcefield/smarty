# Testing slow down when storing ChemicalEnvironments

We were concerned that storing chemical environments would be slower than storing strings. 
Since ChemicalEnvironments can easily be converted to and from SMIRKS strings you could store a list of SMIRKS instead of a list of chemical environments when sampling parameter types (such as smirky). 
The notebook `testing_smirks_speed.ipynb` logs the time to store a list of SMIRKS or environments for a number of steps. It uses `Torsion_0_0.00e+00_results.smarts` as an example of the complext SMIRKS patterns that can be generated during a smirky simulation. 

Below are the results for this test. For each test data is reported in this order:
* Parameter type list
    - generic: starts with only `"[*:1]~[*:2]~[*:3]~[*:4]"`
    - short: starts with first 10 SMIRKS in `*_results.smarts`
    - long: starts with all 82 SMIRKS in `*_results.smarts`
* Time in minutes to do X iterations storing SMIRKS strings
* Time in minutes to do X iterations storing Chemical Environments for each input SMIRKS
* Difference in Chemical Environment and SMIRKS time in minutes 

```
------------------------------  2 Iterations  ------------------------------
               short    1.97e-05    6.54e-05    4.57e-05
                long    1.93e-05    4.58e-04    4.39e-04
             generic    1.34e-05    1.82e-05    4.84e-06


------------------------------  10 Iterations  ------------------------------
               short    7.12e-05    1.16e-04    4.53e-05
                long    8.27e-05    5.40e-04    4.58e-04
             generic    6.60e-05    6.47e-05    -1.23e-06


------------------------------  100 Iterations  ------------------------------
               short    6.19e-04    7.01e-04    8.20e-05
                long    7.44e-04    1.36e-03    6.12e-04
             generic    5.49e-04    6.28e-04    7.92e-05


------------------------------  1000 Iterations  ------------------------------
               short    7.59e-03    1.73e-02    9.76e-03
                long    8.42e-03    2.10e-02    1.26e-02
             generic    6.89e-03    1.61e-02    9.20e-03


------------------------------  10000 Iterations  ------------------------------
               short    8.89e-02    1.09e+00    9.98e-01
                long    9.37e-02    1.17e+00    1.08e+00
             generic    7.18e-02    1.12e+00    1.05e+00


------------------------------  30000 Iterations  ------------------------------
               short    3.61e-01    1.04e+01    1.00e+01
                long    4.51e-01    1.08e+01    1.04e+01
             generic    3.13e-01    1.01e+01    9.76e+00
```

We concluded from this that while the timing difference isn't so significant on the number of iterations typically run with smirky, future move proposal engines would probably benefit from storing SMIRKS patterns rather than Chemical Environments. 
