Each of these folders contains elastic nets run in the same fashion:
- s6.sh (or in general sX.sh, with X a natural number): Schedules sbatch jobs for all the combinations considered in the grid, runs runnerX.sh
- runnerX.sh: specifies the sbatch node requirments and runs the psX.py python code
- psX.py contains the actuall python code that is run

- Run_one_XXX.ipynb: Runs one example for one set of parameter choices. I used this to implement the psX.py code and test it.
- Load_results_XXX.ipynb: Loads the result from the sX.sh-runs and makes plots out of the results.
