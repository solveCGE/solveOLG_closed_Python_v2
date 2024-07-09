# solveOLG closed economy in Python
Solves a simple AK-OLG-model for a closed economy in Python/Numba

## About
Shows how to solve a simple deterministic overlapping-generations model (OLG) of Auerbauch-Kotlikoff type, solving for the transition path between two steady-states. The code is partially optimized for speed by using Numba's just-in-time compiler (@jit) for the most frequently called household solution finding function. The algorithm could be improved by solving the household problems of all finitely lived representative cohorts in parallel.

A model description can be found here: <https://github.com/solveCGE/solveOLG_doc>.

## How to run
Parameters can be set in `calib.py`. Policy shocks are defined in `run.py` (or just uncomment some of the predefined exemplary shocks). The model is then solved by just running `run.py`. 

## Author
Philip Schuster
