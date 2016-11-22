# Create_Planet
This code takes in a MESA (Modules for Experiments in Stellar Astrophysics) gas-envelope profile and generates a core by integrating backwards from the innermost grid point.

## How to use
Use the RUN_COREGEN_DIFF.py script to generate a differentiated core, which solves for the correct mantle/core mass-fractions assuming equations of state from Seager 2007.
Use the RUN_COREGEN.py script to generate a solid core, which solves for the core density at the boundary assuming a polytropic equation of state.

## Notes
Other equations of state are in the process of being included
