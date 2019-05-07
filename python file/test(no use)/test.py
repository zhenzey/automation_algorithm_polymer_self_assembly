import pyswarms as ps
import numpy as np
from pyswarms.utils.functions import single_obj as fx

def main():
    # Set-up hyperparameters
    options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}

    # Call instance of PSO
    optimizer = ps.single.GlobalBestPSO(n_particles=10, dimensions=2, options=options)

    # Perform optimization
    cost, pos = optimizer.optimize(fx.sphere_func, print_step=1, iters=100, verbose=3)
    return cost, post
if __name__ == "__main__":
    main()