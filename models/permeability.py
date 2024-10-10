import numpy as np
import scipy.io

#########################################################
# This class handles permeability matrix initialization #
# and evaluation based on user input.                   #
#########################################################

class Permeability:
    def __init__(
        self, 
        simulation_type : str,
        size_of_grid : int,
        x : np.ndarray,
        y : np.ndarray
    ):
        self.simulation_type = simulation_type
        self.KK = self.KKdef(size_of_grid, x, y)
        
    def KKdef(self, size_of_grid : int, x : np.ndarray, y : np.ndarray):
        # defining permeability matrix
        if self.simulation_type == 'rectilinear homogenous':
            K_max = 100
            return K_max * np.ones(size_of_grid + 1)
        elif self.simulation_type == 'rectilinear heterogenous':
            K_max = 100
            return K_max * (0.5 * (1 - 10 ** (-7)) * (np.sin(6 * np.pi * np.cos(x)) * np.cos(4 * np.pi * np.sin(3 * y)) - 1) + 1)
        elif self.simulation_type == 'quarter-five spot heterogenous':
            # Simulating with Upper Ness permeability from SPE10
            data = scipy.io.loadmat('./data/KK30Ness.mat')
            print('Tarbert formation permeability loaded')
            return data['KK']
        else:
            raise Exception('Inputted simulation type is invalid when performing KKdef')
            

