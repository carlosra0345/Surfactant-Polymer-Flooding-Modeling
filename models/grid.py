###########################################################
# This class defines the grid and mesh operations,        #
# including initialization and evaluation of grid points. #
###########################################################

from config import Config
import numpy as np
from typing import Tuple

class Grid:
    def __init__(self, grid_size : int, simulation_type : str):
        self.grid_size = grid_size
        self.simulation_type = simulation_type
        self.dx = self.calculate_dx()
        self.dy = self.calculate_dy()
        self.x, self.y = self.create_meshgrid()
        self.phi_test = self.get_phi_test()

    def calculate_dx(self) -> float:
        # Calculate dx based on grid size
        return (Config.BOX['right'] - Config.BOX['left']) / self.grid_size

    def calculate_dy(self) -> float:
        # Calculate dy based on grid size
        return (Config.BOX['top'] - Config.BOX['bottom']) / self.grid_size

    def create_meshgrid(self) -> Tuple[np.ndarray, np.ndarray]:
        # Create and return meshgrid based on box limits and size
        x = np.arange(Config.BOX['left'], Config.BOX['right'] + self.dx, self.dx) 
        y = np.arange(Config.BOX['bottom'], Config.BOX['top'] + self.dy, self.dy)
        return np.meshgrid(x, y)
    
    def get_phi_test(self) -> np.ndarray:
        #  Determine which subdomain the grid points belong to
        
        #  Evaluate the level set function on the grid points
        #  Negative output signifies the domain $$\Omega^-$$.
        #  Positive output signifies the domain $$\Omega^+$$.
        ii, jj = np.meshgrid(
            np.arange(1, self.grid_size + 2),
            np.arange(1, self.grid_size + 2)
        )
        return self.z_func_test(
            Config.BOX['left'] + (ii - 1) * self.dx,
            Config.BOX['bottom'] + (jj - 1) * self.dy
        )
    
    def z_func_test(self, x : np.ndarray, y : np.ndarray) -> np.ndarray:
        # Specifying the initial position of the water front
        # A function describing the initial position of 
        # the water front in the shape of a circular arc 
        # $$ z(x,y) = x^2+y^2-0.015 $$ 
        # This can take array input
        init_front_hs = 0.1
        if self.simulation_type == 'rectilinear homogenous':
            return y - init_front_hs + 0.01 * (np.cos((80 * np.pi * x)))
        elif self.simulation_type == 'rectilinear heterogenous':
            return y - init_front_hs
        elif self.simulation_type == 'quarter-five spot heterogenous':
            return np.square(x) + np.square(y) - 0.015
        else:
            raise Exception('Inputted simulation type is invalid when performing z_func_test')
            
