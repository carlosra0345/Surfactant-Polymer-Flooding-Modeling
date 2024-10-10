#################################################
# This file contains configuration settings and #
# constants used across the simulation.         #
#################################################

class SimulationData:
    def __init__(
        self,
        size_of_grid : int,
        shear_thinning_on : bool,
        initial_polymer_concentration : float,
        initial_surfactant_concentration : float,
        plot_type : str,
        polymer_type : str,
        simulation_type : str
    ):
        self.size_of_grid = size_of_grid
        self.shear_thinning_on = shear_thinning_on
        self.initial_polymer_concentration = initial_polymer_concentration
        self.initial_surfactant_concentration = initial_surfactant_concentration
        self.plot_type = plot_type
        self.polymer_type = polymer_type
        self.simulation_type = simulation_type
        
class Config:
    BOX = {
        "left": 0,
        "right": 1,
        "bottom": 0,
        "top": 1,
    }
    DEFAULT_C0 = 0.001
    DEFAULT_G0 = 0.1
    DEFAULT_CFL = 1
    DEFAULT_T_STOP = 500
    DEFAULT_EPSILON = 10
    # Add other constants and configuration options as needed
