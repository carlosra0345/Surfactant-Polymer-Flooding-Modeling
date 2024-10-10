from typing import List
from config import SimulationData
import json

################################################################
# This utility handles user inputs and settings configuration. #
################################################################
        
def get_user_input(file_name : str) -> List[SimulationData]:
    simulations : List[SimulationData] = []
    with open(file_name, 'r') as input_file:
                input_simulations = json.load(input_file)
                for simulation in input_simulations:
                    simulation_data = SimulationData(
                        size_of_grid=simulation['size_of_grid'],
                        shear_thinning_on=simulation['shear_thinning_on'],
                        initial_polymer_concentration=simulation['initial_polymer_concentration'],
                        initial_surfactant_concentration=simulation['initial_surfactant_concentration'],
                        plot_type=simulation['plot_type'],
                        polymer_type=simulation['polymer_type'],
                        simulation_type=simulation['simulation_type']
                    )
                    simulations.append(simulation_data)
    return simulations
