from typing import List
import json

class SimulationData:
    def __init__(
                 self, 
                 size_of_grid: int,
                 shear_thinning_on : bool, 
                 initial_polymer_concentration : float,
                 initial_surfactant_concentration : float,
                 plot_type : str,
                 polymer_type : str,
                 permeability_field_type : str
                ):
        self.size_of_grid = size_of_grid
        self.shear_thinning_on = shear_thinning_on
        self.initial_polymer_concentration = initial_polymer_concentration
        self.initial_surfactant_concentration = initial_surfactant_concentration
        self.plot_type = plot_type
        self.polymer_type = polymer_type
        self.permeability_type = permeability_field_type
        
    def __repr__(self):
        return (f"SimulationData(size_of_grid={self.size_of_grid}, "
                f"shear_thinning_on={self.shear_thinning_on}, "
                f"initial_polymer_concentration={self.initial_polymer_concentration}, "
                f"initial_surfactant_concentration={self.initial_surfactant_concentration}, "
                f"plot_type='{self.plot_type}', "
                f"polymer_type='{self.polymer_type}', "
                f"permeability_field_type='{self.permeability_type}')")

class InputFile:
    def __init__(self, file : str): 
        self.file= file
        self.simulations : List[SimulationData] = []
        self._get_json()
        
    def _get_json(self) -> None:
        with open(self.file, 'r') as input_file:
            input_simulations = json.load(input_file)
            for simulation in input_simulations:
                simulation_data = SimulationData(
                    size_of_grid=simulation['size_of_grid'],
                    shear_thinning_on=simulation['shear_thinning_on'],
                    initial_polymer_concentration=simulation['initial_polymer_concentration'],
                    initial_surfactant_concentration=simulation['initial_surfactant_concentration'],
                    plot_type=simulation['plot_type'],
                    polymer_type=simulation['polymer_type'],
                    permeability_field_type=simulation['permeability_field_type']
                )
                self.simulations.append(simulation_data)