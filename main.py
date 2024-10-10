########################################################
# The entry point of your program where the simulation #
# is initiated and user inputs are handled.            #
########################################################

from models.simulation import Simulation
from utils.input_handler import get_user_input

def main():
    user_inputs = get_user_input('simulation_input.json')
    for simulation_input in user_inputs:
        try:         
            simulation = Simulation(simulation_input)
            simulation.run()
        except Exception as e:
            print(f'An error has occured: {e}')
            
        
if __name__ == "__main__":
    main()
    
# FINISHED ON LINE 202 on master_surf_grid on compvis function
# started working on DIVERGENCE function in compvis function