from config import SimulationData, Config
import numpy as np
from typing import Tuple

###################################################
# This class handles the main simulation loop and # 
# coordinates the various components              #
# (grid, viscosity, permeability, etc.).          #
###################################################

class Simulation:
    def __init__(self, user_inputs : SimulationData):
        self.user_inputs = user_inputs
        self.grid = None
        self.viscosity = None
        self.permeability = None
        
        # main source terms
        self.f = None
        self.MFW = None
        self.src = None # magnitude of mass flow rate at source 
        
        # main rp variables
        self.s0 = None # initial residual water saturation inside the reservoir
        self.c0 = None # initial polymer concentration
        self.g0 = None # initial surfactant concentration
        self.c0_array = None
        self.UU = None # saturation matrix
        self.CC = None # polymer concentration matrix
        self.GG = None # surfactant concentration matrix
        self.interface = None
        self.miuo = None
        self.miuw = None
        self.beta1 = None
        self.miup = None
        self.miup_array = None
        self.swr0 = None 
        self.sor0 = None
        self.TL_mixing_flag = None # flag to enable todd-longstaff mixing parameter
        self.shear_flag = None # flag to enable shear effects with ad-hoc model from eclipse code
        self.t = None
        self.t_cal = None
        self.u = None
        self.v = None
        self.dt = None
        self.COC = None
        self.prod_rate = None
        self.CROIP = None
        self.conc_save = None
        self.sum_UU = None
        

    def run(self):
        self.initialize_grid()
        self.initialize_source_terms()
        self.initialize_permeability()
        self.initialize_rp_variables()
        # self.initialize_viscosity()
        
        # # Main simulation loop
        while self.should_continue():
            self.update_simulation()
        # self.plot_results()

    def initialize_grid(self):
        from models.grid import Grid
        self.grid = Grid(
            self.user_inputs.size_of_grid, 
            self.user_inputs.simulation_type
        )
        
    def initialize_source_terms(self):
        # Defining right hand side of elliptic system - source terms
        self.f = np.zeros((self.user_inputs.size_of_grid + 1, self.user_inputs.size_of_grid + 1))
        self.MFW = np.zeros((self.user_inputs.size_of_grid + 1, self.user_inputs.size_of_grid + 1))
        self.src = 120000
        
        # setting permeability state
        if self.user_inputs.simulation_type in ['rectilinear homogenous', 'rectilinear heterogenous']:
            self.f[ : , 0] = self.src # Intensity of Injection well
            self.f[ : , self.user_inputs.size_of_grid] = -self.src # intensity of production well
        elif self.user_inputs.simulation_type == 'quarter-five spot heterogenous':
            self.f[0,0] = self.src # intensity of injection well
            self.f[self.user_inputs.size_of_grid, self.user_inputs.size_of_grid] = -self.src # intensity of production well
        else:
            raise Exception('Inputted simulation type is invalid when performing initialize_source_terms')

    def s0c0(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        # function to initialize s,c,g in the domain
        s0 = np.zeros((self.user_inputs.size_of_grid + 1, self.user_inputs.size_of_grid + 1))
        c0 = np.zeros_like(s0)

        D = (self.grid.phi_test > 10**(-10)) + (np.abs(self.grid.phi_test) < 10**(-10))
        
        if not self.user_inputs.shear_thinning_on:
            g0 = np.array([])
            s0 = (~D) + D * (1 - self.s0)
            c0 = (~D) * self.c0
        else:
            s0 = (~D) + D * (1 - self.s0)
            c0 = (~D) * self.c0
            g0 = (~D) * self.g0
            
        return s0, c0, g0
    
    def initialize_rp_variables(self):
        self.s0 = 0.79 # initial residual water saturation inside the reservoir
        self.c0 = self.user_inputs.initial_polymer_concentration # initial polymer concentration
        self.g0 = self.user_inputs.initial_surfactant_concentration # initial surfactant concentration
        self.c0_array = self.c0 * np.ones((self.user_inputs.size_of_grid + 1, self.user_inputs.size_of_grid + 1))
        # UU - saturation matrix
        # CC - polymer concentration matrix
        # GG - surfactant concentration matrix
        self.UU, self.CC, self.GG = self.s0c0()
        self.interface = np.zeros((60,1))
        self.miuo = 10
        self.miuw = 1.26
        self.beta1 = 15000
        self.miup = self.miuw * (1 + self.beta1 * self.c0)
        self.miup_array = self.miup * np.ones((self.user_inputs.size_of_grid + 1, self.user_inputs.size_of_grid + 1))
        self.swr0 = 0.1
        self.sor0 = 0.3
        self.TL_mixing_flag = 0
        self.shear_flag = 0
        self.t = 0
        self.t_cal = 0
        self.u = np.zeros((self.user_inputs.size_of_grid + 1, self.user_inputs.size_of_grid + 1))
        self.v = np.zeros_like(self.u)
        self.dt = Config.DEFAULT_CFL * self.grid.dx / self.src
        if (self.user_inputs.simulation_type == 'quarter-five spot heterogenous'):
            self.dt = Config.DEFAULT_CFL * self.grid.dx / self.src * 100
        self.COC = np.zeros((1, 2000))
        self.prod_rate = np.zeros((1, int(np.floor(Config.DEFAULT_T_STOP / self.dt))))
        self.CROIP = np.zeros((1, int(np.floor(Config.DEFAULT_T_STOP / self.dt))))
        self.src_total = 0
        self.sum_UU = 0
    
    def initialize_viscosity(self):
        from models.viscosity import Viscosity
        self.viscosity = Viscosity(self.user_inputs.shear_thinning_on)

    def initialize_permeability(self):
        from models.permeability import Permeability
        self.permeability = Permeability(
            self.user_inputs.simulation_type, 
            self.user_inputs.size_of_grid, 
            self.grid.x, self.grid.y
        )

    def should_continue(self):
        return (self.t < Config.DEFAULT_T_STOP 
                and 
                self.UU[self.user_inputs.size_of_grid + 1,
                        self.user_inputs.size_of_grid + 1] <= 0.70
        )
        
    # Function to compute divergence using finite differences
    def divergence(U, V):
        dU_dx = np.gradient(U, axis=0)  # Derivative of U with respect to x (rows)
        dV_dy = np.gradient(V, axis=1)  # Derivative of V with respect to y (columns)
        return dU_dx + dV_dy
    
    def compvis(self):
        # function to compute viscosity of injected
        # displacing phase containing polymer
        # miuo = Displaced phase viscosity, c = polymer concentration
        gamma_dot = np.zeros
        if not self.user_inputs.shear_thinning_on:
            n, m = self.CC.shape
            if self.c0 == 0:
                miua = self.miuw * np.ones((n, m))
            else:
                miua = self.miuw * (1 + self.beta1 * self.CC)
        else:
            rho_water = 1000
            rho_xanthane = 1500
            w1 = rho_xanthane * self.CC
            w2 = rho_water * (1 - self.CC)
            wppm = ((w1) / (w1 + w2)) * 1e6
            w10 = rho_xanthane * self.c0_array
            w20 = rho_water * (1 - self.c0_array)
            wppm0 = ((w10) / (w10 + w20)) * 1e6
            
            if self.user_inputs.polymer_type == 'xanthane':
                eps_coeff = [1.15410398e-04, 2.04937780e+00]
                n_coeff = [ 3.05428284, -0.27294817]
            elif self.user_inputs == 'schizophyllan':
                eps_coeff = [0.03647214, 1.32175949]
                n_coeff = [ 4.86265534, -0.41570227]
                
            epsilon_0 = eps_coeff[0] * (wppm0 ** eps_coeff[1])
            power_n0 = np.min(n_coeff[0] * (wppm ** n_coeff[1]), 1)
            epsilon = eps_coeff[0] * (wppm ** eps_coeff[1])
            power_n = np.min(n_coeff[0] * (wppm ** n_coeff[0]), 1)
            
            n, m = self.CC.shape
            miua = self.miuw * np.ones((n, m))
            a1 = self.divergence(self.grid.x, self.v)
            a2 = self.divergence(self.grid.y, self.u)
            a3 = self.divergence(self.grid.x, self.u)
            a4 = self.divergence(self.grid.y, self.v)
            pi_D = np.abs(-0.25 * ((a1 + a2) ** 2) + a3 * a4)
            
            for ii in range(miua.shape[0]):
                for jj in range(miua.shape[1]):
                    if (self.CC[ii, jj] > 0):
                        gamma_dot[ii, jj] = 2 * np.sqrt(pi_D[ii, jj])
                        self.miup_array[ii, jj] = epsilon_0[ii, jj] * (gamma_dot[ii, jj] ** (power_n0[ii, jj] - 1))
                        
                        miua[ii, jj] = epsilon[ii, jj] * (gamma_dot[ii, jj] ** (power_n[ii, jj] - 1))
                        if miua[ii, jj] < self.miuw:
                            miua[ii, jj] = self.miuw
                        if miua[ii, jj] > 100:
                            miua[ii, jj] = 100
                        if self.miup_array[ii, jj] < self.miuw:
                            self.miup_array[ii, jj] = self.miuw
                        if self.miup_array[ii, jj] > 100:
                            self.miup_array[ii, jj] = 100
                            
            gamma_max = np.max(gamma_dot)
                            
        return miua, gamma_dot
            
    def TLmixing(self, miua, c):
        # implements the shear effects of polymer on the EOR simulation
        # Todd-Longstaff mixing model
        omega = 0
        n, m = miua.shape
        mu_w_eff = np.zeros((n, m))
        mu_p_eff = np.zeros((n, m))
        mu_w_e = np.zeros((n, m))
        C_bar = np.zeros((n, m))
        m_mu = np.zeros((n, m))
        for ii in range(n):
            for jj in range(m):
                if c[ii, jj] > 0:
                    mu_p_eff[ii, jj] = (miua[ii, jj] ** omega) * (self.miup_array[ii, jj] ** (1 - omega))
                    mu_w_e[ii, jj] = (miua[ii, jj] ** omega) * (self.miuw ** (1 - omega))
                    C_bar[ii, jj] = c[ii, jj] / self.c0
                    mu_w_eff[ii, jj] = mu_w_e[ii, jj] * mu_p_eff[ii, jj] / (
                                        (1 - C_bar[ii, jj]) * mu_p_eff[ii, jj] +
                                        C_bar[ii, jj] * mu_w_e[ii, jj]
                                    )
                    # calculating shear multipliers
                    m_mu[ii, jj] = (1 + self.beta1 * c[ii, jj])
                else:
                    mu_w_eff[ii, jj] = miua[ii, jj]
                    m_mu[ii, jj] = 1
                    
        return mu_w_eff, m_mu
    
    def compres(self, sigma, u, v, miua):
        # function to compute the residual saturations as a function of surfactant
        # concentration via capillary number variation. Hence it varies with change
        # in surf concentration and velocity evolution. Must be recomputed at everytime step.
        Nco0 = 1.44 * 10 ** (-4)
        Nca0 = 1.44 * 10 ** (-4)
        
        # compute capillary number
        nca = np.sqrt((u ** 2) + (v ** 2)) * miua / sigma
        nco = np.sqrt((u ** 2) + (v ** 2)) * self.miuo / sigma
        Nca = np.linalg.norm(nca) # compute 2-norm of Nc Matrix
        Nco = np.linalg.norm(nco)
        
        # define residual saturations as functions of capillary numbers
        if Nco < Nco0:
            sor = self.sor0
        else:
            sor = self.sor0 * (Nco0 / Nco) ** 0.5213
            
        if Nca < Nca0:
            swr = self.swr0
        else:
            swr = self.swr0 * (Nca0 / Nca) ** 0.1534
        
        return swr, sor
        
    def update_simulation(self):
        self.src_total = self.src_total + self.src
        self.t = self.t + self.dt
        inner_iter = 0
        epsilon = Config.DEFAULT_EPSILON
        sigma = 10.001 / (self.GG + 1) - 0.001 # this is surface tension
        miua, shear = self.compvis()
        
        if not self.user_inputs.shear_thinning_on:
            self.miup = np.max(miua[0, : ])
            self.miup_array = self.miup * np.ones(self.user_inputs.size_of_grid + 1, self.user_inputs.size_of_grid + 1)
            
        if self.TL_mixing_flag:
            mu_w_eff, m_mu = self.TLmixing(miua, self.CC)
            miua = mu_w_eff
        
        mShear, nShear = shear.shape
        
        # compute residual saturation s_{ra}, s_{ro}
        swr, sor = self.compres(sigma, self.u, self.v, miua)
        
        # recompute mobilities (with surfactants)
        
        
        
        

    def plot_results(self):
        from utils.plotting import plot_results
        plot_results(self)
