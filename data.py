from dataclasses import dataclass
import numpy as np
@dataclass
class ProblemParameters:
        n = int(input("Enter number of retailers: "))
    # Decay and deterioration rates (tuned for sensitivity)
        delta_1: float = 0.05,   # Demand growth rate (higher = steeper cup sides)
        delta_2: float = 0.04,   # Production decay (higher = sharper peak)
        delta_3: float = 0.06,   # Reprocessing decay
        theta: float =  0.03,     # Deterioration rate

        # Initial inventory (balanced for typical operations)
        P0: float = 120,         # Initial production
        R0: float = 35,          # Initial reworkable stock
        D0: float = 80,          # Initial demand
        lmbd: float = 0.5,      # Liters per packet
        eeta: float = 0.92,      # Reprocessing efficiency
        mu0: float = 250,        # Procurement rate
        T: float = 24,           # Time horizon

        # Waste proportions and rates
        sigma: float = 0.03, 
        gamma: float = 0.04,
        alpha: float = 0.02,
        beta: float = 0.05,
        w_z: float = 0.06,
        w_f: float = 0.08,
        w_ri: float = np.clip(np.random.normal(0, 3, size=n), 0.02, 0.06),# individual retailer waste factors list
        #Vehicle Constraint
        maxQ: float = 1200, 
        minQ: float = 900,
        # Cost parameters (in Rs.)
        CP: float = 35,        # procurement cost per liter
        c_p: float = 2,        # processing cost per liter per hr
        c_r: float = 1,        # reprocessing cost per liter per hr
        c_s: float = 0.15,      # holding cost for serviceable milk per liter per time
        c_raw: float = 0.04,    # holding cost for raw milk per liter per time
        c_rw: float = 0.05,     # holding cost for reworkable milk
        c_w: float = 0.06,      # holding cost for waste
        c_t: float = 5,       # transportation cost per packet per dispatch
        c_pack: float = 2,     # packaging cost per liter of final product
        #Fixed cost
        setup_p: float = 1000, # processing plant setup
        setup_w: float = 2000,  # waste processing unit setup
        setup_cs: float = 1000, # cold storage unit setup
        # Revenue rates
        sp: float = 65,        # selling price per liter serviceable milk
        sp_re: float = 45,     # resale price per liter of returned milk
        omega: float = 0.85,   # fraction of resale returns
        r_f: float = 100,       # revenue per liter of fat by-product
        r_s: float = 50,       # revenue per liter of solids by-product
        # Penalty
        p_w: float = 2,       # penalty on unutilized waste
        p_b: float = 5,        # penalty on backorder
        p_delay: float =0.1,
        # Emissions parameters
        e_p: float = 0.3,      # emissions per liter production
        e_rp: float = 0.07, # emissions multiplier for reprocessed milk
        e_c: float = 0.02,     # emissions per unit time of cold storage
        e_w: float = 0.1,      # emissions per liter waste processing
        e_dump: float = 0.15,  # emissions per liter dumped waste
        # Carbon market parameters
        emission_cap: float = 800,    # emission limit before penalty applies
        credit_price: float = 10,    # price per emission unit over the cap
