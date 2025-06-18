import numpy as np
import pandas as pd
from scipy.integrate import quad
from math import log
from tqdm import tqdm


def safe_log(x):
    if x <= 0:
        return -np.inf
    return np.log(x)

class DairySupplyChainOptimizer:
    def __init__(self, params):
        """Initialize with all parameters"""
        self.params = params
        
        # frequently used parameters
        self.delta_1 = params['delta_1']
        self.delta_2 = params['delta_2']
        self.delta_3 = params['delta_3']
        self.theta   = params['theta']

        self.K =  ((1 - params['sigma']) * (1 - params['gamma']) * (1 - params['alpha']) + 
                   params['gamma'] * params['beta']) * (1 - params['w_z']) * (1 - params['w_f'])

        self.H = ((1 - params['beta']) * params['gamma'] + 
                  params['alpha'] * (1 - params['gamma'])) * (1 - params['w_z']) * (1 - params['w_f'])

        self.P0  = params['P0']
        self.R0  = params['R0']
        self.D0  = params['D0']
        self.lmbd = params['lmbd']
        self.eeta = params['eeta']

        self.w   =  params['w_z'] * (1 - params['w_f']) + params['w_f']
        self.w_r = np.sum(params['w_ri'])
        self.mu0 = params['mu0']
        self.T   = params['T']
        self.n   = params['n']

    # Solve Time Points
    def solve_T5(self, T1, T4):
        term1 = (self.H * self.P0 / self.delta_2) * (np.exp(-self.delta_2 * T4) - np.exp(-self.delta_2 * T1))
        term2 = (self.delta_3 / self.R0) * term1 + np.exp(-self.delta_3 * T4)
        return -safe_log(term2) / self.delta_3
    
    def solve_T2(self, T1):
        numerator   = self.K * self.P0 * (np.exp(-self.delta_2 * T1) - 1)
        denominator = self.delta_2 * (self.D0 * self.lmbd - self.K * self.P0)
        if abs(denominator) < 1e-10:
            raise ValueError("Denominator is zero — no feasible T2 solution")
        return numerator / denominator
    
    def solve_T3(self, T1, T4):
        T5 = self.solve_T5(T1, T4)
        T2 = self.solve_T2(T1)

        part1 = (self.eeta * self.R0) / (self.lmbd * (self.delta_3 - self.theta))
        part2 = self.D0 / (self.delta_1 + self.theta)
        term1 = (part1 / part2) * (np.exp(-(self.delta_3 - self.theta) * T4 - self.theta * T5) - np.exp(-self.delta_3 * T5))

        part3 = (self.K * self.P0) / (self.lmbd * (self.delta_2 - self.theta))
        term2 = (part3 / part2) * (np.exp(-(self.delta_2 - self.theta) * T2 - self.theta * T5) - 
                                   np.exp(-(self.delta_2 - self.theta) * T4 - self.theta * T5))

        term3 = (np.exp((self.delta_1 + self.theta) * T2 - self.theta * T5) - 
                 np.exp((self.delta_1 + self.theta) * self.T - self.theta * T5) +
                 np.exp(self.delta_1 * T5))

        log_arg = term1 + term2 + term3
        return (self.theta * T5 + safe_log(log_arg)) / (self.delta_1 + self.theta)


    #Solution of Differential Equations of the Inventories
    def I_s0(self, t):
        return (self.D0 / self.delta_1) * (1 - np.exp(self.delta_1 * t))
    def I_s1(self, t, T1):
        T2 = self.solve_T2(T1)
        return ((self.K * self.P0) / (self.lmbd * self.delta_2)) * (np.exp(-self.delta_2 * T2) - np.exp(-self.delta_2 * t)) + \
               (self.D0 / self.delta_1) * (np.exp(self.delta_1 * T2) - np.exp(self.delta_1 * t))
    def I_s2(self, t, T1):
        T2 = self.solve_T2(T1)
        return ((self.K * self.P0) / (self.lmbd * (self.delta_2 - self.theta))) * \
               (np.exp(-(self.delta_2 - self.theta) * T2 - self.theta * t) - np.exp(-self.delta_2 * t)) + \
               (self.D0 / (self.delta_1 + self.theta)) * \
               (np.exp((self.delta_1 + self.theta) * T2 - self.theta * t) - np.exp(self.delta_1 * t))
    def I_s3(self, t, T1, T4):
        T2 = self.solve_T2(T1)
        T3 = self.solve_T3(T1, T4)
        return ((self.K * self.P0) / (self.lmbd * (self.delta_2 - self.theta))) * \
               (np.exp(-(self.delta_2 - self.theta) * T2 - self.theta * t) - np.exp(-self.delta_2 * t)) + \
               (self.D0 / (self.delta_1 + self.theta)) * \
               (np.exp((self.delta_1 + self.theta) * T2 - self.theta * t) - np.exp((self.delta_1 + self.theta) * T3 - self.delta_1 * t))
    def I_s4(self, t, T1, T4):
        T2 = self.solve_T2(T1)
        T3 = self.solve_T3(T1, T4)
        return ((self.eeta * self.R0) / (self.lmbd * (self.delta_3 - self.theta))) * \
               (np.exp(-(self.delta_3 - self.theta) * T4 - self.theta * t) - np.exp(-self.delta_3 * t)) + \
               ((self.K * self.P0) / (self.lmbd * (self.delta_2 - self.theta))) * \
               (np.exp(-(self.delta_2 - self.theta) * T2 - self.theta * t) - np.exp(-(self.delta_2 - self.theta) * T4 - self.theta * t)) + \
               (self.D0 / (self.delta_1 + self.theta)) * \
               (np.exp((self.delta_1 + self.theta) * T2 - self.theta * t) - np.exp((self.delta_1 + self.theta) * T3 - self.delta_1 * t))
    def I_s5(self, t):
        return (self.D0 / (self.delta_1 + self.theta)) * \
               (np.exp((self.delta_1 + self.theta) * self.T - self.theta * t) - np.exp(self.delta_1 * t))

    def I_w0(self, t, T1):
        return ((self.w * self.P0) / self.delta_2) * (np.exp(-self.delta_2 * T1) - np.exp(-self.delta_2 * t)) + \
               self.w_r * (t - T1)
    def I_w1(self, t, T1, T4):
        T3 = self.solve_T3(T1, T4)
        return ((self.w * self.P0) / self.delta_2) * (np.exp(-self.delta_2 * T1) - np.exp(-self.delta_2 * t)) + \
               self.w_r * (T3 - T1)
    def I_w2(self, t, T1, T4):
        T3 = self.solve_T3(T1, T4)
        return (((1 - self.eeta) * self.R0) / (self.lmbd * self.delta_3)) * \
               (np.exp(-self.delta_3 * T4) - np.exp(-self.delta_3 * t)) + \
               ((self.w * self.P0) / self.delta_2) * \
               (np.exp(-self.delta_2 * T1) - np.exp(-self.delta_2 * T4)) + \
               self.w_r * (T3 - T1)
    
    def I_r0(self, t, T1):
        return (self.mu0 / 2) * (t**2 - T1**2) + \
               (self.P0 / self.delta_2) * (np.exp(-self.delta_2 * t) - np.exp(-self.delta_2 * T1))
    def I_r1(self, t, T4):
        return (self.P0 / self.delta_2) * (np.exp(-self.delta_2 * t) - np.exp(-self.delta_2 * T4))

    def I_rw0(self, t, T1):
        return (self.H * self.P0 / self.delta_2) * (np.exp(-self.delta_2 * T1) - np.exp(-self.delta_2 * t))
    def I_rw1(self, t, T1, T4):
        T5 = self.solve_T5(T1, T4)
        return (self.R0 / self.delta_3) * (np.exp(-self.delta_3 * t) - np.exp(-self.delta_3 * T5))
    

    def calculate_total_emission(self,e_p,e_rp,e_c,e_w,e_dump,e_t,total_waste,unuse_waste_generate,T1,T4,x):
            T3 = self.solve_T3(T1,T4)
            T5 = self.solve_T5(T1,T4)
            # Process Emissions
            epsilon_p = (e_p * (self.P0 * np.exp(-self.delta_2 * (T4 - T1)) + e_rp * self.R0 * np.exp(-self.delta_3 * (T5 - T4))))
            
            # Cold Storage Emissions
            epsilon_c = e_c * ((T3-T1) + (self.T-T1) + (T4-T1) + (self.T-T1))
            
            # Waste Processing Emissions
            epsilon_w = e_w * total_waste + e_dump * max(0, unuse_waste_generate)
            
            # Transportation Emissions
            epsilon_t = e_t * x

            return epsilon_p + epsilon_c + epsilon_w + epsilon_t

    # Cost Function
    def calculate_total_profit(self, T1, T4, x):
        try:
            T2 = self.solve_T2(T1)
            T3 = self.solve_T3(T1, T4)
            T5 = self.solve_T5(T1, T4)
            p = self.params
    

            # Milk procurement and processing
            total_milk_procured = self.mu0 * (T3 - T1)
            procurement_cost = p['CP'] * total_milk_procured
            processing_cost = p['c_p'] * total_milk_procured

            # Reprocessing
            total_milk_to_reprocess, _ = quad(lambda t: self.I_rw0(t, T1), T1, T4)
            reprocess_cost = p['c_r'] * total_milk_to_reprocess

            # Holding costs for serviceable milk
            number_s1, _ = quad(lambda t: self.I_s1(t, T1), T1, T2)
            number_s2, _ = quad(lambda t: self.I_s2(t, T1), T2, T3)
            number_s3, _ = quad(lambda t: self.I_s3(t, T1, T4), T3, T4)
            number_s4, _ = quad(lambda t: self.I_s4(t, T1, T4), T4, T5)
            number_s5, _ = quad(lambda t: self.I_s5(t), T5, self.T)
            total_serviceable_milk_held = number_s1 + number_s2 + number_s3 + number_s4 + number_s5
            holding_cost_s = p['c_s'] * total_serviceable_milk_held

            # Holding costs for raw milk
            liters_r0, _ = quad(lambda t: self.I_r0(t, T1), T1, T3)
            liters_r1, _ = quad(lambda t: self.I_r1(t, T4), T3, T4)
            holding_cost_r = p['c_raw'] * (liters_r0 + liters_r1)

            # Holding cost for milk to rework
            numbers_to_rework, _ = quad(lambda t: self.I_rw0(t, T1), T1, T3)
            holding_cost_rw = p['c_rw'] * numbers_to_rework

            # Waste costs
            waste_w0, _ = quad(lambda t: self.I_w0(t, T1), T1, T3)
            waste_w1, _ = quad(lambda t: self.I_w1(t, T1, T4), T3, T4)
            waste_w2, _ = quad(lambda t: self.I_w2(t, T1, T4), T4, T5)
            total_waste = waste_w0 + waste_w1 + waste_w2
            holding_cost_w = p['c_w'] * total_waste

            # Transportation & vehicle costs
            transportation_cost = p['c_t'] * p['c'] * x

            #Fixed Cost
            setup_cost = p['setup_p'] + p['setup_w'] + p['setup_cs'] * 4
            packaging_cost = p['c_pack'] * (self.K * self.P0 * np.exp(-self.delta_2 * (T4 - T1)) +
                                        (p['eeta'] * self.R0 * np.exp(-self.delta_3 * (T5 - T4)) / p['lmbd']))

            # Revenues
            r_1 = p['sp'] * self.D0 * (np.exp(self.delta_1 * (T3 - T1)) + np.exp(self.delta_1 * (self.T - T5)))
            revenue_from_solids = p['w_f'] * (self.P0 * np.exp(-self.delta_2 * (T4 - T1)) +
                                            self.R0 * np.exp(-self.delta_2 * (T5 - T4))) * p['r_f']
            revenue_from_cream = p['w_z'] * (self.P0 * np.exp(-self.delta_2 * (T4 - T1)) +
                                            self.R0 * np.exp(-self.delta_2 * (T5 - T4))) * p['r_s']
            revenue_from_return = p['sp_re'] * p['omega']*self.w_r
            total_revenue = r_1 + revenue_from_solids + revenue_from_cream + revenue_from_return
                
            unuse_waste_generate = total_waste + (p['omega']+1)*self.w_r - ((p['w_f'] + p['w_z'])* (self.P0 * np.exp(-self.delta_2 * (T4 - T1)) + self.R0 * np.exp(-self.delta_2 * (T5 - T4))))
            penalty_unuse_waste = p['p_w'] * unuse_waste_generate
            penalty_backorder = p['p_b'] * number_s1
            
            e_p = p['e_p']
            e_rp = p['e_rp']
            e_c = p['e_c']
            e_w = p['e_w']
            e_dump = p['e_dump']
            e_t = p['e_t']
            total_emission = self.calculate_total_emission(e_p,e_rp,e_c,e_w,e_dump,e_t,total_waste,unuse_waste_generate,T1,T4,x)
            excess_emission = total_emission - p['emission_cap']
            carbon_balance = excess_emission * p['credit_price']
                
            # Total Profit
            total_profit = total_revenue-(procurement_cost + processing_cost + reprocess_cost +
                        holding_cost_s + holding_cost_r + holding_cost_rw + holding_cost_w +
                        transportation_cost + penalty_unuse_waste + penalty_backorder + setup_cost + packaging_cost + 
                        carbon_balance)
            
            return {
                'ProcessingStartTime': T1, 'ProcessingEndTime': T4, 'NumberOfDispatches': x,
                'BackOrderReplenishTime': T2,'MorningDispatchEndTime': T3,'EveningDispatchStartTime': T5,
                'Procurement': procurement_cost, 'Processing': processing_cost,
                'HoldingService': holding_cost_s, 'HoldingRaw': holding_cost_r,
                'HoldingRework': holding_cost_rw, 'HoldingWaste': holding_cost_w,
                'Transport': transportation_cost, 'Setup': setup_cost, 'Packaging': packaging_cost,
                'Total Waste':total_waste,'Emissions': total_emission, 'CarbonBalance': carbon_balance,
                'Revenue': total_revenue, 
                'TotalProfit': total_profit,
            }

        except ValueError as e:
            print(f"Error in calculation: {str(e)}")
            return np.inf


if __name__ == "__main__":
    # Dummy parameter 
    n = int(input("Enter the number of retailers you have?"))

    params = {
        'n':n,
        # Decay and deterioration rates
        'delta_1': 0.03, 'delta_2': 0.05, 'delta_3': 0.04, 'theta': 0.01,

        # Initial inventory and production parameters
        'P0': 100, 'R0': 20, 'D0': 50, 'lmbd': 0.8, 'eeta': 0.9, 'mu0': 200, 'T': 24, 'n': n,

        # Waste proportions and rates
        'sigma': 0.02, 'gamma': 0.03, 'alpha': 0.01, 'beta': 0.04,
        'w_z': 0.05, 'w_f': 0.07, 'w_ri':  np.clip(np.random.normal(0.03, 0.005, size=n), 0.01, 0.05),  # individual retailer waste factors list

        # Cost parameters (INR per unit or per operation)
        'CP': 35,        # procurement cost per liter
        'c_p': 4,        # processing cost per liter
        'c_r': 6,        # reprocessing cost per liter
        'c_s': 0.5,      # holding cost for serviceable milk per liter per time
        'c_raw': 0.4,    # holding cost for raw milk per liter per time
        'c_rw': 0.6,     # holding cost for reworkable milk
        'c_w': 0.8,      # holding cost for waste
        'c_t': 10,       # transportation cost per trip
        'c': 1.5,        # cost per vehicle
        'c_pack': 2,     # packaging cost per liter of final product

        # Setup Costs (Fixed Costs)
        'setup_p': 1000, # processing plant setup
        'setup_w': 800,  # waste processing unit setup
        'setup_cs': 500, # cold storage unit setup

        # Revenue rates
        'sp': 55,        # selling price per liter serviceable milk
        'sp_re': 40,     # resale price per liter of returned milk
        'omega': 0.85,   # fraction of resale returns
        'r_f': 30,       # revenue per liter of fat by-product
        'r_s': 25,       # revenue per liter of solids by-product

        # Penalty and incentive rates
        'p_w': 15,       # penalty per liter unutilized waste
        'p_b': 5,
        # Emissions parameters
        'e_p': 0.3,      # emissions per liter production
        'e_rp': 0.6, # emissions multiplier for reprocessed milk
        'e_c': 0.02,     # emissions per unit time of cold storage
        'e_w': 0.1,      # emissions per liter waste processing
        'e_dump': 0.15,  # emissions per liter dumped waste
        'e_t': 0.05,     # emissions per transport trip

        # Carbon market parameters
        'emission_cap': 500,    # emission limit before penalty applies
        'credit_price': 500,    # price per emission unit over the cap
    }

    optimizer = DairySupplyChainOptimizer(params)
    results=[]
    for T1 in tqdm(range(0, 13)):
        for T4 in range(12, 25):
            for x in range(2, 21):
                res = optimizer.calculate_total_profit(T1, T4, x)
                if res:
                    results.append(res)

    df = pd.DataFrame(results)
    print(df.head())
    df.to_csv("ScenarioBackorderResults.csv", index=False)
