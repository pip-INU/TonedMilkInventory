import numpy as np
def safe_exp(x, max_value=700):
    """Prevent overflow in exponential calculations"""
    return np.exp(np.clip(x, -max_value, max_value))

def safe_log(x, epsilon=1e-20):
    """Prevent log(0) errors"""
    return np.log(np.maximum(x, epsilon))