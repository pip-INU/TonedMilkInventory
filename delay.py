import numpy as np
def generate_random_delays(n, max_delay=1.0, small_delay_prob=0.1):
    high_delay_prob = 0.06
    delays = []
    for _ in range(n):
        r = np.random.rand()
        if r < high_delay_prob:
            delays.append(max_delay)
        elif r < high_delay_prob + small_delay_prob:
            delays.append(np.random.uniform(0.2, 0.5))
        else:
            delays.append(0)
    return delays