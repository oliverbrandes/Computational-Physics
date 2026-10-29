import numpy as np
import matplotlib.pyplot as plt
import integrator

def main():
    n = 128
    x_max = 1.0
    x_min = 0.0
    delta_x = (x_max - x_min) / n
    x = np.zeros((n))
    for i in range(n):
        x[i] = x_min + (i + 0.5) * delta_x

    speed = 1.0
    courant = 0.5
    rho_max = 1.0
    v_max = 1.0
    integrate = integrator.Integrator(x, speed, courant, rho_max, v_max)

    rho_floor = 0.3
    amp = 0.1
    width = 0.125
    center = 0.5
    integrate.initialize(rho_floor, amp, width, center)

    final_time = 2.0
    n_plot = 10
    integrate.integrate(final_time, n_plot)


  
if __name__ == '__main__':
    main()
        
