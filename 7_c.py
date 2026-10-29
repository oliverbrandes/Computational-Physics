import numpy as np
import matplotlib.pyplot as plt
import integrator

def main():
    n = 128
    x_min, x_max = 0.0, 1.0
    delta_x = (x_max - x_min) / n
    x = x_min + (np.arange(n) + 0.5) * delta_x

    rho_max, v_max = 1.0, 1.0
    rho_left, rho_right = 0.6 * rho_max, 0.8 * rho_max
    x_c = 0.5 * (x_min + x_max)

    integ = integrator.Integrator(x, 1.0, 0.5, rho_max, v_max, fixed_boundaries=True)
    for i in range(n):
        integ.old.rho[i] = rho_left if x[i] < x_c else rho_right

    s = v_max * (1.0 - (rho_left + rho_right) / rho_max)
    print("analytical shock speed s =", s)

    t = 0.0
    for t_plot in [0.25, 0.5, 0.75]:
        while t < t_plot - 1e-12:
            integ.RK4()
            t += integ.delta_t
        x_shock = x_c + s * t
        rho_exact = np.where(x < x_shock, rho_left, rho_right)
        plt.plot(x, integ.old.rho.data, "o", markersize=3, label="numerical")
        plt.plot(x, rho_exact, "-", label="analytical")
        plt.title("t = {:.2f}".format(t))
        plt.xlabel("x"); plt.ylabel("rho")
        plt.legend(); plt.grid(True)
        plt.show()

if __name__ == '__main__':
    main()