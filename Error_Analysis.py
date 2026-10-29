import numpy as np
import matplotlib.pyplot as plt
import integrator


def run(n):
    x_max = 2.0
    x_min = -2.0
    delta_x = (x_max - x_min) / n
    x = np.zeros((n))
    for i in range(n):
        x[i] = x_min + (i + 0.5) * delta_x

    speed = 1.0
    courant = 0.5
    integrate = integrator.Integrator(x, speed, courant)

    amp = 1.0
    width = x_max / 4
    integrate.initialize(amp, width)

    final_time = 4.0
    n_plot = 10**9
    time, error = integrate.integrate(final_time, n_plot)
    return time, error

def main():
    for n in [64, 128, 256]:
       time, error = run(n)
       plt.plot(time, error, label = "N = " + str(n))

    plt.xlabel("t")
    plt.ylabel("error")
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()

