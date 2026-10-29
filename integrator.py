import numpy as np
import state
import gridfunction

class Integrator:

    def __init__(self, x, v, courant, rho_max, v_max, fixed_boundaries=False):
        self.rho_max = rho_max
        self.v_max = v_max
        self.old = state.State(x, "old")
        self.k1 = state.State(x, "k1")
        self.k2 = state.State(x, "k2")
        self.k3 = state.State(x, "k3")
        self.k4 = state.State(x, "k4")
        self.inter = state.State(x, "intermediate")
        self.delta_x = x[1] - x[0]
        self.delta_t = courant * self.delta_x
        self.n = np.size(x)
        self.v = v
        self.x = x
        self.fixed_boundaries = fixed_boundaries   # for part (c)
        # grid functions
        self.rho_L = gridfunction.GridFunction(x, "rho_L")
        self.rho_R = gridfunction.GridFunction(x, "rho_R")
        self.f = gridfunction.GridFunction(x, "f")

    def initialize(self, rho_floor, amp, width, center):
        for i in range(self.n):
            temp = (self.x[i] - center) / width
            self.old.rho[i] = rho_floor + amp * np.exp(-temp**2)
        self.old.plot(0.0)

    def integrate(self, t_max, n_plot):
        t = 0.0
        n_step = 0
        while (t < t_max):
            self.RK4()
            t += self.delta_t
            n_step += 1
            if (n_step % n_plot == 0):
                self.old.plot(t)

    #HRSC
    def flux(self, rho):
        # traffic flux f = rho * v_max * (1 - rho/rho_max)
        return rho * self.v_max * (1.0 - rho / self.rho_max)

    def minmod(self, a, b):
        # opposite signs (local max/min): flat, so no new extrema
        if a * b <= 0.0:
            return 0.0
        # same sign: take the smaller slope
        elif abs(a) < abs(b):
            return a
        else:
            return b

    def reconstruct(self, rho):
        for i in range(self.n):
            im1 = rho.ind(i - 1)
            ip1 = rho.ind(i + 1)
            # limited slope times delta_x
            slope = self.minmod(rho[i] - rho[im1], rho[ip1] - rho[i])
            # left edge of cell i = right side of interface i - 1/2
            self.rho_R[i] = rho[i] - 0.5 * slope
            # right edge of cell i = left side of interface i + 1/2
            self.rho_L[ip1] = rho[i] + 0.5 * slope

    def compute_fluxes(self):
        # Godunov flux: exact Riemann solution at each interface i - 1/2
        rho_star = 0.5 * self.rho_max  # peak of f (sonic point)
        for i in range(self.n):
            rL = self.rho_L[i]
            rR = self.rho_R[i]
            if rL <= rR:
                # shock case: minimum of f on [rL, rR]
                self.f[i] = min(self.flux(rL), self.flux(rR))
            else:
                # rarefaction case: maximum of f on [rR, rL]
                if rR <= rho_star <= rL:
                    self.f[i] = self.flux(rho_star)
                else:
                    self.f[i] = max(self.flux(rL), self.flux(rR))

    def compute_derivs(self, current, derivs):
        rho = current.rho
        self.reconstruct(rho)
        self.compute_fluxes()
        for i in range(self.n):
            ip1 = rho.ind(i + 1)
            # cars in through left face minus cars out through right face
            derivs.rho[i] = -(self.f[ip1] - self.f[i]) / self.delta_x
        if self.fixed_boundaries:
            # part (c): boundary points are not updated
            derivs.rho[0] = 0.0
            derivs.rho[self.n - 1] = 0.0

    def RK4(self):
        # step 1:
        self.compute_derivs(self.old, self.k1)
        # step 2:
        self.inter.equals(self.old)
        self.inter.plus(self.delta_t / 2., self.k1)
        self.compute_derivs(self.inter, self.k2)
        # step 3:
        self.inter.equals(self.old)
        self.inter.plus(self.delta_t / 2., self.k2)
        self.compute_derivs(self.inter, self.k3)
        # step 4:
        self.inter.equals(self.old)
        self.inter.plus(self.delta_t, self.k3)
        self.compute_derivs(self.inter, self.k4)
        #
        self.old.plus(self.delta_t / 6., self.k1)
        self.old.plus(self.delta_t / 3., self.k2)
        self.old.plus(self.delta_t / 3., self.k3)
        self.old.plus(self.delta_t / 6., self.k4)