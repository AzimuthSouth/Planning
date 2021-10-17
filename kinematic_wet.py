import algebra
import numpy


class plate_wet():
    def set_data(self, hksi, bt):
        """
        set data for plate geometry and depth
        """
        self.h = hksi
        self.beta = bt


    def calc_deadrise_koeff(self, ksi):
        bt = algebra.get_value(self.beta, ksi)
        self.k_beta = numpy.pi / 2 * (numpy.pi / 2 / bt - 1)**2
        self.LL = 2 - numpy.cos(bt)
        return [bt, self.k_beta, self.LL]


    def h_ksi(self, ksi):
        """
        osadka in current point ksi
        ksi - current coordinate along kil line,
        """
        h = algebra.get_value(self.h, ksi)
        return h

    def V_n(self, V, ksi, phi):
        return -V * phi

    def a_n(self, V, ksi, phi):
        pass

    def kil_line_straight(self):
        """
        kil line at straight orientation phi=0
        """
        p1 = [self.xycm[0] - self.ki[0], self.xycm[1] - self.ki[1]]
        p2 = [self.xycm[0] - self.ki[0] + self.l, self.xycm[1] - self.ki[1]]
        return [p1, p2, self.xycm]

    def kil_line_rotate(self, phi):
        """
        kil line rotate at angle phi
        phi - diferent angle in radians
        """
        p1, p2, _ = self.kil_line_straight()
        a = algebra.calc_transform_matrix_XY(phi)
        shift = [self.xycm[0], self.xycm[1], 0.0]
        pp1 = [p1[0] - shift[0], p1[1] - shift[1], 0.0]
        p1_ = algebra.new_coord(pp1, a) + shift
        pp2 = [p2[0] - shift[0], p2[1] - shift[1], 0.0]
        p2_ = algebra.new_coord(pp2, a) + shift
        #print(f"straight p1={p1}, p2={p2}")
        #print(f"rotate p1={p1_}, p2={p2_}")
        return [p1_, p2_, self.xycm]


    def YM_hydrostatic_numerical(self, phi, rho):
        g = 9.81 # m * s^-2
        ksi_0 = self.h[0][0]
        ksi_1 = self.h[-1][0]
        def f(x):
            beta = algebra.get_value(self.beta, x)
            return (self.h_ksi(x))**2 * rho * g / numpy.tan(beta)
        fy = algebra.kvadratura(f, ksi_0, ksi_1)
        def f(x):
            beta = algebra.get_value(self.beta, x)
            return (self.h_ksi(x))**2 * x * rho * g / numpy.tan(beta)
        mz = algebra.kvadratura(f, ksi_0, ksi_1)
        return [fy, mz]

    def khi_y_khi_m(self, phi, ksi):
        beta = algebra.get_value(self.beta, ksi)
        khi_y = 1 - numpy.pi * phi / 4 / numpy.tan(beta)
        khi_m = khi_y * (1 + 1.0 / (1 + (1 + 2 * numpy.tan(beta) / numpy.pi / phi)**0.5))
        return [khi_y, khi_m]

    def YM_hydrodynamic_numerical(self, phi, rho, V):
        ksi_0 = self.h[0][0]
        ksi_1 = self.h[-1][0]
        #for hydrodynamic force
        def f(x):
            khi_y, _ = self.khi_y_khi_m(phi, x)
            _, k_beta, _ = self.calc_deadrise_koeff(x)
            res = 2 * self.h_ksi(x) * (self.V_n(V, x, phi))**2 * rho * k_beta * khi_y
            #return 2 * self.h_ksi(x) * (self.V_n(V, x, phi))**2 * rho * k_beta * khi_y
            return res
        fy = algebra.kvadratura_N(f, ksi_0, ksi_1, 10)
        print("moment")
        #for hydrodynamic moment
        def f(x):
            _, khi_m = self.khi_y_khi_m(phi, x)
            _, k_beta, _ = self.calc_deadrise_koeff(x)
            return 2 * self.h_ksi(x) * (self.V_n(V, x, phi))**2 * rho * k_beta * x * khi_m
        mz = algebra.kvadratura_N(f, ksi_0, ksi_1, 10)
        #for shear moment
        cf = 0.001 #???
        print("shear")
        def f(x):
            beta = algebra.get_value(self.beta, x)
            return (0.14 - numpy.pi / 4 * self.h_ksi(x)) * self.h_ksi(x) * -rho * numpy.pi / 2 * cf  / numpy.sin(beta) * V**2
        mtr = algebra.kvadratura(f, ksi_0, ksi_1)

        return [fy, mz, mtr]
