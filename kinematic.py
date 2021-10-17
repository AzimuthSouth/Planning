import algebra
import wave
import numpy


class plate():
    def set_data(self, cm_xy, ki, len):
        """
        set data for plate cm in global and local coordinate systems,
        cm_xy - center mass
        ki = distance from cm to trance, to kil
        len - plate length along kil line
        """
        self.xycm = cm_xy
        self.ki = ki
        self.l = len

    def set_deadrise_angle(self, bt):
        """
        bt - deadrise angle in radians
        k_beta - Vagner coefficient
        LL - Longinovich coefficient
        """
        self.beta = bt
        self.k_beta = numpy.pi / 2 * (numpy.pi / 2 / bt - 1)**2
        # print(f"beta={self.beta}, k_beta={self.k_beta}")
        self.LL = 2 - numpy.cos(bt)

    def h_ksi(self, ksi, phi, wave=None, t=0):
        """
        osadka in current point ksi
        ksi - current coordinate along kil line,
        phi - different in radians
        """
        y_wave = 0
        if wave:
            x = (ksi - self.ki[0]) * numpy.cos(phi)
            y_wave = wave.y_wave(x, t)
        return self.ki[1] - self.xycm[1] - ksi * phi + y_wave

    def h_ksi_move(self, ksi, phi, y_c=None, wave=None, t=0):
        """
        osadka in current point ksi
        ksi - current coordinate along kil line,
        phi - different in radians
        y_c - current cm position
        """
        y_wave = 0
        if wave:
            x = (ksi - self.ki[0]) * numpy.cos(phi)
            y_wave = wave.y_wave(x, t)
        yc = self.xycm[1]
        if y_c:
            yc = y_c
        return self.ki[1] - yc - ksi * phi + y_wave


    def V_n(self, V, ksi, phi, yc_t = 0, phi_t = 0, wave=None, t=0):
        v_wave = 0
        if wave:
            x = (ksi - self.ki[0]) * numpy.cos(phi)
            v_wave = wave.d_y_wave_dt(x, t)
        return -V * phi - v_wave + yc_t + ksi * phi_t


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

    def profil_rotate(self, phi):
        """
        kil line rotate at angle phi
        phi - diferent angle in radians
        """
        p1, p2, cm = self.kil_line_straight()
        a = algebra.calc_transform_matrix_XY(phi)
        shift = [self.xycm[0], self.xycm[1], 0.0]
        pp1 = [p1[0] - shift[0], p1[1] - shift[1], 0.0]
        p1_ = algebra.new_coord(pp1, a) + shift
        pp2 = [p2[0] - shift[0], p2[1] - shift[1], 0.0]
        p2_ = algebra.new_coord(pp2, a) + shift
        dy = 2 * (cm[1] - p1[1])
        p3 = [p1[0], p1[1] + dy]
        p4 = [p2[0], p2[1] + dy]
        pp3 = [p3[0] - shift[0], p3[1] - shift[1], 0.0]
        p3_ = algebra.new_coord(pp3, a) + shift
        pp4 = [p4[0] - shift[0], p4[1] - shift[1], 0.0]
        p4_ = algebra.new_coord(pp4, a) + shift
        #print(f"straight p1={p1}, p2={p2}")
        #print(f"rotate p1={p1_}, p2={p2_}")
        return [p1_, p2_, p4_, p3_, self.xycm]

    def moisten_length(self, phi):
        p1, p2, _ = self.kil_line_rotate(phi)
        vec = algebra.vec(p1, p2)
        if (p2[1] - p1[1]) == 0:
            t = 1.0
        else:
            t = -p1[1] / (p2[1] - p1[1])
        t = min(t, 1.0)
        xy_front = algebra.line_point(p1, vec, t)
        return algebra.distance(p1, xy_front)

    def front_back_h_ksi(self, phi):
        ksi_back = -self.ki[0]
        h_ksi_back = self.h_ksi(ksi_back, phi)
        ksi_front = -self.ki[0] + self.moisten_length(phi)
        h_ksi_front = self.h_ksi(ksi_front, phi) #should be 0
        return [h_ksi_back, h_ksi_front]


    def YM_hydrostatic(self, phi, rho):
        g = 9.81 # m * s^-2
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.moisten_length(phi)
        # print(f"k0={ksi_0}, k1={ksi_1}")
        Y_hs = rho * g / numpy.tan(self.beta)
        #for hydrostatic force
        s1 = (self.ki[1] - self.xycm[1])**2 * (ksi_1 - ksi_0)
        s2 = -(self.ki[1] - self.xycm[1]) * phi * (ksi_1**2 - ksi_0**2)
        s3 = 1.0 / 3 * phi**2 * (ksi_1**3 - ksi_0**3)
        #for hydrostatic moment
        s1m = 0.5 * (self.ki[1] - self.xycm[1])**2 * (ksi_1**2 - ksi_0**2)
        s2m = -2.0 / 3.0 * (self.ki[1] - self.xycm[1]) * phi * (ksi_1**3 - ksi_0**3)
        s3m = 1.0 / 4 * phi**2 * (ksi_1**4 - ksi_0**4)
        return [Y_hs * (s1 + s2 + s3), Y_hs * (s1m + s2m + s3m)]

    def YM_hydrostatic_numerical(self, phi, rho):
        g = 9.81 # kg * m * s^-2
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.moisten_length(phi)
        def f(x):
            return (self.h_ksi(x, phi))**2
        fy = algebra.kvadratura(f, ksi_0, ksi_1) * rho * g / numpy.tan(self.beta)
        def f(x):
            return (self.h_ksi(x, phi))**2 * x
        mz = algebra.kvadratura(f, ksi_0, ksi_1) * rho * g / numpy.tan(self.beta)
        return [fy, mz]

    def YM_hydrostatic_numerical_full_length(self, phi, rho, wave=None, t=0):
        g = 9.81 # kg * m * s^-2
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        def f(x):
            h = self.h_ksi(x, phi, wave, t)
            if h <=0:
                return 0
            return h**2
        fy = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * rho * g / numpy.tan(self.beta)
        def f(x):
            h = self.h_ksi(x, phi, wave, t)
            if h <= 0:
                return 0
            return h**2 * x
        mz = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * rho * g / numpy.tan(self.beta)
        return [fy, mz]

    def YM_hydrostatic_numerical_full_length_move(self, phi, rho, yc=None, wave=None, t=0):
        g = 9.81 # kg * m * s^-2
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            if h <=0:
                return 0
            return h**2
        fy = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * rho * g / numpy.tan(self.beta)
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            if h <= 0:
                return 0
            return h**2 * x
        mz = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * rho * g / numpy.tan(self.beta)
        return [fy, mz]



    def khi_y_khi_m(self, phi):
        khi_y = 1 - numpy.pi * phi / 4 / numpy.tan(self.beta)
        khi_m = khi_y * (1 + 1.0 / (1 + (1 + 2 * numpy.tan(self.beta) / numpy.pi / phi)**0.5))
        # print(f"khi_y={khi_y}, khi_m={khi_m}")
        return [khi_y, khi_m]


    def YM_hydrodynamic(self, phi, rho, V):
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.moisten_length(phi)
        khi_y, khi_m = self.khi_y_khi_m(phi)
        Y_hd = 2 * rho * self.k_beta
        #for hydrodynamic force
        s1 = (self.ki[1] - self.xycm[1]) * V**2 * phi**2 * (ksi_1 - ksi_0)
        s2 = -0.5 * V**2 * phi**3 * (ksi_1**2 - ksi_0**2)
        #for hydrodynamic moment
        s1m = 0.5 * (self.ki[1] - self.xycm[1]) * V**2 * phi**2 * (ksi_1**2 - ksi_0**2)
        s2m = -1.0 / 3 * V**2 * phi**3 * (ksi_1**3 - ksi_0**3)
        #for shear moment
        cf = 0.001 #???
        M_tr = -rho * numpy.pi / 2 * cf  / numpy.sin(self.beta) * V**2
        s1t = (self.ki[1] * (1 - numpy.pi / 4) + self.xycm[1] * numpy.pi / 4) * (self.ki[1] - self.xycm[1])
        s1t *= (ksi_1 - ksi_0)
        s2t = -0.5 * ( (1 - numpy.pi / 2) * self.ki[1] + numpy.pi / 2 * self.xycm[1]) * phi
        s2t *= (ksi_1**2 - ksi_0**2)
        s3t = -numpy.pi / 12 * phi**2 * (ksi_1**3 - ksi_0**3)
        return [Y_hd * (s1 + s2) * khi_y, Y_hd * (s1m + s2m) * khi_m, M_tr * (s1t + s2t + s3t)]

    def YM_hydrodynamic_numerical(self, phi, rho, V):
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.moisten_length(phi)
        khi_y, khi_m = self.khi_y_khi_m(phi)
        Y_hd = rho * self.k_beta * khi_y
        #for hydrodynamic force
        def f(x):
            return 2 * self.h_ksi(x, phi) * (self.V_n(V, x, phi))**2
        fy = algebra.kvadratura(f, ksi_0, ksi_1) * Y_hd
        #for hydrodynamic moment
        M_hd = rho * self.k_beta * khi_m
        def f(x):
            return 2 * self.h_ksi(x, phi) * (self.V_n(V, x, phi))**2 * x
        mz = algebra.kvadratura(f, ksi_0, ksi_1) * M_hd
        #for shear moment
        cf = 0.001 #???
        M_tr = -rho * numpy.pi / 2 * cf  / numpy.sin(self.beta) * V**2
        def f(x):
            return (self.ki[1] - numpy.pi / 4 * self.h_ksi(x, phi)) * self.h_ksi(x, phi)
        mtr = algebra.kvadratura(f, ksi_0, ksi_1) * M_tr
        return [fy, mz, mtr]


    def YM_hydrodynamic_numerical_full_length(self, phi, rho, V, wave=None, t=0):
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        khi_y, khi_m = self.khi_y_khi_m(phi)
        Y_hd = rho * self.k_beta * khi_y
        #for hydrodynamic force
        def f(x):
            h = self.h_ksi(x, phi, wave, t)
            if h <= 0:
                return 0
            return 2 * h * (self.V_n(V, x, phi))**2
        fy = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * Y_hd
        #for hydrodynamic moment
        M_hd = rho * self.k_beta * khi_m
        def f(x):
            h = self.h_ksi(x, phi, wave, t)
            if h <= 0:
                return 0
            return 2 * h * (self.V_n(V, x, phi))**2 * x
        mz = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * M_hd
        #for shear moment
        cf = 0.001 #???
        M_tr = -rho * numpy.pi / 2 * cf  / numpy.sin(self.beta) * V**2
        def f(x):
            h = self.h_ksi(x, phi, wave, t)
            if h <= 0:
                return 0
            return (self.ki[1] - numpy.pi / 4 * h) * h
        mtr = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * M_tr
        return [fy, mz, mtr]


    def YM_hydrodynamic_numerical_full_length_move(self, phi, rho, V, yc = None, wave=None, t=0):
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        khi_y, khi_m = self.khi_y_khi_m(phi)
        Y_hd = rho * self.k_beta * khi_y
        #for hydrodynamic force
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            if h <= 0:
                return 0
            return 2 * h * (self.V_n(V, x, phi))**2
        fy = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * Y_hd
        #for hydrodynamic moment
        M_hd = rho * self.k_beta * khi_m
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            if h <= 0:
                return 0
            return 2 * h * (self.V_n(V, x, phi))**2 * x
        mz = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * M_hd
        #for shear moment
        cf = 0.001 #???
        M_tr = -rho * numpy.pi / 2 * cf  / numpy.sin(self.beta) * V**2
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            if h <= 0:
                return 0
            return (self.ki[1] - numpy.pi / 4 * h) * h
        mtr = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * M_tr
        return [fy, mz, mtr]


    def YM_total(self, phi, rho, V, yc=None, wave=None, t=0):
        Y_hsn, M_hsn = self.YM_hydrostatic_numerical_full_length(phi, rho, wave, t)
        Y_hdn, M_hdn, M_trn = self.YM_hydrodynamic_numerical_full_length(phi, rho, V, yc, wave, t)
        return [Y_hsn + Y_hdn, M_hsn + M_hdn + M_trn]


    def a1(self, phi, phi_t, yc, yc_t, m, I, rho, V, wave, t):
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        khi_y, _ = self.khi_y_khi_m(phi)
        k1 = rho * self.k_beta * khi_y * (2 - numpy.cos(self.beta))
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            Vn = self.V_n(V, x, phi, yc_t, phi_t, wave, t)
            if h <= 0:
                return 0
            if Vn >= 0:
                return 0
            return h**2
        a1 = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * k1 + m
        print(f"a1={a1}")
        return a1

    def a2(self, phi, phi_t, yc, yc_t, m, I, rho, V, wave, t):
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        khi_y, _ = self.khi_y_khi_m(phi)
        k1 = rho * self.k_beta * khi_y * (2 - numpy.cos(self.beta))
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            Vn = self.V_n(V, x, phi, yc_t, phi_t, wave, t)
            if h <= 0:
                return 0
            if Vn >= 0:
                return 0
            return h**2 * x
        a2 = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * k1
        print(f"a2={a2}")
        return a2

    def a3(self, phi, phi_t, yc, yc_t, m, I, rho, V, wave, t):
        g = 9.81 # m/s^2
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        khi_y, _ = self.khi_y_khi_m(phi)
        k1 = rho * self.k_beta * khi_y
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            Vn = self.V_n(V, x, phi, yc_t, phi_t, wave, t)
            if h <= 0:
                return 0
            if Vn >= 0:
                return 0
            s1 = 2 * h * Vn**2
            s2 = h**2 * (2 - numpy.cos(self.beta)) * (2 * V * phi_t + wave.d2_y_wave_dt2(x, t))
            return s1 + s2
        a3 = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * k1 - m * g
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            Vn = self.V_n(V, x, phi, yc_t, phi_t, wave, t)
            if h <= 0:
                return 0
            if Vn >= 0:
                return 0
            return h**2
        a3 += rho * g / numpy.tan(self.beta) * algebra.kvadratura_N(f, ksi_0, ksi_1, 4)
        print(f"a3={a3}")
        return a3

    def b1(self, phi, phi_t, yc, yc_t, m, I, rho, V, wave, t):
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        _, khi_m = self.khi_y_khi_m(phi)
        k1 = rho * self.k_beta * khi_m * (2 - numpy.cos(self.beta))
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            Vn = self.V_n(V, x, phi, yc_t, phi_t, wave, t)
            if h <= 0:
                return 0
            if Vn >= 0:
                return 0
            return h**2 * x
        b1 = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * k1
        print(f"b1={b1}")
        return b1


    def b2(self, phi, phi_t, yc, yc_t, m, I, rho, V, wave, t):
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        _, khi_m = self.khi_y_khi_m(phi)
        k1 = rho * self.k_beta * khi_m * (2 - numpy.cos(self.beta))
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            Vn = self.V_n(V, x, phi, yc_t, phi_t, wave, t)
            if h <= 0:
                return 0
            if Vn >= 0:
                return 0
            return h**2 * x**2
        b2 = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * k1 + I
        print(f"b2={b2}")
        return b2


    def b3(self, phi, phi_t, yc, yc_t, m, I, rho, V, wave, t):
        g = 9.81 # m/s^2
        ksi_0 = -self.ki[0]
        ksi_1 = -self.ki[0] + self.l
        _, khi_m = self.khi_y_khi_m(phi)
        k1 = rho * self.k_beta * khi_m
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            Vn = self.V_n(V, x, phi, yc_t, phi_t, wave, t)
            if h <= 0:
                return 0
            if Vn >= 0:
                return 0
            s1 = 2 * h * Vn**2 * x
            s2 = h**2 * (2 - numpy.cos(self.beta)) * (2 * V * phi_t + wave.d2_y_wave_dt2(x, t)) * x
            return s1 + s2
        b3 = algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * k1
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            Vn = self.V_n(V, x, phi, yc_t, phi_t, wave, t)
            if h <= 0:
                return 0
            if Vn >= 0:
                return 0
            return h**2 * x
        b3 += rho * g / numpy.tan(self.beta) * algebra.kvadratura_N(f, ksi_0, ksi_1, 4)
        cf = 0.001 #???
        M_tr = -rho * numpy.pi / 2 * cf  / numpy.sin(self.beta) * V**2
        def f(x):
            h = self.h_ksi_move(x, phi, yc, wave, t)
            Vn = self.V_n(V, x, phi, yc_t, phi_t, wave, t)
            if h <= 0:
                return 0
            if Vn >= 0:
                return 0
            return (self.ki[1] - numpy.pi / 4 * h) * h
        b3 += algebra.kvadratura_N(f, ksi_0, ksi_1, 4) * M_tr
        print(f"b3={b3}")
        return b3


    def system2(self, y, t, phi, phi_t, m, I, rho, V, wave):
        y1, y2 = y
        print(f"y={y}, phi={phi}, phi_t={phi_t}, m={m}, rho={rho}, V={V}, t={t}")
        dydt = [y2, self.a3(phi, phi_t, y1, y2, m, I, rho, V, wave, t) / self.a1(phi, phi_t, y1, y2, m, I, rho, V, wave, t)]
        return dydt

    def system2f(self, phi, t, yc, yc_t, m, I, rho, V, wave):
        y3, y4 = phi
        print(f"y={yc}, phi={phi}, m={m}, I={I}, rho={rho}, V={V}, t={t}")
        dydt = [y4, self.b3(y3, y4, yc, yc_t, m, I, rho, V, wave, t) / self.b2(y3, y4, yc, yc_t, m, I, rho, V, wave, t)]
        return dydt

    def system4(self, yy, t, m, I, rho, V, wave):
        y1, y2, y3, y4 = yy
        print(f"yy={yy}, m={m}, I={I}, rho={rho}, V={V}, t={t} wave={wave}")
        a1_ = self.a1(y3, y4, y1, y2, m, I, rho, V, wave, t)
        a2_ = self.a2(y3, y4, y1, y2, m, I, rho, V, wave, t)
        a3_ = self.a3(y3, y4, y1, y2, m, I, rho, V, wave, t)
        b1_ = self.b1(y3, y4, y1, y2, m, I, rho, V, wave, t)
        b2_ = self.b2(y3, y4, y1, y2, m, I, rho, V, wave, t)
        b3_ = self.b3(y3, y4, y1, y2, m, I, rho, V, wave, t)
        dydt = [y2, (a3_ * b2_ - a2_ * b3_) / (a1_ * b2_ - a2_ * b1_),
                y4, (a3_ * b1_ - a1_ * b3_) / (a2_ * b1_ - a1_ * b2_)]
        return dydt

    def system4_1(self, yy, t, m, I, rho, V, wave):
        y1, y2, y3, y4 = yy
        print(f"yy={yy}, m={m}, I={I}, rho={rho}, V={V}, t={t} wave={wave}")
        a1_ = self.a1(0, y1, m, rho, wave, t)
        a2_ = self.a2(0, y1, m, rho, wave, t)
        a3_ = self.a3(0, 0, y1, y2, m, rho, V, wave, t)
        b1_ = self.b1(0, y1, m, I, rho, wave, t)
        b2_ = self.b2(0, y1, m, I, rho, wave, t)
        b3_ = self.b3(0, 0, y1, y2, m, I, rho, V, wave, t)
        dydt = [y2, (a3_ * b2_ - a2_ * b3_) / (a1_ * b2_ - a2_ * b1_),
                y4, (a3_ * b1_ - a1_ * b3_) / (a2_ * b1_ - a1_ * b2_)]
        return dydt


    def swimmingy(self, yc, yc_t, phi, phi_t, m, I, rho, V, wave):
        y0 = [yc, yc_t]
        T, _ = wave.get_T_c()

        t = numpy.linspace(0, 5 * T, 501)
        from scipy.integrate import odeint
        sol = odeint(self.system2, y0, t, args=(phi, phi_t, m, I, rho, V, wave))
        print(f"yc={sol[:, 0]}")
        print(f"yc_t={sol[:, 1]}")
        import matplotlib.pyplot as plt
        plt.plot(t, sol[:, 0], 'b', label='y_c')
        plt.plot(t, sol[:, 1], 'g', label='dy_c/dt')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

        return [t, sol]


    def swimmingf(self, yc, yc_t, phi, phi_t, m, I, rho, V, wave):
        ph = [9 * i for i in range(0, 10)]
        phi = [phi_i * numpy.pi / 180.0 for phi_i in ph ]
        ang = phi[1] / numpy.pi * 180.0
        ang_rad = phi[1]

        y0 = [ang_rad, 0.0]
        T, _ = wave.get_T_c()

        t = numpy.linspace(0, 5 * T, 501)
        from scipy.integrate import odeint
        sol = odeint(self.system2f, y0, t, args=(yc, yc_t, m, I, rho, V, wave))
        print(type(sol))
        import matplotlib.pyplot as plt
        plt.plot(t, sol[:, 0], 'b', label='phi')
        plt.plot(t, sol[:, 1], 'g', label='dphi/dt')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

        return [t, sol]


    def swimming(self, yc, yc_t, phi, phi_t, m, I, rho, V, wave):
        y0 = [yc, yc_t, phi, phi_t]
        print(f"y0={y0}")
        T, _ = wave.get_T_c()

        t = numpy.linspace(0, 10 * T, 501)
        from scipy.integrate import odeint
        sol = odeint(self.system4, y0, t, args=(m, I, rho, V, wave))

        import matplotlib.pyplot as plt
        plt.plot(t, sol[:, 0], label='yc')
        plt.plot(t, sol[:, 1], label='dyc/dt')
        plt.plot(t, sol[:, 2], label='phi')
        plt.plot(t, sol[:, 3], label='dphi/dt')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

        return [t, sol]

    def postprocessing(self, t, sol, m, I, rho, V, wave):
        Y = []
        M = []
        for i in range(len(sol)):
            resi = numpy.zeros((1,3))
            Yhs, Mhs = self.YM_hydrostatic_numerical_full_length_move(sol[i][2], rho, sol[i][0], wave, t[i])
            Yhd, Mhd, Mtr = self.YM_hydrodynamic_numerical_full_length_move(sol[i][2], rho, V, sol[i][0], wave, t[i])
            Y.append(Yhs + Yhd)
            M.append(Mhs + Mhd + Mtr)
        import matplotlib.pyplot as plt
        plt.plot(t, Y, label='Y')
        plt.plot(t, M, label='M')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

        return [t, Y, M]


    def postprocessingf(self, t, sol, yc, m, I, rho, V, wave):
        Y = []
        M = []
        for i in range(len(sol)):
            Yhs, Mhs = self.YM_hydrostatic_numerical_full_length_move(sol[i][0], rho, yc, wave, t[i])
            Yhd, Mhd, Mtr = self.YM_hydrodynamic_numerical_full_length_move(sol[i][0], rho, V, yc, wave, t[i])
            Y.append(Yhs + Yhd)
            M.append(Mhs + Mhd + Mtr)
        import matplotlib.pyplot as plt
        plt.plot(t, Y, label='Y')
        plt.plot(t, M, label='M')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

        return [t, Y, M]

    def postprocessingy(self, t, sol, phi, m, I, rho, V, wave):
        Y = []
        M = []
        for i in range(len(sol)):
            Yhs, Mhs = self.YM_hydrostatic_numerical_full_length_move(phi, rho, sol[i][0], wave, t[i])
            Yhd, Mhd, Mtr = self.YM_hydrodynamic_numerical_full_length_move(phi, rho, V, sol[i][0], wave, t[i])
            Y.append(Yhs + Yhd)
            M.append(Mhs + Mhd + Mtr)
        print(f"phi={phi}")
        print(f"Y={Y}")
        print(f"M={M}")

        import matplotlib.pyplot as plt
        plt.plot(t, Y, label='Y')
        plt.plot(t, M, label='M')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.grid()
        plt.show()

        return [t, Y, M]
