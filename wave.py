import numpy

class wave():
    def set_data(self, hwave, lwave, L_=0):
        self.h = hwave
        self.l = lwave
        self.L = L_
        _, c = self.get_T_c()
        self.shift_t = -self.L / c

    def y_wave(self, x, t):
        g = 9.81
        k = 2 * numpy.pi / self.l
        c = (g * self.l / 2 / numpy.pi)**0.5
        return self.h / 2 * numpy.cos(k * (x + c * (t + self.shift_t)))

    def d_y_wave_dt(self, x, t):
        g = 9.81
        k = 2 * numpy.pi / self.l
        c = (g * self.l / 2 / numpy.pi)**0.5
        return -k * c * self.h / 2 * numpy.sin(k * (x + c * (t + self.shift_t)))

    def d2_y_wave_dt2(self, x, t):
        g = 9.81
        k = 2 * numpy.pi / self.l
        c = (g * self.l / 2 / numpy.pi)**0.5
        return (k * c)**2 * self.h / 2 * numpy.cos(k * (x + c * (t + self.shift_t)))

    def get_T_c(self):
        g = 9.81
        c = (g * self.l / 2 / numpy.pi)**0.5
        T = (2 * numpy.pi * self.l / g )**0.5
        return [T, c]
