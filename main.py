'''
Glissirovanie ploskokilevatoj plastini po volne bez smachivanija skul
'''
import kinematic
import kinematic_wet
import kinematic
import graphics
import wave
import numpy

def test_poplavok():
    p = kinematic_wet.plate_wet()
    depth = [(0.48, 0.0), (4.09, 0.0), (4.8, 0.2737), (5.66, 0.290),
              (5.67, 0.0), (9.04, 0.0), (9.05, 0.028), (9.55, 0.048),
              (9.91, 0.009), (9.92, 0.0), (10.7, 0.0)]
    deadrise = [(0.48, 35 * numpy.pi / 180), (5.66, 35 * numpy.pi / 180),
                (9.04, 26 * numpy.pi / 180), (10.7, 26 * numpy.pi / 180)]
    cmx = 5.747
    depth1 = []
    deadrise1 = []
    for i in range(len(depth)):
        depth1.append((depth[i][0] - cmx, depth[i][1]))
    for i in range(len(deadrise)):
        deadrise1.append((deadrise[i][0] - cmx, deadrise[i][1]))
    p.set_data(depth1, deadrise1)
    ang = 9.37
    ang_rad = ang * numpy.pi / 180.0
    Y_hsn, M_hsn = p.YM_hydrostatic_numerical(ang_rad, 998.2)
    print(f"at {ang} degrees Y_hs = {Y_hsn}, M_hs={M_hsn}")
    Y_hdn, M_hdn, Mtr = p.YM_hydrodynamic_numerical(ang_rad, 998.2, 21.4)
    print(f"at {ang} degrees Y_hd = {Y_hdn}, M_hd={M_hdn}, M_tr={Mtr}")
    print(f"at {ang} degrees Y = {Y_hsn + Y_hdn}, M={M_hsn + M_hdn + Mtr}")

def test_plate_wet():
    p = kinematic_wet.plate_wet()
    depth = [(-0.45, 0.13669), (0.4097, 0.00164)]
    deadrise = [(-0.45, 30 * numpy.pi / 180), (0.4097, 30 * numpy.pi / 180)]
    p.set_data(depth, deadrise)
    ang = 9.0
    ang_rad = ang * numpy.pi / 180.0
    Y_hsn, M_hsn = p.YM_hydrostatic_numerical(ang_rad, 998.2)
    print(f"at {ang} degrees Y_hs = {Y_hsn}, M_hs={M_hsn}")
    Y_hdn, M_hdn, Mtr = p.YM_hydrodynamic_numerical(ang_rad, 998.2, 8.0)
    print(f"at {ang} degrees Y_hd = {Y_hdn}, M_hd={M_hdn}")


def test_plate():
    p = kinematic.plate()
    p.set_data([0.0, 0.15], [0.745, 0.4], 5.0)
    p.set_deadrise_angle(30.0 * numpy.pi / 180)
    p1, p2, cm = p.kil_line_straight()
    dy = 2 * (cm[1] - p1[1])
    '''
    traces = [[[p1[0], p2[0], p2[0], p1[0], p1[0]],
               [p1[1], p2[1], p2[1] + dy, p1[1] + dy, p1[1]]]]
    lbl = ['straight']
    '''
    traces = []
    lbl = []
    ph = [9 * i for i in range(0, 10)]
    phi = [phi_i * numpy.pi / 180.0 for phi_i in ph ]
    '''
    for i in range(len(phi)):
        p1, p2, cm = p.kil_line_rotate(phi[i])
        rotate = [[p1[0], p2[0], cm[0], p1[0]], [p1[1], p2[1], cm[1], p1[1]]]
        traces.append(rotate)
        lbl.append('rotate' + str(ph[i]))
        #graphics.plot_data(traces, 'x', lbl)
    '''
    ang = phi[1] / numpy.pi * 180.0
    ang_rad = phi[1]
    l = p.moisten_length(ang_rad)
    print(f"l at {ang} degrees = {l}")
    h1, h2 = p.front_back_h_ksi(ang_rad)
    print(f"at {ang} degrees h1 = {h1}, h2={h2}")
    Y_hs, M_hs = p.YM_hydrostatic(ang_rad, 998.2)
    print(f"at {ang} degrees Y_hs = {Y_hs}, M_hs={M_hs}")
    Y_hsn, M_hsn = p.YM_hydrostatic_numerical_full_length(ang_rad, 998.2)
    print(f"at {ang} degrees Y_hsn = {Y_hsn}, M_hsn={M_hsn}")
    Y_hsn, M_hsn = p.YM_hydrostatic_numerical_full_length_move(ang_rad, 998.2)
    print(f"at {ang} degrees Y_hsn_move = {Y_hsn}, M_hsn_move={M_hsn}")


    Y_hd, M_hd, M_tr = p.YM_hydrodynamic(ang_rad, 998.2, 8.0)
    print(f"at {ang} degrees Y_hd = {Y_hd}, M_hd = {M_hd}, M_tr = {M_tr}")
    Y_hdn, M_hdn, M_trn = p.YM_hydrodynamic_numerical_full_length(ang_rad, 998.2, 8.0)
    print(f"at {ang} degrees Y_hdn = {Y_hdn}, M_hdn = {M_hdn}, M_trn = {M_trn}")
    print ("Summary:")
    print(f"at {ang} degrees Y = {Y_hsn + Y_hdn}, M={M_hsn + M_hdn + M_trn}")

    print("Wave")
    w = wave.wave()
    w.set_data(0.0, 10.0, L_=5.0-0.745)
    time = 0.0
    Y_hsn, M_hsn = p.YM_hydrostatic_numerical_full_length(ang_rad, 998.2, w, time)
    Y_hdn, M_hdn, M_trn = p.YM_hydrodynamic_numerical_full_length(ang_rad, 998.2, 8.0, w, time)
    print ("Summary wave:")
    print(f"at {ang} degrees Y = {Y_hsn + Y_hdn}, M={M_hsn + M_hdn + M_trn}")

    p1, p2, p3, p4, _ = p.profil_rotate(ang_rad)
    rotate = [[p1[0], p2[0], p3[0], p4[0], p1[0]],
              [p1[1], p2[1], p3[1], p4[1], p1[1]]]
    traces.append(rotate)
    lbl.append('rotate' + str(ang))
    test_wave(w, traces, lbl, -0.745, 5, time)
    graphics.plot_data(traces, 'x', lbl)


def test_run_wave():
    p = kinematic.plate()
    p.set_data([0.0, 0.089], [5.45, 0.15], 10.0)
    p.set_deadrise_angle(30.0 * numpy.pi / 180)
    ph = [9 * i for i in range(0, 10)]
    phi = [phi_i * numpy.pi / 180.0 for phi_i in ph ]
    ang = phi[1] / numpy.pi * 180.0
    ang_rad = phi[1]
    w = wave.wave()
    w.set_data(0.5, 10.0, L_=0)
    T, _ = w.get_T_c()
    t = numpy.linspace(0, T, 20)
    Y = []
    M = []
    for time in t:
        yt, mt = p.YM_total(ang_rad, 998.2, 8.0, w, time)
        Y.append(yt)
        M.append(mt)
    traces = [[t, Y]]
    graphics.plot_data(traces, 't', [['fy']], 1)
    traces = [[t, M]]
    graphics.plot_data(traces, 't', [['mz']], 1)


def solve():
    p = kinematic.plate()
    p.set_data([0.0, 0.15], [0.745, 0.4], 5.0)
    p.set_deadrise_angle(30.0 * numpy.pi / 180)
    p1, p2, cm = p.kil_line_straight()
    dy = 2 * (cm[1] - p1[1])
    traces = []
    lbl = []
    #traces = [[[p1[0], p2[0], p2[0], p1[0], p1[0]],
    #           [p1[1], p2[1], p2[1] + dy, p1[1] + dy, p1[1]]]]
    #lbl = ['straight']
    ph = [9 * i for i in range(0, 10)]
    phi = [phi_i * numpy.pi / 180.0 for phi_i in ph ]
    ang = phi[1] / numpy.pi * 180.0
    ang_rad = phi[1]
    Y_hs, M_hs = p.YM_hydrostatic(ang_rad, 998.2)
    print(f"at {ang} degrees Y_hs = {Y_hs}, M_hs={M_hs}")
    Y_hsn, M_hsn = p.YM_hydrostatic_numerical_full_length(ang_rad, 998.2)
    print(f"at {ang} degrees Y_hsn = {Y_hsn}, M_hsn={M_hsn}")
    Y_hsn, M_hsn = p.YM_hydrostatic_numerical_full_length_move(ang_rad, 998.2)
    print(f"at {ang} degrees Y_hsn_move = {Y_hsn}, M_hsn_move={M_hsn}")


    Y_hd, M_hd, M_tr = p.YM_hydrodynamic(ang_rad, 998.2, 8.0)
    print(f"at {ang} degrees Y_hd = {Y_hd}, M_hd = {M_hd}, M_tr = {M_tr}")
    Y_hdn, M_hdn, M_trn = p.YM_hydrodynamic_numerical_full_length(ang_rad, 998.2, 8.0)
    print(f"at {ang} degrees Y_hdn = {Y_hdn}, M_hdn = {M_hdn}, M_trn = {M_trn}")
    Y_hdn, M_hdn, M_trn = p.YM_hydrodynamic_numerical_full_length_move(ang_rad, 998.2, 8.0)
    print(f"at {ang} degrees Y_hdn_move = {Y_hdn}, M_hdn_move = {M_hdn}, M_trn_move = {M_trn}")
    print ("Summary:")
    print(f"at {ang} degrees Y = {Y_hsn + Y_hdn}, M={M_hsn + M_hdn + M_trn}")
    m = (Y_hsn + Y_hdn) / 9.81
    p1, p2, p3, p4, _ = p.profil_rotate(ang_rad)
    rotate = [[p1[0], p2[0], p3[0], p4[0], p1[0]],
              [p1[1], p2[1], p3[1], p4[1], p1[1]]]
    traces.append(rotate)
    lbl.append('rotate' + str(ang))

    w = wave.wave()
    xl = 5-0.745
    w.set_data(0.0, 3.0, L_=xl)
    time = 0.0
    test_wave(w, traces, lbl, -0.745, 5, time)

    print("Wave")
    Y_hsn, M_hsn = p.YM_hydrostatic_numerical_full_length(ang_rad, 998.2, w, time)
    Y_hdn, M_hdn, M_trn = p.YM_hydrodynamic_numerical_full_length(ang_rad, 998.2, 8.0, w, time)
    print ("Summary wave:")
    print(f"at {ang} degrees Y = {Y_hsn + Y_hdn}, M={M_hsn + M_hdn + M_trn}")

    Y_hsn, M_hsn = p.YM_hydrostatic_numerical_full_length_move(ang_rad, 998.2, 0.15, w, time)
    Y_hdn, M_hdn, M_trn = p.YM_hydrodynamic_numerical_full_length_move(ang_rad, 998.2, 8.0, 0.15, w, time)
    print ("Summary wave move:")
    print(f"at {ang} degrees Y = {Y_hsn + Y_hdn}, M={M_hsn + M_hdn + M_trn}")


    graphics.plot_data(traces, 'x', lbl)
    #t, sol = p.swimmingy(0.15, 0, ang_rad, 0, m,1000.0, 998.2, 8.0, w)
    #p.postprocessingy(t, sol, ang_rad, m, 1000.0, 998.2, 8.0, w)
    #t, sol = p.swimmingf(0.15, 0, 0, 0, m, 1000.0, 998.2, 8.0, w)
    #p.postprocessingf(t, sol, 0.15, m, 1000.0, 998.2, 8.0, w)
    t, sol = p.swimming(0.17, 0.0, 0.9*ang_rad, 0.0, m, 10000.0, 998.2, 8.0, w)
    p.postprocessing(t, sol, m, 10000.0, 998.2, 8.0, w)

def test_wave(w, traces, lbl, ksi0, L, time):
    g = 9.81
    T, c = w.get_T_c()
    #L-ksi0 - front side of the ship
    #yv(x, t) = h/2 * cos(2pi / lamda*(x + c * t))
    #yv(L + ksi0, t0) = h/2:
    #cos(2pi / lamda* (L + ksi0 + c * t0)) = 1
    #2pi / lamda * (L + ksi0 + c * t0) = 2pi * k
    #L + ksi0 + c * t0 = lamda * k
    #t0 = -(L + ksi0) / c
    t = numpy.linspace(0, T, 10)
    print(f"T={T}")
    x = numpy.linspace(ksi0, L + ksi0, 201)
    hh = w.y_wave(x, time)
    vv = w.d_y_wave_dt(x, time)
    aa = w.d2_y_wave_dt2(x, time)
    traces.append([x, hh])
    #traces.append([x, vv])
    #traces.append([x, aa])
    lbl.append(f"wave_{time}")
    #lbl.append(f"vv_{time}")
    #lbl.append(f"aa_{time}")


def calc_wave():
    w = wave.wave()
    w.set_data(0.5, 10.0, x_0=0, L_=0)
    time = 1.0
    xl = 5
    traces = []
    lbl = []
    test_wave(w, traces, lbl, -0.745, xl, time)
    graphics.plot_data(traces, 'x', lbl)


#test_poplavok()
#test_plate()
#test_plate_wet()
#test_wave()
#test_run_wave()
solve()
#calc_wave()
