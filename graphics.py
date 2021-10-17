import matplotlib.pyplot as plt


def plot_data(traces, xname, lbl, code=0):
    for i in range(len(traces)):
        plt.plot(traces[i][0], traces[i][1], marker=".", markersize=3,
                label=lbl[i])
    plt.xlabel(xname)
    if code==0:
        plt.axis('square')
        plt.axis('equal')
    plt.grid()
    #plt.ylim((0.25, 0.32))
    plt.show()
