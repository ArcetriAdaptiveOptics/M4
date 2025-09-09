import numpy as np
def spiralize(ott, npos, step):
    p = spiral_pos(npos, step)
    npos = int(len(p) / 2)
    p0 = ott.parabola.getPosition()
    for i in range(npos):
        #p0 = ott.parabola.getPosition()
        p1 = p0 + np.array([0, 0, 0, p[i, 0], p[i, 1], 0])
        print("New Par command:")
        print(p1)
        ott.parabola.setPosition(p1)


def spiral_pos(nstep, step):
    p = np.array([0,0])
    pp = []
    for i in range(nstep):
            direct = (-1)**i*step
            mov = i+1
            for j in range(mov):
                    p = p+np.array([direct,0])
                    pp.append(p)
            for j in range(mov):
                    p = p+np.array([0,direct])
                    pp.append(p)
    pp = np.array(pp)
    plt.plot(pp[:,0],pp[:,1],'-x')
    return pp

