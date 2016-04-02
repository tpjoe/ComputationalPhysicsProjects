import numpy as np
from numpy import linalg as LA
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager


# Set the font dictionaries (for plot title and axis titles)
title_font = {'fontname':'Arial', 'size':'16', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Arial', 'size':'16'}

# Set the font properties (for use in legend)
font_path = 'C:\Windows\Fonts\Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=14)

class solarsystem:

    count = 0

    def __init__(self):
        self.planets = []

    def add(self, a_planet):
        self.planets.append(a_planet)
        solarsystem.count += 1

    def showPlanet(self):
        for a_planet in self.planets:
            print a_planet.name

    def compilePlanet(self):
        y_i = np.zeros([solarsystem.count,7])
        ii = 0
        for a_planet in self.planets:
            y_i[ii,6] = a_planet.mass
            y_i[ii,5] = a_planet.vz
            y_i[ii,4] = a_planet.vy
            y_i[ii,3] = a_planet.vx
            y_i[ii,2] = a_planet.z
            y_i[ii,1] = a_planet.y
            y_i[ii,0] = a_planet.x
            ii += 1
        return y_i

    def force(self, x, y, z, Mothers):
        G = 4*np.pi**2
        distance = x**2 + y**2 + z**2
        F = G*Mothers/distance**1.5
        return F

    def kinetic(self, vec, i):
        KE = 0.5*vec[i,6]*(vec[i,3]**2+vec[i,4]**2+vec[i,5]**2)
        return KE

    def potential(self, vec, i):
        G = 4*np.pi**2
        PE = -2*G*vec[i,6]/(vec[i,0]**2+vec[i,1]**2+vec[i,2]**2)**0.5
        return PE

    def angMomentum(self, vec, i):
        L = 2*np.pi*vec[i,6]*(vec[i,0]**2+vec[i,1]**2+vec[i,2]**2)/(vec[i,3]**2+vec[i,4]**2+vec[i,2]**5)**0.5
        return L


    def positionPlot(self, position):
        plt.figure()
        for i in range(solarsystem.count):
            plt.plot(position[i,0], position[i,1])
        plt.xlabel("X", **axis_font)
        plt.ylabel("Y", **axis_font)
        plt.title("Trajectory", **title_font)
        plt.tick_params(axis = 'x', labelsize = 15)
        plt.tick_params(axis = 'y', labelsize = 15)
        plt.show()
        plt.show()

    def distance(self, vec, i):
        d = (vec[i,0]**2 + vec[i,1]**2 + vec[i,2]**2)**0.5
        return d


    def KEPlot(selfself, KE, h, tmax):
        plt.figure()
        m = range(int(tmax/h)+1)
        for i in range(solarsystem.count):
            if i != 0:
                plt.plot(m, KE[i]+KE[0])
        plt.xlabel("Time steps", **axis_font)
        plt.ylabel("Kinetic Energy", **axis_font)
        plt.title("Kinetic Energy", **title_font)
        plt.tick_params(axis = 'x', labelsize = 15)
        plt.tick_params(axis = 'y', labelsize = 15)
        plt.show()

    def PEPlot(self, PE, h, tmax):
        plt.figure()
        m = range(int(tmax/h)+1)
        for i in range(solarsystem.count):
            if i != 0:
                plt.plot(m, PE[i]+PE[0])
        plt.xlabel("Time steps", **axis_font)
        plt.ylabel("Potential Energy", **axis_font)
        plt.title("Potential Energy", **title_font)
        plt.tick_params(axis = 'x', labelsize = 15)
        plt.tick_params(axis = 'y', labelsize = 15)
        plt.show()


    def distancePlot(self, distance, h, tmax):
        plt.figure()
        m = range(int(tmax/h)+1)
        for i in range(solarsystem.count):
            if i != 0:
                plt.plot(m, distance[i])
        plt.xlabel("Time steps", **axis_font)
        plt.ylabel("Distance from Sun (A.U.)", **axis_font)
        plt.title("Distance from Sun", **title_font)
        plt.tick_params(axis = 'x', labelsize = 15)
        plt.tick_params(axis = 'y', labelsize = 15)
        plt.show()

    def LPlot(self, L, h, tmax):
        plt.figure()
        m = range(int(tmax/h)+1)
        for i in range(solarsystem.count):
            if i != 0:
                plt.plot(m, L[i])
        plt.xlabel("Time steps", **axis_font)
        plt.ylabel("Angular Momentum", **axis_font)
        plt.title("Angular Momentum", **title_font)
        plt.tick_params(axis = 'x', labelsize = 15)
        plt.tick_params(axis = 'y', labelsize = 15)
        plt.show()

def derivative(dat, de, n):
    for i in range(n):
        acc_x  = 0.0
        acc_y = 0.0
        acc_z = 0.0
        for j in range(n):
            if i != j:
                mod_force = SS.force(dat[j,0]-dat[i,0], dat[j,1]-dat[i,1], dat[j,2]-dat[i,2], dat[j,6])
                acc_x += mod_force*(dat[j,0]-dat[i,0])
                acc_y += mod_force*(dat[j,1]-dat[i,1])
                acc_z += mod_force*(dat[j,2]-dat[i,2])
        de[i,3] = acc_x
        de[i,4] = acc_y
        de[i,5] = acc_z
    for i in range(n):
        de[i,0] = dat[i,3]
        de[i,1] = dat[i,4]
        de[i,2] = dat[i,5]


def sum_matrix(results, coeff_one, first, coeff_two, second, n):
    for i in range(n):
        for j in range(6):
            results[i,j] = coeff_one*first[i,j] + coeff_two*second[i,j]
        results[i,6] = first[i,6]


def solverRK4(h, tmax):
    y_temp = np.zeros([solarsystem.count, 7])
    position = np.zeros((solarsystem.count,7,int(tmax/h)+1))
    KE = np.zeros((solarsystem.count, int(tmax/h)+1))
    PE = np.zeros((solarsystem.count, int(tmax/h)+1))
    L = np.zeros((solarsystem.count, int(tmax/h)+1))
    distance = np.zeros((solarsystem.count, int(tmax/h)+1))
    k1 = np.zeros([solarsystem.count, 7])
    k2 = np.zeros([solarsystem.count, 7])
    k3 = np.zeros([solarsystem.count, 7])
    k4 = np.zeros([solarsystem.count, 7])
    y_i = SS.compilePlanet()
    t = 0.0
    n=0

    while t < tmax:
        derivative(y_i, k1, solarsystem.count)
        sum_matrix(y_temp, 1, y_i, 0.5*h, k1, solarsystem.count)
        derivative(y_temp, k2, solarsystem.count)
        sum_matrix(y_temp, 1, y_i, 0.5*h, k2, solarsystem.count)
        derivative(y_temp, k3, solarsystem.count)
        sum_matrix(y_temp, 1, y_i, h, k3, solarsystem.count)
        derivative(y_temp, k4, solarsystem.count)

        for i in range(solarsystem.count):
            if i!=0:
                for j in range(6):
                    y_i[i,j] = y_i[i,j] + h*(k1[i,j] + 2*k2[i,j] + 2*k3[i,j] + k4[i,j])/6
                    position[i,j,n] = y_i[i,j]
                    distance[i,n] = SS.distance(y_i,i)
                    KE[i,n] = SS.kinetic(y_i,i)
                    PE[i,n] = SS.potential(y_i,i)
                    PE[0,n] = SS.potential(y_i,i)/y_i[i,6]*y_i[0,6]
                    if i != 0:
                        L[i,n] = SS.angMomentum(y_i,i)
        t += h
        n += 1

# fixing minor bugs
    position[:,:,int(tmax/h)] = position[:,:,int(tmax/h)-1]
    KE[:,int(tmax/h)] = KE[:,int(tmax/h)-1]
    PE[:,int(tmax/h)] = PE[:,int(tmax/h)-1]
    distance[:,int(tmax/h)] = distance[:,int(tmax/h)-1]

    return [position, KE, PE, distance, L]






class planet:

    def __init__(self, name, mass, x, y, z, vx, vy, vz):
        self.name = name
        self.mass = mass
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

Sun = planet('Sun', 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
Mercury = planet('Mercury', 1.2e-7, 0.39, 0.0, 0.0, 0.0, 9.96, 0.0)
Venus = planet('Venus', 2.4e-6, 0.72, 0.0, 0.0, 0.0, 7.36, 0.0)
Earth = planet('Earth', 1.5e-6, 1.0, 0.0, 0.0, 0.0, 6.28, 0.0)
Mars = planet('Mars', 3.3e-7, 1.52, 0.0, 0.0, 0.0, 5.06, 0.0)
Jupiter = planet('Jupiter', 9.5e-4, 5.2, 0.0, 0.0, 0.0, 2.75, 0.0)
Saturn = planet('Saturn', 2.75e-4, 9.54, 0.0, 0.0, 0.0, 2.04, 0.0)
Uranus = planet('Uranus', 4.4e-5, 19.19, 0.0, 0.0, 0.0, 1.43, 0.0)
Neptune = planet('Neptune',5.1e-5, 30.06, 0.0, 0.0, 0.0, 1.14, 0.0)
Pluto = planet('Pluto', 5.6e-9, 39.53, 0.0, 0.0, 0.0, 0.99, 0.0)

SS = solarsystem()

SS.add(Sun)
SS.add(Earth)
#SS.add(Mercury)
#SS.add(Venus)
#SS.add(Mars)
#SS.add(Jupiter)
#SS.add(Saturn)
#SS.add(Uranus)
#SS.add(Neptune)
#SS.add(Pluto)


[h,tmax] = [0.01,10]

[position, KE, PE, distance, L] = solverRK4(h, tmax)
SS.positionPlot(position)
#SS.KEPlot(KE, h, tmax)
#SS.PEPlot(PE, h, tmax)
SS.distancePlot(distance, h, tmax)
#SS.LPlot(L,h, tmax)
