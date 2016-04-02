import numpy as np


class solarsystem:
    def __init__(self):
        self.planets = []

    def add(self, a_planet):
        self.planets.append(a_planet)

    def showPlanet(self):
        for a_planet in self.planets:
            print a_planet.name

    def compilePlanet(self):
        y_i = np.zeros([planet.count,7])
        i = 0
        for a_planet in self.planets:
            y_i[i,6] = a_planet.mass
            y_i[i,5] = a_planet.vz
            y_i[i,4] = a_planet.vy
            y_i[i,3] = a_planet.vx
            y_i[i,2] = a_planet.z
            y_i[i,1] = a_planet.y
            y_i[i,0] = a_planet.x
            i += 1
        return y_i

    def force(self, x, y, z, Mothers):
        G = 4*np.pi**2
        distance = x**2 + y**2 + z**2
        F = G*Mothers/distance**1.5
        return F

def derivative(dat, de, n):
    for i in range(planet.count):
        acc_x = acc_y = acc_z = 0.0
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
        results[i,6] = first [i,6]

def solverRK4(h, tmax):
    y_i_temp = k1 = k2 = k3 = k4 = np.zeros([planet.count, 7])
    y_i = SS.compilePlanet()
    t = 0.0

    while t<tmax:
        derivative(y_i, k1, planet.count)
        sum_matrix(y_i_temp, 1, y_i, 0.5*h, k1, planet.count)
        derivative(y_i, k2, planet.count)
        sum_matrix(y_i_temp, 1, y_i, 0.5*h, k2, planet.count)
        derivative(y_i, k3, planet.count)
        sum_matrix(y_i_temp, 1, y_i, 0.5*h, k3, planet.count)
        derivative(y_i, k4, planet.count)

        for i in range(planet.count):
            for j in range(6):
                y_i[i,j] = y_i[i,j] + h*(k1[i,j] + 2*k2[i,j] + 2*k3[i,j] + k4[i,j])/6
        t += h


class planet:

    count = 0

    def __init__(self, name, mass, x, y, z, vx, vy, vz):
        self.name = name
        self.mass = mass
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

        planet.count += 1


Sun = planet('Sun', 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
Mercury = planet('Mercury', 1.2e-7, 0.39, 0.0, 0.0, 0.0, 9.96, 0.0)
Venus = planet('Venus', 2.4e-6, 0.72, 0.0, 0.0, 0.0, 7.36, 0.0)
Earth = planet('Earth', 1.5e-6, 1.0, 0.0, 0.0, 0.0, 6.26, 0.0)
Mars = planet('Mars', 3.3e-7, 1.52, 0.0, 0.0, 0.0, 5.06, 0.0)
Jupiter = planet('Jupiter', 9.5e-4, 5.2, 0.0, 0.0, 0.0, 2.75, 0.0)
Saturn = planet('Saturn', 2.75e-4, 9.54, 0.0, 0.0, 0.0, 2.04, 0.0)
Uranus = planet('Uranus', 4.4e-5, 19.19, 0.0, 0.0, 0.0, 1.43, 0.0)
Neptune = planet('Neptune',5.1e-5, 30.06, 0.0, 0.0, 0.0, 1.14, 0.0)
Pluto = planet('Pluto', 5.6e-9, 39.53, 0.0, 0.0, 0.0, 0.99, 0.0)

SS.add(Sun)
SS.add(Earth)
SS.add(Mercury)
SS.add(Venus)
SS.add(Mars)
SS.add(Jupiter)
SS.add(Saturn)
SS.add(Uranus)
SS.add(Neptune)
SS.add(Pluto)
