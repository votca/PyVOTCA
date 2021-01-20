import numpy as np

class Molecule:
    def __init__(self):
        self.__elements = []
        self.__coordinates = []

    def add_atom(self,element : str , x : float, y : float, z : float ):
        self.__elements.append(element)
        self.__coordinates.append(np.array([x, y, z]))

    def print(self):
        for (element, coordinates) in zip(self.__elements, self.__coordinates):
            print(element, coordinates)

    def readXYZfile(self, filename):
        xyzfile=open(filename)
        natoms = int(xyzfile.readline())
        title = xyzfile.readline()[:-1]
        for i in range(natoms):
            line = xyzfile.readline().split()
            self.__elements.append(line[0])
            self.__coordinates.append(np.array([float(line[1]), float(line[2]), float(line[3])]))
        xyzfile.close()
    
    def writeXYZfile(self, filename):
        xyzfile=open(filename,"w")
        natoms=len(self.__elements)
        xyzfile.write("%d\n%s\n" % (natoms, "Created by PyVOTCA"))
        for i in range(natoms):
            xyzfile.write("%s  %f  %f  %f\n" % (self.__elements[i], float(self.__coordinates[i][0]), float(self.__coordinates[i][1]),float(self.__coordinates[i][2]) ))
        xyzfile.close()



