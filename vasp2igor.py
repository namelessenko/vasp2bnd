from curses.ascii import isalpha
from vasprun import Vasprun
import numpy as np
import math
import func

HARTREE = 0.036749308136649
BOHR = 0.5291

class Vasp2Igor(Vasprun):
    """
    Conevert Vasp out to Igor input
    """
    def __init__(self, vaspfile):
        super().__init__(vaspfile)
        super().get_structure()
        super().get_efermi()
        super().get_eigenvalues()
        super().get_kpoints()
        super().get_electronicparam()
        super().get_subdivision()
        self.reducebasis()
        self.kpath()
        
    
    def reducebasis(self):
        self._reduced_a = self._a
        self._reduced_b = self._b / self._a
        self._reduced_c = self._c / self._a

        self.__reduced_lattiece_constant = self._reduced_a/BOHR, self._reduced_b, self._reduced_c #reduced constant *BOHR
        self.__reduced_basis = np.divide(self._basis, (self._reduced_a))
        self.__reduced_rec_basis = np.multiply(self._recbasis, (2 * np.pi * self._reduced_a))

    
    def bandstransform(self):
        """Reading and replace eigenvalue to eig*=HARTEE
        """
        transformedbands = []

        for i in super().eigenvalues:
            transformedbands.append(i*HARTREE)
        transformedbands = list(func.list_split(transformedbands, int(self.get_electronicparam()["NBANDS"])))
        
        return transformedbands

    
    def transformkpoints(self):
        """Reading and replace kpoint to kpoint*=2/1.73250 (squere area of Hexagonal BZ is 3*((3)^1/2)/2) it works only for HEXAGONAL BZ
        """
        transformedkpoints = []
        for kpoints in super().kpoints:
            for kpoint in kpoints:
                transformedkpoints.append(kpoint*2/math.sqrt(3))
        transformedkpoints = list(func.list_split(transformedkpoints, 3))

        return transformedkpoints
        
    def kpath(self):
        with open('ignore/Cu_KPOINTS') as f:
            self._kpath = f.readlines()[4:]
            for i in range(len(self._kpath)):
                self._kpath[i] = self._kpath[i].split()
            self._kpath = list(filter(None, self._kpath))
            for i in self._kpath:
                for j in range(len(i)):
                    if j == 3:
                        i[j] = str(i[j])
                    else:
                        i[j] = float(i[j])
            print(self._kpath)

                
                        
    def writefile(self, name):
        with open(name, 'w') as f:
            #Firsr line lat. const
            for i in self.__reduced_lattiece_constant:
                f.write(str(i)+"\t")
            f.write("\n\n")

            #Second line efermi, kpoints info, num of band 
            f.write("8  ")
            f.write(str(self.efermi*HARTREE)+"  ")
            f.write(str(np.prod(self.subdivision))+"    ")
            for i in self.subdivision:
                f.write(str(i))
            f.write(str(self.get_electronicparam()["NBANDS"]))
            f.write("   1	0	0	0	0")
            f.write("\n\n")

            #K-Path from KPOINTS file
            for i in self._kpath[::2]:
                f.write(str(i)[1:-1]+"\n")
            f.write("\n\n")

            #Something fatuous numbers
            f.write("1 /	nopused \n\
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 \n\
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 \n\
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 \n\
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 \n")

            #Reduced basis
            for i in self.__reduced_basis:
                for j in i:
                    f.write(str(j)+"\t")
                f.write("\n")
            f.write("\n")

            #Reduced rec_basis
            for i in self.__reduced_rec_basis:
                for j in i:
                    f.write(str(j)+"\t")
                f.write("\n")
            f.write("\n")

            #Kpoints and bands
            for i, (kpoint, band) in enumerate(zip(self.transformkpoints(), self.bandstransform()),start=1):
                f.write('  '.join(map(str, kpoint))+ '        /' + str(i) + '\n')
                f.write('  '.join(map(str, band))+'\n')