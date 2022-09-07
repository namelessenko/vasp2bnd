from array import array
from numba import njit
from numba import jit


from itertools import islice
import xml.etree.ElementTree as ET
import func
import time
import numpy as np
from math import ceil


HARTREE = 0.036749308136649
BOHR = 0.5291

class Vasprun:
    """parser for vasprun.xml"""

    def __init__(self, vaspfile):
        self.root_node = ET.parse(vaspfile).getroot()
        self._basis = []
        self._recbasis = []
        self.__eigenvalues = []
        self.__kpoints = []

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def b1(self):
        return self._b1

    @property
    def b2(self):
        return self._b2

    @property
    def b3(self):
        return self._b3

    @property
    def alpha(self):
        return self._alpha

    @property
    def beta(self):
        return self._beta

    @property
    def gamma(self):
        return self._gamma


    def get_structure(self):
        """
        Parsing of structural constants 
        """
        
        for tag in self.root_node.findall('structure'):
            if tag.get('name') == 'finalpos':
                value = tag.findall('crystal/varray')
                for i in value:
                    if i.get('name') == 'basis':
                        for j in i:
                            self._basis.append(list(map(float ,j.text.split())))
                    elif i.get('name') == 'rec_basis':
                        for j in i:
                            self._recbasis.append(list(map(float ,j.text.split())))  

        # BASIS, LATTIECE CONSTANT 
        self._basis = np.asanyarray(self._basis)
        self._recbasis = np.asanyarray(self._recbasis)    
        #LAT CONST
        self._a=np.sqrt(np.sum(self._basis[0,:]**2))
        self._b=np.sqrt(np.sum(self._basis[1,:]**2))
        self._c=np.sqrt(np.sum(self._basis[2,:]**2))
        #REC LAT CONST
        self._b1 = np.sqrt(np.sum(self.recbasis[0,:]**2))
        self._b2 = np.sqrt(np.sum(self.recbasis[1,:]**2))
        self._b3 = np.sqrt(np.sum(self.recbasis[2,:]**2))
        #Angl of lattiece
        self._alpha = np.arccos(np.dot(self._basis[1,:],self._basis[2,:])/self._b/self._c)*360/2/np.pi
        self._beta  = np.arccos(np.dot(self._basis[0,:],self._basis[2,:])/self._a/self._c)*360/2/np.pi
        self._gamma = np.arccos(np.dot(self._basis[0,:],self._basis[1,:])/self._a/self._b)*360/2/np.pi


    @property
    def eigenvalues(self):
        return self.__eigenvalues
    
    @property
    def basis(self):
        return self._basis

    @property
    def recbasis(self):
        return self._recbasis

    def get_eigenvalues(self):
        """
        Parsing eigenvalues for band structure in list format
        """

        for num, tag in enumerate(self.root_node.findall('calculation/eigenvalues/array/set/set/set'), start=1):
            value = tag.get('comment')
            if value == 'kpoint {}'.format(num):
                for i in tag.findall('r'):
                    self.__eigenvalues.append(float(i.text[0:-11]))

    @property
    def kpoints(self):
        return self.__kpoints

    def get_kpoints(self):
        """
        Parsing kpoints
        """

        self.__kpoints = []

        for tag in self.root_node.findall('kpoints/varray'):
            value = tag.get('name')
            if value == 'kpointlist':
                for item in tag.findall('v'):
                    self.__kpoints.append(item.text.split())
                    for i in self.__kpoints:
                        for j,item in enumerate(i):
                            i[j] = float(item)


    def get_electronicparam(self):
        """_parsing of electronic parameters for bandtransform function_

        Returns:
            _dict_: _string_
        """        

        electronicparam = {}

        for tag in self.root_node.findall('parameters/separator'):
            key = tag.get('name')
            if key == 'electronic':
                for value in tag.findall('i'):
                    key = value.get('name')
                    electronicparam[key] = value.text

        return electronicparam

    @property
    def efermi(self):
        return self.__efermi

    @efermi.setter
    def efermi(self, value):
        print("Be careful!\nYou're changing the efermi\n")
        self.__efermi = value

    def get_efermi(self):
        self.__efermi = self.root_node.findall('calculation/dos/i')
        self.__efermi = [float(i.text) for i in self.__efermi]
        self.__efermi = self.__efermi[0]
        
    def get_subdivision(self):
        self.subdivision = []
        for tag in self.root_node.findall('kpoints/generation/v'):
            key = tag.get('name')
            if key == 'genvec1':
                self.subdivision.append(list(map(float, tag.text.split())))
            elif key == 'genvec2':
                self.subdivision.append(list(map(float, tag.text.split())))
            elif key == 'genvec3':
                self.subdivision.append(list(map(float, tag.text.split())))

        self.subdivision = np.asanyarray(self.subdivision)
        self.h1=np.sqrt(np.sum(self.subdivision[0,:]**2))
        self.h2=np.sqrt(np.sum(self.subdivision[1,:]**2))
        self.h3=np.sqrt(np.sum(self.subdivision[2,:]**2))

        self.N1 = int(ceil((self.b1*3)/self.h1))
        self.N2 = int(ceil((self.b2*3*(self.b1/self.b2))/self.h2))
        self.N3 = int(ceil((self.b3*3*(self.b1/self.b3))/self.h3))
       
       
        # tmp1 = 0
        # tmp2 = 0
        # tmp3 = 0
        # for tag in self.root_node.findall('kpoints/varray'):
        #     value = tag.get('name')
        #     if value == 'kpointlist':
        #         for item in tag.findall('v'):
        #             self.__kpoints.append(item.text.split())

        # for i in range(len(self.__kpoints)):
        #     if self.__kpoints[i][0] > 0:
        #         tmp1 = self.__kpoints[i][0]
        #         break
        # for i in range(len(self.__kpoints)):
        #     if self.__kpoints[i][1] > 0:
        #         tmp2 = self.__kpoints[i][1]
        #         break
        # for i in range(len(self.__kpoints)):
        #     if self.__kpoints[i][2] > 0:
        #         tmp3 = self.__kpoints[i][2]
        #         break
                
        # print(tmp1)
        # print(tmp2)
        # print(tmp3)
                    
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
        self.hkpoints()
        
    
    
    def reducebasis(self):
        self._reduced_a = self._a
        self._reduced_b = self._b / self._a
        self._reduced_c = self._c / self._a


        self.__reduced_lattiece_constant = self._reduced_a/BOHR, self._reduced_b, self._reduced_c #reduced constant *BOHR
        self.__reduced_basis = np.divide(self._basis, (self._reduced_a))
        self.__reduced_rec_basis = np.multiply(self._recbasis, (2 * np.pi * self._reduced_a))

    
    def bandstransform(self):
        """Reading and replace eigevalues to eig*=HARTEE"""

        transformedbands = []

        for i in super().eigenvalues:
            transformedbands.append(i*HARTREE)#and round just for macros

        transformedbands = list(func.list_split(transformedbands, int(self.get_electronicparam()["NBANDS"])))

        return transformedbands

    
    def transformkpoints(self):

        transformedkpoints = []

        for kpoints in super().kpoints:
            for kpoint in kpoints:
                transformedkpoints.append(kpoint*(2/1.73250))#and round just for macros
        transformedkpoints = list(func.list_split(transformedkpoints, 3))
        return transformedkpoints
        
    def hkpoints(self):
        with open('ignore/KPOINTS') as f:
            self._txt = f.readlines()[4:]
            

    def writefile(self, name):
        with open(name, 'w') as f:
            for i in self.__reduced_lattiece_constant:
                f.write(str(i)+"\t")
            f.write("\n\n")
            
            for i in self._txt:
                f.write(i)
            f.write("\n\n")

            f.write(str(self.efermi*HARTREE))
            f.write("\n\n")

            for i in self.__reduced_basis:
                for j in i:
                    f.write(str(j)+"\t")
                f.write("\n")
            f.write("\n")

            for i in self.__reduced_rec_basis:
                for j in i:
                    f.write(str(j)+"\t")
                f.write("\n")
            f.write("\n")


            for i, (kpoint, band) in enumerate(zip(self.transformkpoints(), self.bandstransform()),start=1):
                # f.writelines(str(kpoint)[1:-1]+ '        /' + str(i) + '\n')
                # f.writelines(str(band)[1:-1]+'\n')
                f.write('  '.join(map(str, kpoint))+ '        /' + str(i) + '\n')
                f.write('  '.join(map(str, band))+'\n')
           


if __name__ == "__main__":
    start_time = time.time()
    v = Vasp2Igor('ignore/vasprunISYM.xml') 
    print(v.N1)
    print(v.N2)
    print(v.N3) 
    v.writefile("test.BND")  
    print("--- %s seconds ---" % (time.time() - start_time))
 