import xml.etree.ElementTree as ET
import numpy as np


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
        # self.subdivision = []

        for tag in self.root_node.findall('kpoints/generation/v'):
            key=tag.get('name')
            if key == 'divisions':
                self.subdivision = (list(map(int, tag.text.split())))
        # for tag in self.root_node.findall('kpoints/generation/v'):
        #     key = tag.get('name')
        #     if key == 'genvec1':
        #         self.subdivision.append(list(map(float, tag.text.split())))
        #     elif key == 'genvec2':
        #         self.subdivision.append(list(map(float, tag.text.split())))
        #     elif key == 'genvec3':
        #         self.subdivision.append(list(map(float, tag.text.split())))


        # ____CALC SUBDIVISION____
        # self.subdivision = np.asanyarray(self.subdivision)
        # self.h1=np.sqrt(np.sum(self.subdivision[0,:]**2))
        # self.h2=np.sqrt(np.sum(self.subdivision[1,:]**2))
        # self.h3=np.sqrt(np.sum(self.subdivision[2,:]**2))

        # self.N1 = int(ceil((self.b1*3)/self.h1))
        # self.N2 = int(ceil((self.b2*3*(self.b1/self.b2))/self.h2))
        # self.N3 = int(ceil((self.b3*3*(self.b1/self.b3))/self.h3))
       
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