
import numpy
import scipy.integrate
import copy

PINOCCHIO_CORE_COSMO_DEFAULT_SEED=0
PINOCCHIO_CORE_COSMO_DEFAULT_NAME='test'
PINOCCHIO_CORE_COSMO_DEFAULT_OM=0.29
PINOCCHIO_CORE_COSMO_DEFAULT_OK=0.0
PINOCCHIO_CORE_COSMO_DEFAULT_OB=0.046
PINOCCHIO_CORE_COSMO_DEFAULT_H=0.69
PINOCCHIO_CORE_COSMO_DEFAULT_S8=0.82
PINOCCHIO_CORE_COSMO_DEFAULT_NS=0.96
PINOCCHIO_CORE_COSMO_DEFAULT_W0=-1.0
PINOCCHIO_CORE_COSMO_DEFAULT_W1=0.0

PINOCCHIO_CORE_COSMO_DEFAULT_C=3.0e5 # [km/s]
PINOCCHIO_CORE_COSMO_DEFAULT_G=4.302e-9 # [Mpc Msolar^-1 (km/s)^2]



class Cosmology(object):

    def __DE_EOS_W0_W1__(self, z):
        return (1.0+z)**(3.0*(1.0+self._w0+self._w1)) * numpy.exp(-3.0*self._w1*z/(1.0+z))

    def __DE_EOS_W0__(self, z):
        return (1.0+z)**(3.0*(1.0+self._w0))

    def __DE_EOS__(self, z):
        return 1.0

    def __SET_DE_EOS__(self):
        if self._w1 != 0.0:
            self.DE_EOS = self.__DE_EOS_W0_W1__
        elif self._w0 != -1.0:
            self.DE_EOS = self.__DE_EOS_W0__
        else:
            self.DE_EOS = self.__DE_EOS__
            
    
    def __init__(self, name=PINOCCHIO_CORE_COSMO_DEFAULT_NAME, seed=PINOCCHIO_CORE_COSMO_DEFAULT_SEED, om=PINOCCHIO_CORE_COSMO_DEFAULT_OM, ok=PINOCCHIO_CORE_COSMO_DEFAULT_OK, ob=PINOCCHIO_CORE_COSMO_DEFAULT_OB, h=PINOCCHIO_CORE_COSMO_DEFAULT_H, s8=PINOCCHIO_CORE_COSMO_DEFAULT_S8, ns=PINOCCHIO_CORE_COSMO_DEFAULT_NS, w0=PINOCCHIO_CORE_COSMO_DEFAULT_W0, w1=PINOCCHIO_CORE_COSMO_DEFAULT_W1, c=PINOCCHIO_CORE_COSMO_DEFAULT_C, G=PINOCCHIO_CORE_COSMO_DEFAULT_G):
        
        self.name = name
        self.seed = seed
        self.om = om
        self.ok = ok
        self.ob = ob
        self.h = h
        self.s8 = s8
        self.ns = ns

        self.c = c
        self.G = G

        self._w0 = w0
        self._w1 = w1

        self.__SET_DE_EOS__()

    @property
    def ol(self):
        return (1.0 - self.om - self.ok)

    @property
    def w0(self):
        return self._w0

    @w0.setter
    def w0(self, value):
        self._w0 = value
        self.__SET_DE_EOS__()

    @property
    def w1(self):
        return self._w1

    @w1.setter
    def w1(self, value):
        self._w1 = value
        self.__SET_DE_EOS__()

    def relative_cosmology(self, z):
        c = Cosmology()

        c.name = 'relative'
        c.seed = self.seed
        c.om = self.om_z(z)
        c.ok = self.ok_z(z)
        c.ob = self.ob_z(z)
        c.h = self.h_z(z)

        c.s8 = 0.0
        c.ns = self.ns

        c._w0 = (self.w0 + self.w1 - self.w1 / (1.0 + z))
        c._w1 = self.w1 / (1.0 + z)

        c.__SET_DE_EOS__()

        return c


    def om_z(self, z):
        
        z1 = (1.0+z)
        z2 = z1 * z1
        z3 = z2 * z1

        omz = self.om * z3
        okz = self.ok * z2
        olz = self.ol * self.DE_EOS(z)

        result = omz / (omz + okz + olz)

        return result

    def ob_z(self, z):

        z1 = (1.0+z)
        z2 = z1 * z1
        z3 = z2 * z1

        omz = self.om * z3
        okz = self.ok * z2
        olz = self.ol * self.DE_EOS(z)

        obz = self.ob * z3

        result = obz / (omz + okz + olz)

        return result

    def ok_z(self, z):
        
        z1 = (1.0+z)
        z2 = z1 * z1
        z3 = z2 * z1

        omz = self.om * z3
        okz = self.ok * z2
        olz = self.ol * self.DE_EOS(z)

        result = okz / (omz + okz + olz)

        return result

    def ol_z(self, z):
        
        z1 = (1.0+z)
        z2 = z1 * z1
        z3 = z2 * z1

        omz = self.om * z3
        okz = self.ok * z2
        olz = self.ol * self.DE_EOS(z)

        result = olz / (omz + okz + olz)

        return result

    def h_z(self, z):
        
        z1 = (1.0+z)
        z2 = z1 * z1
        z3 = z2 * z1

        omz = self.om * z3
        okz = self.ok * z2
        olz = self.ol * self.DE_EOS(z)
        
        result = self.h * numpy.sqrt(omz+okz+olz)

        return result

    def int_z(self, z):

        z1 = (1.0+z)
        z2 = z1 * z1
        z3 = z2 * z1

        omz = self.om * z3
        okz = self.ok * z2
        olz = self.ol * self.DE_EOS(z)
        
        result = 1.0/numpy.sqrt(omz+okz+olz)

        return result

    def Di(self, z):

        f = lambda x,c: c.int_z(x)

        result,error = scipy.integrate.quad(f,0.0,z,args=(self))

        return result

    def Dc(self, z):

        di = self.Di(z)

        result = 3000.0/self.h * di

        return result

    def Dt(self, z):

        di = self.Di(z)

        if self.ok > 0.0:
            sqrok = numpy.sqrt(self.ok)
            return 3000.0/(self.h * sqrok) * numpy.sinh(sqrok * di)
        elif self.ok < 0.0:
            sqrok = numpy.sqrt(-1.0 * self.ok)
            return 3000.0/(self.h * sqrok) * numpy.sin(sqrok * di)
        else:
            return 3000.0/self.h * di

    def Da(self, z):
        
        dt = self.Dt(z)

        result = dt / (1.0 + z)

        return result


    def Da_rel(self, z1, z2):
        
        c = self.relative_cosmology(z1)

        #print c.h, c.om, c.ok, c.ol, c.w0, c.w1

        z21 = (1.0 + z2) / (1.0 + z1) - 1.0

        result = c.Da(z21)

        return result
