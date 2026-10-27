import astropy.constants as c
import astropy.units as u
import numpy as np
from scipy.integrate import quad as num_int

def kappa_nu (nu, beta=2): #eq 10
	nu_0 = 1e12*u.Hz
	if nu <= nu_0: 
		rescaled = (nu/nu_0)**beta
		k = 0.1 * rescaled * u.cm**2/u.g
		return k



def M_dust (M_gas, DGR): #eq 11
	return (M_gas * DGR).cgs


def M_gas (M_str, z): #eq 12
	coeff = 3.87 * 1e9
	term1 = (1+z)**(1.35)
	rescaled = 1e-10*(M_str/c.M_sun)
	term2 = (rescaled)**(0.49)

	return c.M_sun*(coeff*term1*term2).cgs


def Z_logOH (SFR, M_str): #eq 13
	term1 = 0.14*np.log10(SFR/(c.M_sun/u.yr))
	term2 = 0.37*np.log10(M_str/(c.M_sun))

	return term2 - term1 + 4.82


def DGR (Z, Z_sun = 8.69, logDGRsun = -2.21): #eq 14, Z_sun from https://www.aanda.org/articles/aa/pdf/2006/43/aa5216-06.pdf
	if Z > 0.26*Z_sun:
		logDGR = logDGRsun + np.log10(Z/Z_sun)
	else:
		logDGR = logDGRsun + 3.15*np.log10(Z/Z_sun) + 1.25

	return 10**(logDGR)


def Planck (T, nu):
	term1 = (2*c.h*nu**3)/(c.c**2) 
	exp = (c.h*nu)/(c.k_B*T)
	term2 = (-1+np.e**exp)
	return (term1/term2).cgs



def I_z (z, M=0.27, L=0.73): #from 29.168 in Carroll & Ostlie
	term1 = M*(1+z)**3
	term2 = L
	final = (term1 + term2)**(-0.5)
	return final



def d_Lum (z, Ho=70*u.km/(u.s*1e6*u.pc)): #from 29.184 in Carroll & Ostlie
	coeff = (c.c/Ho)*(1+z)
	S_z = num_int(I_z, 0, z)[0]
	return (coeff*S_z).cgs



def S_nu (nu, T_d, M_str, z, SFR): #eq 17
	kappa = kappa_nu(nu)
	B_nu = Planck(T_d, nu)
	d_L = d_Lum(z)
	Z_oh = Z_logOH(SFR, M_str)
	dtog = DGR(Z_oh)
	M_g = M_gas(M_str, z)
	M_d = M_dust(M_g, dtog)

	num = kappa * B_nu * M_d * (1+z)
	den = d_L**2

	return (num/den).cgs





def T_from_Snu (nu,Snu,M_str,z,SFR):#rearranging eq17
	kappa = kappa_nu(nu)
	d_L = d_Lum(z)
	Z_oh = Z_logOH(SFR, M_str)
	dtog = DGR(Z_oh)
	M_g = M_gas(M_str, z)
	M_d = M_dust(M_g, dtog)

	num = (1+z)*kappa*2*c.h*M_d*nu**3
	den = Snu*(c.c*d_L)**2

	print(num, den)

	ln = np.log(num/den + 1)
	coeff = c.k_B/(c.h*nu)

	Tinv = coeff*ln

	return 1/Tinv



