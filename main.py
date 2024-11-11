import numpy
import astropy.units as u
import astropy.cosmology
import matplotlib.pyplot as plt
import matplotlib
import scipy.integrate as integrate
from scipy.fft import fft2, rfft2, fftshift
from scipy.special import jv
matplotlib.use("TkAgg")

cosmo = astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

tmax=18 # days
vmax=1e4 # km/s
r = vmax * tmax * 3600 * 24

z=1200 # km/s
d_z = z/ 70 # Mpc
d_z = d_z * 3.09e13  # Mpc to km

theta = r/d_z

d = 176 / (theta * 1e6)

r_snIa = 2.43e-5 / u.yr / (u.Mpc)**3 # SNe/yr/Mpc^3/h^3_70
r_cc = 9.1e-5 / u.yr / (u.Mpc)**3 # SNe/yr/Mpc^3/h^3_70
r_II = r_cc* 0.649
r_IIL = r_cc* 0.079
r_IIP = r_cc* 0.57
r_IIb = r_cc * 0.109
r_IIn = r_cc * 0.047
r_Ic = r_cc * 0.075
r_Ib = r_cc * 0.108

M_snIa = (-19.46 + 5 * numpy.log10(70/60))
# M_II= (-18 + 5 * numpy.log10(70/60)) 
M_II =  -15.97
M_IIL =  -18.28
M_IIP =  -16.67
M_IIb =  -16.69
M_IIn =  -17.66
M_Ic =  -17.44
M_Ib =  -18.26

v_snIa = 1e4 / 3.09e19 # Mpc/s
v_cc = 1e4 / 3.09e19 # Mpc/s

t_snIa = 18*3600*24 #[s]
t_cc  = 25*3600*24 #[s]

# r=[r_snIa,r_cc]
# M=[M_snIa, M_cc]
# v=numpy.array([v_snIa, v_cc])
# t=numpy.array([t_snIa, t_cc])

r=[r_snIa,r_II, r_IIL, r_IIP, r_IIb, r_IIn, r_Ic, r_Ib]
M=[M_snIa,M_II, M_IIL, M_IIP, M_IIb, M_IIn, M_Ic, M_Ib]
v=numpy.array([v_snIa, v_cc, v_cc, v_cc, v_cc, v_cc, v_cc, v_cc])
t=numpy.array([t_snIa, t_cc, t_cc, t_cc, t_cc, t_cc, t_cc, t_cc])
type_str = ["Ia","II", "IIL", "IIP", "IIb", "IIn", "Ic", "Ib"]

def dNdm(m, M, r):
	z = astropy.cosmology.z_at_value(cosmo.distmod, (m-M)*u.mag,0.000001,0.5)
	dVdz = cosmo.differential_comoving_volume(z)
	dmdz = 5/numpy.log(10)*(astropy.constants.c.to('km/s')*(1+z)/cosmo.H(z)/cosmo.luminosity_distance(z)+ 1/(1+z))
	return (r*dVdz/dmdz*4*numpy.pi*u.sr/(1+z)).value

def snRate():
	limmag = numpy.linspace(8,12,40)
	rates = []
	zs = []
	for _r, _M in zip(r, M):
		_rates=[]
		_zs=[]
		for l in limmag:
			_rates.append( integrate.quad(dNdm, 0, l, args=(_M,_r))[0])
			_zs.append(astropy.cosmology.z_at_value(cosmo.distmod, (l-_M)*u.mag,0.0000001,0.5).value)
		rates.append(_rates)
		zs.append(_zs)

	fig, ax1 = plt.subplots(constrained_layout = True)
	color = 'tab:red'
	ax1.set_xlabel(r'$m_\text{lim}$')
	ax1.set_ylabel(r'$z_\text{max}$ [dotted]')# , color=color)
	ax1.tick_params(axis='y') #, labelcolor=color)
	for z,ts in zip(zs,type_str):
		ax1.plot(limmag, z, ls='dotted')#, label=r"${} [z_\text{{max}}]$".format(ts))
		# a1$.lot(limmag, zs[1], ls='-.', label=r'CCSN [$z_\text{max}$]')

	ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

	color = 'tab:blue'
	ax2.set_ylabel(r'$N_\text{cum}$ [$\text{yr}^{-1}$] [solid]') #, color=color)  # we already handled the x-label with ax1
	for rate, ts in zip(rates,type_str):
		ax2.plot(limmag, rate, label=r"{0}".format(ts))
	# ax2.plot(limmag, rates[1], label=r'CCSN [$N_\text{cum}$]')
	ax2.tick_params(axis='y') #, labelcolor=color)

	# fig.tight_layouts()  # otherwise the right y-label is slightly clipped
	ax2.legend(loc=2)
	fig.tight_layout() 
	plt.savefig('rates.pdf')
	# plt.show()

snRate()
wfe
def angularSize():

	zs = numpy.linspace(0.0001,0.004,40)
	dA = cosmo.angular_diameter_distance(zs)


	# theta2=(t*v)[:, None]/dA[None,:] * 206265 * 1e6
	# print(176/theta2)
	theta=2*(t*v)[:, None]/dA[None,:]

	fig, ax1 = plt.subplots(constrained_layout = True)
	color = 'tab:red'
	ax1.set_xlabel(r'$z$')
	ax1.set_ylabel(r'$\theta$ [nrad] [dashed]') #, color=color)
	ax1.tick_params(axis='y') #, labelcolor=color)
	ax1.plot(zs, theta[0]*1e9) #, ls='-.', label=r'SN Ia [$\theta$]')
	ax1.plot(zs, theta[1]*1e9) #, ls='-.', label=r'CC [$\theta$]')
	# plt.plot(limmag, zs[1], ls='-.', label=r'CCSN [$z_\text{max}$]')
	ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

	color = 'tab:blue'
	ax2.set_ylabel(r'$d$ [km] [solid]') #, color=color)  # we already handled the x-label with ax1
	# ax2.plot(zs, 176*(550/700)/(theta[0]*1e6), label=r'SN Ia [$r$]')
	# ax2.plot(zs, 176*(550/700)/(theta[1]*1e6), label=r'CCSN [$r$]')
	ax2.plot(zs, 1.22*440e-9/(theta[0])*1e-3, label=r'SN Ia')
	ax2.plot(zs, 1.22*440e-9/(theta[1])*1e-3, label=r'CCSN')
	ax2.tick_params(axis='y') # , labelcolor=color)

	# fig.tight_layouts()  # otherwise the right y-label is slightly clipped
	ax2.legend(loc=2)
	plt.savefig('angle.pdf')
	# plt.show()

angularSize()
wef

def gamma():
	def Pz(p):
		rmax = 2.25
		y2 = rmax**2 - p**2
		cos2  = y2/rmax**2
		return (1-cos2)/(1+cos2)


	u = numpy.linspace(-1.5,1.5,101)
	# plt.plot(u,Pz(u))
	# plt.xlabel('p')
	# plt.ylabel('Pz')
	# plt.savefig('Pz.pdf')
	# plt.clf()
	# wfe
	def intensity(t1, t2, disk=False):
		rho = numpy.sqrt(t1**2+t2**2)
		if (rho>1):
			return 0
		else:
			if disk:
				return 1
			theta = numpy.arctan2(t1,t2)
			return 0.5*(1- Pz(rho)) + Pz(rho) * numpy.cos(theta)**2

	nbin=1001
	u = numpy.linspace(-10.,10.,nbin)
	v = u

	I = numpy.zeros((nbin,nbin))
	for i,_u in enumerate(u):
		for j,_v in enumerate(v):
			I[i,j] = intensity(_u,_v)

	totalI = I.sum()

	# plt.plot(I[nbin//2,450:550])
	# plt.plot(I[450:550,nbin//2])
	# plt.show()
	# wef

	nrange = 16
	plt.imshow(I[(nrange//2-1)*nbin//nrange:(nrange//2+1)*nbin//nrange,(nrange//2-1)*nbin//nrange:(nrange//2+1)*nbin//nrange])
	plt.savefig('intensity.pdf')
	plt.clf()

	gamma = fft2(I)
	gamma2 = numpy.abs(gamma)**2
	# print(gamma2)

	I = numpy.zeros((nbin,nbin))
	for i,_u in enumerate(u):
		for j,_v in enumerate(v):
			I[i,j] = intensity(_u,_v,disk=True)

	I = I/I.sum()*totalI

	dum = fft2(I)
	dum2 = numpy.abs(dum)**2

	nrange = 40
	plt.plot(fftshift(gamma2)[(nrange//2-1)*nbin//nrange:(nrange//2+1)*nbin//nrange,nbin//2],label='u',color='blue'); 
	plt.plot(fftshift(gamma2)[nbin//2,(nrange//2-1)*nbin//nrange:(nrange//2+1)*nbin//nrange],label='y',color='brown'); 
	plt.plot(fftshift(dum2)[nbin//2,(nrange//2-1)*nbin//nrange:(nrange//2+1)*nbin//nrange],label='Airy',color='red'); 
	plt.legend()
	plt.savefig('gamma.pdf')
	plt.clf()

	# plt.imshow(fftshift(gamma2)[(nrange//2-1)*nbin//nrange:(nrange//2+1)*nbin//nrange,(nrange//2-1)*nbin//nrange:(nrange//2+1)*nbin//nrange])
	# plt.savefig('gamma_im.pdf')
	# plt.clf()


gamma()