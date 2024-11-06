import numpy
import astropy.units as u
import astropy.cosmology
import matplotlib.pyplot as plt
import matplotlib
import scipy.integrate as integrate
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

M_snIa = (-19.46 + 5 * numpy.log10(70/60))
M_cc = (-18 + 5 * numpy.log10(70/60)) 

v_snIa = 1e4 / 3.09e19 # Mpc/s
v_cc = 1e4 / 3.09e19 # Mpc/s

t_snIa = 18*3600*24 #[s]
t_cc  = 25*3600*24 #[s]

r=[r_snIa,r_cc]
M=[M_snIa, M_cc]
v=numpy.array([v_snIa, v_cc])
t=numpy.array([t_snIa, t_cc])

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
			_zs.append(astropy.cosmology.z_at_value(cosmo.distmod, (l-_M)*u.mag,0.000001,0.5).value)
		rates.append(_rates)
		zs.append(_zs)

	fig, ax1 = plt.subplots(constrained_layout = True)
	color = 'tab:red'
	ax1.set_xlabel(r'$m_\text{lim}$')
	ax1.set_ylabel(r'$z_\text{max}$', color=color)
	ax1.tick_params(axis='y', labelcolor=color)
	ax1.plot(limmag, zs[0], ls='-.', label=r'SN Ia [$z_\text{max}$]')
	ax1.plot(limmag, zs[1], ls='-.', label=r'CCSN [$z_\text{max}$]')

	ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

	color = 'tab:blue'
	ax2.set_ylabel(r'$N_\text{cum}$ [$\text{yr}^{-1}$]', color=color)  # we already handled the x-label with ax1
	ax2.plot(limmag, rates[0], label=r'SN Ia [$N_\text{cum}$]')
	ax2.plot(limmag, rates[1], label=r'CCSN [$N_\text{cum}$]')
	ax2.tick_params(axis='y', labelcolor=color)

	# fig.tight_layouts()  # otherwise the right y-label is slightly clipped
	fig.legend(loc=2)
	plt.savefig('rates.pdf')
	# plt.show()

snRate()

def angularSize():

	zs = numpy.linspace(0.0001,0.004,40)
	dA = cosmo.angular_diameter_distance(zs)


	# theta2=(t*v)[:, None]/dA[None,:] * 206265 * 1e6
	# print(176/theta2)
	theta=2*(t*v)[:, None]/dA[None,:]

	fig, ax1 = plt.subplots(constrained_layout = True)
	color = 'tab:red'
	ax1.set_xlabel(r'$z$')
	ax1.set_ylabel(r'$\theta$ [nrad]', color=color)
	ax1.tick_params(axis='y', labelcolor=color)
	ax1.plot(zs, theta[0]*1e9, ls='-.', label=r'SN Ia [$\theta$]')
	ax1.plot(zs, theta[1]*1e9, ls='-.', label=r'CC [$\theta$]')
	# plt.plot(limmag, zs[1], ls='-.', label=r'CCSN [$z_\text{max}$]')
	ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

	color = 'tab:blue'
	ax2.set_ylabel(r'$d$ [km]', color=color)  # we already handled the x-label with ax1
	# ax2.plot(zs, 176*(550/700)/(theta[0]*1e6), label=r'SN Ia [$r$]')
	# ax2.plot(zs, 176*(550/700)/(theta[1]*1e6), label=r'CCSN [$r$]')
	ax2.plot(zs, 1.22*440e-9/(theta[0])*1e-3, label=r'SN Ia [$d$]')
	ax2.plot(zs, 1.22*440e-9/(theta[1])*1e-3, label=r'CCSN [$d$]')
	ax2.tick_params(axis='y', labelcolor=color)

	# fig.tight_layouts()  # otherwise the right y-label is slightly clipped
	fig.legend(loc=2)
	plt.savefig('angle.pdf')
	# plt.show()

angularSize()