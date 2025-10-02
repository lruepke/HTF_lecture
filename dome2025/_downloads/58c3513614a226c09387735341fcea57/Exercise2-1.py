xmin,xmax=-1700,1700
Tmin,Tmax=20,1000
wGauss,x0=400,0
xigma = wGauss/2.355
# see https://en.wikipedia.org/wiki/Full_width_at_half_maximum
fig=plt.figure(figsize=(6,3))
ax=plt.gca()
x=np.linspace(xmin,xmax,500)
y=Tmin+(Tmax-Tmin)*np.exp(-(x-x0)*(x-x0)/(2*xigma*xigma))
ax.plot(x,y)
ax.axvspan(x0-wGauss/2,x0+wGauss/2,ymax=0.5, alpha=0.5, color='lightgreen')
ax.text(0.5,0.25,'FWHM\n%.0f m'%(wGauss), va='center',ha='center',transform=ax.transAxes)
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.set_xlabel('Distance (m)')
ax.set_title('Gaussian shape temperature boundary condition')
ax.set_ylim(y.min(),y.max())
plt.tight_layout()
plt.show()