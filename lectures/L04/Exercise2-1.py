xmin,xmax=-1700,1700
Tmin,Tmax=20,1000
wGauss,x0=10,0
fig=plt.figure(figsize=(6,3))
ax=plt.gca()
x=np.linspace(xmin,xmax,500)
y=Tmin+(Tmax-Tmin)*np.exp(-(x-x0)*(x-x0)/(2*wGauss*wGauss))
ax.plot(x,y)
ax.set_ylabel('Temperature ($^{\circ}$C)')
ax.set_xlabel('Distance (m)')
ax.set_title('Gaussian shape temperature boundary condition')
ax.set_ylim(y.min(),y.max())
plt.tight_layout()
plt.show()