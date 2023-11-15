xmin,xmax=0,3000
width_half,q_total = 500, 86000
c=width_half/np.sqrt(2*np.log(2))
qmax,qmin,x0=q_total/(c*np.sqrt(2*np.pi)),0,1500
fig=plt.figure(figsize=(6,3))
ax=plt.gca()
x=np.linspace(xmin,xmax,1000)
q=qmin+(qmax-qmin)*np.exp(-(x-x0)*(x-x0)/(2*c*c))
ax.plot(x,q)
ax.fill_between(x,q,0,color='lightgreen',alpha=0.5)
ax.hlines(0.5, 0,0.5,transform=ax.transAxes, ls='dashed',lw=0.5, color='gray')
ax.vlines(0.5, 0.5,1,transform=ax.transAxes, ls='dashed',lw=0.5, color='gray')
ax.plot([x0-width_half,x0],[qmax/2, qmax/2],color='r',lw=2)
ax.text(x0-width_half/2, qmax/2, 'half width\n%.0f'%(width_half),va='top',ha='center')
q_total_cal=0
dx=(x.max()-x.min())/(len(x)-1)
# integrate heat flux along profile
for i in range(0,len(x)):
    q_total_cal = q_total_cal + dx*q[i]
ax.text(x0,qmin*1.2,'Integrated total heat flux\n%.1f kw'%(q_total_cal/1000),va='bottom',ha='center',fontweight='bold')
text_parms='$x_{0}$ = %.0f\n$q_{min}$ = %.0f\nc = %.1f m\n$q_{max}$ = %.1f w/m$^2$'%(x0,qmin,c,qmax)
ax.text(0.98,0.98,text_parms,va='top',ha='right',transform=ax.transAxes)
ax.set_ylabel('Heat flux (W/m$^{\mathregular{2}}$)')
ax.set_xlabel('Distance (m)')
ax.set_title('Gaussian shape heat flux profile')
ax.set_ylim(q.min(),q.max())
plt.tight_layout()
plt.show()