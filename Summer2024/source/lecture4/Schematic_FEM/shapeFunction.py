#!/Users/zguo/.pyenv/shims/python
# -*-coding:utf-8-*-
# Plot schematic diagram of shape functions of FEM
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2021/05/25, GEOMAR
# ===============================================================

import sys
import argparse
import sys
import os
from colored import fg, bg, attr
C_GREEN = fg('green')
C_RED = fg('red')
C_BLUE = fg('blue')
C_DEFAULT = attr('reset')
#===============================================================
import linecache
import numpy as np
import meshio
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import LightSource
from nice import niceAxis,text3d
from matplotlib.collections import PolyCollection

mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
figpath='.'
fmt_figs=['pdf','svg']

def usage(argv):
    basename = argv[0].split('/')
    basename = basename[len(basename)-1]
    description='Plot schematic diagram of shape functions of FEM'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2021/05/25, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' example usage'+C_DEFAULT)
    print('='*len(head))
def createMesh_1D(xmin=0,xmax=10,dx=2):
    x=np.arange(xmin,xmax,dx)
    el2nod=[]
    for i in range(0,len(x)-1):
        el2nod.append([i, i+1])
    return {'GCOORD':{'x':x}, 'EL2NOD':el2nod}
def createMesh_2D(xmin=0,dx=1,nx=5,ymin=0,dy=1,ny=4):
    x=np.arange(xmin,xmin+(nx)*dx,dx)
    y=np.arange(ymin,ymin+(ny)*dy,dy)
    xx,yy=np.meshgrid(x,y)
    xx,yy=xx.reshape(-1,1),yy.reshape(-1,1)
    GCOORD=np.column_stack((xx,yy))
    # el2nod
    EL2NOD=[]
    for j in range(0,ny-1):
        for i in range(0,nx-1):
            LL=i+(nx)*j
            EL2NOD.append(np.array([LL, LL+1, LL+nx+1, LL+nx],dtype=int))
    EL2NOD=np.array(EL2NOD,dtype=int)
    # print(GCOORD.shape,EL2NOD)
    return {'GCOORD':GCOORD, 'EL2NOD':EL2NOD}
def Linear1D():
    colors=mpl.cm.get_cmap('tab10').colors
    dx=2
    mesh=createMesh_1D(dx=2)
    x=mesh['GCOORD']['x']
    el2nod=mesh['EL2NOD']
    # figure
    fig=plt.figure(figsize=(8,2.5))
    ax=plt.gca()
    # shape function definition
    Ni_x  = lambda xi, dx, x : 1 - (x-xi)/dx
    Ni1_x = lambda xi, dx, x : (x-xi)/dx
    # plot mesh and shape function
    for i, el in enumerate(el2nod):
        colors_node=[colors[el[0]%len(colors)], colors[el[1]%len(colors)]]
        # element
        ax.plot(x[el],x[el]*0-0.1,lw=4,clip_on=False,color='k')
        ax.text(x[el].mean(),-0.15,'Element %d'%(i),ha='center',va='top')
        # node
        ax.scatter(x[el],x[el]*0-0.1,marker='o',fc=colors_node,ec='w',s=100,clip_on=False,zorder=10)
        ax.text(x[el][0],-0.1, '%d'%(el[0]),ha='center',va='center',zorder=11,fontsize=9,color='w')
        if(i==(len(el2nod)-1)): 
            ax.text(x[el][1],-0.1, '%d'%(el[1]),ha='center',va='center',zorder=11,fontsize=9,color='w')
        # calculate shape function value of each node
        tmp_x=np.linspace(x[el][0],x[el][1],2)
        Ni = Ni_x(x[el][0], dx, tmp_x)
        Ni1 = Ni1_x(x[el][0], dx, tmp_x)
        ax.plot(tmp_x, Ni,   color=colors_node[0])
        ax.plot(tmp_x, Ni1,  color=colors_node[1])
    # axis settings
    ax.set_ylim(0,1)
    ax.set_xlim(x.min(),x.max())
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.xaxis.set_ticks([])
    ax.spines['left'].set_position(("axes", -0.02))
    ax.set_ylabel('Shape function value')
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    for fmt in fmt_figs:
        figurename=str('%s/shapeFunction_Linear_1D.%s'%(figpath,fmt))
        plt.savefig(figurename, bbox_inches='tight')

def loadMesh2D(gmshfile):
    """Node connection of triangle of gmsh is clockwise, we need to make it counterclockwise

    Args:
        gmshfile ([type]): gmsh(.msh) file name

    Returns:
        [type]: {'GCOORD':GCOORD, 'EL2NOD':EL2NOD}
    """
    mesh=meshio.read(gmshfile)
    GCOORD=mesh.points[:,0:2]
    EL2NOD=[]
    for cells in mesh.cells:
        if((cells.type=='quad') | (cells.type=='triangle')):
            EL2NOD=cells.data
    EL2NOD_old=EL2NOD.copy()
    for i in range(1,EL2NOD.shape[1]):
        EL2NOD[:,i]=EL2NOD_old[:,-i]
    return {'GCOORD':GCOORD, 'EL2NOD':EL2NOD}
def plotMesh_2D(mesh,colorsName='tab20',figname='mesh2D_structured',extend_x=0):
    colors=mpl.cm.get_cmap(colorsName).colors
    x,y=mesh['GCOORD'][:,0],mesh['GCOORD'][:,1]
    el2nod=mesh['EL2NOD']
    len_x,len_y=x.max()-x.min(),y.max()-y.min()
    figwidth=12
    figheight=len_y/len_x*figwidth
    fig=plt.figure(figsize=(figwidth,figheight))
    ax=plt.gca()
    # plot element
    for i,element in enumerate(el2nod):
        ax.fill(x[element],y[element],color=colors[i%len(colors)],closed=True,ec='k',clip_on=False)
        ax.text(x[element].mean(),y[element].mean(),'Element %d'%(i),ha='center',va='center',bbox={'fc':'lightgray','ec':'None','boxstyle':'round'})
    # plot node
    nodes=[]
    for el in el2nod:
        for node in el:
            nodes.append(node)
    nodes=np.array(nodes,dtype=int)
    for nod in np.unique(nodes):
        ax.plot(x[nod],y[nod],marker='o',ms=15,mfc='w',mec='k',clip_on=False)
        ax.text(x[nod],y[nod],'%d'%(nod),ha='center',va='center')
    ax.axis('off')
    ax.axis('scaled')
    ax.set_xlim(x.min()-len_x*extend_x,x.max()+len_x*extend_x)
    ax.set_ylim(y.min(),y.max())
    for fmt in fmt_figs:
        figurename=str('%s/%s.%s'%(figpath,figname,fmt))
        plt.savefig(figurename, bbox_inches='tight')
def plotConnectivityMatrix_2D(mesh,colorsName='tab20',figname='Matrix2D_structured'):
    colors=mpl.cm.get_cmap(colorsName).colors
    x,y=mesh['GCOORD'][:,0],mesh['GCOORD'][:,1]
    el2nod=mesh['EL2NOD']
    len_x,len_y=x.max()-x.min(),y.max()-y.min()
    nnel,nnod=len(el2nod),len(x)
    nnel1=int(np.sqrt(nnel))
    nnel2=int(nnel/nnel1)
    rows=np.min([nnel1,nnel2])
    cols=int(nnel/rows)
    width_ratios=[1]*(cols*2+1)
    width_ratios[cols]=0.1
    figheight=12
    figwidth=figheight*2.1
    fig,axes=plt.subplots(rows,cols*2+1, sharex=True, sharey=True,figsize=(figwidth,figheight),
    gridspec_kw={"hspace":0.15,'wspace':0.15,"width_ratios":width_ratios})
    # plot matrix of each element
    gs=axes[0][0].get_gridspec()
    ax_global = fig.add_subplot(gs[:, cols+1:])
    ax_equal = fig.add_subplot(gs[:,cols])
    for ax in axes[:,cols:]:
        for a in ax:
            a.remove()
    for i in range(0,rows):
        for j in range(0,cols):
            ind_el = j+i*cols
            element=el2nod[ind_el]
            ax_local=axes[i][j]
            for ax in [ax_local,ax_global]:
                ax.axis('scaled')
                ax.set_xlim(-0.5,nnod-0.5)
                ax.set_ylim(-0.5,nnod-0.5)
                ax.xaxis.set_ticks(np.arange(0,nnod+1,1)-0.5)
                ax.yaxis.set_ticks(np.arange(0,nnod+1,1)-0.5)
                ax.grid(axis='both',which='major',clip_on=False)
                ax.tick_params(axis='both',which='both',color='None')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                for spine in ax.spines:
                    ax.spines[spine].set_visible(False)
                ax.invert_yaxis()
                ax_local.set_title('Element %d'%(ind_el))
                for nod_col in element:
                    y=[nod_col-0.5, nod_col-0.5, nod_col+0.5, nod_col+0.5]
                    for nod_row in element:
                        x=[nod_row-0.5, nod_row+0.5, nod_row+0.5, nod_row-0.5]
                        ax.fill(x,y,color=colors[ind_el%len(colors)],alpha=0.8)
                if(j<(cols-1)):
                    ax_local.plot(1.06,0.5,'P',color='k',transform=ax_local.transAxes,clip_on=False)
                if((i>0) & (j==0)):
                    ax_local.plot(-0.06,0.5,'P',color='k',transform=ax_local.transAxes,clip_on=False)
    ax_equal.plot([0,1],[0.505,0.505],lw=3,color='k',transform=ax_equal.transAxes,clip_on=False)
    ax_equal.plot([0,1],[0.495,0.495],lw=3,color='k',transform=ax_equal.transAxes,clip_on=False)
    ax_equal.axis('off')
    # for i,element in enumerate(el2nod):
    #     ind=np.append(element,element[0])
    #     axes[0][0].fill(x[ind],y[ind],color=colors[i%len(colors)])
        # ax.text(x[element].mean(),y[element].mean(),'Element %d'%(i),ha='center',va='center',bbox={'fc':'lightgray','ec':'None','boxstyle':'round'})
    # axes[0].remove()
    for fmt in fmt_figs:
        figurename=str('%s/%s.%s'%(figpath,figname,fmt))
        plt.savefig(figurename, bbox_inches='tight')
# shape function
def shapes_Quad(xi, eta, fem='Q1'):
    if(fem=='Q1'):
        N1, f1, x1, y1 = 0.25*(1-xi)*(1-eta), '$\\frac{1}{4}(1-\\xi)(1-\\eta)$', -1, -1
        N2, f2, x2, y2 = 0.25*(1+xi)*(1-eta), '$\\frac{1}{4}(1+\\xi)(1-\\eta)$', 1, -1
        N3, f3, x3, y3 = 0.25*(1+xi)*(1+eta), '$\\frac{1}{4}(1+\\xi)(1+\\eta)$', 1, 1
        N4, f4, x4, y4 = 0.25*(1-xi)*(1+eta), '$\\frac{1}{4}(1-\\xi)(1+\\eta)$', -1, 1
        return {'name':'4-node quadrilateral', 'shapes':[N1, N2, N3, N4],'formula':[f1,f2,f3,f4], 'x':[x1,x2,x3,x4],'y':[y1,y2,y3,y4]}
    elif(fem=='Q2'):
        N1, f1, x1, y1 = 0.25*(xi-1)*(eta-1)*xi*eta, '$\\frac{1}{4}(\\xi-1)(\\eta-1)\\xi\\eta$', -1, -1
        N2, f2, x2, y2 = 0.25*(xi+1)*(eta-1)*xi*eta, '$\\frac{1}{4}(\\xi+1)(\\eta-1)\\xi\\eta$', 1, -1
        N3, f3, x3, y3 = 0.25*(xi+1)*(eta+1)*xi*eta, '$\\frac{1}{4}(\\xi+1)(\\eta+1)\\xi\\eta$', 1, 1
        N4, f4, x4, y4 = 0.25*(xi-1)*(eta+1)*xi*eta, '$\\frac{1}{4}(\\xi-1)(\\eta+1)\\xi\\eta$', -1, 1
        N5, f5, x5, y5 = 0.5*(1-xi**2)*eta*(eta-1), '$\\frac{1}{2}(1-\\xi^2)\\eta(\\eta-1)$', 0, -1
        N6, f6, x6, y6 = 0.5*(1-eta**2)*xi*(xi+1),  '$\\frac{1}{2}(1-\\eta^2)\\xi(\\xi+1)$',  1, 0
        N7, f7, x7, y7 = 0.5*(1-xi**2)*eta*(eta+1), '$\\frac{1}{2}(1-\\xi^2)\\eta(\\eta+1)$', 0, 1
        N8, f8, x8, y8 = 0.5*(1-eta**2)*xi*(xi-1),  '$\\frac{1}{2}(1-\\eta^2)\\xi(\\xi-1)$', -1, 0
        N9, f9, x9, y9 = (1-xi**2)*(1-eta**2), '$(1-\\xi^2)(1-\\eta^2)$', 0, 0
        return {'name':'9-node quadrilateral', 
                'shapes': [N1, N2, N3, N4, N5, N6, N7, N8, N9],
                'formula':[f1, f2, f3, f4, f5, f6, f7, f8, f9], 
                'x':      [x1, x2, x3, x4, x5, x6, x7, x8, x9],
                'y':      [y1, y2, y3, y4, y5, y6, y7, y8, y9]}
    else:
        print('Undefined fem')
        exit()
def triangulate_polygon(x,y,num=40):
    # using gmsh to triangulate a polygon
    geofile,mshfile='tmp.geo','tmp.msh'
    fpout=open(geofile,'w')
    lc_x,lc_y=(x.max()-x.min())/num,(y.max()-y.min())/num
    lc=np.max([lc_x, lc_y])
    nnod=len(x)
    fpout.write('lc=%f;\n'%(lc))
    for i in range(0,nnod):
        fpout.write('Point(%d)={%f, %f, 0, lc};\n'%(i,x[i],y[i]))
    for i in range(0,nnod-1):
        fpout.write('Line(%d)={%d, %d};\n'%(i,i,i+1))
    fpout.write('Line(%d)={%d, %d};\n'%(nnod-1,nnod-1,0))
    # fpout.write('Line(%d)={%d, %d};\n'%(i,i,i+1))
    fpout.write('Curve Loop(1) = {%d:%d};\n'%(0,nnod-1))
    fpout.write('Plane Surface(1) = {1};\n')
    fpout.close()
    # meshing
    cmd='gmsh %s -2 -o %s'%(geofile,mshfile)
    print(cmd)
    os.system(cmd)
    # read mesh
    mesh=meshio.read(mshfile)
    triangles=[]
    for cells in mesh.cells:
        if(cells.type=='triangle'):
            triangles=cells.data
    # remove tmp files
    os.system('rm %s %s'%(geofile,mshfile))
    return {'GCOORD':mesh.points[:,0:2],'EL2NOD':triangles}
def shapes_Triangle(x,y, fem='Q1'):
    # see also https://www.iue.tuwien.ac.at/phd/nentchev/node25.html
    # see also http://www.sd.ruhr-uni-bochum.de/downloads/Shape_funct.pdf
    J=np.array([[1, x[0], y[0]],[1, x[1], y[1]],[1, x[2], y[2]]],dtype=float)
    detJ=np.linalg.det(J)
    mesh=triangulate_polygon(x,y)
    xi,eta=mesh['GCOORD'][:,0],mesh['GCOORD'][:,1]
    tris=mesh['EL2NOD']
    # shape function
    A = np.abs(detJ/2) # area of the triangle element
    xi0, xi1, xi2 = np.zeros_like(xi),np.zeros_like(xi),np.zeros_like(xi)
    for i in range(0,len(xi)):
        x0,y0=xi[i],eta[i]
        J0=np.array([[1, x0, y0],[1, x[1], y[1]],[1, x[2], y[2]]],dtype=float)
        J1=np.array([[1, x0, y0],[1, x[0], y[0]],[1, x[2], y[2]]],dtype=float)
        J2=np.array([[1, x0, y0],[1, x[1], y[1]],[1, x[0], y[0]]],dtype=float)
        A0=np.abs(np.linalg.det(J0))/2.0
        A1=np.abs(np.linalg.det(J1))/2.0
        A2=np.abs(np.linalg.det(J2))/2.0
        xi0[i]=A0/A 
        xi1[i]=A1/A 
        xi2[i]=A2/A
    N,f,x_nods,y_nods,ind_axis,name=[],[],[],[],[],''
    if(fem=='Q1'):
        N1, N2, N3 = 1-xi1-xi2, xi1, xi2
        f1,f2,f3='$1-\\xi_1-\\xi_2$','$\\xi_1$','$\\xi_2$'
        N,f=[N1,N2,N3],[f1,f2,f3]
        x_nods,y_nods=x,y
        ind_axis=[2,3,4]
        name='3-node triangle'
    elif(fem=='Q2'):
        N1, N2, N3 = (1-xi1-xi2)*(1-2*xi1-2*xi2), xi1*(2*xi1-1), xi2*(2*xi2-1)
        f1,f2,f3='$(1-\\xi_1-\\xi_2)(1-2\\xi_1-2\\xi_2)$', '$\\xi_1(2\\xi1-1)$', '$\\xi_2(2\\xi_2-1)$'
        ci=4
        N4, N5, N6 = ci*xi1*(1-xi1-xi2), ci*xi1*xi2, ci*xi2*(1-xi1-xi2)
        f4,f5,f6 = '$4\\xi_1(1-\\xi_1-\\xi_2$)', '$4\\xi_1\\xi_2$', '$4\\xi_2(1-\\xi_1-\\xi_2)$'
        N,f=[N1,N2,N3,N4, N5, N6],[f1,f2,f3,f4,f5,f6]
        x_nods=[x[0], x[1], x[2], (x[0]+x[1])/2.0, (x[1]+x[2])/2.0, (x[2]+x[0])/2.0]
        y_nods=[y[0], y[1], y[2], (y[0]+y[1])/2.0, (y[1]+y[2])/2.0, (y[2]+y[0])/2.0]
        ind_axis=[2,3,4,6,7,8]
        name='6-node triangle'
    elif(fem=='Q3'):
        c1=0.5
        N1 = c1*(1-3*xi1-3*xi2)*(2-3*xi1-3*xi2)*(1-xi1-xi2)
        c2=4.5
        N2 = c2*xi1*(xi1-1/3.0)*(xi1-2.0/3.0)
        N3 = c2*xi2*(xi2-1/3.0)*(xi2-2.0/3.0)
        N4 = c2*xi1*(2-3*xi1-3*xi2)*(1-xi1-xi2)
        N5 = c2*xi1*xi2*(3*xi1-1)
        N6 = c2*xi2*(3*xi2-1)*(1-xi1-xi2)
        N7 = c2*xi1*(3*xi1-1)*(1-xi1-xi2)
        N8 = c2*xi1*xi2*(3*xi2-1)
        N9 = c2*xi2*(2-3*xi1-3*xi2)*(1-xi1-xi2)
        N10= 27*xi1*xi2*(1-xi1-xi2)
        N=[N1,N2,N3,N4,N5,N6,N7,N8,N9,N10]
        f=['']*10

        x_nods=[x[0],x[1],x[2], x[0]+(x[1]-x[0])/3.0, x[1]+(x[2]-x[1])/3.0, x[2]+(x[0]-x[2])/3.0, 
                x[0]+2*(x[1]-x[0])/3.0, x[1]+2*(x[2]-x[1])/3.0, x[2]+2*(x[0]-x[2])/3.0, x.mean()]
        y_nods=[y[0],y[1],y[2], y[0]+(y[1]-y[0])/3.0, y[1]+(y[2]-y[1])/3.0, y[2]+(y[0]-y[2])/3.0, 
                y[0]+2*(y[1]-y[0])/3.0, y[1]+2*(y[2]-y[1])/3.0, y[2]+2*(y[0]-y[2])/3.0, y.mean()]
        ind_axis=[2,3,4,6,7,8,10,11,12,9]
        name='10-node triangle'
    else:
        print('Undefined fem')
        exit()
    return {'name':name, 'shapes':N,'formula':f, 
            'x':x_nods,'y':y_nods, 'ind_axis':ind_axis,
            'triangles':mpl.tri.Triangulation(xi,eta,tris)}
def write2VTU(vtufile, xx,yy,zz):
    ncols=xx.shape[0]
    nrows=xx.shape[1]
    x=xx.reshape(-1)
    y=yy.reshape(-1)
    z=zz.reshape(-1)
    npoints=len(x)
    nCells=(ncols-1)*(nrows-1)
    VTK_CELLTYPE=9 #四边形
    np_per_cell=4
    fpout=open(vtufile,'w')
    fpout.write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n')
    fpout.write('  <UnstructuredGrid>\n')
    fpout.write('    <Piece NumberOfPoints="%.0f" NumberOfCells="%.0f">\n'%(npoints,nCells))
    fpout.write('      <PointData Scalars="ShapeFunction">\n')
    fpout.write('        <DataArray type="Float64" Name="ShapeFunction" format="ascii">\n')
    fpout.write('          ')
    # Info('Writing ShapeFunction field ...')
    for i in range(0,len(z)):
        fpout.write('%f '%(z[i]))
    fpout.write('\n        </DataArray>\n')
    fpout.write('      </PointData>\n')
    fpout.write('      <CellData>\n')
    fpout.write('      </CellData>\n')
    fpout.write('      <Points>\n')
    fpout.write('        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">\n')
    # Info('Writing xyz data ...')
    for i in range(0,len(x)):
        fpout.write('          %f %f %f\n'% (x[i],y[i],z[i]))
    fpout.write('        </DataArray>\n')
    fpout.write('      </Points>\n')
    fpout.write('      <Cells>\n')
    fpout.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
    # Info('Writing connectivity ...')
    for nrow in range(0,nrows-1):
        for ncol in range(0,ncols-1):
            LL=ncol + nrow*ncols
            fpout.write('          %.0f %.0f %.0f %.0f\n'%(LL, LL+1, LL+1+ncols, LL+ncols))
    fpout.write('        </DataArray>\n')
    fpout.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
    fpout.write('          ')
    # Info('Writing offsets ...')
    for i in range(0,nCells):
        fpout.write('%.0f '%(i*np_per_cell))
    fpout.write('        </DataArray>\n')
    fpout.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
    fpout.write('          ')
    # Info('Writing cell type ...')
    for i in range(0,nCells):
        fpout.write('%.0f '%(VTK_CELLTYPE))
    fpout.write('        </DataArray>\n')
    fpout.write('      </Cells>\n')
    fpout.write('    </Piece>\n')
    fpout.write('  </UnstructuredGrid>\n')
    fpout.write('</VTKFile>\n')
    # Info('xyz to vtu Done')
    fpout.close()
    # Info('Converting ASCII to binary')
    os.system('meshio-binary '+vtufile)
def FE_Quad(fem='Q1',cmap=cm.GnBu):
    xi=np.linspace(-1,1,50)
    eta=np.linspace(-1,1,50)
    xi,eta=np.meshgrid(xi,eta)
    # get shape function of each node
    shapefunc=shapes_Quad(xi,eta,fem=fem)
    nip=len(shapefunc['shapes'])
    cols=5
    rows=int(nip/5)+1
    fig = plt.figure(figsize=(16,8))
    # plot element
    ax_el=fig.add_subplot(rows,cols,1, facecolor='None')
    ax_el.axis('scaled')
    ax_el.set_xlim(-1,1)
    ax_el.set_ylim(-1,1)
    ax_el.fill([-1,1,1,-1],[-1,-1,1,1],color='lightgray',ec='k',lw=3,closed=True)
    ax_el.text(0.5,0.25,'%s\nElement'%(shapefunc['name']),ha='center',va='center',transform=ax_el.transAxes,fontweight='bold',fontsize=14)
    ax_el.axis('off')
    # plot shape function of each node
    for i, N, formula, x0,y0 in zip(range(0,nip),shapefunc['shapes'],shapefunc['formula'],shapefunc['x'],shapefunc['y']):
        # write to VTU
        # write2VTU('N%d.vtu'%(i),xx,yy,N)
        # node on element
        l_node,=ax_el.plot(x0,y0,marker='o',mec='w',clip_on=False,ms=15)
        ax_el.text(x0,y0,'%d'%(i),color='w',ha='center',va='center',fontweight='bold')
        # shape function
        ax = fig.add_subplot(rows,cols,i+2, projection='3d',facecolor='None')
        ls = LightSource(270, 45)
        norm = mpl.colors.TwoSlopeNorm(vmin=N.min()-1E-5, vcenter=0, vmax=N.max())
        rgb = ls.shade(N, cmap=cmap, norm=norm, vert_exag=0.1, blend_mode='soft')
        CS=ax.plot_surface(xi,eta, N, rstride=1, cstride=1, facecolors=rgb, alpha=0.9, linewidth=0)
        # CS=ax.plot_surface(xi,eta, N, rstride=1, cstride=1, cmap=cmap, norm=norm,linewidth=0.01,edgecolor='k', shade=False,vmin=-1,vmax=1)
        ax.set_xlabel('$\\xi$',labelpad=0)
        ax.set_ylabel('$\\eta$',labelpad=0)
        ax.set_zlabel('Shape function',labelpad=-3)
        ax.zaxis.set_minor_locator(MultipleLocator(0.1))
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_zlim(0,1)
        ax.zaxis.set_ticks([0,1])
        ax.set_title('N$_{\mathregular{%d}}=$%s'%(i,formula),color=l_node.get_color(),fontweight='bold',fontsize=14)
        # 重新自定义坐标轴属性
        niceAxis(ax,fill_pane=False,label3D=True,fs_label=0.1, scaled=False)
    plt.subplots_adjust(wspace=0)
    for fmt in fmt_figs:
        figurename=str('%s/shapeFunction_2D_%s.%s'%(figpath,fem,fmt))
        plt.savefig(figurename, bbox_inches='tight')
def FE_Tri(fem='Q1',cmap=cm.Spectral_r):
    mesh=loadMesh2D('mesh/unstructure.msh')
    tri=mesh['EL2NOD'][0]
    x,y=mesh['GCOORD'][:,0],mesh['GCOORD'][:,1]
    x_tri,y_tri=x[tri],y[tri]
    shapefunc=shapes_Triangle(x_tri,y_tri,fem=fem)

    nip=len(shapefunc['shapes'])
    triangles=shapefunc['triangles']
    # plot 
    cols=4
    rows=int(nip/cols)+1
    fig = plt.figure(figsize=(16,4*rows))
    # plot element
    ax_el=fig.add_subplot(rows,cols,1, facecolor='None')
    ax_el.axis('scaled')
    ax_el.set_xlim(x_tri.min(),x_tri.max())
    ax_el.set_ylim(y_tri.min(),y_tri.max())
    ax_el.fill(x_tri,y_tri,color='lightgray',ec='k',lw=3,closed=True,clip_on=False)
    ax_el.text(0,1,'%s\nElement'%(shapefunc['name']),ha='left',va='top',transform=ax_el.transAxes,fontweight='bold',fontsize=14)
    ax_el.axis('off')
    for i, ind_axis,N, formula, x0,y0 in zip(range(0,nip),shapefunc['ind_axis'],shapefunc['shapes'],shapefunc['formula'],shapefunc['x'],shapefunc['y']):
        # write to VTU
        # write2VTU('N%d.vtu'%(i),xx,yy,N)
        # node on element
        l_node,=ax_el.plot(x0,y0,marker='o',mec='w',clip_on=False,ms=15)
        ax_el.text(x0,y0,'%d'%(i),color='w',ha='center',va='center',fontweight='bold')
        # shape function
        ax = fig.add_subplot(rows,cols,ind_axis, projection='3d',facecolor='None')
        ax.view_init(elev=45, azim=-70)
        norm = mpl.colors.TwoSlopeNorm(vmin=N.min()-1E-5, vcenter=0, vmax=N.max())
        ax.plot_trisurf(triangles.x, triangles.y, N, triangles=triangles.triangles, 
            cmap=cmap,edgecolor=(0,0,0,0.5),linewidth=0,alpha=0.9,zorder=11)
        # plot element
        ax.plot(x_tri[[0,1,2]],y_tri[[0,1,2]],[0]*3,color='k',zorder=11)
        ax.plot(x_tri[[2,0]],y_tri[[2,0]],[0,0],color='k',ls='dashed')
        ax.plot(x0,y0,0,marker='o',mec='w',clip_on=False,ms=15,color=l_node.get_color(),zorder=12,alpha=0.8)
        # poly = PolyCollection([[(x_tri[0],y_tri[0]), (x_tri[1],y_tri[1]), (x_tri[2],y_tri[2])]], 
        #     facecolors='lightgray', alpha=0.7,edgecolors='k')
        # ax.add_collection3d(poly, zs=0, zdir='z')
        # make 3d cylinder
        for tmpx,tmpy in zip(x_tri,y_tri):
            ax.plot([tmpx, tmpx],[tmpy,tmpy],[0,1],ls='solid',color='g',zorder=10)
        ax.plot([x0,x0],[y0,y0],[0,1],color='k',lw=2,zorder=10)
        ax.plot(x_tri[[0,1,2,0]],y_tri[[0,1,2,0]],[1]*4,color='g',zorder=10)
        # set axis
        ax.set_xlim(x_tri.min(),x_tri.max())
        ax.set_ylim(y_tri.min(),y_tri.max())
        ax.set_zlim(0,1)
        ax.axis('off')
        ax.set_title('N$_{%d}$ = %s'%(i,formula),color=l_node.get_color(),fontweight='bold',fontsize=14)
        
    plt.subplots_adjust(wspace=0)
    for fmt in fmt_figs:
        figurename=str('%s/shapeFunction_2D_triangle_%s.%s'%(figpath,fem,fmt))
        plt.savefig(figurename, bbox_inches='tight')
def connectivity():
    # 1. structured mesh with structured indexing
    mesh=createMesh_2D()
    plotMesh_2D(mesh,extend_x=1/3.0)
    plotConnectivityMatrix_2D(mesh)
    # structured mesh with unstructured indexing
    mesh=loadMesh2D('mesh/structure.msh')
    plotMesh_2D(mesh,figname='mesh2D_structured_usi',extend_x=1/3.0)
    plotConnectivityMatrix_2D(mesh,figname='Matrix2D_structured_usi')
    # unstructured mesh
    mesh=loadMesh2D('mesh/unstructure.msh')
    plotMesh_2D(mesh,figname='mesh2D_unstructured')
    plotConnectivityMatrix_2D(mesh,figname='Matrix2D_unstructured')

def main(argv):
    # # 1. 1D element
    # Linear1D()
    # # 2. 2D Quad element
    # FE_Quad(fem='Q1',cmap=cm.plasma)
    # FE_Quad(fem='Q2',cmap=cm.Spectral_r)
    # # 3. 2D triangle element
    # FE_Tri(fem='Q1')
    # FE_Tri(fem='Q2')
    # FE_Tri(fem='Q3')
    # # 4. matrix connectivity of 2D mesh 
    connectivity()
if __name__ == '__main__':
    sys.exit(main(sys.argv))