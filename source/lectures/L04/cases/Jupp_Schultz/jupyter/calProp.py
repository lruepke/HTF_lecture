#!/Users/zguo/.pyenv/shims/python
# -*-coding:utf-8-*-
# Calculate water properties using iapws package
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2020/12/01, GEOMAR
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
import numpy as np 
import iapws
from console_progressbar import ProgressBar

def usage(argv):
    basename = argv[0].split('/')
    basename = basename[len(basename)-1]
    description='Calculate water properties using iapws package'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2020/12/01, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' '+C_DEFAULT)
    print('='*len(head))

def calProps():
    T=np.linspace(0,1000, 1000)
    P=np.linspace(1, 1000, 1000) #bar
    P=P*1E5
    maxnum=len(T)
    pb = ProgressBar(total=maxnum,prefix=C_BLUE+'Progress: '+C_DEFAULT, suffix=' Completed'+C_DEFAULT, decimals=3, length=50, fill=C_GREEN+'#', zfill=C_DEFAULT+'-')
    f_T=open('data_iapws/TT.dat','w')
    f_P=open('data_iapws/PP.dat','w')
    f_rho=open('data_iapws/RHO.dat','w')
    f_h=open('data_iapws/H.dat','w')
    f_mu=open('data_iapws/MU.dat','w')
    for i in range(0,len(P)):
        for j in range(0,len(T)):
            steam=iapws.IAPWS97(P=P[i]/1E6, T=t[j]+273.15)
            f_T.write('%.6E '%(T[j]))
            f_P.write('%.6E '%(P[i]))
            f_rho.write('%.6E '%(steam.rho))
            f_h.write('%.6E '%(steam.h))
            f_mu.write('%.6E '%(steam.mu))
        f_T.write('\n')
        f_P.write('\n')
        f_rho.write('\n')
        f_h.write('\n')
        f_mu.write('\n')
        pb.print_progress_bar(i+1)
    f_rho.close()
    f_mu.close()
    f_h.close()
def main(argv):

    calProps()

if __name__ == '__main__':
    sys.exit(main(sys.argv))