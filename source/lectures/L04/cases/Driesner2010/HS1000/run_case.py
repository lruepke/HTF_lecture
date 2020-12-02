#!/Users/zguo/.pyenv/shims/python
# -*-coding:utf-8-*-
# Run a case in background
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2020/10/24, GEOMAR
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
import psutil
import os
import os.path
import numpy
import sys
def usage(argv):
    basename = argv[0].split('/')
    basename = basename[len(basename)-1]
    description='Run a case in background'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2020/10/24, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' caseDir'+C_DEFAULT)
    print('='*len(head))


def main(argv):
    argc=len(argv)
    if(not argc==2):
        usage(argv)
        exit(0)
    caseDir=argv[1]
    # run
    pidfile = str('%s/pid'%caseDir)
    # 1. check
    if(os.path.exists(pidfile)):
        pid = int(numpy.loadtxt(pidfile))
        if(psutil.pid_exists(pid)):
            print(C_BLUE+"Solver is still running in %s%s%s : "%(C_RED,caseDir,C_DEFAULT)+C_GREEN+str(pid)+C_DEFAULT)
            print(C_RED+"you can not run a second same model in the same folder"+C_DEFAULT)
            exit(0)
    cmd=str('cd %s && nohup ./run.sh >log.log & echo $! > %s'%(caseDir,pidfile))
    os.system(cmd)
    print('Case is running in %s%s%s : %s%s%s'%(C_BLUE,caseDir,C_DEFAULT, C_GREEN,str(int(numpy.loadtxt(pidfile))),C_DEFAULT))

if __name__ == '__main__':
    sys.exit(main(sys.argv))