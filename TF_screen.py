#%matplotlib inline
import os
import sys
import tifffile as tiff
import pandas as pd
import matplotlib
# For visualising
#matplotlib.use('MacOSX') 
# For saving
#matplotlib.use('Agg') 
matplotlib.use('MacOSX') #TkAgg
import matplotlib.pyplot as plt

import numpy as np
import itertools as IT
import time
import glob
import skimage
from skimage import io as skio
from skimage import color

sys.path.append("./src")
from utilities import *

#plt.interactive(True)
matplotlib.rcParams.update({'font.size': 12})
np.set_printoptions(precision=3)



#os.chdir("/Users/Povilas/Dropbox/Projects/Metformin_TF_acs-2/")
#sourceloc = "Users/Povilas/Dropbox/Projects/Metformin_TF_acs-2/"

os.chdir("/home/pnorv/Dropbox/Projects/Metformin_TF_acs-2/")
sourceloc = "home/pnorv/Dropbox/Projects/Metformin_TF_acs-2/"

odir="."


structure=pd.read_csv('./Allfiles_summary.csv')
fileinfo=readTF('./Allfiles_annotated.csv')

genes=fileinfo.keys()
#genes=['OP50-C']

repflds=['Rep1_23-8-17','Rep2_24-8-17']


for gene in genes:
    if gene=='OP50-C':
        nflds=6*2
    else:
        nflds=4*2
    gval=fileinfo[gene]
    gmax = gval['Max']
    print gmax
    
    fwidth=6*nflds
    fheight=5*gmax
    
    print nflds, gmax, fwidth, fheight
    
    
    fig, axes = plt.subplots(nrows=gmax, ncols=nflds, figsize=(fwidth,fheight), dpi=300)
    
    ofile="{}_summary.pdf".format(gene)
    
    print ofile
    col=0
    for rep in [1,2]:
        repval = gval[rep]
        repfld=repflds[rep-1]
        for tp in ['C','T']:
            tpval=repval[tp]
            
            for fld, fldvals in tpval.iteritems():
                print "Gene: {} Replicate: {} Type: {} Folder: {}".format(gene, rep, tp, fld)
                fldcontsel=fldvals['ByFIN']
                print fldcontsel.keys()
                
                
                row=0
                for fli in fldcontsel.keys():
                    fl=fldcontsel[fli]
                    
                    flfull="./{}/{}/{}".format(repfld,fld,fl)
                    print flfull #, os.path.isfile(flfull)

                    image = tiff.imread(flfull)
                    imghsv = color.rgb2hsv(image)
                    imgrgb = img_as_float(image)
                    v=imghsv[:,:,2]
                    comb,labeled_worms,contours,worms=labeller4(v)
                    plotcontours(comb,axes[row,col],ylabel=fl,title="{} - {}".format(rep,fld),plotcontour=False,plotlabels=False,white=False)
                    plotcontours(imgrgb,axes[row,col+1],labeled_worms,plotlabels=False)
                    row+=1
                col+=2

    #fig.tight_layout()
    fig.savefig(ofile, bbox_inches='tight')
    plt.close(fig)
