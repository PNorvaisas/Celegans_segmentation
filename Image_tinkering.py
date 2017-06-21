import os
import sys

import skimage
import matplotlib.pyplot as plt
#%matplotlib inline
import numpy as np

import matplotlib

#plt.interactive(True)

from skimage import io as skio

from skimage import color
from skimage import data
from skimage import img_as_float
from skimage import filters
from skimage.util.dtype import dtype_range
from skimage.util import img_as_ubyte


from skimage.feature import canny
from scipy import ndimage as ndi
from skimage.filters import sobel
from skimage.morphology import watershed


from skimage.morphology import disk, opening, dilation, square
from skimage.morphology import erosion, white_tophat, black_tophat, closing
from skimage import exposure

from skimage.filters.rank import median, mean
from scipy import misc, ndimage


import tifffile as tiff

import textwrap as tw

import csv

import pandas

def filename(ifile):
	if ifile.split('.')[0]=='':
		ipat=''
		iname=''
		itype=ifile.split('.')[1]
	else:
		if "\\" in ifile.split('.')[0]:
			sep="\\"
		elif "/" in ifile.split('.')[0]:
			sep="/"
		else:
			ipat=''
			iname=ifile.split('.')[0]
			itype=ifile.split('.')[1]
			return ipat, iname, itype
		allpath=ifile.split('.')[0]
		iname=allpath.split(sep)[-1]
		ipath=allpath.split(sep)[:-1]
		ipat='/'.join(ipath)
		itype=ifile.split('.')[1]
	return ipat, iname, itype


def thresholding(file,thresholds,sizes):
    #file = '20170503/image_011.tif'
    #Changed from skio.imread, which could not handle 16bit TIFF
    image = tiff.imread(file)

    vo_size,vd_size,ho_size,hd_size=sizes
    
    vmin,vmax,hmin,hmax=thresholds
    #image.shape
    #image.dtype
    imghsv=color.rgb2hsv(image)
    
    imgrgb=color.hsv2rgb(imghsv)

    h =imghsv[:,:,0]
    #s =imghsv[:,:,1]
    v =imghsv[:,:,2]

    vmask = (v > vmin).astype(np.uint8);  # Thresholding in the Brightness
    vmask = (vmask< vmax).astype(np.uint8);

    disk_oelem = disk(vo_size) # Remove small regions
    disk_delem = disk(vd_size) # Join regions
    vopened = opening(vmask, selem=disk_oelem)
    vdilated = dilation(vopened, selem=disk_delem)

    #plt.imshow(vopened); plt.title('Brightness mask')
    #plt.imshow(vdilated); plt.title('Brightness mask')
    hmask = (h > hmin).astype(np.uint8);  # Thresholding in the Hue-channel
    hmask = (hmask< hmax).astype(np.uint8);

    disk_oelem = disk(ho_size) # Remove small regions
    disk_delem = disk(hd_size) # Remove small regions
    hopened = opening(hmask, selem=disk_oelem)
    hdilated = dilation(hopened, selem=disk_delem)

    #plt.imshow(hopened); plt.title('Hue mask')
    #plt.imshow(hdilated); plt.title('Hue mask')

    img2 = imghsv.copy()
    img2[hdilated.astype(bool), :] = 0; # Set the pixels to zero by Hue mask
    img2[vdilated.astype(bool), :] = 0; # Set the pixels to zero by Lightness mask

    img2rgb=color.hsv2rgb(img2)

    return image,imgrgb,vdilated,hdilated,img2,img2rgb



def thresholding2(imghsv,thresholds,sizes,delworms):
    #file = '20170503/image_011.tif'
    #Changed from skio.imread, which could not handle 16bit TIFF
    
    if len(sizes)==2:
        vc_size,hc_size=sizes
    else:
        print 'Fix your sizes parameters!'
        vc_size=sizes[0]
        hc_size=sizes[1]
    
    vmin,vmax,hmin,hmax=thresholds
    #image.shape
    #image.dtype


    h =imghsv[:,:,0]
    #s =imghsv[:,:,1]
    v =imghsv[:,:,2]
    
    vmask = (v > vmin).astype(np.float64);  # Thresholding in Brightness
    vmask = (vmask< vmax).astype(np.float64)
    
    hmask = (h > hmin).astype(np.float64);  # Thresholding in Hue
    hmask = (hmask< hmax).astype(np.float64)
    
    #Works reallt well
    #Morphological filtering
    v_closing = closing(vmask, selem=disk(vc_size))
    h_closing = closing(hmask, selem=disk(hc_size))
    #plt.imshow(hopened); plt.title('Hue mask')
    #plt.imshow(hdilated); plt.title('Hue mask')

    img2 = imghsv.copy()
    img2[h_closing.astype(bool), :] = 0; # Set the pixels to zero by Hue mask
    img2[v_closing.astype(bool), :] = 0; # Set the pixels to zero by Lightness mask
    for worm in delworms:
        x1,y1,x2,y2=worm
        xs=[x1,x2]
        ys=[y1,y2]
        
        img2[min(ys):max(ys),min(xs):max(xs), :] = 0

    img2rgb=color.hsv2rgb(img2)

    return v_closing,h_closing,img2,img2rgb


def thresholding3(file,thresholds,sizes,alpha,delworms):
    #file = '20170503/image_011.tif'
    #Changed from skio.imread, which could not handle 16bit TIFF
    image = tiff.imread(file)
    if len(sizes)==2:
        vc_size,hc_size=sizes
    else:
        print 'Fix your sizes parameters!'
        vc_size=sizes[0]
        hc_size=sizes[1]
    
    vmin,vmax,hmin,hmax=thresholds
    #image.shape
    #image.dtype
    imghsv=color.rgb2hsv(image)
    imgrgb=color.hsv2rgb(imghsv)

    h =imghsv[:,:,0]
    #s =imghsv[:,:,1]
    v =imghsv[:,:,2]
    
    v_filter_blurred_l = ndimage.gaussian_filter(v, 1)
    vsharp = v + alpha * (v - v_filter_blurred_l)

    
    h_filter_blurred_l = ndimage.gaussian_filter(h, 1)
    hsharp = h + alpha * (h - h_filter_blurred_l)
    
    
    
    vmask = (vsharp > vmin).astype(np.float64);  # Thresholding in Brightness
    vmask = (vmask< vmax).astype(np.float64)
    
    hmask = (hsharp > hmin).astype(np.float64);  # Thresholding in Hue
    hmask = (hmask< hmax).astype(np.float64)
    
    #Works reallt well
    #Morphological filtering
    v_closing = closing(vmask, selem=disk(vc_size))
    h_closing = closing(hmask, selem=disk(hc_size))
    #plt.imshow(hopened); plt.title('Hue mask')
    #plt.imshow(hdilated); plt.title('Hue mask')

    img2 = imghsv.copy()
    img2[h_closing.astype(bool), :] = 0; # Set the pixels to zero by Hue mask
    img2[v_closing.astype(bool), :] = 0; # Set the pixels to zero by Lightness mask
    for worm in delworms:
        x1,y1,x2,y2=worm
        xs=[x1,x2]
        ys=[y1,y2]
        
        img2[min(ys):max(ys),min(xs):max(xs), :] = 0

    img2rgb=color.hsv2rgb(img2)

    return image,imgrgb,v_closing,h_closing,img2,img2rgb




def slice_list(input, size):
    input_size = len(input)
    slice_size = input_size / size
    remain = input_size % size
    result = []
    iterator = iter(input)
    for i in range(size):
        result.append([])
        for j in range(slice_size):
            result[i].append(iterator.next())
        if remain:
            result[i].append(iterator.next())
            remain -= 1
    return result

class NestedDict(dict):
	def __getitem__(self, key):         
		if key in self: return self.get(key)
		return self.setdefault(key, NestedDict())

def numerize(s):
	try:
		if s=='NAN':
			return s
		float(s)
		if float(s).is_integer():
			return int(float(s))
		elif float(s)==0:
			return float(s)
		else:
			return float(s)

	except ValueError:
		return s
    
def indx2well(ind,start=0,rowln=12):
    
    indt=ind-start
    rows=['A','B','C','D','E','F','G','H']
    row=indt//rowln
    col=indt-row*rowln+1
    well='{}{}'.format(rows[row],col)
    return well


def readmet(ifile):
	print ifile
	nutrients=NestedDict()
	rdr=csv.reader(open(ifile,'r'), delimiter=',')
	data=[ln for ln in rdr]
	headers=data[0]
	headin={ hd : headers.index(hd) for hd in headers}
	nec=['Metabolite','EcoCycID','Plate','Well','Index']
	if all(n in headers for n in nec):
		print 'Necessary headers found!'
	else:
		print 'Missing essential headers in metabolites file!'
		print headers
		sys.exit(0)
	for ln in data[1:]:
		pl=str(numerize(ln[headin['Plate']])).strip().encode('ascii','ignore')
		wl=str(numerize(ln[headin['Well']])).strip().encode('ascii','ignore')
		for hd in headin.keys():
			nutrients[pl][wl][hd]=str(numerize(ln[headin[hd]])).strip().encode('ascii','ignore')
	return nutrients




def readdel(ifile):
    print ifile
    
    deletions=NestedDict()
    rdr=csv.reader(open(ifile,'r'), delimiter=',')
    data=[ln for ln in rdr]
    headers=data[0]
    headin={ hd : headers.index(hd) for hd in headers}
    nec=['Replicate','Plate','Index','File','X1','Y1','X2','Y2']#,'Well'
    if all(n in headers for n in nec):
        print 'Necessary headers found!'
    else:
        print 'Missing essential headers in deletions file!'
        print headers
        sys.exit(0)
    for ln in data[1:]:
        #print ln
        #Reading
        rep=ln[headin['Replicate']].strip().encode('ascii','ignore')
        pl=ln[headin['Plate']].strip().encode('ascii','ignore')
        fl=ln[headin['File']].strip().encode('ascii','ignore')
        indx=numerize(ln[headin['Index']])
        coords=[numerize(ln[headin[hdr]]) for hdr in ['X1','Y1','X2','Y2']]
        #Storing
        deletions[rep][pl][indx]['File']=fl
        
        if 'Worms' in deletions[rep][pl][indx].keys():
            worms=deletions[rep][pl][indx]['Worms']
            worms.append(coords)
            deletions[rep][pl][indx]['Worms']=worms
        else:
            deletions[rep][pl][indx]['Worms']=[coords]
        
    return deletions




def plot_comparison(figs, labels):
    nfig=len(figs)
    nlab=len(labels)
    fig, axes = plt.subplots(ncols=nfig, figsize=(4*nfig, 4), sharex=True,
                                   sharey=True)
    
    for fid,fplot in enumerate(figs):
        ax = axes[fid]
        ax.imshow(fplot) #, cmap=plt.cm.gray
        ax.set_title(labels[fid])
        ax.axis('off')
        ax.set_adjustable('box-forced')
    
#    ax2.imshow(filtered, cmap=plt.cm.gray)
#    ax2.set_title(filter_name)
#    ax2.axis('off')
#    ax2.set_adjustable('box-forced')
#    
#    ax3.imshow(extract) #cmap=plt.cm.gray
#    ax3.set_title('Extract')
#    ax3.axis('off')
#    ax3.set_adjustable('box-forced')
    
    return fig,axes
    

def freqtable(vector):
    my_series=pandas.Series(vector)
    counts = my_series.value_counts()
    ftable=[[key,value] for key,value in dict(counts).iteritems() ]
    return ftable


def writecsv(data,ofile,sep=','):
    f=open(ofile,"wb")
    ofile=csv.writer(f, delimiter=sep) # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
    for row in data:
    	#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
    	ofile.writerow(row)
    f.close()
    
def labeling(v,thr_1,thr_2,elevation='sobel'):
    
    markers = np.zeros_like(v)
    markers[v < thr_1] = 1
    markers[v > thr_2] = 2
    
    if elevation=='canny':
        elevation_map = canny(v)
        
    else:
        elevation_map = sobel(v)
    
    
    
    segmentation = watershed(elevation_map, markers)
    segmentation = ndi.binary_fill_holes(segmentation - 1)
    labeled_worms, _ = ndi.label(segmentation)
    
    return labeled_worms




os.chdir("/users/povilas/Dropbox/Projects/2015-Metformin/Biolog_Met_NGM/Celegans")
matplotlib.rcParams.update({'font.size': 12})

sourceloc="/Users/Povilas/Dropbox/Projects/Biolog_Celegans_Metformin"
replicate="Rep1"
plate="PM1"

odir='.'

nutrients=readmet('/Users/Povilas/Dropbox/Projects/2015-Metformin/Biolog/Biolog_metabolites_EcoCyc.csv')

alldeletions=readdel('{}/Deletions.csv'.format(odir))


settings={'Rep1': [[0.02,1,0.345,0.4], [7,7]],
          'Rep2': [[0.02,1,0.25,0.5], [7,7]],
          'Rep3': [[0.02,1,0.33,0.5], [7,7]]}

#Analyse multiple images
#thresholds=[0.05,1,0.355,0.4]
#sizes=[3,1,3,1]



#==============================================================================
# #@#Test parameters
#==============================================================================

#
#odir='2017-05-17_Trial_analysis2'
#
#sourceloc="2017-05-17_Trial"


replicate="Rep2"
plate="PM1"

plate='Control_1930'




thresholds,sizes=settings[replicate]

#hstep=0.01
#hmin=0.3
#hmax=0.4
#
#bstep=0.005
#bmin=0.01
#bmax=0.05



hstep=0.01
hmin=0.05
hmax=0.15

bstep=0.01
bmin=0.01
bmax=0.08



hues=np.arange(hmin,hmax+hstep,hstep)
brights=np.arange(bmin,bmax+bstep,bstep)

totalf=len(brights)*len(hues)


#files=["image_{}.tif".format(str(im).zfill(3)) for im in range(minimage,maximage+1)]

#imid=33


#files=["{}_{}.tif".format(plate,str(im).zfill(3)) for im in range(minimage,maximage+1)]

#filenum=[1,2,9,15,18,28,33,35,61]
filenum=[1,2,9,28,33]

files=["{}_{}.tif".format(plate,str(imid).zfill(3)) for imid in filenum]



data=[]

sizei=sizes

size_v_closing=[7]
size_h_closing=[7]


overall=len(files)*totalf*len(size_v_closing)*len(size_h_closing)
print overall

indo=0.0
for svc in size_v_closing:
    sizei[0]=svc
    for shc in size_h_closing:
        sizei[1]=shc
        sizestr='|'.join([str(sz) for sz in sizei])
        for file in files:
            fpat, fname, ftype = filename(file)
            print file
            location = "{}/{}/{}/{}".format(sourceloc,replicate,plate,file)
            
            well=indx2well(int(fname.split('_')[1]),start=1)
            fig, axes = plt.subplots(nrows=len(brights), ncols=len(hues),sharex=True, sharey=True, figsize=(40,24), dpi=300)
            fig.suptitle('{} - {} | {}'.format(fname,well,nutrients[plate][well]['Metabolite']))
            #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.5, wspace=0, hspace=0)
            ind=0.0
            
            image = tiff.imread(location)
            imghsv=color.rgb2hsv(image)
            #imgrgb=color.hsv2rgb(imghsv)
            
            v=imghsv[:,:,2]
            
            for bi, bval in enumerate(brights):
                for hi,hval in enumerate(hues):
                    
                    ind=ind+1.0
                    indo=indo+1.0
                    thresi=thresholds
                    thresi[0]=bval
                    thresi[2]=hval
                    thrstr='|'.join([str(th) for th in thresi])
                    
                    indx=int(fname.split('_')[1])
                    if indx in alldeletions[replicate][plate].keys():
                        print 'Deleting some worms!'
                        delworms=alldeletions[replicate][plate][indx]['Worms']
                    else:
                        delworms=[]
                        
                    #v_closing,h_closing,img2hsv,img2rgb=thresholding2(imghsv,thresi,sizei,delworms)
                    labeled_worms=labeling(v,bval,hval)
                    
                    plt.sca(axes[bi, hi]);ax = axes[bi, hi]
                    #plt.imshow(img2hsv)
                    plt.imshow(labeled_worms)
                    plt.setp(ax.get_yticklabels(), visible=False)
                    plt.setp(ax.get_xticklabels(), visible=False)
        
                    if hi == 0:
                        ax.yaxis.set_label_position("left")
                        plt.ylabel(bval, rotation='vertical')
                    if bi == 0:
                        plt.title(hval)
        
        
                    prcf=ind*100.0/totalf
                    prco=indo*100.0/overall
                    print '{:}:{:6.1f}% |{:6.1f}%'.format(fname,prcf,prco)
        
            #fig.tight_layout()
            fig.savefig('{}/{}_{}_{}_labeling.pdf'.format(odir,replicate,fname,sizestr), bbox_inches='tight')
            plt.close(fig)
            
            


#==============================================================================
# #@#96 well previews
#==============================================================================
#rows=['A','B','C','D','E','F','G','H']

ttllen=24
ttlsize=24

#Rep1
#thresholds=[0.02,1,0.345,0.4]
#sizes=[7,7]


minimage=1
maximage=96


replicate='Rep3'
plates=['PM1','PM2A','PM3B','PM4A']

plate=plates[0]

figures=['Original','Mask','Labeling']

figures=['Labeling']

thresi,sizei=settings[replicate]
thrstr='|'.join([str(th) for th in thresi])
sizestr='|'.join([str(sz) for sz in sizei])


thr_1,thr_2=[0.01,0.07]

lthrstr='{}|{}'.format(thr_1,thr_2)

data=[]

maxpix=0


#files=["{}_{}.tif".format('NGM_NoMetf',str(imid).zfill(3)) for imid in range(1,9)]+["{}_{}.tif".format('NGM_Metf',str(imid).zfill(3)) for imid in range(1,9)]
#plates=['Control_1930']

#plates=['PM1']
indo=0.0
for plate in plates:
    files=["{}_{}.tif".format(plate,str(im).zfill(3)) for im in range(1,96+1)]
    total=len(plates)*len(files)
    figall=[]
    axall=[]
    for fi,fign in enumerate(figures):
        fig, axes = plt.subplots(nrows=8, ncols=12, figsize=(72,47), dpi=300)
        fig.suptitle('{}-{}'.format(replicate,plate),fontsize=40)
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=2, wspace=0.5, hspace=0.1)
        figall.append(fig)
        axall.append(axes)
    row=0
    col=0
    
    for flin,fln in enumerate(files):
        fpat, fname, ftype = filename(fln)
        #print plate, fln, row, col
        location = "{}/{}/{}/{}".format(sourceloc,replicate,plate,fln)
        indx=flin+1
        if indx in alldeletions[replicate][plate].keys():
            print 'Deleting some worms!'
            delworms=alldeletions[replicate][plate][indx]['Worms']
        else:
            delworms=[]
            
        image = tiff.imread(location)
        imghsv = color.rgb2hsv(image)
        
        well=indx2well(flin)
        ttl='{}-{} | {}'.format(flin+1,well,nutrients[plate][well]['Metabolite'])
        ttlw='\n'.join(tw.wrap(ttl,ttllen))
        
        for fi,fign in enumerate(figures):
            fig=figall[fi]
            axes=axall[fi]
            #plt.figure(fi)
            #print row,col
            try:
                plt.sca(axes[row,col])
            except Exception as e:
                print e
                sys.exit(1)
            try:
                ax = axes[row,col]
            except Exception as e:
                print e
                sys.exit(1)
            
            if fign=='Original':
                imgrgb=color.hsv2rgb(imghsv)
                plt.imshow(imgrgb)
            elif fign=='Mask':
                v_closing,h_closing,img2hsv,img2rgb=thresholding2(imghsv,thresi,sizei,delworms)
                plt.imshow(img2hsv)
            elif fign=='Labeling':
                v=imghsv[:,:,2]
                labeled_worms=labeling(v,thr_1,thr_2)
                plt.imshow(labeled_worms)
            else:
                plt.imshow(imgrgb)
            
            plt.title(ttlw,fontsize=ttlsize)
            plt.setp(ax.get_yticklabels(), visible=False)
            plt.setp(ax.get_xticklabels(), visible=False)
            
        if fign=='Mask':
            v2 =img2hsv[:,:,2]
            bright1D=[el for el in list(v2.ravel()) if el>0]
            ftable=freqtable(bright1D)
            rowhead=[replicate,plate,well,fln]#+thresi+sizei
            for fv in ftable:
                data.append(rowhead+fv)
        
        col=col+1
        if col==12:
            col=0
            row=row+1
        
        indo=indo+1
        prc=(flin+1)*100.0/len(files)        
        prco=(indo)*100.0/total

        print '{:} {:} {:}:{:6.1f}% | {:6.1f}%'.format(replicate,plate,fname,prc,prco)
        
    for fi,fign in enumerate(figures):
        fig=figall[fi]
        axes=axall[fi]
        #plt.figure(fi)
        fig.tight_layout()
        
        if fign=='Original':
            ofignm='{}/{}_{}_{}.pdf'.format(odir,replicate,plate,fign)
        elif fign=='Mask':
            ofignm='{}/{}_{}_{}_{}_{}.pdf'.format(odir,replicate,plate,thrstr,sizestr,fign)
        elif fign=='Labeling':
            ofignm='{}/{}_{}_{}_{}.pdf'.format(odir,replicate,plate,lthrstr,fign)
        else:
            ofignm='{}/{}_{}_{}.pdf'.format(odir,replicate,plate,fign)
        fig.savefig(ofignm, bbox_inches='tight')
        plt.close(fig)



header=['Replicate','Plate','Well','File','Brightness_value','Frequency']
data.insert(0,header)
f=open('{}/Summary_{}_{}_{}.csv'.format(odir,replicate,thrstr,sizestr),"wb")
ofile=csv.writer(f, delimiter='\t') # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
for row in data:
	#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
	ofile.writerow(row)
f.close()



##########
#Just get the data
##########

replicate='Rep2'
thresi,sizei=settings[replicate]


files=["{}_{}.tif".format('NGM_NoMetf',str(imid).zfill(3)) for imid in range(1,9)]+["{}_{}.tif".format('NGM_Metf',str(imid).zfill(3)) for imid in range(1,9)]
folders=['Controls_1930','Controls_2240']

data=[]

#plates=['PM1']
indo=0.0

fig, axes = plt.subplots(nrows=4, ncols=8, figsize=(72,47), dpi=300)
fig.suptitle('Negative controls',fontsize=40)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=2, wspace=0.5, hspace=0.1)


row=0
col=0

for folder in folders:
    #files=["{}_{}.tif".format(plate,str(im).zfill(3)) for im in range(1,96+1)]
    total=len(folders)*len(files)

    for flin,fln in enumerate(files):
        fpat, fname, ftype = filename(fln)
        location = "{}/{}/{}/{}".format(sourceloc,replicate,folder,fln)
        if os.path.exists(location):
            image = tiff.imread(location)
            imghsv = color.rgb2hsv(image)
            imgrgb = color.hsv2rgb(imghsv)
            v_closing,h_closing,img2hsv,img2rgb=thresholding2(imghsv,thresi,sizei,[])

    
            v2 =img2hsv[:,:,2]
            bright1D=[el for el in list(v2.ravel()) if el>0]
            ftable=freqtable(bright1D)
            rowhead=[replicate,folder,fln]#+thresi+sizei
            for fv in ftable:
                data.append(rowhead+fv)
                
        try:
            plt.sca(axes[row,col])
        except Exception as e:
            print e
            sys.exit(1)
        try:
            ax = axes[row,col]
        except Exception as e:
            print e
            sys.exit(1)
        
        plt.imshow(imgrgb)
        if col == 0:
            ax.yaxis.set_label_position("left")
            plt.ylabel('{}_{}'.format(folder,'-'.join(fname.split('_')[0:1])), rotation='vertical')
        
            
        col=col+1
        if col==8:
            col=0
            row=row+1
            
        indo=indo+1
        prc=(flin+1)*100.0/len(files)        
        prco=(indo)*100.0/total
        print '{:} {:} {:}:{:6.1f}% | {:6.1f}%'.format(replicate,plate,fname,prc,prco)


fig.tight_layout() 
ofignm='{}/Negative_Controls_comparison.pdf'.format(odir)
fig.savefig(ofignm, bbox_inches='tight')
plt.close(fig)       
        
        


header=['Replicate','Folder','File','Brightness_value','Frequency']
data.insert(0,header)
f=open('{}/Summary_{}_{}_{}_Negative_controls.csv'.format(odir,replicate,thrstr,sizestr),"wb")
ofile=csv.writer(f, delimiter='\t') # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
for row in data:
	#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
	ofile.writerow(row)
f.close()




#Check controls


replicate='Rep2'
folder='Controls_1930'
indx=2
#rfile="{}/{}/NGM_Control_Rep1.tif".format(sourceloc,replicate)
tfile="{}/{}/{}/{}_{}.tif".format(sourceloc,replicate,folder,'NGM_NoMetf',str(indx).zfill(3))


thresi,sizei=settings[replicate]

    
image = tiff.imread(tfile)
imghsv=color.rgb2hsv(image)
imgrgb=color.hsv2rgb(imghsv)

v_closing,h_closing,img2hsv,img2rgb=thresholding2(imghsv,thresi,sizei,[])


#plt.imshow(imgrgb)


plt.imshow(img2hsv)






#alldeletions=readdel('{}/Deletions.csv'.format(odir))
from skimage.filters import threshold_otsu, threshold_adaptive
from skimage.filters import rank
from skimage.feature import peak_local_max

replicate='Rep3'
plate='PM1'
indx=3
#rfile="{}/{}/NGM_Control_Rep1.tif".format(sourceloc,replicate)
tfile="{}/{}/{}/{}_{}.tif".format(sourceloc,replicate,plate,plate,str(indx).zfill(3))


thresi=[0.02,1,0.33,0.4]
sizei=[7,7]
alpha=2


if indx in alldeletions[replicate][plate].keys():
    delworms=alldeletions[replicate][plate][33]['Worms']
else:
    delworms=[]
    
image = tiff.imread(tfile)
imghsv=color.rgb2hsv(image)
imgrgb=color.hsv2rgb(imghsv)

#v_closing,h_closing,img2hsv,img2rgb=thresholding2(imghsv,thresi,sizei,delworms)

#s_image,s_imgrgb,s_v_closing,s_h_closing,s_img2hsv,s_img2rgb=thresholding3(tfile,thresi,sizei,alpha,delworms)

#Watershed labelling


plt.imshow(imgrgb)


plt.imshow(imghsv)


h=imghsv[:,:,0]
s=imghsv[:,:,1]
v=imghsv[:,:,2]


dh = rank.median(h, disk(2))
dv = rank.median(v, disk(2))

#v=imghsv[:,:,0]

plt.imshow(h)




thr_1,thr_2=[0.27,0.33]

thr_1,thr_2=[0.01,0.07]


markers = np.zeros_like(v)

markers[v < thr_1] = 1
markers[v > thr_2] = 2
plt.imshow(markers)

#markers[np.logical_not(binary_adaptive)] = 1
#markers[binary_adaptive] = 2

#elevation_map = canny(h)
elevation_map = sobel(v)
plt.imshow(elevation_map)

segmentation = watershed(elevation_map, markers)
segmentation = ndi.binary_fill_holes(segmentation - 1)
plt.imshow(segmentation)
labeled_worms, _ = ndi.label(segmentation)

plt.imshow(labeled_worms)



imgsel=labeled_worms[400:700,700:1200]
plt.imshow(imgsel)


distance = ndi.distance_transform_edt(imgsel)
plt.imshow(distance)



local_maxi = peak_local_max(distance, indices=False,
                            footprint=np.ones((3, 3)),
                            labels=imgsel)
plt.imshow(local_maxi)
markers = ndi.label(local_maxi)[0]
labels = watershed(-distance, markers, mask=imgsel)

plt.imshow(labels)


#labeled_worms=labeling(v,0.02,0.08)

x, y = np.indices((80, 80))
x1, y1, x2, y2 = 28, 28, 44, 52
r1, r2 = 16, 20
mask_circle1 = (x - x1)**2 + (y - y1)**2 < r1**2
mask_circle2 = (x - x2)**2 + (y - y2)**2 < r2**2
image = np.logical_or(mask_circle1, mask_circle2)

    

#Canny filter does not work
edges = canny(v)
plt.imshow(edges)

fill_worms = ndi.binary_fill_holes(edges)
plt.imshow(fill_worms)


label_objects, nb_labels = ndi.label(fill_worms)
sizes = np.bincount(label_objects.ravel())
mask_sizes = sizes > 20
mask_sizes[0] = 0
worms_cleaned = mask_sizes[label_objects]

plt.imshow(worms_cleaned)


#Other testing
v=img2hsv[:,:,2]
#s_v=s_img2hsv[:,:,2]

figs=[imgrgb,v_closing,h_closing,img2hsv,v]
labels=['Original','V_Mask','H_Mask','Mask','Extract']
fig=plot_comparison(figs,labels)


#Adaptive thresholding

#block_size = 25
#binary_adaptive = threshold_adaptive(h, block_size,offset=0.01)#, offset=10
#plt.imshow(binary_adaptive)
#imgrgb_t=imgrgb.copy()
#
#
#imgrgb_t[binary_adaptive,:]=0;
#plt.imshow(imgrgb_t)


#fig = plt.figure()
#ax = fig.add_subplot(111)
#plt.imshow(s_v)

#coords = []

def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d'%(ix, iy)

    global coords
    coords.append([ix, iy])

    if len(coords) == 2:
        fig.canvas.mpl_disconnect(cid)

    return coords

cid = fig.canvas.mpl_connect('button_press_event', onclick)




#Rectangle
from matplotlib.patches import Rectangle

class Annotate(object):
    def __init__(self):
        self.ax = plt.gca()
        self.rect = Rectangle((0,0), 1, 1)
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)

    def on_press(self, event):
        print 'press'
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_release(self, event):
        print 'release'
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.ax.figure.canvas.draw()

a = Annotate()
plt.show()




#refhsv=color.rgb2hsv(ref)

#image.shape
#image.dtype

#plt.imshow(image)
#image = tiff.imread(file)
#lena = misc.lena()
#blurred_l = ndimage.gaussian_filter(lena, 3)
##increase the weight of edges by adding an approximation of the Laplacian:
    
blur =imghsv[:,:,2]

filter_blurred_l = ndimage.gaussian_filter(blur, 1)
alpha = 30
sharpened = blur + alpha * (blur - filter_blurred_l)
plt.imshow(sharpened)



vc_size,hc_size=[7,7]

vmin,vmax,hmin,hmax=[0.05,1,0.35,0.4]
#image.shape
#image.dtype
imghsv=color.rgb2hsv(image)

imgrgb=color.hsv2rgb(imghsv)

h =imghsv[:,:,0]
#s =imghsv[:,:,1]
v =imghsv[:,:,2]


vmask = (v > vmin).astype(np.float64);  # Thresholding in the Brightness
vmask = (vmask< vmax).astype(np.float64);

plt.imshow(vmask)


v_closing = closing(vmask, selem=disk(vc_size))
plt.imshow(v_closing)


image=h_closing

filter_blurred_l = ndimage.gaussian_filter(image, 1)
alpha = 30
sharpened = image + alpha * (image - filter_blurred_l)
plt.imshow(sharpened)





#testing

thresholds=[0.02,1,0.355,0.4]
sizes=[7,7]

minimage=1
maximage=96
maxpix=0
data=[]

plot=True
saveimg=True


files=["PM1_{}.tif".format(str(im).zfill(3)) for im in range(minimage,maximage+1)]




rowslices=slice_list(files,8) #Slice list to 8 row slices


for file in files:
	fpat, fname, ftype = filename(file)
	print file
	location = "{}/{}/{}/{}".format(sourceloc,replicate,plate,file)

	image,v_opening,h_opening,img2hsv,img2rgb=thresholding(location,thresholds,sizes)

	histinfo_o = exposure.histogram(image)
	histinfo_t = exposure.histogram(img2rgb)
	histinfo_t[0][0] = 0

	h =image[:,:,0]
	s =image[:,:,1]
	v =image[:,:,2]

	#h2 =img2[:,:,0]
	#s2 =img2[:,:,1]
	v2 =img2hsv[:,:,2]

	if plot:
		fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(8.27,11.69), dpi=100)
		fig.suptitle(fname)
		plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.2, hspace=0.2)

		plt.sca(axes[0,0]);ax = axes[0,0]
		plt.imshow(image);plt.title('Original')
		plt.sca(axes[0,1]);ax = axes[0,1]
		plt.imshow(v*10);plt.title('Brightness original')
		plt.sca(axes[1,0]);ax = axes[1,0]
		plt.imshow(img2hsv);plt.title('Mask based on hue and brightness')
		plt.sca(axes[1, 1]);ax = axes[1, 1]
		plt.imshow(v2);plt.title('Extract')

		plt.sca(axes[2,0]);ax=axes[2,0]
		plt.plot(histinfo_o[1]/255.0,histinfo_o[0]);plt.title("Brightness histogram - Original");plt.xlabel("Brightness");plt.ylabel("Pixels");plt.xlim([0,0.1])
		plt.sca(axes[2,1]);ax=axes[2,1]
		plt.plot(histinfo_t[1],histinfo_t[0]);plt.title("Brightness histogram - Thresholding");plt.xlabel("Brightness");plt.ylabel("Pixels");plt.xlim([0,1])


		fig.tight_layout()
		if saveimg:
			fig.savefig('{}/{}_{}_analysis.pdf'.format(odir,replicate,fname), bbox_inches='tight')
		plt.close(fig)

	bright1D=[el for el in list(v2.ravel()) if el>0]
	if len(bright1D)>maxpix:
		maxpix=len(bright1D)
	row=[fname]+thresholds+sizes+bright1D
	data.append(row)


maxpix


header=['File','Brightness_min','Brightness_max','Hue_min','Hue_max','Bright_opening','Bright_dilation','Hue_opening','Hue_dilation']+range(1,maxpix+1)
data.insert(0,header)



ofile='{}/Summary.csv'.format(odir)

writecsv(data,ofile)



#==============================================================================
# Mark sick worms
#==============================================================================

#def onclick(event):
#    global ix, iy
#    global coords
#    #global worms
#    #coords=[]
#
#    ix, iy = event.xdata, event.ydata
#    print 'x = %d, y = %d'%(ix, iy)
#    coords.append([ix, iy])
#
#    #if len(coords) == 2:
#    fig.canvas.mpl_disconnect(cid)
#                    
#    return coords



deletions={'Rep1-PM1':[3,29,33,41,44,57,62,93],
           'Rep1-PM2A':[30,40,41,45,47,48,75],
           'Rep1-PM3B':[10,31,33,41,90],
           'Rep1-PM4A':[32,58,82]}

dfiles=[ "{}/{}/{}/{}_{}.tif".format(sourceloc,replicate,plate,plate,str(im).zfill(3)) for im in  deletions.iteritems()]

thresi=[0.02,1,0.345,0.4]
sizei=[7,7]



#from matplotlib.patches import Rectangle

deldata=[]
plt.ion() ## Note this correction
for dkey in deletions.keys():
    rep,plate=dkey.split('-')
    dfiles=deletions[dkey]
    
    thresi,sizei=settings[rep]
    for indx in dfiles:
        dfile="{}/{}/{}/{}_{}.tif".format(sourceloc,rep,plate,plate,str(indx).zfill(3))
        dname=os.path.basename(dfile)
        well=indx2well(indx)
        
        
        print rep,plate,well,dname
        
        image = tiff.imread(location)
        imghsv=color.rgb2hsv(image)
        imgrgb=color.hsv2rgb(imghsv)
        
        
        
        v_closing,h_closing,img2hsv,img2rgb=thresholding2(imghsv,thresi,sizei,delworms)
        
        #v=img2hsv[:,:,2]
        #plt.imshow(s_img2hsv[822:964,749:1014,:])
        worms=[]
        coords=[]
        fig,axes=plot_comparison(imgrgb,img2hsv,'Pick worms: {}'.format(dname))

        while True:
            #while len(coords)<2:
            coords=plt.ginput(2)
            print coords
            
            x1,y1=coords[0]
            x2,y2=coords[1]
            
            xs=[x1,x2]
            ys=[y1,y2]
            
            excordt=[x1,y1,x2,y2]
            excord=[int(c) for c in excordt]
            
            #ax2=axes[1]
            #ax2.add_patch(Rectangle(( max(xs)- min(xs), max(ys) -min(ys)), 0.2, 0.2,
            #          alpha=0.5, facecolor='none'))
            #fig.canvas.draw()
            #cid = fig.canvas.mpl_connect('button_press_event', onclick)
            #plt.waitforbuttonpress()
            #print coords,cid
            
            #plt.pause(0.05)
            save=raw_input("Save coordinates (y/n)?: ")
            if save=='y':
                worms.append(excord)
                coords=[]
                multi=raw_input("Mark other worms (y/n)?: ")
                if multi=='y':
                    continue
                else:
                    break
            else:
                coords=[]
                multi=raw_input("Mark other worm (y/n)?: ")
                if multi=='y':
                    continue
                else:
                    break
            
        headr=[rep,plate,indx,dname]
        for worm in worms:
            deldata.append(headr+worm)
        plt.close(fig)


header=['Replicate','Plate','Index','File','X1','Y1','X2','Y2']
deldata.insert(0,header)

#ofile='{}/Deletions.csv'.format(odir)
#writecsv(deldata,ofile)









