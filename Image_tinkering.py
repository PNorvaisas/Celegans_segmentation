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


from skimage.morphology import disk, opening, dilation, square
from skimage import exposure


import tifffile as tiff

import textwrap as tw

import csv

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







os.chdir("/users/povilas/Dropbox/Projects/2015-Metformin/Worm_imaging")
matplotlib.rcParams.update({'font.size': 12})

sourceloc="/Users/Povilas/Dropbox/Projects/Biolog_Celegans_Metformin"
replicate="Rep1"
plate="PM1"

odir='Metformin_Biolog_Celegans_Screen_results'

nutrients=readmet('/Users/Povilas/Dropbox/Projects/2015-Metformin/Biolog/Biolog_metabolites_EcoCyc.csv')
#Analyse multiple images
#thresholds=[0.05,1,0.355,0.4]
#sizes=[3,1,3,1]


thresholds=[0.02,1,0.355,0.4]
sizes=[3,1,3,1]




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

	image,vdilated,hdilated,img2hsv,img2rgb=thresholding(location,thresholds,sizes)

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



f=open('{}/Summary.csv'.format(odir),"wb")
ofile=csv.writer(f, delimiter='\t') # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
for row in data:
	#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
	ofile.writerow(row)
f.close()







#==============================================================================
# #@#Test parameters
#==============================================================================

#
#odir='2017-05-17_Trial_analysis2'
#
#sourceloc="2017-05-17_Trial"

thresholds=[0.05,1,0.355,0.4]
sizes=[3,1,3,1]

hstep=0.005
hmin=0.34
hmax=0.375

bstep=0.0025
bmin=0.01
bmax=0.03

hues=np.arange(hmin,hmax+hstep,hstep)
brights=np.arange(bmin,bmax+bstep,bstep)

totalf=len(brights)*len(hues)


#files=["image_{}.tif".format(str(im).zfill(3)) for im in range(minimage,maximage+1)]

replicate="Rep1"
plate="PM1"
#imid=33


#files=["{}_{}.tif".format(plate,str(im).zfill(3)) for im in range(minimage,maximage+1)]

files=["{}_{}.tif".format(plate,str(imid).zfill(3)) for imid in [1,2,9,15,18,28,33,35,61]]

data=[]

sizei=sizes

size_v_opening=[3]
size_v_dilation=[1]
size_h_opening=[3]
size_h_dilation=[1]


overall=len(files)*totalf*len(size_v_opening)*len(size_v_dilation)*len(size_h_opening)*len(size_h_dilation)
print overall

indo=0.0
for svo in size_v_opening:
    sizei[0]=svo
    for svd in size_v_dilation:
        sizei[1]=svd
        for sho in size_h_opening:
            sizei[2]=sho
            for shd in size_h_dilation:
                sizei[3]=shd
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
                    for bi, bval in enumerate(brights):
                        for hi,hval in enumerate(hues):
                            
                            ind=ind+1.0
                            indo=indo+1.0
                            thresi=thresholds
                            thresi[0]=bval
                            thresi[2]=hval
                            thrstr='|'.join([str(th) for th in thresi])
                            
                            image,imgrgb,vdilated,hdilated,img2hsv,img2rgb=thresholding(location,thresi,sizes)
                            plt.sca(axes[bi, hi]);ax = axes[bi, hi]
                            plt.imshow(img2hsv)
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
                    fig.savefig('{}/{}_{}_{}_parameters.pdf'.format(odir,replicate,fname,sizestr), bbox_inches='tight')
                    plt.close(fig)
                

#==============================================================================
# End Block
#==============================================================================



#==============================================================================
# #@#96 well previews
#==============================================================================



#rows=['A','B','C','D','E','F','G','H']

ttllen=24
ttlsize=24

thresholds=[0.03,1,0.365,0.4]
sizes=[3,1,3,1]


minimage=1
maximage=96


thresi=thresholds
sizei=sizes

thrstr='|'.join([str(th) for th in thresi])
sizestr='|'.join([str(sz) for sz in sizei])


#
#odir='2017-05-17_Trial_analysis2'
#sourceloc="2017-05-17_Trial"

plates=['PM1','PM2A','PM3B','PM4A']

figures=['Original','Mask']




#plates=['PM1']
for plate in plates:
    files=["{}_{}.tif".format(plate,str(im).zfill(3)) for im in range(1,96+1)]
    
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
        
        image,imgrgb,vdilated,hdilated,img2hsv,img2rgb=thresholding(location,thresi,sizei)
        #image = tiff.imread(location)
        
        #imghsv=color.rgb2hsv(image)
        #imgrgb=color.hsv2rgb(imghsv)
        #h =imghsv[:,:,0]
        #s =imghsv[:,:,1]
        #v =imghsv[:,:,2]
        well=indx2well(flin)
        ttl='{}-{} | {}'.format(flin+1,well,nutrients[plate][well]['Metabolite'])
        ttlw='\n'.join(tw.wrap(ttl,ttllen))
        
        for fi,fign in enumerate(figures):
            fig=figall[fi]
            axes=axall[fi]
            #plt.figure(fi)
            plt.sca(axes[row,col])
            ax = axes[row,col]
            
            if fign=='Original':
                plt.imshow(imgrgb)
            elif fign=='Mask':
                plt.imshow(img2hsv)
            else:
                plt.imshow(image)
            
            plt.title(ttlw,fontsize=ttlsize)
            plt.setp(ax.get_yticklabels(), visible=False)
            plt.setp(ax.get_xticklabels(), visible=False)
        
        col=col+1
        if col==12:
            col=0
            row=row+1
        prc=(flin+1)*100.0/len(files)
        print '{:} {:} {:}:{:6.1f}%'.format(replicate,plate,fname,prc)
        
    for fi,fign in enumerate(figures):
        fig=figall[fi]
        axes=axall[fi]
        #plt.figure(fi)
        fig.tight_layout()
        
        if fign=='Original':
            ofignm='{}/{}_{}_{}.pdf'.format(odir,replicate,plate,fign)
        elif fign=='Mask':
            ofignm='{}/{}_{}_{}_{}_{}.pdf'.format(odir,replicate,plate,thrstr,sizestr,fign)
        else:
            ofignm='{}/{}_{}_{}.pdf'.format(odir,replicate,plate,fign)
        fig.savefig(ofignm, bbox_inches='tight')
        plt.close(fig)
    




rfile="{}/{}/NGM_Control_Rep1.tif".format(sourceloc,replicate)

tfile="{}/{}/{}/PM1_001.tif".format(sourceloc,replicate,'PM1')

ref=tiff.imread(rfile)
image = tiff.imread(tfile)


refhsv=color.rgb2hsv(ref)

image.shape
image.dtype

plt.imshow(image)


#import opencv as cv2















#Get brightness channel


plt.imshow(imgrgb)

bright=imghsv[:,:,2]
bright1D=bright.ravel()
bright1D=bright1D*255.0
len(bright1D)


bright1D_filt=[el for el in list(bright1D) if el>0]

len(bright1D_filt)


np.mean(bright1D_filt)

histinfo_t = exposure.histogram(imgrgb)
histinfo_t[0][0] = 0


np.sum(histinfo_t[1]*histinfo_t[0])/np.sum(histinfo_t[0])




#image 6
0.064363794254385251
#image 5
0.044202210383164464




plt.plot(histinfo[1],histinfo[0])
plt.title("Gaussian Histogram");plt.xlabel("Value");plt.ylabel("Frequency")



def plot_img_and_hist(image, axes, bins=256):
    """Plot an image along with its histogram and cumulative histogram.

    """
    ax_img, ax_hist = axes
    ax_cdf = ax_hist.twinx()

    # Display image
    ax_img.imshow(image, cmap=plt.cm.gray)
    ax_img.set_axis_off()

    # Display histogram
    ax_hist.hist(image.ravel(), bins=bins)
    ax_hist.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    ax_hist.set_xlabel('Pixel intensity')

    xmin, xmax = dtype_range[image.dtype.type]
    ax_hist.set_xlim(xmin, xmax)

    # Display cumulative distribution
    img_cdf, bins = exposure.cumulative_distribution(image, bins)
    ax_cdf.plot(bins, img_cdf, 'r')

    return ax_img, ax_hist, ax_cdf


# Load an example image
img = img_as_ubyte(color.rgb2gray(imgthres))

#img = img_as_ubyte(color.rgb2gray(img2rgb))

# Global equalize
img_rescale = exposure.equalize_hist(img)

# Equalization
selem = disk(30)
img_eq = rank.equalize(img, selem=selem)


# Display results
fig = plt.figure(figsize=(8, 5))
axes = np.zeros((2, 3), dtype=np.object)
axes[0, 0] = plt.subplot(2, 3, 1, adjustable='box-forced')
axes[0, 1] = plt.subplot(2, 3, 2, sharex=axes[0, 0], sharey=axes[0, 0],
                         adjustable='box-forced')
axes[0, 2] = plt.subplot(2, 3, 3, sharex=axes[0, 0], sharey=axes[0, 0],
                         adjustable='box-forced')
axes[1, 0] = plt.subplot(2, 3, 4)
axes[1, 1] = plt.subplot(2, 3, 5)
axes[1, 2] = plt.subplot(2, 3, 6)

ax_img, ax_hist, ax_cdf = plot_img_and_hist(img, axes[:, 0])
ax_img.set_title('Low contrast image')
ax_hist.set_ylabel('Number of pixels')

ax_img, ax_hist, ax_cdf = plot_img_and_hist(img_rescale, axes[:, 1])
ax_img.set_title('Global equalise')

ax_img, ax_hist, ax_cdf = plot_img_and_hist(img_eq, axes[:, 2])
ax_img.set_title('Local equalize')
ax_cdf.set_ylabel('Fraction of total intensity')


# prevent overlap of y-axis labels
fig.tight_layout()
plt.show()













#Colorize

def colorize(image, hue, saturation=1):
    """ Add color of the given hue to an RGB image.

    By default, set the saturation to 1 so that the colors pop!
    """
    hsv = color.rgb2hsv(image)
    hsv[:, :, 1] = saturation
    hsv[:, :, 0] = hue
    return color.hsv2rgb(hsv)

hue_rotations = np.linspace(0, 1, 6)

fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True)

for ax, hue in zip(axes.flat, hue_rotations):
    # Turn down the saturation to give it that vintage look.
    tinted_image = colorize(image, hue, saturation=0.3)
    ax.imshow(tinted_image, vmin=0, vmax=1)
    ax.set_axis_off()
    ax.set_adjustable('box-forced')
fig.tight_layout()


#Rank
from skimage.filters import rank


image=imggray.copy()
# Square regions defined as slices over the first two dimensions.
top_left = (slice(100),) * 2
bottom_right = (slice(-100, None),) * 2

sliced_image = image.copy()
sliced_image[top_left] = colorize(image[top_left], 0.82, saturation=0.5)
sliced_image[bottom_right] = colorize(image[bottom_right], 0.5, saturation=0.5)

# Create a mask selecting regions with interesting texture.
noisy = rank.entropy(grayscale_image, np.ones((9, 9)))
textured_regions = noisy > 4
# Note that using `colorize` here is a bit more difficult, since `rgb2hsv`
# expects an RGB image (height x width x channel), but fancy-indexing returns
# a set of RGB pixels (# pixels x channel).
masked_image = image.copy()
masked_image[textured_regions, :] *= red_multiplier

fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(8, 4), sharex=True, sharey=True)
ax1.imshow(sliced_image)
ax2.imshow(masked_image)
ax1.set_adjustable('box-forced')
ax2.set_adjustable('box-forced')

plt.show()



############



file = '20170503/image_011.tif'
img = skio.imread(file)

img.shape
img.dtype

plt.imshow(img)

imghsv=color.rgb2hsv(img)
imggray=color.rgb2gray(img)

plt.imshow(imghsv)

imgfromgray=color.gray2rgb(imggray)

sobel = filters.sobel(imgfromgray)


plt.imshow(sobel)

plt.imshow(imgfromgray)

grayscale_image = img_as_float(img[::2, ::2])


image = color.gray2rgb(grayscale_image)


plt.imshow(image*[10,10,10])




#Hue gradients
hue_gradient = np.linspace(0, 1)
hsv = np.ones(shape=(1, len(hue_gradient), 3), dtype=float)
hsv[:, :, 0] = hue_gradient

all_hues = color.hsv2rgb(hsv)

fig, ax = plt.subplots(figsize=(5, 2))
# Set image extent so hues go from 0 to 1 and the image is a nice aspect ratio.
ax.imshow(all_hues, extent=(0, 1, 0, 0.2))
ax.set_axis_off()





#Split components
h =imghsv[:,:,0]
s =imghsv[:,:,1]
v =imghsv[:,:,2]


plt.figure(1, figsize=(15, 15))

plt.subplot(4,2,1); plt.imshow(h, cmap='gray'); plt.title('Hue')
plt.subplot(4,2,2); plt.imshow(s, cmap='gray'); plt.title('Saturation')
plt.subplot(4,2,3); plt.imshow(v, cmap='gray'); plt.title('Value')
plt.tight_layout()



#Set hue range




#Split components





