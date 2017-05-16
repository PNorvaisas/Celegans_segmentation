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
	image = skio.imread(file)

	vo_size,vd_size,ho_size,hd_size=sizes
	vmin,vmax,hmin,hmax=thresholds
	#image.shape
	#image.dtype
	imghsv=color.rgb2hsv(image)

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


	return image,vdilated,hdilated,img2,img2rgb




os.chdir("/users/povilas/Dropbox/Projects/2015-Metformin/Worm_imaging/Trial images")
matplotlib.rcParams.update({'font.size': 12})



#Analyse multiple images


thresholds=[0.05,1,0.355,0.4]
sizes=[3,1,3,1]


odir='20170503_analysis'

minimage=9
maximage=11
maxpix=0
data=[]

plot=True
saveimg=True


files=["image_{}.tif".format(str(im).zfill(3)) for im in range(minimage,maximage+1)]

for file in files:
	fpat, fname, ftype = filename(file)
	print file
	location = "20170503/{}".format(file)

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
			fig.savefig('{}/{}_analysis.pdf'.format(odir,fname), bbox_inches='tight')
		plt.close(fig)

	bright1D=[el for el in list(v2.ravel()) if el>0]
	if len(bright1D)>maxpix:
		maxpix=len(bright1D)
	row=[fname]+thresholds+sizes+bright1D
	data.append(row)


maxpix


header=['File','Brightness_min','Brightness_max','Hue_min','Hue_max','Bright_opening','Bright_dilation','Hue_opening','Hue_dilation']+range(1,maxpix+1)
data.insert(0,header)



f=open('20170503_analysis/Summary.csv',"wb")
ofile=csv.writer(f, delimiter='\t') # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
for row in data:
	#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
	ofile.writerow(row)
f.close()






#Test parameters
odir='20170503_parameters'

thresholds=[0.05,1,0.355,0.4]
sizes=[3,1,3,1]


hstep=0.01
hmin=0.3
hmax=0.4

bstep=0.01
bmin=0.04
bmax=0.07


hues=np.arange(hmin,hmax+hstep,hstep)
brights=np.arange(bmin,bmax+bstep,bstep)




minimage=9
maximage=9
maxpix=0


total=len(brights)*len(hues)
files=["image_{}.tif".format(str(im).zfill(3)) for im in range(minimage,maximage+1)]

data=[]
for file in files:
	fpat, fname, ftype = filename(file)
	print file
	location = "20170503/{}".format(file)

	fig, axes = plt.subplots(nrows=len(brights), ncols=len(hues),sharex=True, sharey=True, figsize=(24,12), dpi=100)
	fig.suptitle(fname)

	#plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.5, wspace=0, hspace=0)

	ind=0.0
	for bi, bval in enumerate(brights):
		for hi,hval in enumerate(hues):
			ind=ind+1.0
			thresi=thresholds
			thresi[0]=bval
			thresi[2]=hval
			image,vdilated,hdilated,img2hsv,img2rgb=thresholding(location,thresi,sizes)

			plt.sca(axes[bi, hi]);ax = axes[bi, hi]
			plt.imshow(img2hsv)
			plt.setp(ax.get_yticklabels(), visible=False)
			plt.setp(ax.get_xticklabels(), visible=False)

			if hi == 0:
				ax.yaxis.set_label_position("left")
				plt.ylabel(bval, rotation='horizontal')
			if bi == 0:
				plt.title(hval)


			prc=ind*100.0/total
			print '{:}:{:6.1f}%'.format(fname,prc)

	#fig.tight_layout()
	fig.savefig('{}/{}_parameters.pdf'.format(odir,fname), bbox_inches='tight')
	plt.close(fig)










header=['File','Brightness_min','Brightness_max','Hue_min','Hue_max','Bright_opening','Bright_dilation','Hue_opening','Hue_dilation']+range(1,maxpix+1)
data.insert(0,header)



f=open('20170503_analysis/Summary.csv',"wb")
ofile=csv.writer(f, delimiter='\t') # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
for row in data:
	#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
	ofile.writerow(row)
f.close()










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





