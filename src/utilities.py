import os
import sys
import tifffile as tiff
import textwrap as tw
import csv
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


#plt.interactive(False)


import numpy as np
import itertools as IT
import time


import skimage
from skimage import io as skio
from skimage import color
from skimage import data
from skimage import img_as_float

from skimage.util.dtype import dtype_range
from skimage.util import img_as_ubyte

from skimage.feature import canny
from skimage.feature import peak_local_max

from skimage import filters
from skimage.filters import sobel

from skimage.morphology import disk, opening, dilation, square, watershed

from skimage.morphology import erosion, white_tophat, black_tophat, closing

from skimage import exposure




from skimage.filters import rank
from skimage.filters.rank import median, mean
from skimage import measure
from skimage.measure import label, regionprops



import scipy.optimize as opt
from scipy import ndimage as ndi



def filename(ifile):
	if ifile.split('.')[0] == '':
		ipat = ''
		iname = ''
		itype = ifile.split('.')[1]
	else:
		if "\\" in ifile.split('.')[0]:
			sep = "\\"
		elif "/" in ifile.split('.')[0]:
			sep = "/"
		else:
			ipat = ''
			iname = ifile.split('.')[0]
			itype = ifile.split('.')[1]
			return ipat, iname, itype
		allpath = ifile.split('.')[0]
		iname = allpath.split(sep)[-1]
		ipath = allpath.split(sep)[:-1]
		ipat = '/'.join(ipath)
		itype = ifile.split('.')[1]
	return ipat, iname, itype


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
		if s == 'NAN':
			return s
		float(s)
		if float(s).is_integer():
			return int(float(s))
		elif float(s) == 0:
			return float(s)
		else:
			return float(s)

	except ValueError:
		return s


def indx2well(ind, start=0, rowln=12):
	indt = ind - start
	rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	row = indt // rowln
	col = indt - row * rowln + 1
	well = '{}{}'.format(rows[row], col)
	return well

def readTF(ifile):
    print ifile
    TFs = NestedDict()
    allfiles=pd.read_csv(ifile)
    headers = allfiles.columns
    print headers
    nec = ['Replicate', 'Gene', 'Folder','File','FileNo','FileInd','Max']
    if all(n in headers for n in nec):
        print 'Necessary headers found!'
    else:
        print 'Missing essential headers in metabolites file!'
        print headers
        sys.exit(0)
    for index, row in allfiles.iterrows():
        rep=row['Replicate']
        gene=row['Gene']
        tp=row['Type']
        fld=row['Folder']
        flno=row['FileNo']
        flin=row['FileInd']
        fl=row['File']
        gmax=row['Max']
        TFs[gene][rep][tp][fld]['ByFNO'][flno] = fl
        TFs[gene][rep][tp][fld]['ByFIN'][flin] = fl
        if not 'Max' in TFs[gene].keys():
            TFs[gene]['Max']=gmax
        
    return TFs

def readmet(ifile):
	print ifile
	nutrients = NestedDict()
	rdr = csv.reader(open(ifile, 'r'), delimiter=',')
	data = [ln for ln in rdr]
	headers = data[0]
	headin = {hd: headers.index(hd) for hd in headers}
	nec = ['Metabolite', 'EcoCycID', 'Plate', 'Well', 'Index']
	if all(n in headers for n in nec):
		print 'Necessary headers found!'
	else:
		print 'Missing essential headers in metabolites file!'
		print headers
		sys.exit(0)
	for ln in data[1:]:
		pl = str(numerize(ln[headin['Plate']])).strip().encode('ascii', 'ignore')
		wl = str(numerize(ln[headin['Well']])).strip().encode('ascii', 'ignore')
		for hd in headin.keys():
			nutrients[pl][wl][hd] = str(numerize(ln[headin[hd]])).strip().encode('ascii', 'ignore')
	return nutrients


def readdel(ifile):
	print ifile

	deletions = NestedDict()
	rdr = csv.reader(open(ifile, 'r'), delimiter=',')
	data = [ln for ln in rdr]
	headers = data[0]
	headin = {hd: headers.index(hd) for hd in headers}
	nec = ['Replicate', 'Plate', 'Index', 'File', 'X1', 'Y1', 'X2', 'Y2']  # ,'Well'
	if all(n in headers for n in nec):
		print 'Necessary headers found!'
	else:
		print 'Missing essential headers in deletions file!'
		print headers
		sys.exit(0)
	for ln in data[1:]:
		# print ln
		# Reading
		rep = ln[headin['Replicate']].strip().encode('ascii', 'ignore')
		pl = ln[headin['Plate']].strip().encode('ascii', 'ignore')
		fl = ln[headin['File']].strip().encode('ascii', 'ignore')
		indx = numerize(ln[headin['Index']])
		coords = [numerize(ln[headin[hdr]]) for hdr in ['X1', 'Y1', 'X2', 'Y2']]
		# Storing
		deletions[rep][pl][indx]['File'] = fl

		if 'Worms' in deletions[rep][pl][indx].keys():
			worms = deletions[rep][pl][indx]['Worms']
			worms.append(coords)
			deletions[rep][pl][indx]['Worms'] = worms
		else:
			deletions[rep][pl][indx]['Worms'] = [coords]

	return deletions



def wormdel(img, worms):
	for worm in delworms:
		x1, y1, x2, y2 = worm
		xs = [x1, x2]
		ys = [y1, y2]
		img[min(ys):max(ys), min(xs):max(xs), :] = 0
	return img


def readman(ifile):
	print ifile
	manual = NestedDict()
	rdr = csv.reader(open(ifile, 'r'), delimiter=',')
	data = [ln for ln in rdr]
	headers = data[0]
	headin = {hd: headers.index(hd) for hd in headers}
	nec = ['Replicate', 'Plate', 'Index', 'Threshold']  # ,'Well'
	if all(n in headers for n in nec):
		print 'Necessary headers found!'
	else:
		print 'Missing essential headers in deletions file!'
		print headers
		sys.exit(0)

	for ln in data[1:]:
		# print ln
		# Reading
		rep = ln[headin['Replicate']].strip().encode('ascii', 'ignore')
		pl = ln[headin['Plate']].strip().encode('ascii', 'ignore')
		indx = numerize(ln[headin['Index']])
		thrs = numerize(ln[headin['Threshold']])
		thrs = thrs if thrs != "" else 0
		manual[rep][pl][indx] = thrs if thrs != 0.0 else 0.0
	return manual




def freqtable(vector):
	my_series = pd.Series(vector)
	counts = my_series.value_counts()
	ftable = [[key, value] for key, value in dict(counts).iteritems()]
	return ftable


def writecsv(data, ofile, sep=','):
	f = open(ofile, "wb")
	ofile = csv.writer(f, delimiter=sep)  # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
	for row in data:
		# row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
		ofile.writerow(row)
	f.close()
    
    
def labeling3(imghsv, hthres,cthres,size):
	h = imghsv[:, :, 0]
	# s=imghsv[:,:,1]
	v = imghsv[:, :, 2]

	# Filtering
	dh = rank.median(h, disk(3)) / 255.0
	dv = rank.median(v, disk(3)) / 255.0
	# np.max(dh)
	# hmax=np.percentile(h,95)

	hf = dh.copy()
	hf[(hf < hthres).astype(bool)] = 0
	hf = hf / hthres
	# plt.imshow(hf)

	comb_r = hf * dv
	comb = opening(comb_r, selem=disk(3))

	markers = np.zeros_like(comb)
	# Mark background
	markers[comb == 0] = 1
	markers[comb > cthres] = 2

	elevation_map = sobel(v)
	segmentation = watershed(elevation_map, markers)
	segmentation = ndi.binary_fill_holes(segmentation - 1)
	labeled_worms, _ = ndi.label(segmentation)

	for w in list(np.unique(labeled_worms)):
		# print labeled_worms[labeled_worms==w].shape[0]
		if labeled_worms[labeled_worms == w].shape[0] < size:
			labeled_worms[labeled_worms == w] = 0

	return labeled_worms


def labeller4(v,vmin=0.04,vmax=0.06,step=0.001,smin=40000,smax=60000,pmin=3,pmax=5,blobsize=1000,fd=3,od=5):
    #v = imghsv[:, :, 2]

    # Filtering noise
    dv = rank.median(v, disk(fd)) / 255.0
    #Expanding 
    comb = opening(dv, selem=disk(od))
    imagef=skimage.img_as_ubyte(comb)


    umarks=[0]
    wsize=0
    imgprc=100

    iter=0
    total=np.sum(comb >0)

    hthri=vmin
    step=step
    #wsize < smin or wsize > smax or 
    while (imgprc<pmin or imgprc>pmax) and iter<100 and hthri<vmax:
        markers = np.zeros_like(comb)
        # Mark background
        imgprc=np.float(np.sum(comb > hthri)*100)/total

        markers[comb > hthri] = 2
        umarks=np.unique(markers)
        wsize=np.count_nonzero(markers[markers==2])
        #print("Threshold: {:.4f}, Worm size: {:d} Image covered: {:.2f}%".format(hthri,wsize,imgprc))
        hthri+=step
        iter+=1

    markers[comb == 0] = 1

    print("Threshold: {:.4f}, Worm size: {:d} Image covered: {:.2f}%".format(hthri,wsize,imgprc))
    distance = ndi.distance_transform_edt(markers)
    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3)),
                                labels=markers)
    markersmaxi = ndi.label(local_maxi)[0]
    segmentation = watershed(distance, markersmaxi, mask=markers)
    segmentation_fill = ndi.binary_fill_holes(segmentation)
    labeled_worms, _ = ndi.label(segmentation_fill)

    for w in list(np.unique(labeled_worms)):
        # print labeled_worms[labeled_worms==w].shape[0]
        if labeled_worms[labeled_worms == w].shape[0] < blobsize:
            labeled_worms[labeled_worms == w] = 0

    #Smoother worms
    labeled_worms = opening(labeled_worms, selem=disk(10)).astype('uint8')

    wormind = list(np.unique(labeled_worms))
    worms = {w: labeled_worms[labeled_worms == w].shape[0] for w in wormind}
    worms = {wk: wp for wk, wp in worms.items() if wp < 1000000}
    
    contours = measure.find_contours(labeled_worms, 0.8)
    
    return [comb,labeled_worms,contours,worms]


def plotcontours(image,ax,labeled_worms=False,title="",plotcontour=True,plotlabels=True,white=True,ylabel="",xlabel=""):

    plt.sca(ax)
    if white:
        extract = image.copy()
        extract[labeled_worms == 0, :] = [1, 1, 1]
    else:
        extract = image



    ax.set_xlim(0, 1344)
    ax.set_ylim(1024,0)
    
    ax.set_label(None)
    ax.set_axes(None)

    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklines(), visible=False)
    plt.setp(ax.get_xticklines(), visible=False)
    
    ax.set_title(title)
    
    if ylabel!="":
        plt.ylabel(ylabel)
    if xlabel!="":
        plt.xlabel(xlabel)
    
    

    plt.imshow(extract)

    if plotcontour:
        contours = measure.find_contours(labeled_worms, 0.8)
        for n, contour in enumerate(contours):
            plt.plot(contour[:, 1], contour[:, 0], linewidth=1)
    if plotlabels:
        for region in regionprops(labeled_worms):
            minr, minc, maxr, maxc = region.bbox
            plt.text(maxc, maxr, region.label)

def jetimage(image,typeadjust=True):
    image_rescale=exposure.rescale_intensity(image)
    cmap =plt.get_cmap('jet')
    image_rgb = cmap(image_rescale)
    
    if typeadjust:
        image_rgb = skimage.img_as_ubyte(np.delete(image_rgb, 3, 2))
    
    return image_rgb

def map2table(data):
	table = []
	for a in data.keys():
		for b in data[a].keys():
			for c in data[a][b].keys():
				row = [a, b, c, data[a][b][c]]
				table.append(row)
	return table


def gauss(x, p):  # p[0]==mean, p[1]==stdev
	return 1.0 / (p[1] * np.sqrt(2 * np.pi)) * np.exp(-(x - p[0]) ** 2 / (2 * p[1] ** 2))

def fitgauss(xr, yr):
	if xr[0]==0.0:
		xr=xr[1:]
		yr=yr[1:]

	diffs = [j - i for i, j in zip(xr[:-1], xr[1:])]
	# Renormalize to a proper PDF
	Y = yr / np.sum(yr[1:] * diffs)
	X = xr
	# Fit a guassian
	p0 = [0.3, 0.05]  # Inital guess is a normal distribution
	errfunc = lambda p, x, y: gauss(x, p) - y  # Distance to the target function
	p1, success = opt.leastsq(errfunc, p0[:], args=(X, Y))
	fit_mu, fit_stdev = p1
	FWHM = 2 * np.sqrt(2 * np.log(2)) * fit_stdev
	print "FWHM: {}, mu: {}, SD: {}".format(FWHM, fit_mu, fit_stdev)
	return X, Y, fit_mu, fit_stdev

