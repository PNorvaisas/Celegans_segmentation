import os
import sys
import tifffile as tiff
import textwrap as tw
import csv
import pandas
import matplotlib

#matplotlib.use('MacOSX') # For visualising
#matplotlib.use('Agg') # For saving
#matplotlib.use('MacOSX') #TkAgg

import matplotlib.pyplot as plt


plt.interactive(False)


import numpy as np
import itertools as IT
import time


import skimage
from skimage import io as skio
from skimage import color
from skimage import data
from skimage import img_as_float
from skimage import filters
from skimage.util.dtype import dtype_range
from skimage.util import img_as_ubyte
from skimage.feature import canny
from skimage.filters import sobel
from skimage.morphology import disk, opening, dilation, square, watershed
from skimage.morphology import erosion, white_tophat, black_tophat, closing
from skimage import exposure


from skimage.filters import rank
from skimage.filters.rank import median, mean
from skimage import measure
from skimage.measure import label, regionprops

#from skimage.filters import threshold_otsu, threshold_adaptive

#from skimage.feature import peak_local_max
#from scipy import misc


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
	my_series = pandas.Series(vector)
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

