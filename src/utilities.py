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

from skimage.morphology import disk, opening, dilation, square, watershed, skeletonize

from skimage.morphology import erosion, white_tophat, black_tophat, closing

from skimage import exposure

import datetime
from multiprocessing import Pool, Manager

import cPickle



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



def wormdel(img, delworms):
    for worm in delworms:
        x1, y1, x2, y2 = worm
        xs = [x1, x2]
        ys = [y1, y2]
        if np.ndim(img)==2:
            img[min(ys):max(ys), min(xs):max(xs)] = 0
        else:
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
        wsize=np.count_nonzero(markers[markers>1])
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

def multiprocessor(cores,func,args,nargs=0):

    #if total==0:
    if nargs==0:
        alen=len(args)
    else:
        alen=nins
    
    m = Manager()
    q = m.Queue()
    p=Pool(cores)
    argsq=[ arg + (q,) for arg in args]
    result = p.map_async(func, argsq)
    start=time.time()
    prcprev=0
    while True:
        if result.ready():
            break
        else:
            lcomp=q.qsize()
            #print lcomp
            prc = float(lcomp)*100/float(alen)
            if prc>prcprev:
                timepassed=float(time.time()-start)
                if prc!=0:
                    total=int((timepassed*100/prc))
                    remaining=int((timepassed*100/prc)-timepassed)
                else:
                    total=int('inf')
                    remaining=int('inf')
                print "Done {0:3d}%, {1:>5}/{2:<5} remaining: {3:<5} total: {4:<5}".format(int(prc),str(lcomp),str(alen),str(datetime.timedelta(seconds=remaining)),str(datetime.timedelta(seconds=total)))
                prcprev=prc
        time.sleep(5)
    print 'Collecting results....'
    results=result.get()
    print 'Done!'
    #print results
    return results

def collector_M((fln,settings,q),savetype='BW',colordepth='uint8'):
    
    for key,val in settings.items():
        exec(key + '=val')
    
    image = tiff.imread(fln)
    #Save only one layer, as they are redundant
    if np.ndim(image)>3:
        raise ValueError('Tiff image has more than 3 dimensions!')
    if np.ndim(image)==3:
        if image.shape[2]>3:
            raise ValueError('Tiff image has more than 3 layers!')
            

    if savetype=='RGB':
        if np.ndim(image)==3:
            imagesave=image
        else:
            imagesave=color.gray2rgb(image)
            
    elif savetype=='BW':
        if np.ndim(image)==3:
            if np.sum(image[:,:,0]-image[:,:,1])==0 and np.sum(image[:,:,1]-image[:,:,2])==0:
                imagesave=image[:,:,0]
            else:
                imagesave=color.rgb2hsv(image)[:,:,2]
        else:
            imagesave=image
            
    elif savetype=='HV':
        if np.ndim(image)==3:
            imagesave=color.rgb2hsv(image)[:,:,(0,2)].astype('float16')
        else:
            v=color.rgb2hsv(color.gray2rgb(image))[:,:,2].astype('float16')
            h=np.zeros_like(v).astype('float16')
            imagesave=np.stack([h,v],axis=2)
     
    elif savetype=='HSV':
        if np.ndim(image)==3:
            imagesave=color.rgb2hsv(image).astype('float16')
        else:
            v=color.rgb2hsv(color.gray2rgb(image))[:,:,2].astype('float16')
            h=np.zeros_like(v).astype('float16')
            imagesave=np.stack([h,h,v],axis=2)
    else:
        imagesave=image
            

    if imagesave.dtype!=colordepth:
        if colordepth=='uint8':
            imagesave=skimage.img_as_ubyte(imagesave)
        elif colordepth=='uint16':
            imagesave=skimage.img_as_uint(imagesave)
        elif colordepth=='int':
            imagesave=skimage.img_as_int(imagesave)
        elif colordepth=='float':   
            imagesave=skimage.img_as_float(imagesave)
        else:
            print "Unknown color depth: {}".format(colordepth)
            
    q.put(fln)
    
    return [fln,imagesave]

def collectsubset(ofile,imagesets,imagesindex,size=100):
    
    oimagest="{}.pkl".format(ofile)
    rndsel=np.random.choice(range(0,imagesets[0].shape[0]-1), size, replace=False)

    imagesets_ar=[]
    imagesindex_f={ ik:ival for ik,ival in imagesindex.items() if ival in rndsel }
    
    
    #print len(rndsel), len(imagesindex_f)

    imagesindex_t={ik:ii for ii,(ik,ival) in enumerate(imagesindex_f.items())}
    
    for iset in imagesets:
        isetdata=[]
        for ik,ival in imagesindex_f.items():
            isetdata.append(iset[ival])
        imagesets_ar.append(isetdata)


    print "Stacking..."
    imagesets_s=tuple([np.stack(iset, axis=0) for iset in imagesets_ar])
    iwrite=(imagesindex_t,)+imagesets_s
    f = open(oimagest, "w")
    print "Writing..."
    cPickle.dump(iwrite, f)
    f.close()
    print "Zipping..."
    os.system("gzip -f {}".format(oimagest))
    print("Saving images complete!")






def labeller4_M((fln, image_raw, settings,q), vmin=0.04, vmax=0.06, step=0.001, smin=40000, smax=60000, pmin=3, pmax=5, blobsize=1000, maxsize=1000000, fd=3 ,od=5, sd=10):
    
    for key,val in settings.items():
        exec(key + '=val')
    
    if np.ndim(image_raw)==2:
        #Stays uint8
        image=color.gray2rgb(image_raw)
        v=image_raw
    else:
        image=image_raw
        if np.sum(image[:,:,0]-image[:,:,1])==0 and np.sum(image[:,:,1]-image[:,:,2])==0:
            v=image_raw
        else:
            #Will cause change of type
            v=skimage.img_as_ubyte(color.rgb2hsv(image_raw))[:,:,2]
    
    # Filtering noise
    #Gets back to 255 scale!
    dv = rank.median(v, disk(fd))
    #Expanding 
    comb = opening(dv, selem=disk(od))
    
    comb=skimage.img_as_float(comb)
    
    imagef=skimage.img_as_ubyte(comb)

    wsize=0
    imgprc=100

    iter=0
    isize=comb.shape[0]*comb.shape[1]
    
    
    vthri=vmin
    step=step
    #wsize < smin or wsize > smax or 
    while (imgprc<pmin or imgprc>pmax) and iter<100 and vthri<vmax:
        vthrio=vthri
        vthri+=step
        # Mark background
        imgprc=np.float(np.sum(comb > vthri)*100)/float(isize)
        
        wsize=np.count_nonzero(comb[comb>vthri])
        #print("Threshold: {:.4f}, Worm size: {:d} Image covered: {:.2f}%".format(vthri,wsize,imgprc))
        iter+=1
        
    markers = np.zeros_like(comb)
    markers[comb > vthrio] = 2
    markers[comb == 0] = 1

    labeled_worms=waterseg(markers,sd,blobsize,maxsize)
    
    #Data collection step
    #In Velocity background is summed in only one channel!
    results_data=getresults(image,labeled_worms)
    
    q.put(fln)
    
    return [fln,imagef,labeled_worms,results_data]





def labeller5_M((fln,image_raw,settings,q), wgoal=60000, vmin=0.04, vmax=0.06, step=0.001, smin=40000, smax=60000, pmin=3, pmax=5, blobsize=1000, maxsize=1000000, fd=20, od=5,sd=10):
    
    for key,val in settings.items():
        exec(key + '=val')
    
    if np.ndim(image_raw)==2:
        image=color.gray2rgb(image_raw)
        v=image_raw
    elif np.ndim(image_raw)>3:
        raise ValueError('Tiff image has more than 3 dimensions!')
    else:
        if image_raw.shape[2]==3:
            image=image_raw
            v=skimage.img_as_ubyte(color.rgb2hsv(image_raw))[:,:,2]
        elif image_raw.shape[2]==2:
            image=image_raw[:,:,1]
            v=skimage.img_as_ubyte(image)
        else:
            raise ValueError('Tiff image has more than 3 layers!')
        
    
    # Filtering noise
    #Will change type!
    dv = rank.median(v, disk(fd))
    #Expanding 
    #comb = opening(dv, selem=disk(od))
    
    comb=skimage.img_as_float(dv)
    
    imagef=skimage.img_as_ubyte(comb)
    
    isize=comb.shape[0]*comb.shape[1]

    wsize=isize
    wsizeo=isize
    vthri=vmin
    step=0.001
    wdiff=abs(wgoal-wsize)
    wdiffo=abs(wgoal-wsize)
    
    iter=0

    while wdiff<=wdiffo and iter<100 and vthri<0.06:
        vthrio=vthri
        vthri+=step

        wdiffo=wdiff
        wsize=np.count_nonzero(comb[comb > vthri])

        wdiff=abs(wsize-wgoal)

        imgprc=np.float(np.sum(comb > vthri)*100)/isize
        #print("Threshold: {:.4f}, Worm size: {:d}, Wdiff: {:d} Image covered: {:.2f}%".format(vthri,wsize,wdiff,imgprc))
        iter += 1


    #umarks=np.unique(markers)
    markers = np.zeros_like(comb)
    markers[comb > vthrio] = 2
    markers[comb ==0] = 1

    labeled_worms=waterseg(markers,sd,blobsize)
    
    results_data=getresults(image,labeled_worms)

    #Data collection step
    #In Velocity background is summed in only one channel!
    q.put(fln)
    
    return [fln,imagef,labeled_worms,results_data]


#def waterseg(markers,sd,blobsize,maxsize):
#    
#    distance = ndi.distance_transform_edt(markers)
#    local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3)),
#                                labels=markers)
#    markersmaxi = ndi.label(local_maxi)[0]
#    segmentation = watershed(distance, markersmaxi, mask=markers)
#    segmentation_fill = ndi.binary_fill_holes(segmentation)
#    labeled_worms, _ = ndi.label(segmentation_fill)
#
#    #Filter out noise and patches
#    for w in list(np.unique(labeled_worms)):
#        # print labeled_worms[labeled_worms==w].shape[0]
#        if labeled_worms[labeled_worms == w].shape[0] < blobsize:
#            labeled_worms[labeled_worms == w] = 0
#            
#    #Smoother worms
#    labeled_worms = opening(labeled_worms, selem=disk(sd)).astype('uint8')
#    
#    return (labeled_worms)

def waterseg(markers,sd,blobsize,approach="TF"):

    if approach=='Biolog':
        #Old approach is wrong
        elevation_map = sobel(v)
        segmentation = watershed(elevation_map, markers)
        segmentation = ndi.binary_fill_holes(segmentation - 1)
        labeled_worms, _ = ndi.label(segmentation)
    else:
        distance = ndi.distance_transform_edt(markers)
        local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3)),
                                    labels=markers)
        markersmaxi = ndi.label(local_maxi)[0]
        segmentation = watershed(distance, markersmaxi, mask=markers)
        segmentation_fill = ndi.binary_fill_holes(segmentation)
        labeled_worms, _ = ndi.label(segmentation_fill)
    
        
    #Filter out noise and patches
    for w in list(np.unique(labeled_worms)):
        # print labeled_worms[labeled_worms==w].shape[0]
        if labeled_worms[labeled_worms == w].shape[0] < blobsize:
            labeled_worms[labeled_worms == w] = 0
        
    #Smoother worms
    labeled_worms = skimage.img_as_ubyte(opening(labeled_worms, selem=disk(sd)))

    return labeled_worms

def getdistributions(layer,labeled_worms,levels):
    
    wormind=[win for win in list(np.unique(labeled_worms)) if win!=0]
    if len(wormind)>0:
        results_data=[]
        for w in wormind:
            bright1D = np.around(layer[labeled_worms == w].ravel(), 3)
            ftable = np.array(freqtable(bright1D))
            ftable = ftable[np.argsort(ftable[:, 0])]
            fdict = {"{0:.3f}".format(freq[0]): int(freq[1]) for freq in ftable}
            data.append([w]+[fdict[key] if key in fdict.keys() else 0 for key in levels])
    else:
        results_data=[[]]

    return results_data


def getresults(image,labeled_worms):
    b_sz=np.sum(np.isin(labeled_worms,0))
    b_sum=np.sum(image[labeled_worms==0,0])
    b_mean=b_sum/b_sz

    wormind=[win for win in list(np.unique(labeled_worms)) if win!=0]

    if len(wormind)>0:
        results_data=[]
        for w in wormind:
            w_sz=np.sum(np.isin(labeled_worms, w))
            #In Velocity worm outline is summed in all 3 channels!
            w_sum=np.sum(image[labeled_worms==w,:])
            w_mean=w_sum/w_sz
            results_data.append([w,w_sz,w_sum,w_mean,b_sz,b_sum,b_mean])
    else:
        results_data=[[]]

    return results_data



def plotcontours(image,ax,labeled_worms=False,title="",plotcontour=True,plotlabels=True,white=True,ylabel="",xlabel=""):

    plt.sca(ax)
    if image.dtype=='uint8':
        vmax=255
    else:
        vmax=1
        
    if white:
        extract = image.copy()
        if np.ndim(extract)==2:
            extract[labeled_worms == 0] = vmax
        else:
            extract[labeled_worms == 0, :] = [vmax,vmax,vmax]
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

