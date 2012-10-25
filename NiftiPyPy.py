# This is the main program for processing of postmortem brain scans. 
# It has no external library dependencies, allowing it to be run under the 
# PyPy interpreter, which produces a >100x performance increase, decreasing
# run times from hours to a few minuets. Nibabel (and by extension Numpy) 
# library dependencies are offloaded into NiftiLoad.py and NiftiSave.py. 
# You must run NiftiLoad and NiftiSave before and after, respectively.
# temp files (-tmp) are generated but not deleted during this process.
# This decision is designed to allow one to test variable settings more 
# quickly.

##############################################################################
# Variables for user control. 

directory = r"C:/NiftiSubjects/B-5807_B_01/"
file_in = "B-5807_B_01.nii"

# How far from each border of the image should the program remove crud, as a 
# percentage? This does not apply in the rostral-caudal axis, where there is 
# typically no space. 
# This option is present because many postmortem images have a bar of bright
# crud on at least one side that makes by-section normalization difficult.
# It would probably be most ideal to manually remove this bar from the image
# prior to sectioning.
pct_border = .05

# This setting normalizes data with inconsistent scan brightness within a scan.
normalize_slices = True 

# This setting normalizes the image intensity to a standard value across images
# of differing brightness. If set to false, you will need to play with your 
# gray and white matter cutoffs below to get good results
normalize_image = True

# Tell the program whether to invert the slice values around a threshold 
# (wm_max, see below). This is the main step that inverts gray and white 
# intensities.
intensity_correct = True

# This setting removes the bright "rind" artifact produced by inverting 
# dim voxels (half brain/half background) on the surface of the brain. 
remove_wm_rind = True

# Force pixels? This is to produce better freesurfer output. Everything that's
# above wm_min gets one value, and everything below it gets another.
make_wm_mask = False

# If making a mask, we need to a series of blurs and sharpens to remove 
# occasional misidentified, isolated pixels? How many times should this be 
# performed?
cleanup_reps = 5

##############################################################################
# It is strongly suggested that these variables not be altered unless you 
# know exactly what you are doing!

# Set your intensity values manually here.
# Do not alter if you are using normalization!
# Normalization assumes wm_max & gm_min = 1700, wm_min = 300, gm_max = 3000, 
# background = 0
wm_max = 1700
wm_min = 300
gm_max = 3000
gm_min = 1700
background = 0

# Normalization intensity
norm_intens = 1785

# The intensity we will set our mask at. Don't overlap this with voxel values!
mask_set = 10000

# How far away should I check for isolation? If a pixel below the minimum wm 
# threshhold encounters nothing but space in both directions in all 3 
# dimensions, it will be deleted.
# This step is performed because otherwise lots of background becomes bright 
# during normalization.
search_dist = 5

# How many pixels need to be seen, total, in all 3 axes, with a blur of one 
# pixel from the search destination, for the pixel to be included? 
inclusion_req = 3

# How many slices away in a rostral-caudal axis should sections be averaged?
# If you see an occasional section that sticks out in intensity, try increasing
# this distance.
slice_dist = 3

# How thick is the typical rind, in mm?
rind_thickness = 0.7

# How far beyond the rind should we look for open space, in mm?
# Note that exploration distance is thickness + exploration divided by pixel 
# size in each dimension, then rounded *down* (because that is how Python 2 
# handles integers.
rind_explore = 1.3

# How many bg voxels need to be found in this space for it to count as a 
# border cell?
rind_cutoff = 4

# How far from cortical surfacet to place protective mask, in mm?
mask_dist = 2.0

# Blur settings: this is how far away voxels draw information from while 
# fixing the rind, in mm
blur_range = 2.5

# the first file is unchanged, so you can try multiple processing approaches
# without reloading from the raw file.
raw_data_suffix = "-tmprawdata"
proc_data_suffix = "-tmpprocdata"
zoom_suffix = "-tmpzoom"

##############################################################################
# Program begins here
# Note: I'm using marshal instead of the more standard pickle to move my files 
# around because pickle appears to be incompatible with the current release of 
# PyPy.
import time
import marshal

# For background blur at end of processing set
# Experimental as of 10/25/12
import random


def raw_name(file_in):
    file_cutoff = file_in.find(".")
    raw_name = file_in[:file_cutoff]
    return raw_name


def make_copy(img_data):
    """This is the equivalent of the copy.deepcopy function, but much faster
    because it's tailored to this data type. Using deepcopy instead costs 
    an extra ~5 seconds/copy on my machine with a data set this large"""
    img_out = ([[[int(e) for e in sag_row] for sag_row in ax_row] 
                for ax_row in img_data])
    return img_out


def remove_boundaries(img_data, pct_border, background):
    for ax_index, ax_row in enumerate(img_data):
        for sag_index, sag_row in enumerate(ax_row):
            for cor_index, voxel in enumerate(sag_row):
                if (ax_index < pct_border * len(img_data) or ax_index > 
                    (len(img_data) - pct_border * len(img_data))):
                    img_data[ax_index][sag_index][cor_index] = background
                if (sag_index < pct_border * len(img_data) or sag_index > 
                    (len(img_data) - pct_border * len(img_data))):
                    img_data[ax_index][sag_index][cor_index] = background
        print ("border_strip row " + str(ax_index) + " done" + " time: " + 
               str(time.clock()))
    return img_data


def isolation_check(search_dist, inclusion_req, ax_index, sag_index, 
                    cor_index, wm_min):
    for search_loc in range(search_dist - 1, search_dist + 2):
        inclusion_count = 0
        if inclusion_count < inclusion_req:
            try:
                if (img_data[ax_index + search_loc][sag_index][cor_index] > 
                    wm_min):
                    inclusion_count +=1
            except:
                pass
        if inclusion_count < inclusion_req:
            try:
                if (img_data[ax_index][sag_index + search_loc][cor_index] > 
                    wm_min):
                    inclusion_count +=1
            except:
                pass
        if inclusion_count < inclusion_req:
            try:
                if (img_data[ax_index][sag_index][cor_index + search_loc] > 
                    wm_min):
                    inclusion_count +=1
            except:
                pass
    return inclusion_count


def clean_bg(img_data, search_dist, inclusion_req, wm_min, background):
    for ax_index, ax_row in enumerate(img_data):
        for sag_index, sag_row in enumerate(ax_row):
            for cor_index, voxel in enumerate(sag_row):
                if voxel < wm_min:
                    inclusion_count = isolation_check(search_dist, 
                                                      inclusion_req, ax_index, 
                                                      sag_index, cor_index, 
                                                      wm_min)
                    if inclusion_count < inclusion_req:
                        img_data[ax_index][sag_index][cor_index] = background
                        
    for ax_index, ax_row in enumerate(img_data):
        for sag_index, sag_row in enumerate(ax_row):
            for cor_index, voxel in enumerate(sag_row): 
                iso_count = 0
                for ax_loc in range(ax_index - 1, ax_index + 2):
                    for sag_loc in range(sag_index - 1, sag_index + 2):
                        for cor_loc in range(cor_index - 1, cor_index + 2):
                            try:
                                if (img_data[ax_loc][sag_loc][cor_loc] == 
                                    background):
                                    iso_count += 1
                            except:
                                pass
                if iso_count >= 26:
                    img_data[ax_index][sag_index][cor_index] = background
    return img_data


def intensity_norm_calc(intensity_list, cor_loc, slice_dist, gm_max):
    intensity_num = 1
    intensity_sum = intensity_list[cor_loc]
    for cor_slice in range(-(slice_dist), slice_dist + 1):
        try:
            if intensity_list[cor_loc + cor_slice] < gm_max:
                intensity_sum = intensity_sum + intensity_list[cor_loc + 
                                                               cor_slice]
                intensity_num += 1
        except:
            pass
    intensity_avg = intensity_sum / intensity_num
    return intensity_avg


def slice_normalization(img_data, slice_dist, wm_min, wm_max, gm_min, gm_max):
    all_intensity_list = []
    rep_intensity_list = []
    for cor_loc in range(len(img_data[0][0])):
        cor_list = []
        for ax_index, ax_row in enumerate(img_data):
            for sag_index, sag_row in enumerate(ax_row):
                #if img_mask_com[ax_index][sag_index][cor_loc] == mask_set:
                if img_data[ax_index][sag_index][cor_loc] > wm_min:
                    cor_list.append(img_data[ax_index][sag_index][cor_loc])
        
        if cor_list and len(cor_list) > 10:
            cor_list.sort()
            max_intensity = cor_list[-10]
            rep_intensity_list.append(max_intensity)
        else:
            max_intensity = gm_max
        all_intensity_list.append(max_intensity) 
    avg_intensity_loc = (len(rep_intensity_list) / 2)
    avg_high_intensity = rep_intensity_list[avg_intensity_loc]
    #avg_high_intensity = 2560
    for cor_loc in range(len(img_data[0][0])):
        slice_intensity_avg = intensity_norm_calc(all_intensity_list, cor_loc, 
                                                  slice_dist, gm_max)
        for ax_index, ax_row in enumerate(img_data):
            for sag_index, sag_row in enumerate(ax_row):

                img_data[ax_index][sag_index][cor_loc] = (
                    int(img_data[ax_index][sag_index][cor_loc] 
                    * float(avg_high_intensity) / float(slice_intensity_avg)))
        print ("norm row " + str(cor_loc) + " done" + " time: " + 
               str(time.clock()))
    return img_data  


def image_normalization(img_data):
    img_raw = []
    # still trying to figure out a sharp list comprehension for this
    for sag_row in img_data:
        for ax_row in sag_row:
            for voxel in ax_row:
                img_raw.append(voxel)
    # bin the raw data and look for the mode            
    bin_size = 5
    bin_min = 0
    bin_max = 10000
    bin_list = [0] * ((bin_max - bin_min) / bin_size)
    for voxel in img_raw:
        binselect = int(float(voxel) / bin_size - bin_min)
        bin_list[binselect] = bin_list[binselect] + 1
    max_bin_val, max_bin_loc = 0, 0
    for binlist_bin in range(5, bin_max, bin_size):
        if bin_list[binlist_bin / bin_size] > max_bin_val:
            max_bin_val = bin_list[binlist_bin / bin_size]
            max_bin_loc = binlist_bin
    peak_diff = max_bin_loc - 1785
    for ax_index, ax_row in enumerate(img_data):
        for sag_index, sag_row in enumerate(ax_row):
            for cor_index, voxel in enumerate(sag_row):
                if voxel > 0:
                    img_data[ax_index][sag_index][cor_index] = voxel - peak_diff
    return img_data


def smart_voxel_flip(img_data, wm_min, wm_max, gm_min, gm_max):
    wm_range = float(wm_max - wm_min)
    gm_range = float(gm_max - gm_min)
    for ax_index, ax_row in enumerate(img_data):
        for sag_row in ax_row:
            for cor_index, voxel in enumerate(sag_row):
                if voxel >= wm_min and voxel < wm_max:
                    intensity_rating = (wm_max - voxel) / wm_range 
                    sag_row[cor_index] = int(gm_min + gm_range * 
                                             (intensity_rating))
                elif voxel >= gm_min and voxel < gm_max:
                    intensity_rating = (voxel-gm_min) / gm_range 
                    sag_row[cor_index] = int(wm_max - wm_range * 
                                             (intensity_rating))
                else: 
                    sag_row[cor_index] = int(background)
        print ("flip row " + str(ax_index) + " done" + " time: " + 
               str(time.clock()))
    return img_data


def row_search(row, mask_dist, img_zoom, axis):
    start_mask_track, stop_mask_track = 0, 0
    start_mask, stop_mask = False, False
    for row_index, voxel in enumerate(row):
        if not start_mask and voxel > gm_min:
            start_mask_track += 1
        if not start_mask and start_mask_track >= mask_dist/img_zoom[axis]:
            start_mask = row_index
    row.reverse()
    for row_index, voxel in enumerate(row):
        if not stop_mask and voxel > gm_min:
            stop_mask_track += 1
        if not stop_mask and stop_mask_track >= mask_dist/img_zoom[axis]:
            stop_mask = len(row)- row_index
    return start_mask, stop_mask


def brain_mask_cor(img_mask, img_data, img_zoom, mask_set, gm_min, background, 
                   mask_dist=mask_dist):
    for ax_index, ax_row in enumerate(img_mask):
        for sag_index, sag_row in enumerate(ax_row):
            start_mask, stop_mask = row_search(sag_row, mask_dist, img_zoom, 
                                               axis=2)
            if start_mask and stop_mask:
                for cor_index, voxel in enumerate(sag_row):
                    if cor_index > start_mask and cor_index < stop_mask:
                        sag_row[cor_index] = mask_set                    
        print "cor mask row: " + str(ax_index) + " time: " + str(time.clock())
    return img_mask


def brain_mask_ax(img_mask, img_data, img_zoom, mask_set, gm_min, background, 
                  mask_dist=mask_dist):
    for cor_loc in range(len(img_data[0][0])):
        cor_grid = [[img_data[ax_entry][sag_entry][cor_loc] 
                     for sag_entry in xrange(len(img_data[0]))] 
                    for ax_entry in xrange(len(img_data))]
        for ax_index, ax_row in enumerate(cor_grid):
            start_mask, stop_mask = row_search(ax_row, mask_dist, img_zoom, 
                                               axis=0)
            if start_mask and stop_mask:
                for sag_index, voxel in enumerate(ax_row):
                    if (sag_index > start_mask and 
                        sag_index < stop_mask and 
                        img_mask[ax_index][sag_index][cor_loc] == mask_set):
                        pass
                    else: 
                        img_mask[ax_index][sag_index][cor_loc] = img_data[ax_index][sag_index][cor_loc]                  
        print "ax mask row: " + str(cor_loc) + " time: " + str(time.clock())
    return img_mask
 
 
def brain_mask_sag(img_mask, img_data, img_zoom, mask_set, gm_min, background, 
                   mask_dist=mask_dist):         
    for cor_loc in range(len(img_data[0][0])):
        #make an array containing the contents of each coronal slice
        cor_grid = [[img_data[ax_entry][sag_entry][cor_loc] 
                     for ax_entry in xrange(len(img_data))] 
                    for sag_entry in xrange(len(img_data[0]))]    
        for sag_index, sag_row in enumerate(cor_grid):
            start_mask, stop_mask = row_search(sag_row, mask_dist, img_zoom, 
                                               axis=1)
            if start_mask and stop_mask:
                for ax_index, voxel in enumerate(sag_row):
                    if (ax_index > start_mask and 
                        ax_index < stop_mask and 
                        img_mask[ax_index][sag_index][cor_loc] == mask_set):
                        pass
                    else:
                        img_mask[ax_index][sag_index][cor_loc] = img_data[ax_index][sag_index][cor_loc]                    
        print "sag mask row: " + str(cor_loc) + " time: " + str(time.clock())   
    return img_mask


def bound_check(img_data, ax_index, r_dist_ax, sag_index, r_dist_sag, 
                cor_index, r_dist_cor, background):
    boundary_load = 0
    for ax_loc in range(ax_index - r_dist_ax, 
                         ax_index + r_dist_ax + 1):
        for sag_loc in range(sag_index - r_dist_sag, 
                            sag_index + r_dist_sag + 1):
            for cor_loc in range(cor_index - r_dist_cor, 
                                 cor_index + r_dist_cor + 1):
                try:
                    if (img_data[ax_loc][sag_loc][cor_loc] == 
                        background):
                        boundary_load += 1         
                except:
                    pass
    return boundary_load


def remove_rind(img_data, img_zoom, gm_min, gm_max, wm_min, wm_max, 
                background, img_mask, mask_set, rind_thickness, 
                rind_explore, rind_cutoff):
    img_chg = make_copy(img_data)
    # pre-calculate distance search settings here for efficiency
    r_dist_ax = int((rind_thickness + rind_explore)/img_zoom[0])
    r_dist_sag = int((rind_thickness + rind_explore)/img_zoom[1])
    r_dist_cor = int((rind_thickness + rind_explore)/img_zoom[2])
    r_blur_ax = int(blur_range/img_zoom[0])
    r_blur_sag = int(blur_range/img_zoom[1])
    r_blur_cor = int(blur_range/img_zoom[2])

    for ax_index, ax_row in enumerate(img_data):
        for sag_index, sag_row in enumerate(ax_row):
            for cor_index, voxel in enumerate(sag_row):
                if voxel > wm_min:
                    boundary_load = bound_check(img_data, ax_index, r_dist_ax, 
                                                sag_index, r_dist_sag, 
                                                cor_index, r_dist_cor, 
                                                background)
                    if boundary_load >= rind_cutoff:
                        voxel_val, voxel_add = 0, 0
                        voxel_mask = img_mask[ax_index][sag_index][cor_index]
                        gm_max_temp = gm_max
                        while not voxel_add:
                            for ax_loc in range(ax_index - r_blur_ax, 
                                                ax_index + r_blur_ax + 1):
                                for sag_loc in range(sag_index - r_blur_sag, 
                                                     sag_index + 
                                                     r_blur_sag + 1):
                                    for cor_loc in range(cor_index - 
                                                         r_blur_cor, 
                                                         cor_index + 
                                                         r_blur_cor + 1):
                                        try:
                                            voxel_iter = img_data[ax_loc][sag_loc][cor_loc]
                                            if voxel_mask == mask_set and voxel_iter > gm_min:
                                                voxel_val = voxel_val + voxel_iter
                                                voxel_add += 1
                                            elif voxel_iter > gm_min and voxel_iter < gm_max_temp:
                                                voxel_val = voxel_val + voxel_iter
                                                voxel_add += 1
                                        except:
                                            pass
                            gm_max_temp = gm_max_temp + (wm_max - wm_min) * .10
                        img_chg[ax_index][sag_index][cor_index] = int(float(voxel_val) / float(voxel_add))
        print "rind row: " + str(ax_index) + " time: " + str(time.clock())
    return img_chg

def force_mask(img_data, gm_min, gm_max, wm_min, wm_max, background):
    for ax_index, ax_row in enumerate(img_data):
        for sag_index, sag_row in enumerate(ax_row):
            for cor_index, voxel in enumerate(sag_row):
                if voxel > wm_min:
                    img_data[ax_index][sag_index][cor_index] = wm_max
                elif voxel > background:
                    img_data[ax_index][sag_index][cor_index] = gm_min
                else:
                    img_data[ax_index][sag_index][cor_index] = background
    return img_data


def pixel_cleanup(img_data, img_zoom, gm_min, gm_max, wm_min, wm_max, background):
#    r_blur_ax = int(blur_range/img_zoom[0])
#    r_blur_sag = int(blur_range/img_zoom[1])
#    r_blur_cor = int(blur_range/img_zoom[2])
    
    r_blur_ax = 3
    r_blur_sag = 3
    r_blur_cor = 1
    
    for ax_index, ax_row in enumerate(img_data):
        for sag_index, sag_row in enumerate(ax_row):
            for cor_index, voxel in enumerate(sag_row):
                if voxel > background:
                    voxel_sum, voxel_num = 0, 0
                    for ax_loc in range(ax_index - r_blur_ax, ax_index + r_blur_ax + 1):
                        for sag_loc in range(sag_index - r_blur_sag, sag_index + r_blur_sag + 1):
                            for cor_loc in range(cor_index - r_blur_cor, 
                                                         cor_index + r_blur_cor + 1):
                                voxel_sum = voxel_sum + img_data[ax_loc][sag_loc][cor_loc]
                                voxel_num += 1 
    return img_data


def bright_image(img_data):
    for ax_index, ax_row in enumerate(img_data):
        for sag_index, sag_row in enumerate(ax_row):
            for cor_index, voxel in enumerate(sag_row):
                if voxel > background:
                    img_data[ax_index][sag_index][cor_index] = voxel * 1.5
    return img_data


def bg_natural(img_data):
    for ax_index, ax_row in enumerate(img_data):
        for sag_index, sag_row in enumerate(ax_row):
            for cor_index, voxel in enumerate(sag_row):
                if voxel == background:
                    img_data[ax_index][sag_index][cor_index] = random.randint(1,100)
    return img_data


print time.clock()
print "loading data"
raw_name = raw_name(file_in)
datafile = open((directory + raw_name + raw_data_suffix), 'rb')
img_data = marshal.load(datafile)
zoomfile = open((directory + raw_name + zoom_suffix), 'rb')
img_zoom = marshal.load(zoomfile)

##############################################################################
# functions that remove crud go here
print "removing junk from the borders of the image"
img_data = remove_boundaries(img_data, pct_border, background)

print "removing bright speckles from background prior to normalization"
img_data = clean_bg(img_data, search_dist, inclusion_req, wm_min, background)

##############################################################################
# functions that flip pixel intensity go here
# normalization is NOT ACTIVE yet
# if intensity_correct and not use_manual_values:

if intensity_correct and normalize_slices: 
    img_data = slice_normalization(img_data, slice_dist, wm_min, wm_max, 
                                   gm_min, gm_max)

if normalize_image:
    img_data = image_normalization(img_data)

if intensity_correct:
    img_data = smart_voxel_flip(img_data, wm_min, wm_max, gm_min, gm_max)

##############################################################################
# functions that maintain pixel intensity and clean up image go below here

wm_min, gm_min = gm_min, wm_min
wm_max, gm_max = gm_max, wm_max

if remove_wm_rind:  
    img_mask = make_copy(img_data)
    img_mask = brain_mask_cor(img_mask, img_data, img_zoom, mask_set, gm_min, 
                              background, mask_dist=mask_dist)
    img_mask = brain_mask_ax(img_mask, img_data, img_zoom, mask_set, gm_min, 
                             background, mask_dist=mask_dist)
    img_mask = brain_mask_sag(img_mask, img_data, img_zoom, mask_set, gm_min, 
                              background, mask_dist=mask_dist)
        
    print "removing rind"

    img_data = remove_rind(img_data, img_zoom, gm_min, gm_max, wm_min, 
                           wm_max, background, img_mask, mask_set, 
                           rind_thickness, rind_explore, rind_cutoff)
    
if make_wm_mask:
    img_data = pixel_cleanup(img_data, img_zoom, gm_min, gm_max, wm_min, 
                             wm_max, background)
    
    for _ in xrange(cleanup_reps):
        img_data = force_mask(img_data, gm_min, gm_max, wm_min, wm_max, 
                              background)
        
        img_data = pixel_cleanup(img_data, img_zoom, gm_min, gm_max, wm_min, 
                                 wm_max, background)
    

    img_data = force_mask(img_data, gm_min, gm_max, wm_min, wm_max, 
                          background)

# Brighten the image up from its natural levels to improve freesurfer's 
# ability to read it    
img_data = bright_image(img_data)

# Boost the background to random 0-100 values in the hopes of improving seg.
# with GM by making things more natural looking
# Experimental as of 10/25/12
img_data = bg_natural(img_data)
    
print time.clock()
print "saving raw data, use NiftiSave.py to convert to .nii"
datafile = open((directory + raw_name + proc_data_suffix), 'w+b')
marshal.dump(img_data, datafile)
print "done"

# I've left these lines in should you want to save the mask and
# examine it.
#print "saving mask data"
#imagemask = open((directory + raw_name + "-imgmask"), 'w+b')
#marshal.dump(img_mask, imagemask)
