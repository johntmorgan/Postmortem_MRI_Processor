# This is the main program for processing of postmortem brain scans. 
# It has no external library dependencies, allowing it to be run under the 
# PyPy interpreter, which produces a >100x performance increase, decreasing
# run times from hours to a few minutes. NiBabel (and by extension Numpy) 
# library dependencies are offloaded into NiftiLoad.py and NiftiSave.py. 
# You must run NiftiLoad and NiftiSave before and after, respectively.
# temp files (-tmp) are generated but not deleted during this process.
# This decision is designed to allow one to test variable settings more 
# quickly.

##############################################################################
# Variables for user control. 

directory = r"C:/NiftiSubjects/MIND Cases/KC Project Cases/"
file_in = "UMB-4226_L_01.nii"

# the first file is unchanged, so you can try multiple processing approaches
# without reloading from the raw file.
raw_data_suffix = "-tmprawdata"
proc_data_suffix = "-tmpprocdata"
zoom_suffix = "-tmpzoom"

# Blank out everything on the postive side of the mirror. If this is 
# on the right side of the image, use make left and make right as you would 
# expect. If this is on the left, *flip* which command you use.
# If the blanking is on the top or bottom half of the image, turn on 
# flip_vertical.
# run with mirror_line = 0
check_orient = False

# Where you want to mirror, how far from midline of image?
# You will probably need to try different mirror_line values to make a nice
# mirror. This setting is *in voxels* so you should adjust based on the 
# coordinate output.
mirror_line = 0
make_left_mirror = False
make_right_mirror = True

# push the image to one side or the other
move_size = -137
move_image_on = True

# flip horizontal and vertical axes
flip_vertical = False

##############################################################################
# Program begins here
# Note: I'm using marshal instead of the more standard pickle to move my files 
# around because pickle appears to be incompatible with the current release of 
# PyPy.
import time
import marshal


def raw_name(file_in):
    """
    Takes an input filename and cuts off the extension so that there is a 
    string that can be used for various input and output files.
    """
    file_cutoff = file_in.find(".")
    raw_name = file_in[:file_cutoff]
    return raw_name


def make_copy(img_data):
    """
    This is the equivalent of the copy.deepcopy function for the MRI images, 
    but much faster because it's tailored to this data type. Using deepcopy 
    instead costs an extra ~5 seconds/copy on my machine with a data set 
    this large.
    """
    img_out = ([[[int(e) for e in sag_row] for sag_row in ax_row] 
                for ax_row in img_data])
    return img_out


def find_params(img_data):
    return (len(img_data), len(img_data[0]), len(img_data[0][0]))

def expand_coronal(img_data):

    return img_data


def move_image(coords, move_size, img_data):
    img_copy = make_copy(img_data)
    for ax_ind, ax_row in enumerate(img_data):
        for sag_ind, sag_row in enumerate(ax_row):
            for cor_ind, voxel in enumerate(sag_row):
                try:
                    sag_row[cor_ind] = img_copy[ax_ind][sag_ind - move_size][cor_ind]
                except:
                    pass
    return img_data


def move_image_90(coords, move_size, img_data):
    img_copy = make_copy(img_data)
    for sag_ind, sag_row in enumerate(img_data):
        for ax_ind, ax_row in enumerate(sag_row):
            for cor_ind, voxel in enumerate(ax_row):
                try:
                    ax_row[cor_ind] = img_copy[sag_ind - move_size][ax_ind][cor_ind]
                except:
                    pass
    return img_data


def mirror_data(mirror_line, coords, img_data):
    mirror_point = (coords[1] / 2) + mirror_line
    for ax_ind, ax_row in enumerate(img_data):
        for sag_ind, sag_row in enumerate(ax_row):
            for cor_ind, voxel in enumerate(sag_row):
                if make_left_mirror and sag_ind < mirror_point:
                    try:
                        sag_row[cor_ind] = img_data[ax_ind][mirror_point + (mirror_point - sag_ind)][cor_ind]
                    except:
                        pass
                elif make_right_mirror and sag_ind > mirror_point:
                    try:
                        sag_row[cor_ind] = img_data[ax_ind][mirror_point - (sag_ind - mirror_point)][cor_ind]
                    except:
                        pass
                elif check_orient and sag_ind > mirror_point:
                    try:        
                        sag_row[cor_ind] = 0
                    except:
                        pass  
    return img_data


def mirror_data_90(mirror_line, coords, img_data):
    mirror_point = (coords[1] / 2) + mirror_line
    for sag_ind, sag_row in enumerate(img_data):
        for ax_ind, ax_row in enumerate(sag_row):
            for cor_ind, voxel in enumerate(ax_row):
                if make_left_mirror and sag_ind < mirror_point:
                    try:
                        ax_row[cor_ind] = img_data[mirror_point + (mirror_point - sag_ind)][ax_ind][cor_ind]
                    except:
                        pass
                elif make_right_mirror and sag_ind > mirror_point:
                    try:
                        ax_row[cor_ind] = img_data[mirror_point - (sag_ind - mirror_point)][ax_ind][cor_ind]
                    except:
                        pass
                elif check_orient and sag_ind > mirror_point:
                    try:        
                        ax_row[cor_ind] = 0
                    except:
                        pass   
    return img_data

print time.clock()
print "Loading data."
raw_name = raw_name(file_in)
datafile = open((directory + raw_name + raw_data_suffix), 'rb')
img_data = marshal.load(datafile)
zoomfile = open((directory + raw_name + zoom_suffix), 'rb')
img_zoom = marshal.load(zoomfile)

coords = find_params(img_data)
print coords

if not flip_vertical:
    if move_image_on:
        img_data = move_image(coords, move_size, img_data)
    img_data = mirror_data(mirror_line, coords, img_data)
if flip_vertical: 
    if move_image_on:
        img_data = move_image_90(coords, move_size, img_data)
    img_data = mirror_data_90(mirror_line, coords, img_data)

print "Saving raw data, use NiftiSave.py to convert to .nii."
datafile = open((directory + raw_name + proc_data_suffix), 'w+b')
marshal.dump(img_data, datafile)
print "Done, final time " + str(time.clock()) + " seconds"
