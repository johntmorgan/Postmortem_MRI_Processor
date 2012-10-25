directory = r"C:/NiftiSubjects/B-5807_B_01/"
file_in = "B-5807_B_01.nii"
file_out = "B-5807_B_01-bright50randbg.nii"

# add a temp_file name if you *don't* want to use the input name plus 
# "-tmpdata" and "-tmpaffine"
# you will need to edit the save and processing files as well
proc_data_suffix = "-tmpprocdata"
affine_suffix= "-tmpaffine"
zoom_suffix = "-tmpzoom"


import marshal
import nibabel as nib
import numpy as np


def save_file(directory, file_out, img_data, img_affine):
    save_image = nib.Nifti1Image(img_data, img_affine)
    path = directory + file_out
    nib.save(save_image, path)


def raw_name(file_in):
    file_cutoff = file_in.find(".")
    raw_name = file_in[:file_cutoff]
    return raw_name


raw_in = raw_name(file_in)
raw_out= raw_name(file_out)

print "loading raw data"
datafile = open((directory + raw_in + proc_data_suffix), 'rb')
affinefile = open((directory + raw_in + affine_suffix), 'rb')

img_data = marshal.load(datafile)
img_affine = marshal.load(affinefile)

print "converting lists into arrays"
img_data = np.array(img_data)
img_affine = np.array(img_affine)

print "saving " + file_out
save_file(directory, file_out, img_data, img_affine)

# save mask data if you are testing the mask settings
#print "doing mask stuff"
#maskfile = open((directory + raw_in + "-imgmask"), 'rb')
#img_mask = marshal.load(maskfile)
#img_mask = np.array(img_mask)
#mask_name = raw_out + "-mask.nii"
#save_file(directory, mask_name, img_mask, img_affine)

print "done"