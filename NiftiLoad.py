import marshal
import nibabel as nib
import numpy as np

directory = r"C:/NiftiSubjects/BTB-3714_B_01/"
file_in = "BTB-3714_B_01.nii"

#add a temp_file name if you *don't* want to use the input name plus "-tmpdata" and
#"-tmpaffine"
#you will need to edit the save and processing files as well
raw_data_suffix = "-tmprawdata"
affine_suffix = "-tmpaffine"
zoom_suffix = "-tmpzoom"

def load_file(directory, file_in):
    path = directory + file_in
    img = nib.load(path)
    img_header = img.get_header()
    return img.get_data(), img.get_affine(), img_header.get_zooms()


def raw_name(file_in):
    file_cutoff = file_in.find(".")
    raw_name = file_in[:file_cutoff]
    return raw_name

print "reading files"
img_data, img_affine, img_zoom = load_file(directory, file_in)
raw_name = raw_name(file_in)

print "converting arrays and tuples to lists"
img_data = img_data.tolist()
img_affine = img_affine.tolist()
img_zoom = np.array(img_zoom)
img_zoom = img_zoom.tolist()
                
print "storing data using marshal format"
datafile = open((directory + raw_name + raw_data_suffix), 'w+b')
marshal.dump(img_data, datafile)

affinefile = open((directory + raw_name + affine_suffix), 'w+b')
marshal.dump(img_affine, affinefile)

zoomfile = open((directory + raw_name + zoom_suffix), 'w+b')
marshal.dump(img_zoom, zoomfile)

print "done"