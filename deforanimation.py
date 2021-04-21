###############################################################################
# deforanimation.py is created for creating an animation for dynamic simulation    #
# h5test.py should be run at first.                                           #
# Created by Zheng Zhou 'zhougeohhu@foxmail.com'                              #
# Please follow the explanation and edit it for your simulation result        #
###############################################################################

import os,re,random
import cv2
# Python libraries above should be checked bedore your beginning
# Python 3.X is suggested
# If you use Spyder, you could have an interface similar to MatLAB

########################### deformation animation #############################
picPath = '/root/msys_sandbox/BM/Deformation/pics/'
def get_images(path):
    file_list = []
    for root,dirs,files in os.walk(path):
        if not files:
            continue
        for file in files:
            if file.endswith('.png'):
                file_list.append(os.path.join(root,file))
    return file_list

def key_sort(image_path):
    pattern = re.compile("\d+")
    image_name = os.path.basename(image_path)
    return int(pattern.findall(image_name)[0])

# Attention please, pictures are suggested to name by number for 'VideoWriter'
# function could read pictures files in the actually orders.
path        =   picPath
file_list   =   get_images(path)
file_list.sort(key=key_sort)
# video basic setting
fps             =   24

_files = os.listdir(path)
number = random.randint(1, len(_files))

for number in file_list:
    simg = cv2.imread(number)
    sp = simg.shape
    # shape [height, width, color]
img_size        =   (sp[1],sp[0])
# be careful, (width, height)
# folder and name for video
save_path       =   '/root/msys_sandbox/BM/Deformation/video/test.mp4'
fourcc          =   cv2.VideoWriter_fourcc('M','J','P','G')
video_writer    =   cv2.VideoWriter(save_path,fourcc,fps,img_size)
for file_name in file_list:
    print (file_name)
    img = cv2.imread(file_name)
    video_writer.write(img)

video_writer.release()
print ('finish')
