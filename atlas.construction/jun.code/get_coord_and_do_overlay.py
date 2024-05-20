import sys
import os
import glob
import cv2
from IPython.display import Image, display
import numpy as np
import copy
import json
from PIL import Image, ImageDraw, ImageFilter, ImageEnhance, ImageOps

do_overlay_offset = 231   # directive to do overlay with offset x  (0: don't)

if do_overlay_offset > 0:
    im_ol = Image.open("../outline.png")
    os.makedirs("res", exist_ok=True)

for tgtfile in glob.glob('*.png'):
    gray_img = cv2.imread(tgtfile, cv2.IMREAD_GRAYSCALE)
    contours, hierarchy = cv2.findContours(gray_img, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    
    for i in range(len(contours)):
        
        ca = cv2.contourArea(contours[i])
        
        if ca < 10000:
            continue  # remove small items. expected: 17698.5 
        
        if ca > 30000:
            continue  # remove large items, backgrounds, main plot
        
        x, y, w, h = cv2.boundingRect(contours[i])
        if x > 1000 or y > 1000 or w < 1000 or h < 1000:
            continue
        
        print (x, y, w, h, ca, tgtfile, sep="\t")
        
        if w > 1926 and w < 1930 and h == 1935 and do_overlay_offset > 0:
            offset = x - do_overlay_offset
            im = Image.open(tgtfile)
            im.paste(im_ol, (offset, 0), im_ol)
            print ("offset found, writing to res/"+tgtfile)
            im.save("res/"+tgtfile, "PNG")

