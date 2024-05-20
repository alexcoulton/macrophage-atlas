### Build 'Blank Map' and overlay on feature expression plots
### Last modified: 14 Jun 2022 by Jun Murai
### Created: 10 Jun 2022 by Jun Murai

#! python3

import sys
import os
from PIL import Image, ImageDraw, ImageFilter, ImageEnhance, ImageOps

## function processImage: build a 'blank map' outline from Cluster DimPlot
# input file:  [basename].png    : DimPlot with cluster labels
#              [basename]_NL.png : DimPlot without cluster labels
# output:      [basename]_NL_outline.png : The 'Blank Map', outlined cluster maps
#              [basename]_lined.png      : DimPlot overlaid with outlined cluster maps
# intermediate files:   proc/[basename]_NL_MF25.png         : processed image with mode filter (25) 
#                       proc/[basename]_NL_CONTOUR.png      : processed image with mode filter (25) and contour filter
#                       proc/[basename]_NL_CONTOUR_Gray.png : processed image with mode filter (25), contour filter, grayscale conversion and contrast enhancement
def processImage(basename):
    os.makedirs("proc", exist_ok=True)
    nolabel = basename+"_NL"
    print ("input files: " + basename + ".png " + nolabel + ".png ")
    im = Image.open(nolabel+".png")
    im2 = im.filter(filter=ImageFilter.ModeFilter(size=25))
    im2.save("proc/"+nolabel+"_MF25.png", "PNG")
    im2 = im2.filter(filter=ImageFilter.CONTOUR)
    im2.save("proc/"+nolabel+"_MF_CONTOUR.png", "PNG")
    im3 = im2.convert("L")
    enhancer = ImageEnhance.Contrast(im3)  
    im3 = enhancer.enhance(20)
    im3.save("proc/"+nolabel+"_MF_CONTOUR_Gray.png", "PNG")
    im4 = im3.convert("1")
    im4r = im4.convert("RGB")
    im4g = im4.convert("L")
    im4g = ImageOps.invert(im4g)
    im4r.putalpha(im4g)
    print ("result files: " + nolabel+"_outline.png " + basename+"_lined.png")
    im4r.save(nolabel+"_outline.png", "PNG")
    im_l = Image.open(basename+".png")
    im_l.paste(im4r, (0, 0), im4r)
    im_l.save(basename+"_lined.png", "PNG")

## example: 00_macrophage_20.png, 00_macrophage_20_NL.png etc. will be processed 
processImage("00_macrophage_24")
processImage("00_macrophage_20")
processImage("00_macrophage_27")
processImage("00_macrophage_49")
processImage("00_macrophage_80")

## end of script
