from PIL import Image,ImageFont,ImageDraw
import colorsys
import numpy as np

def draw_one_box(img,bbox,id,cat,color):
    '''
    img: PIL format
    bbox: (top,left,bottom,right)
    id: identity of the object
    cat: catagory,string
    color: (r,g,b)
    '''
    thickness = (img.size[0] + img.size[1]) // 600
    left,top,right,bottom = bbox
    top = max(0,np.floor(top+0.5).astype('int32'))
    left = max(0,np.floor(left+0.5).astype('int32'))
    bottom = min(img.size[1], np.floor(bottom + 0.5).astype('int32'))
    right = min(img.size[0], np.floor(right + 0.5).astype('int32'))

    label = 'No.{}:{}'.format(id,cat)
    draw = ImageDraw.Draw(img)
    label_size = draw.textsize(label)

    if top-label_size[1]>=0:
        text_origin = np.array([left, top - label_size[1]])
    else:
        text_origin = np.array([left, top + 1])

    for j in range(thickness):
        draw.rectangle([left + j, top + j, right - j, bottom - j],
                        outline=color)
    draw.rectangle([tuple(text_origin), tuple(text_origin + label_size)],
          fill=tuple(color))
    draw.text(list(text_origin),label,fill=(0,0,0))

    return img
