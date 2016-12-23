import numpy as np
import os
import sys
import array

filename = str(sys.argv[1])
gsx = int(sys.argv[2])
gsy = int(sys.argv[3])
gsz = int(sys.argv[4])

lsx = int(sys.argv[5])
lsy = int(sys.argv[6])
lsz = int(sys.argv[7])

input_file = open(filename,'r')
parent_folder = os.path.dirname(input_file.name)
A = np.fromfile(input_file, dtype='uint8')
B = np.reshape(A, (gsx, gsy, gsz))

S = 'parallel ::: '

count = 1

for z in range(0,(gsz-lsz+1),lsz):
    for y in range(0,(gsy-lsy+1),lsy):
        for x in range(0,(gsx-lsx+1),lsx):
            if (x==(gsx-lsx)):
                cx = 0
            else:
                cx = 1
            if (y==(gsy-lsy)):
                cy = 0
            else:
                cy = 1
            if (z==(gsz-lsz)):
                cz = 0
            else:
                cz = 1
            stc = str(count)
            stx = str(lsx+cx)
            sty = str(lsy+cy)
            stz = str(lsz+cz)
            filestr = parent_folder + os.sep + 'fo_' + stc + '_' + stx + '_' + sty + '_' + stz + '.raw'
            S = S + '"./serialct '+ filestr + ' ' + stx + ' ' + sty + ' ' + stz + ' f' + stc + ' ' + stc + '"' + ' '
            C = B[x: x+lsx+cx, y: y+lsy+cy, z: z+lsz+cz]
            C.tofile(filestr)
            count += 1
with open(parent_folder + os.sep +'tree.sh', "w") as tree:
    tree.write(S)
with open(parent_folder + os.sep +'params.txt', "w") as params:
    params.write(str(gsx/lsx) + ' ' + str(gsy/lsy) + ' ' + str(gsz/lsz) + ' '+ str(lsx) + ' ' + str(lsy) + ' ' + str(lsz))
input_file.close()
