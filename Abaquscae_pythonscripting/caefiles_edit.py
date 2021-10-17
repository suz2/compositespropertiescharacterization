import os
import sys
import numpy as np

os.system('cp composites_cae.py composites_cae_forrun.py')
node_resol = 0.01 # node resolution to get nodes
f = open('fibercentroid_r5vf30domain200/fibercentroid_info.dat')
lines = f.readlines()
size = float(lines[0].strip('\n').lstrip().split(' ')[-1])
print(size)

radius = float(lines[1].lstrip().strip('\n').split(' ')[-1])
print(radius)

fibernum = int(lines[3].strip('\n').lstrip().split(' ')[-1])
print(fibernum)

fibercoord = np.zeros((fibernum,2))

row = 0
for line in lines[4:]:
    list = line.strip('\n').split(',')
    fibercoord[row,:] = list[::]
    row+=1

print(fibercoord)


f=open('composites_cae_forrun.py','r+')
flist=f.readlines()
flist[11] = 's.rectangle(point1=(0.0,0.0), point2=(' + str(size) + ',' + str(size) + '))'

for i in range(fibernum):
    circle_centroid_x = fibercoord[i,0]
    circle_centroid_y = fibercoord[i,1]
    circlep1_centroid_x = fibercoord[i,0] + radius
    circlep1_centroid_y = fibercoord[i,1]

    flist.append( 's2.CircleByCenterPerimeter(center=(' + str(circle_centroid_x-100.0) + ', ' + str(circle_centroid_y-100.0) + \
                    '), point1=(' + str(circlep1_centroid_x-100.0) + ', ' + str(circlep1_centroid_y-100.0) + '))\n')
#    flist.append( 's2.Line(point1=(' + str(circle_centroid_x) + ', ' + str(circle_centroid_y) + \
#                    '), point2=(' + str(circlep1_centroid_x) + ', ' + str(circlep1_centroid_y) + '))\n')


middle_line = '''part1.PartitionFaceBySketch(faces=f, sketch=s2)
s2.unsetPrimaryObject()
del modelObject.sketches['__profile__']
###create mesh in gui

###create nodes set
allNodes = mdb.models[model].parts[part].nodes
nodesset = []\n'''

flist.append(middle_line)


for i in range(fibernum):
    circle_centroid_x = fibercoord[i,0]
    circle_centroid_y = fibercoord[i,1]

    if (circle_centroid_x > 0) and (circle_centroid_x < size) and (circle_centroid_y > 0) and (circle_centroid_y < size):
        circle_centroid_x_min = circle_centroid_x - node_resol
        circle_centroid_x_max = circle_centroid_x + node_resol
        circle_centroid_y_min = circle_centroid_y - node_resol
        circle_centroid_y_max = circle_centroid_y + node_resol

        flist.append('nodesset.append(allNodes.getByBoundingBox(' + \
                    str(circle_centroid_x_min) + ' ,' + str(circle_centroid_y_min) + ', 0.0,' + str(circle_centroid_x_max) + ' ,' + str(circle_centroid_y_max) + ', 0.0))\n')
    
flist.append( 'mdb.models[model].parts[part].Set(name=\'fibercentroid\', nodes=nodesset)\n')

f=open('composites_cae_forrun.py','w+')
f.writelines(flist)
f.close()
