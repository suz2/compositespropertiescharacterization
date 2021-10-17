from abaqus import *
from abaqusConstants import *

model = 'test'
part = 'composite'

mdb.Model(name=model, modelType=STANDARD_EXPLICIT)
modelObject = mdb.models[model]

s1 = modelObject.ConstrainedSketch(name='__profile__', sheetSize=500.0) 
s1.setPrimaryObject(option=STANDALONE)
s1.rectangle(point1=(0.0,0.0), point2=(176.0,176.0))

modelObject.Part(name=part, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
part1 = modelObject.parts[part]
part1.BaseShell(sketch=s1)

f=mdb.models[model].parts[part].faces[0]

s2 = modelObject.ConstrainedSketch(name='__profile__', sheetSize=500.0) 
s2.setPrimaryObject(option=STANDALONE)
