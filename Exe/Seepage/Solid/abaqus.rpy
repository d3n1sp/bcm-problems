# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.11-1 replay file
# Internal Version: 2011_04_21-13.08.34 111392
# Run by Dima on Fri Apr 12 23:11:27 2013
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
#: Abaqus Error: 
#: This error may have occurred due to a change to the Abaqus Scripting
#: Interface. Please see the Abaqus Scripting Manual for the details of
#: these changes. Also see the "Example environment files" section of 
#: the Abaqus Site Guide for up-to-date examples of common tasks in the
#: environment file.
#: Execution of "onCaeGraphicsStartup()" in the site directory failed.
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=179.739581644535, 
    height=132.73149061203)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('Seepage.cae')
#: The model database "c:\Users\Dima\Programs\Bcm_all\Exe\Seepage\Solid\Seepage.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Filtration'].parts['Filtr_R30_octa']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.14439, 
    farPlane=3.70622, width=1.21504, height=1.00628, cameraPosition=(-1.19548, 
    0.726579, 2.48376), cameraUpVector=(-0.113568, 0.730046, -0.673896), 
    cameraTarget=(-0.0832165, -0.133566, -0.083217))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.16351, 
    farPlane=3.63694, width=1.22588, height=1.01525, cameraPosition=(-0.483826, 
    2.52183, 1.07943), cameraUpVector=(0.069214, 0.0781725, -0.994534), 
    cameraTarget=(-0.0835882, -0.134504, -0.0824835))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.11743, 
    farPlane=3.73451, width=1.19977, height=0.993625, cameraPosition=(1.83399, 
    0.608879, 1.97818), cameraUpVector=(-0.242589, 0.812935, -0.529422), 
    cameraTarget=(-0.104852, -0.116955, -0.0907286))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.12021, 
    farPlane=3.75341, width=1.20135, height=0.99493, cameraPosition=(1.48582, 
    -1.08618, 2.16666), cameraUpVector=(0.305885, 0.917704, -0.253483), 
    cameraTarget=(-0.10475, -0.116458, -0.0907839))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.11172, 
    farPlane=3.73898, width=1.19654, height=0.990946, cameraPosition=(0.730719, 
    2.00876, 1.74956), cameraUpVector=(-0.139381, 0.387534, -0.911258), 
    cameraTarget=(-0.107316, -0.105939, -0.0922015))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.14789, 
    farPlane=3.7151, width=1.21704, height=1.00793, cameraPosition=(-0.742719, 
    1.42519, 2.31961), cameraUpVector=(-0.187008, 0.558435, -0.808194), 
    cameraTarget=(-0.106571, -0.105644, -0.0924897))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.13553, 
    farPlane=3.72355, width=1.21004, height=1.00213, cameraPosition=(1.48591, 
    0.579791, 2.26751), cameraUpVector=(-0.529272, 0.797657, -0.289161), 
    cameraTarget=(-0.103028, -0.106988, -0.0925725))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.11775, 
    farPlane=3.72753, width=1.19997, height=0.993787, cameraPosition=(2.01183, 
    1.78772, 0.620122), cameraUpVector=(-0.869686, 0.486763, -0.0819019), 
    cameraTarget=(-0.102543, -0.105874, -0.0940922))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.06885, 
    farPlane=3.76507, width=1.17226, height=0.970842, cameraPosition=(1.41457, 
    1.78309, -1.73363), cameraUpVector=(-0.781983, 0.468639, 0.41095), 
    cameraTarget=(-0.101685, -0.105867, -0.0907122))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.12583, 
    farPlane=3.71456, width=1.20455, height=0.99758, cameraPosition=(1.49755, 
    2.22685, 0.6491), cameraUpVector=(-0.63987, 0.253006, -0.72564), 
    cameraTarget=(-0.101966, -0.107369, -0.098779))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.11585, 
    farPlane=3.72735, width=1.1989, height=0.992899, cameraPosition=(0.666716, 
    1.64283, 2.11595), cameraUpVector=(-0.447826, 0.541186, -0.711737), 
    cameraTarget=(-0.100077, -0.106041, -0.102113))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.13455, 
    farPlane=3.7064, width=1.2095, height=1.00167, cameraPosition=(0.549784, 
    2.29214, 1.44422), cameraUpVector=(-0.376619, 0.266671, -0.887156), 
    cameraTarget=(-0.0998677, -0.107203, -0.100911))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.16814, 
    farPlane=3.67514, width=1.22853, height=1.01743, cameraPosition=(0.235158, 
    1.45071, 2.35336), cameraUpVector=(-0.19105, 0.617937, -0.762663), 
    cameraTarget=(-0.0991832, -0.105372, -0.102889))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.23105, 
    farPlane=3.61506, width=1.26418, height=1.04695, cameraPosition=(0.589689, 
    0.265089, 2.71697), cameraUpVector=(-0.175795, 0.89069, -0.419246), 
    cameraTarget=(-0.099813, -0.103266, -0.103535))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.09809, 
    farPlane=3.74591, width=1.18884, height=0.984558, cameraPosition=(0.868604, 
    1.45914, 2.17354), cameraUpVector=(-0.252008, 0.617194, -0.745361), 
    cameraTarget=(-0.100173, -0.104809, -0.102833))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.09381, 
    farPlane=3.75051, width=1.18642, height=0.982551, cameraPosition=(1.12111, 
    1.38358, 2.10157), cameraUpVector=(-0.303579, 0.639629, -0.706197), 
    cameraTarget=(-0.10059, -0.104684, -0.102714))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.10116, 
    farPlane=3.74449, width=1.19059, height=0.986001, cameraPosition=(2.07276, 
    0.917886, 1.56911), cameraUpVector=(-0.496128, 0.766902, -0.407086), 
    cameraTarget=(-0.102112, -0.103939, -0.101863))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.08686, 
    farPlane=3.75791, width=1.18249, height=0.97929, cameraPosition=(2.03225, 
    1.35779, 1.267), cameraUpVector=(-0.60364, 0.648217, -0.464148), 
    cameraTarget=(-0.102056, -0.104542, -0.101449))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.08574, 
    farPlane=3.75922, width=1.18186, height=0.978767, cameraPosition=(1.98934, 
    1.25334, 1.43103), cameraUpVector=(-0.573847, 0.680214, -0.45608), 
    cameraTarget=(-0.101991, -0.104383, -0.101698))
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.10995, 
    farPlane=3.7349, width=1.19558, height=0.99013, cameraPosition=(2.22919, 
    1.32152, 0.946296), cameraUpVector=(-0.657408, 0.659036, -0.365357), 
    cameraTarget=(-0.102348, -0.104484, -0.100977))
p1 = mdb.models['Seepage'].parts['Box']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p1 = mdb.models['Seepage'].parts['Filtr_R65_full']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.17494, 
    farPlane=10.6543, width=1.77144, height=1.46709, viewOffsetX=0.0247739, 
    viewOffsetY=0.0165511)
p1 = mdb.models['Seepage'].parts['Filtr_R65_half']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.45686, 
    farPlane=9.19353, width=1.3303, height=1.10173, viewOffsetX=0.239434, 
    viewOffsetY=0.0815698)
p1 = mdb.models['Seepage'].parts['Filtr_R65_octa']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p1 = mdb.models['Seepage'].parts['Filtr_R65_quarter']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=4.69015, 
    farPlane=7.95589, width=1.75749, height=1.45553, viewOffsetX=0.0141683, 
    viewOffsetY=0.00564077)
p1 = mdb.models['Seepage'].parts['Filtr_R65_octa']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.47114, 
    farPlane=6.04106, width=1.28064, height=1.06061, viewOffsetX=-0.0217646, 
    viewOffsetY=0.0462363)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.61119, 
    farPlane=5.87083, width=1.33231, height=1.1034, cameraPosition=(3.12746, 
    2.53142, -0.135218), cameraUpVector=(-0.803206, 0.465241, -0.372036), 
    cameraTarget=(-0.473956, -0.555679, -0.48177), viewOffsetX=-0.0226427, 
    viewOffsetY=0.0481017)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.46875, 
    farPlane=6.03009, width=1.27976, height=1.05988, cameraPosition=(2.91325, 
    2.02797, 1.62599), cameraUpVector=(-0.663232, 0.610548, -0.432844), 
    cameraTarget=(-0.47436, -0.554828, -0.489132), viewOffsetX=-0.0217496, 
    viewOffsetY=0.0462044)
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.36748, 
    farPlane=6.11943, width=1.86511, height=1.54466, cameraPosition=(2.63502, 
    2.12792, 1.90238), cameraUpVector=(-0.666207, 0.587698, -0.459107))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.39338, 
    farPlane=6.10384, width=1.87945, height=1.55654, cameraPosition=(2.89747, 
    1.90529, 1.78614), cameraUpVector=(-0.660438, 0.632241, -0.405084), 
    cameraTarget=(-0.456182, -0.55874, -0.498477))
p1 = mdb.models['Seepage'].parts['Filtr_R65_half']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.24356, 
    farPlane=9.44633, width=3.03978, height=2.5175, cameraPosition=(5.02493, 
    3.37701, -3.65941), cameraUpVector=(-0.919247, 0.258886, -0.296586), 
    cameraTarget=(0.0621104, -0.0555703, 0.49346))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.65018, 
    farPlane=9.10803, width=3.27551, height=2.71273, cameraPosition=(5.82304, 
    0.0793298, -4.03239), cameraUpVector=(-0.830037, 0.226686, -0.50956), 
    cameraTarget=(0.0642562, -0.0644366, 0.492457))
session.viewports['Viewport: 1'].view.setValues(nearPlane=6.07107, 
    farPlane=8.63322, width=3.51951, height=2.9148, cameraPosition=(0.887017, 
    0.384421, -6.78863), cameraUpVector=(-0.945389, 0.222539, 0.238151), 
    cameraTarget=(0.0281958, -0.0622077, 0.472321))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.59021, 
    farPlane=9.14672, width=3.24074, height=2.68393, cameraPosition=(5.26184, 
    0.188325, -4.65522), cameraUpVector=(-0.89496, 0.127336, -0.427589), 
    cameraTarget=(0.0442331, -0.0629265, 0.480142))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.2157, 
    farPlane=9.48838, width=3.02363, height=2.50412, cameraPosition=(5.82736, 
    3.48276, 3.32294), cameraUpVector=(-0.0672189, 0.193768, -0.978742), 
    cameraTarget=(0.0475542, -0.0435806, 0.526992))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.35562, 
    farPlane=9.32922, width=3.10474, height=2.5713, cameraPosition=(5.24611, 
    2.25983, -4.11375), cameraUpVector=(-0.830995, -0.21283, -0.513956), 
    cameraTarget=(0.0454323, -0.0480451, 0.499843))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.54588, 
    farPlane=9.18648, width=3.21504, height=2.66265, cameraPosition=(5.55423, 
    -0.841805, -4.26509), cameraUpVector=(-0.871892, -0.0513025, -0.487003), 
    cameraTarget=(0.0461549, -0.0553188, 0.499488))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.37735, 
    farPlane=9.39479, width=3.11734, height=2.58174, cameraPosition=(6.60919, 
    -2.2108, -1.94706), cameraUpVector=(-0.653945, 0.00666184, -0.756513), 
    cameraTarget=(0.0520228, -0.0629334, 0.512381))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.53191, 
    farPlane=9.20549, width=3.20694, height=2.65594, cameraPosition=(5.65721, 
    -0.90494, -4.13439), cameraUpVector=(-0.835988, 0.178076, -0.51905), 
    cameraTarget=(0.044178, -0.0521724, 0.494356))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.5085, 
    farPlane=9.22333, width=3.19337, height=2.6447, cameraPosition=(6.93599, 
    0.901966, -1.81043), cameraUpVector=(-0.622882, 0.0365953, -0.78146), 
    cameraTarget=(0.0517263, -0.0415068, 0.508074))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.61261, 
    farPlane=9.11506, width=3.25372, height=2.69468, cameraPosition=(5.98163, 
    0.0660695, -3.79471), cameraUpVector=(-0.822437, 0.0415009, -0.56734), 
    cameraTarget=(0.0464523, -0.0461262, 0.497108))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.84538, 
    farPlane=8.88228, width=1.51597, height=1.25551, viewOffsetX=-0.00945469, 
    viewOffsetY=0.194019)
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.4883, 
    farPlane=8.8513, width=1.42336, height=1.17881, cameraPosition=(6.3093, 
    1.06193, 3.73912), cameraUpVector=(0.0825718, 0.141685, -0.986462), 
    cameraTarget=(-0.188925, -0.0784199, 0.556156), viewOffsetX=-0.00887713, 
    viewOffsetY=0.182167)
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.71331, 
    farPlane=8.75089, width=1.48172, height=1.22714, cameraPosition=(7.02405, 
    -0.387174, -1.18614), cameraUpVector=(-0.551749, 0.170577, -0.81638), 
    cameraTarget=(-0.0513469, -0.0371986, 0.677954), viewOffsetX=-0.00924108, 
    viewOffsetY=0.189636)
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.52125, 
    farPlane=8.9646, width=1.43191, height=1.18589, cameraPosition=(6.51742, 
    -1.53342, -2.26749), cameraUpVector=(-0.664266, 0.167187, -0.728559), 
    cameraTarget=(-0.0111084, -0.00761902, 0.683588), viewOffsetX=-0.00893043, 
    viewOffsetY=0.183261)
