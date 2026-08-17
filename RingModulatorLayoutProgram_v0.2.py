import tkinter as tk
import tkinter.filedialog as tkf
#from tkinter import ttk
import ttkbootstrap as ttk         # shell: pip install ttkbootstrap
from ctypes import windll
import math

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import random
#import kfactory as kf
import gdsfactory as gf

import time

import numpy as np
from PIL import ImageTk, Image
#import shapely
import os
import pathlib


windll.shcore.SetProcessDpiAwareness(2)


ProgramVersionString = "v0.2"
#ProgramVersionDateString = "2025-12-16"    # start of programming
#ProgramVersionDateString = "2026-02-18"    # First full geometry version
#ProgramVersionDateString = "2026-02-20"     # First version saving full GDS file
#ProgramVersionDateString = "2026-03-13"     # Demo loader program added, fixes
ProgramVersionDateString = "2026-08-17"     # GIT version

GDS_DBU_ini = 1e-09   # GDS database units (1nm)
GDS_UU_ini = 1e-06    # GDS user units (1µm)

#scale=GDS_DBU_ini/GDS_UU_ini

CladdingWidth_ini = 2000   # cladding width around structures in nm
RingDiscretization_ini = 256   # number of vertices in full circle of ring
MinAngle_ini = 85  # minimum 'steep' angle for structure in degrees
BevelLength_ini = 300  # length of bevels to mitigate sharp angels (nm)

RingDiameter_ini = 15000  # diameter of ring center in nm
RingWGWidth_ini = 400     # width of ring waveguide in nm
BusGap_ini = 180          # gap between ring and bus waveguide in nm
DropGap_ini = 180         # gap between ring and drop waveguide in nm
BusWGWidth_ini = 400      # width of bus waveguide in nm
DropWGWidth_ini = 400     # width of drop waveguide in nm

OSlabWidth_ini = 1250       # width of optical slabs (thin etched regions) around rib waveguides in nm
BusWGLength_ini = 45000     # length of bus waveguide
DropWGLength_ini = 45000     # length of drop waveguide
TaperLength_ini = 10000     # length of tapers
SlabBorderExtensionOnTapers_ini=250   # tapers of optical slabs get a border of this size (required by some design rules)
PortWGWidth_ini = 450       # width of port waveguides

ESlabWidth_ini = 2250       # width of electrical contact slabs (non-etched regions) around optical slabs in nm
OuterESlabLength_ini = 18000    # length of outer electrical slabs

JunctionOffset_ini = 0      # Offset position of circular pn-junction, positive values shift the junction line to the outside, negative to the inside in nm
HighDopingDistanceInside_ini = 450    # distance between WG border and start of high doping region for ring inside in nm
HighDopingDistanceOutside_ini = 300   # distance between WG border and start of high doping region for ring outside and bus and drop waveguides in nm
ContactPlugRegionBorderWidth_ini = 200    # width of border between contact doping and region into which contact plugs to first metal layer are placed in nm
ContactPlugSize_ini = 350             # width and height of silicon contact plug square in nm
ContactPlugPitch_ini = 600             # horizontal and vertical pitch (center to center distance) between silicon contact plugs in nm
ContactPlugHex_ini = 1                  # Use hexagonal contact plug pattern
ContactPlugPatternYOffset_ini = 0       # offset contact plug pattern in y-direction in nm


ElectricalLeadProtrusion_ini = 20000    # distance between center and outer edge of metal leads to contacts
ElectricalLeadWidth_ini = 3000         # width of RF metal leads
ElectricalLeadGap_ini = 1100           # gap between inner electrical lead and metal2-ring of outer electrical contact
InnerElectricalRingCutawayOffset_ini = -600    # horizontal offset of cut of top metal layer for inner ring contact, can be used to make room to add some more vias

MetalViaBorderWidth_ini = 140        # width of border between metalization and region into which vias between metal layers are placed
ViaSize_ini = 600                     # width and height of via square
ViaPitch_ini = 1200                    # horizontal and vertical pitch (center to center distance) between vias
ViaHex_ini = 1                  # Use hexagonal via pattern
ViaPatternYOffset_ini = 0        # offset via pattern in y-direction in nm


HeaterWidth_ini = 500         # width of heater structure
HeaterLegSpacing_ini = 3850    # spacing between the heater lines at the narrow part of the Omega shape
HeaterContactExtension_ini = 3500  # distance between contacting leads of heater and the region of OuterESlabLength
HeaterMetal1Overlap_ini = 6500 # Overlap between heater structure and Metal1 for contacting
HeaterMetal1Protrusion_ini = 3100  # Protrusion of Metal1 over heater structure
HeaterMetal2Protrusion_ini = 6000  # Protrusion of Metal1 over heater structure


ShowCladding_ini = 1
ShowSilicon_ini = 1
ShowWGRib_ini = 1
ShowOSlabs_ini = 1
ShowLowDoping_ini = 1
ShowHighDoping_ini = 1
ShowContactDoping_ini = 1
ShowContactPlugRegion_ini = 0
ShowContactPlugs_ini = 1
ShowMetal1_ini = 1
ShowMetal2_ini = 1
ShowMetalViaRegion_ini = 0
ShowVias_ini = 1
ShowHeater_ini = 1
ShowConstructionLines_ini = 0

layer_Cladd_ini = (1,0)
layer_Si_ini = (2,0)
layer_OuterOSlabCore_ini = (3,0)
layer_OuterOSlabCladd_ini = (4,0)
layer_InnerOSlabCladd_ini = (5,0)
layer_OuterLowDope_ini = (6,0)
layer_InnerLowDope_ini = (7,0)
layer_OuterHiDope_ini = (8,0)
layer_InnerHiDope_ini = (9,0)
layer_OuterContactDope_ini = (10,0)
layer_InnerContactDope_ini = (11,0)
layer_M1_ini = (12,0)
layer_M1PerfBlock_ini = (13,0)
layer_M2_ini = (14,0)
layer_M2PerfBlock_ini = (15,0)
layer_Heater_ini = (16,0)
layer_HeaterContact_ini = (17,0)
layer_Vias_ini = (18,0)
layer_ContactPlugs_ini = (19,0)


InnerElectricalMetal2RingMinimumGap=0.900     # minimum gap in µm in inner top metal ring for easier programming
#SlabBorderExtensionOnTapers=0.0             # taperes optical slabs get a border of this size


# Dictionary to hold all your variables
global_vars = {}
# Dictionary to hold all your layers
global_layers = {}


#fastFig: any
#fastAxes: any

#GraphWindowOpen = False

#ConfigString=""
savepath = ""
configpath = ""


gf.CONF.max_cellname_length = 128
gf.config.rich_output()
gf.clear_cache()
gf.gpdk.PDK.activate()

#maincell = gf.Component("RM")





def CalculateAngle(p1, p2, p3):
    #print(p1, '  ',p2, '  ',p3)
    CAdir1 = [p1[0]-p2[0], p1[1]-p2[1]] 
    CAdir2 = [p3[0]-p2[0], p3[1]-p2[1]]
    #print(CAdir1)
    #print(CAdir2)
    #Formula from https://de.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
    angle = np.degrees(np.arctan2(CAdir1[0]*CAdir2[1]-CAdir1[1]*CAdir2[0],CAdir1[0]*CAdir2[0]+CAdir1[1]*CAdir2[1]))
    #print(angle)
    return angle





def findFirstProblematicAngle(A, MinAngle):
# A must not have the same vertices at start and at end
    #print(A)
    fieldindex=-1
    currentindex=0
    while currentindex<len(A):
        #print(currentindex)
        beforeindex=currentindex-1  # vertex number before current vertex
        afterindex=currentindex+1   # vertex number after current vertex
        if beforeindex<0:           # check if out of Range
            beforeindex=len(A)-1
    
        if afterindex>(len(A)-1):     # check if out of Range
            afterindex=0
        
        #print(A, '   ', beforeindex, '   ', A[beforeindex])
        angle=abs(CalculateAngle(A[beforeindex], A[currentindex], A[afterindex]))
        if angle<MinAngle:
            fieldindex=currentindex
            currentindex=len(A)-1     # stop loop
        
        currentindex=currentindex+1
    
    return fieldindex



def dist(p1, p2) -> float:
    dx=p2[0]-p1[0]
    dy=p2[1]-p1[1]
    return (math.sqrt(dx**2+dy**2))



def LineSegmentCircleIntersection(P1, P2, PC, r):
    # find intersection of line segment between P1 and P2 and circle around PC with radius r
    # one of P1 or P2 must be inside circle

    P12=[P2[0]-P1[0], P2[1]-P1[1]]   #P2-P1;
    P1C=[PC[0]-P1[0], PC[1]-P1[1]]   #PC-P1;
    P2C=[PC[0]-P2[0], PC[1]-P2[1]]   #PC-P2;
    nP1C=dist([0,0], P1C)
    nP2C=dist([0,0], P2C)
    if nP1C==r:
        Pintersect=P1
        return Pintersect
    
    if nP2C==r:
        Pintersect=P2
        return Pintersect

    # sort points
    if ((nP1C>r) and (nP2C<r)):
        P3=P1
        P1=P2
        P2=P3
    
    P12=[P2[0]-P1[0], P2[1]-P1[1]]   #P2-P1;
    P1C=[PC[0]-P1[0], PC[1]-P1[1]]   #PC-P1;
    P2C=[PC[0]-P2[0], PC[1]-P2[1]]   #PC-P2;

    nP1C=dist([0,0], P1C)
    nP2C=dist([0,0], P2C)
    
    if ((nP1C<r) and (nP2C>r)): # okay, we should have a single intersection
        
        P3 = [P1[0]+0.5*P12[0], P1[1]+0.5*P12[1]]   #P3=P1+0.5*P12;
        P3C = [PC[0]-P3[0], PC[1]-P3[1]]   #P3C=PC-P3;
        nP3C=dist([0,0], P3C)

        while abs(nP3C-r)>1e-6:
            if nP3C>r:
                P2=P3
            else:
                P1=P3
            
            P3 = [P1[0]+0.5*(P2[0]-P1[0]), P1[1]+0.5*(P2[1]-P1[1])]   #P3=P1+0.5*(P2-P1);
            P3C= [PC[0]-P3[0], PC[1]-P3[1]]   #PC-P3;
            nP3C=dist([0,0], P3C)

        Pintersect=P3

    
    #Pintersect=[];
    return Pintersect






def BevelCorner(A, currentindex, BevelLength):
# A must not have the same vertices at start and at end
    
    beforeindex=currentindex-1  # vertex number before current vertex
    afterindex=currentindex+1   # vertex number after current vertex
    if beforeindex<0:            # check if out of Range
        beforeindex=len(A)-1
    
    beforesegment=[currentindex, beforeindex]

    if afterindex>(len(A)-1):     # check if out of Range
        afterindex=0
    
    aftersegment=[currentindex, afterindex]

    #angle=abs(CalculateAngle(A(beforeindex,:), A(currentindex,:), A(afterindex,:)));

    bevelfound=0
    searchswitch=0  # =0: undefined, =-1: before-direction was choosen, =+1: after-direction was choosen

    #distances between currentindex point (sharp tip) and adjacent points
    distancebefore=dist(A[currentindex], A[beforesegment[1]])
    distanceafter=dist(A[currentindex], A[aftersegment[1]])

    while (bevelfound<1):

        if distancebefore<distanceafter:
            searchradius=distancebefore
            searchswitch=-1

            PHelp=A[beforesegment[1]]
            #Pintersectafter=LineSegmentCircleIntersection(A(aftersegment(1),:), A(aftersegment(2),:), A(currentindex,:), searchradius);
            Pintersectafter=LineSegmentCircleIntersection(A[aftersegment[0]], A[aftersegment[1]], A[currentindex], searchradius)
            #Pintersectafter
            distancebetween=dist(Pintersectafter, PHelp)
        else:
            searchradius=distanceafter
            searchswitch=1

            PHelp=A[aftersegment[1]]
            #Pintersectbefore=app.LineSegmentCircleIntersection(A(beforesegment(1),:), A(beforesegment(2),:), A(currentindex,:), searchradius);
            Pintersectbefore=LineSegmentCircleIntersection(A[beforesegment[0]], A[beforesegment[1]], A[currentindex], searchradius)
            #Pintersectbefore
            distancebetween=dist(Pintersectbefore, PHelp)
        

        if distancebetween<BevelLength:  # we'll have to check further
            #set new indices and check their distances to currentindex-point
            if searchswitch==-1: #search further in before-direction
                newbeforeindex=beforesegment[1]-1
                if newbeforeindex<0:            # check if out of Range
                    newbeforeindex=len(A)-1
                
                beforesegment[0]=beforesegment[1]
                beforesegment[1]=newbeforeindex
                
                distancebefore=dist(A[currentindex], A[beforesegment[1]])

            else:    #search further in after-direction
                newafterindex=aftersegment[1]+1
                if newafterindex>(len(A)-1):     # check if out of Range
                    newafterindex=0
                
                aftersegment[0]=aftersegment[1]
                aftersegment[1]=newafterindex

                distanceafter=dist(A[currentindex], A[aftersegment[1]])
            
                        
        else:    # yay, the bevel must lie in the before and after-segments
            bevelfound=1
            #Bevel einfügen
            #show segments
            # P1=A(beforesegment(1),:);
            # P2=A(beforesegment(2),:);
            # P3=A(aftersegment(1),:);
            # P4=A(aftersegment(2),:);
            # plot([P1(1) P2(1)], [P1(2) P2(2)], 'LineStyle','-','LineWidth',3,'Color',[1 0 0]);
            # plot([P3(1) P4(1)], [P3(2) P4(2)], 'LineStyle','-','LineWidth',3,'Color',[1 0 0]);

            dbefore=dist(A[currentindex], A[beforesegment[0]])
            dafter=dist(A[currentindex], A[aftersegment[0]])
            
            if (dbefore>dafter):
                minsearchradius=dbefore
            else:
                minsearchradius=dafter
            

            #distancebefore
            #distanceafter

            if (distancebefore<distanceafter):
                maxsearchradius=distancebefore
            else:
                maxsearchradius=distanceafter
            

            # minsearchradius
            # maxsearchradius
            
            stopsearch=0
            while (stopsearch==0):
                #calculate radius in between and refine further
                middlesearchradius=(maxsearchradius+minsearchradius)/2.0
                # calculate segment length
                #viscircles(A(currentindex,:),middlesearchradius);
                # plot(A(aftersegment(1),1),A(aftersegment(1),2),'ro');
                # plot(A(aftersegment(2),1),A(aftersegment(2),2),'ro');
                # plot(A(currentindex,1),A(currentindex,2),'go');

                beforePoint=LineSegmentCircleIntersection(A[beforesegment[0]], A[beforesegment[1]], A[currentindex], middlesearchradius)
                afterPoint =LineSegmentCircleIntersection(A[aftersegment[0]], A[aftersegment[1]], A[currentindex], middlesearchradius)

                #print('Beforepoint: ', beforePoint, '   Afterpoint: ', afterPoint)
    
                # beforePoint
                # afterPoint
                # 
                # A(aftersegment(1),:)
                # A(aftersegment(2),:)
                # A(currentindex,:)
                # middlesearchradius
                # 
                # 
                # plot(beforePoint(1),beforePoint(2),'go');
                # plot(afterPoint(1),afterPoint(2),'go');

            #stopsearch=1;

                middledistancebetween=dist(afterPoint, beforePoint)
                #print(middledistancebetween)

                if (abs(middledistancebetween-BevelLength)<0.0001):
                    stopsearch=1
                    #print(" ")
                    #plot(beforePoint(1),beforePoint(2),'ro');
                    #plot(afterPoint(1),afterPoint(2),'ro');
                    #now insert bevel
                    #beforesegment
                    #aftersegment

                    #check for table wrap in bevel region
                    if (beforesegment[1]>aftersegment[1]):
                        #print('table wrap branch')
                        #special case with table wrap around in between bevel: build new table
                        #NewA=[A(aftersegment(2):beforesegment(2),:) ; [beforePoint(1) beforePoint(2)] ; [afterPoint(1) afterPoint(2)]];
                        #NewA=A[aftersegment[1]:beforesegment[1]+1] + [beforePoint[0], beforePoint[1]] + [afterPoint[0], afterPoint[1]]
                        NewA=A[aftersegment[1]:beforesegment[1]+1]
                        NewA=np.concatenate(( NewA, np.array([beforePoint]) ), axis=0)
                        NewA=np.concatenate(( NewA, np.array([afterPoint]) ), axis=0)

                    else:
                        #NewA=[A(1:beforesegment(2),:) ; [beforePoint(1) beforePoint(2)] ; [afterPoint(1) afterPoint(2)] ; A(aftersegment(2):end,:)];
                        #NewA=A[0:beforesegment[1]+1] + [beforePoint[0], beforePoint[1]] + [afterPoint[0], afterPoint[1]] + A[aftersegment[1]:]
                        #print('standarad branch')
                        NewA=A[0:beforesegment[1]+1]
                        #print(len(NewA))
                        NewA=np.concatenate(( NewA, np.array([[beforePoint[0], beforePoint[1]]]) ), axis=0)
                        #print(len(NewA))
                        NewA=np.concatenate(( NewA, np.array([afterPoint]) ), axis=0)
                        #print(len(NewA))
                        NewA=np.concatenate(( NewA, A[aftersegment[1]:]), axis=0)
                        #print(len(NewA))
                    
                    #print('BevelCorner A:')
                    #print(len(A))
                    #print('BevelCorner NewA:')
                    #print(len(NewA))
                    #NewA=[A(1:beforesegment(2),:) ; [beforePoint(1) beforePoint(2)] ; [afterPoint(1) afterPoint(2)] ; A(aftersegment(2):end,:)];

                else:
                    if (middledistancebetween>BevelLength): #look further
                        maxsearchradius=middlesearchradius
                    else:
                        minsearchradius=middlesearchradius
    
    return NewA





def BevelContour(FullA, minAngle, BevelLength) -> np.ndarray:
    A=FullA
    #A
    
    AngleIndex=findFirstProblematicAngle(A, minAngle)
    #print(AngleIndex)
    while (AngleIndex>-1):
        #print('Angle index: ', AngleIndex)
        #print('Old polygon vertices: ', len(A))
        NewA=BevelCorner(A, AngleIndex, BevelLength)
        #print('New polygon vertices: ', len(NewA))
        A=NewA
        AngleIndex=findFirstProblematicAngle(A, minAngle)
    #end
    #NewFullA=[A;A(1,:)]
    #return NewFullA
    return A



def BevelComponent(c: gf.Component, MinAngle, BevelLength, layer)->gf.Component:
    c2=gf.Component()
    allpolygons=c.get_polygons_points(by='tuple', layers=[layer])
    for key in allpolygons:
        #print(key)
        for polygon in allpolygons[key]:
        #if 1:
            #polygon = test[key][1]
            #print(polygon)
            #fastAxes.plot(polygon[:,0], polygon[:,1], 'r.')
            newPolygon = BevelContour(polygon, MinAngle, BevelLength)
            #fastAxes.plot(newPolygon[:,0], newPolygon[:,1], 'bo')
            #print(newPolygon)
            #print("-----")
            c2.add_polygon(newPolygon, layer=key)
    return c2




def GenerateSquarePattern(SSize, SPitch, PatternWidth, PatternHeight, Hexagonal, yOffset, layer)-> gf.Component:
    # generates a pattern of sqares with size SSize and horizontal and vertical pitch of SPitch
    # if 'Hexagonal' =0: square grid, =1: hexagonal grid
    # output is a gdsfactory component

    c = gf.Component()
    cLine = gf.Component()

    if yOffset>0:
        yOffset=yOffset%(2*SPitch)
    else:
        yOffset=-(-yOffset%(2*SPitch))

    cBaseSquare=OffsetTaperComp(width1=SSize, width2=SSize, length=SSize, offset=0, layers=[layer])
    cBaseSquare.move((-SSize/2.0,0))
    #BaseShape=polyshape(BaseSquare(:,1), BaseSquare(:,2));

    #how many rows we need?

    


    height=PatternHeight    #2*(RingDiameter/2.0+app.WaveguideWidth/2.0+app.BusRingGap+app.WaveguideWidth+app.OSlabWidth+app.ESlabWidth+app.CladdingWidth+app.HeaterMetal2Protrusion);
    #height2=2*(app.RingDiameter/2.0+app.WaveguideWidth/2.0+app.BusRingGap+app.WaveguideWidth/2.0+app.HeaterWidth/2.0+app.HeaterMetal2Protrusion);
    #height=max([height height2]);
    width=PatternWidth      #2*(app.RingDiameter/2.0+app.WaveguideWidth/2.0+app.OSlabWidth+app.ESlabWidth+app.CladdingWidth);
    #width2=app.OuterESlabLength;
    #width3=app.ElectricalLeadProtrusion*2;
    #width4=2*(app.OuterESlabLength/2.0+app.HeaterContactExtension+app.HeaterMetal1Overlap);
    #width=max([width width2 width3 width4]);

    NoXhalf=math.floor((width/SPitch)/2.0)
    NoYhalf=math.floor((height/SPitch)/2.0)

    #Pattern=polyshape();

    xOffset=0

    #LinePattern=polyshape();
    cLine.add_ref(gf.components.shapes.rectangle(size=(SSize, SSize), layer=layer, centered=True, port_type=None, port_orientations=None), columns=NoXhalf*2+1, rows=1, column_pitch=SPitch, row_pitch=SPitch)
    cLine.move((-NoXhalf*SPitch, 0))

    #for i in range(-NoXhalf,NoXhalf+1,1):
    #    cLine << cBaseSquare.move((i*SPitch,0))  # addboundary(LinePattern, BaseSquare(:,1)+i*SPitch, BaseSquare(:,2) );
    
    PatternRows = [None]*(NoYhalf*2+1) #len(layers)
    for j in range(-NoYhalf,NoYhalf+1,1):
        y=j*SPitch
        if (Hexagonal>0):
            xOffset=SPitch/2.0*(abs(j)%2)   # if j is even, then xOffset becomes =SPitch/2, and =0 otherwise
        PatternRows[j+NoYhalf] = c << cLine
        PatternRows[j+NoYhalf].move((xOffset,y+yOffset))

        #Pattern=union(Pattern, translate(LinePattern, xOffset, y+yOffset));
    
    
    #plot(Pattern);

    SquarePattern=c#Line#Pattern
    return SquarePattern




#def PolygonArea(polygon):



def RemoveDegenerateSquares(component, size, layer)-> gf.Component:
    c=gf.Component()
    errormargin = 1e-6
    newsize=size-errormargin
    pointdict=component.get_polygons_points(merge=False)
    for layer, polygons in pointdict.items():
        for i, polygon in enumerate(polygons):
            #print(i)
            #print(polygon)
            if len(polygon)==4:     # if the polygon has more or less than 4 corners, than it's out for shure
                #now check side lengths
                miss=False
                for i in range(4):
                    l=dist(polygon[i], polygon[(i+1)%4])
                    #print(l)
                    if l<newsize:
                        miss=True
                if miss==False:
                    c.add_polygon(polygon, layer)
    return c





def PlotComponent(
    component,
    facecolor=((1,0,0.4),0.2),
    edgecolor=('green', 1.0),
    hatch='O.'
    ):

    pointdict=component.get_polygons_points(merge=False)
    for layer, polygons in pointdict.items():
        for i, polygon in enumerate(polygons):
            fastAxes.fill(polygon[:,0], polygon[:,1], facecolor=facecolor, edgecolor=edgecolor, hatch=hatch)
            #plt.plot(polygon[:,0], polygon[:,1], color=edgecolor, marker='o', markersize=3)
# --- End of PlotComponent ---



def OffsetTaperComp(
    width1: float = 3,
    width2: float = 10,
    length: float = 20,
    offset: float = 5,
    layers=[(1, 0)]
) -> gf.Component:
    c = gf.Component()
    P=[(0,width1/2), (0,-width1/2), (length,-width2/2+offset), (length,width2/2+offset)]
    for layer in layers:
        #gf.routing.route_quad(c, port1, port2, layer=layer)
        c.add_polygon(P, layer=layer)
    c.add_port(name="in0",  center=(0,0), width=width1, orientation=0, layer=layers[0])
    c.add_port(name="out0", center=(length,offset), width=width2, orientation=0, layer=layers[0])
    return c



def copyLayer(comp: gf.Component, fromlayer, tolayer)->gf.Component:
    return gf.boolean(comp, comp, 'or', layer1=fromlayer, layer2=fromlayer, layer=tolayer)



def CalculateGeometry(maincellname):    # returns a tuple of (gf.Component, ConfigString)
    print("CalculateGeometry function invoked")
    #cWGRingInnerCircle = gf.Component()
    #WGRingInnerCircle = cWGRingInnerCircle << gf.components.circle(radius = global_vars['RingDiameter'].get()/2.0-global_vars['RingWGWidth'].get()/2.0, angle_resolution = 360.0/global_vars['RingDiscretization'].get(), layer=(global_layers['layer_Si'][0].get(), global_layers['layer_Si'][1].get()) )

    gf.clear_cache()    # Clear the GDSFactory high-level cache: remove any cells from gdsfactory
    gf.kcl.clear()      # Clear the KLayout backend library: remove any cells from underlying KLayout Library (KCL)

    # Remove a specific named cell (https://github.com/gdsfactory/gdsfactory/discussions/4387)
    #rm_cell = gf.kcl.cell("RM")
    #if rm_cell is None:
    #    rm = gf.Component("RM")
    #else:
    #    rm = gf.Component(base=gf.kcl[rm_cell.cell_index()].base)
    #    rm.clear()

    GDS_DBU = global_vars['GDS_DBU'].get()
    GDS_UU = global_vars['GDS_UU'].get()
    global scale
    scale = GDS_DBU/GDS_UU

    # !!! all dimensions are transformed to µm from here on !!!! (or better from database units (nm as Standard) to user units (micrometers as standard))
    CladdingWidth = scale*global_vars['CladdingWidth'].get()
    RingDiscretization = global_vars['RingDiscretization'].get()
    MinAngle = global_vars['MinAngle'].get()
    BevelLength = scale*global_vars['BevelLength'].get()

    RingDiameter = scale*global_vars['RingDiameter'].get()
    RingWGWidth = scale*global_vars['RingWGWidth'].get()
    BusGap = scale*global_vars['BusGap'].get()
    DropGap = scale*global_vars['DropGap'].get()
    BusWGWidth = scale*global_vars['BusWGWidth'].get()
    DropWGWidth = scale*global_vars['DropWGWidth'].get()

    OSlabWidth = scale*global_vars['OSlabWidth'].get()
    BusWGLength = scale*global_vars['BusWGLength'].get()
    DropWGLength = scale*global_vars['DropWGLength'].get()
    TaperLength = scale*global_vars['TaperLength'].get()
    SlabBorderExtensionOnTapers = scale*global_vars['SlabBorderExtensionOnTapers'].get()
    PortWGWidth = scale*global_vars['PortWGWidth'].get()

    ESlabWidth = scale*global_vars['ESlabWidth'].get()
    OuterESlabLength = scale*global_vars['OuterESlabLength'].get()

    JunctionOffset = scale*global_vars['JunctionOffset'].get()
    HighDopingDistanceInside = scale*global_vars['HighDopingDistanceInside'].get()
    HighDopingDistanceOutside = scale*global_vars['HighDopingDistanceOutside'].get()
    ContactPlugRegionBorderWidth = scale*global_vars['ContactPlugRegionBorderWidth'].get()
    ContactPlugSize = scale*global_vars['ContactPlugSize'].get()
    ContactPlugPitch = scale*global_vars['ContactPlugPitch'].get()
    ContactPlugHex = global_vars['ContactPlugHex'].get()
    ContactPlugPatternYOffset = scale*global_vars['ContactPlugPatternYOffset'].get()

    ElectricalLeadProtrusion = scale*global_vars['ElectricalLeadProtrusion'].get()
    ElectricalLeadWidth = scale*global_vars['ElectricalLeadWidth'].get()
    ElectricalLeadGap = scale*global_vars['ElectricalLeadGap'].get()
    InnerElectricalRingCutawayOffset = scale*global_vars['InnerElectricalRingCutawayOffset'].get()

    MetalViaBorderWidth = scale*global_vars['MetalViaBorderWidth'].get()
    ViaSize = scale*global_vars['ViaSize'].get()
    ViaPitch = scale*global_vars['ViaPitch'].get()
    ViaHex = global_vars['ViaHex'].get()
    ViaPatternYOffset = scale*global_vars['ViaPatternYOffset'].get()

    HeaterWidth = scale*global_vars['HeaterWidth'].get()
    HeaterLegSpacing = scale*global_vars['HeaterLegSpacing'].get()
    HeaterContactExtension = scale*global_vars['HeaterContactExtension'].get()
    HeaterMetal1Overlap = scale*global_vars['HeaterMetal1Overlap'].get()
    HeaterMetal1Protrusion = scale*global_vars['HeaterMetal1Protrusion'].get()
    HeaterMetal2Protrusion = scale*global_vars['HeaterMetal2Protrusion'].get()


    layer_Cladd = (global_layers['layer_Cladd'][0].get(), global_layers['layer_Cladd'][1].get())
    layer_Si = (global_layers['layer_Si'][0].get(), global_layers['layer_Si'][1].get())
    layer_OuterOSlabCore = (global_layers['layer_OuterOSlabCore'][0].get(), global_layers['layer_OuterOSlabCore'][1].get())
    layer_OuterOSlabCladd = (global_layers['layer_OuterOSlabCladd'][0].get(), global_layers['layer_OuterOSlabCladd'][1].get())
    layer_InnerOSlabCladd = (global_layers['layer_InnerOSlabCladd'][0].get(), global_layers['layer_InnerOSlabCladd'][1].get())
    layer_OuterLowDope = (global_layers['layer_OuterLowDope'][0].get(), global_layers['layer_OuterLowDope'][1].get())
    layer_InnerLowDope = (global_layers['layer_InnerLowDope'][0].get(), global_layers['layer_InnerLowDope'][1].get())
    layer_OuterHiDope = (global_layers['layer_OuterHiDope'][0].get(), global_layers['layer_OuterHiDope'][1].get())
    layer_InnerHiDope = (global_layers['layer_InnerHiDope'][0].get(), global_layers['layer_InnerHiDope'][1].get())
    layer_OuterContactDope = (global_layers['layer_OuterContactDope'][0].get(), global_layers['layer_OuterContactDope'][1].get())
    layer_InnerContactDope = (global_layers['layer_InnerContactDope'][0].get(), global_layers['layer_InnerContactDope'][1].get())
    layer_M1 = (global_layers['layer_M1'][0].get(), global_layers['layer_M1'][1].get())
    layer_M1PerfBlock = (global_layers['layer_M1PerfBlock'][0].get(), global_layers['layer_M1PerfBlock'][1].get())
    layer_M2 = (global_layers['layer_M2'][0].get(), global_layers['layer_M2'][1].get())
    layer_M2PerfBlock = (global_layers['layer_M2PerfBlock'][0].get(), global_layers['layer_M2PerfBlock'][1].get())
    layer_Heater = (global_layers['layer_Heater'][0].get(), global_layers['layer_Heater'][1].get())
    layer_HeaterContact = (global_layers['layer_HeaterContact'][0].get(), global_layers['layer_HeaterContact'][1].get())
    layer_Vias = (global_layers['layer_Vias'][0].get(), global_layers['layer_Vias'][1].get())
    layer_ContactPlugs = (global_layers['layer_ContactPlugs'][0].get(), global_layers['layer_ContactPlugs'][1].get())


    #global ConfigString
    ConfigString="Ring Modulator Layout\n"
    ConfigString=ConfigString+"generated by Ring Modulator Layout Program 2 "
    ConfigString=ConfigString+ProgramVersionString+" from "+ProgramVersionDateString+"\n"
    ConfigString=ConfigString+"by Dr.-Ing. Marc Schneider (KIT)\n\n"

    ConfigString=ConfigString+"[Variables]\n"
    for key, value in global_vars.items():
        #print(key, "=", str( global_vars[key].get()))
        ConfigString=ConfigString+key+"="+str( global_vars[key].get())+"\n"

    ConfigString=ConfigString+"\n[Layers]\n"
    for key, value in global_layers.items():
        #print(key, "=", str( global_vars[key].get()))
        ConfigString=ConfigString+key+"=("+str( global_layers[key][0].get())+","+str( global_layers[key][1].get())+")\n"
    #print(ConfigString)

    circleAngleResolution = 360.0/RingDiscretization


    cRingRibInnerCircle = gf.components.circle(radius = RingDiameter/2.0-RingWGWidth/2.0, angle_resolution = circleAngleResolution, layer=layer_Si )
    cRingRibCenterCircle = gf.components.circle(radius = RingDiameter/2.0, angle_resolution = circleAngleResolution, layer=layer_Si )
    cRingRibOuterCircle = gf.components.circle(radius = RingDiameter/2.0+RingWGWidth/2.0, angle_resolution = circleAngleResolution, layer=layer_Si )
    cRingRibJunctionCircle = gf.components.circle(radius = RingDiameter/2.0+JunctionOffset, angle_resolution = circleAngleResolution, layer=layer_Si )
    cRingInnerOSlabCircle = gf.components.circle(radius = RingDiameter/2.0-RingWGWidth/2.0-OSlabWidth, angle_resolution = circleAngleResolution, layer=layer_Si )
    cRingOuterOSlabCircle = gf.components.circle(radius = RingDiameter/2.0+RingWGWidth/2.0+OSlabWidth, angle_resolution = circleAngleResolution, layer=layer_Si )

    cRingInnerESlabCircle = gf.components.circle(radius = RingDiameter/2.0-RingWGWidth/2.0-OSlabWidth-ESlabWidth, angle_resolution = circleAngleResolution, layer=layer_Si )
    cRingOuterESlabCircle = gf.components.circle(radius = RingDiameter/2.0+RingWGWidth/2.0+OSlabWidth+ESlabWidth, angle_resolution = circleAngleResolution, layer=layer_Si )

    cRingInnerLoHiDopingCircle = gf.components.circle(radius = RingDiameter/2.0-RingWGWidth/2.0-HighDopingDistanceInside, angle_resolution = circleAngleResolution, layer=layer_Si )
    cRingOuterLoHiDopingCircle = gf.components.circle(radius = RingDiameter/2.0+RingWGWidth/2.0+HighDopingDistanceOutside, angle_resolution = circleAngleResolution, layer=layer_Si )

    cRingInnerCladdCircle = gf.components.circle(radius = RingDiameter/2.0-RingWGWidth/2.0-OSlabWidth-ESlabWidth-CladdingWidth, angle_resolution = circleAngleResolution, layer=layer_Si )
    cRingOuterCladdCircle = gf.components.circle(radius = RingDiameter/2.0+RingWGWidth/2.0+OSlabWidth+ESlabWidth+CladdingWidth, angle_resolution = circleAngleResolution, layer=layer_Si )

    BusRibCenterY = -(RingDiameter/2.0+RingWGWidth/2.0+BusGap+BusWGWidth/2.0)
    if round(BusRibCenterY/scale)%2 == 0:
        BusPortCenterY=BusRibCenterY
    else:
        BusPortCenterY=BusRibCenterY-GDS_DBU/GDS_UU #/2.0

    DropRibCenterY = RingDiameter/2.0+RingWGWidth/2.0+DropGap+DropWGWidth/2.0
    if round(DropRibCenterY/scale)%2 == 0:
        DropPortCenterY=DropRibCenterY
    else:
        DropPortCenterY=DropRibCenterY+GDS_DBU/GDS_UU #/2.0


    #BusPortCenterY=BusRibCenterY-1000
    #DropPortCenterY=DropRibCenterY+1000
    PortPosBusLeft = (-BusWGLength/2.0, BusPortCenterY)
    PortPosBusRight = (BusWGLength/2.0, BusPortCenterY)
    PortPosDropLeft = (-DropWGLength/2.0, DropPortCenterY)
    PortPosDropRight = (DropWGLength/2.0, DropPortCenterY)
    #plt.plot(PortPosBusLeft[0], PortPosBusLeft[1], 'bo')
    #plt.plot(PortPosBusRight[0], PortPosBusRight[1], 'bo')
    #plt.plot(PortPosDropLeft[0], PortPosDropLeft[1], 'bo')
    #plt.plot(PortPosDropRight[0], PortPosDropRight[1], 'bo')
    cBusRibCenter = OffsetTaperComp(width1=BusWGWidth, width2=BusWGWidth, length=BusWGLength-2*TaperLength, offset=0, layers=[layer_Si])
    cBusRibCenter.move((-BusWGLength/2.0+TaperLength, BusRibCenterY ))
    cBusRibLeftTaper = OffsetTaperComp(width1=PortWGWidth, width2=BusWGWidth, length=TaperLength, offset=BusRibCenterY-BusPortCenterY, layers=[layer_Si])
    cBusRibLeftTaper.move(cBusRibLeftTaper['out0'].center, cBusRibCenter['in0'].center)
    cBusRibRightTaper = OffsetTaperComp(width1=BusWGWidth, width2=PortWGWidth, length=TaperLength, offset=BusPortCenterY-BusRibCenterY, layers=[layer_Si])
    cBusRibRightTaper.move(cBusRibRightTaper['in0'].center, cBusRibCenter['out0'].center)
    cBusOSlabCenter = OffsetTaperComp(width1=BusWGWidth+2*OSlabWidth, width2=BusWGWidth+2*OSlabWidth, length=BusWGLength-2*TaperLength, offset=0, layers=[layer_Si])
    cBusOSlabCenter.move(cBusOSlabCenter['in0'].center, cBusRibCenter['in0'].center)
    cBusOSlabLeftTaper = OffsetTaperComp(width1=PortWGWidth, width2=BusWGWidth+2*OSlabWidth, length=TaperLength, offset=BusRibCenterY-BusPortCenterY, layers=[layer_Si])
    cBusOSlabLeftTaper.move(cBusOSlabLeftTaper['in0'].center, cBusRibLeftTaper['in0'].center)
    cBusOSlabRightTaper = OffsetTaperComp(width1=BusWGWidth+2*OSlabWidth, width2=PortWGWidth, length=TaperLength, offset=BusPortCenterY-BusRibCenterY, layers=[layer_Si])
    cBusOSlabRightTaper.move(cBusOSlabRightTaper['in0'].center, cBusRibRightTaper['in0'].center)

    cBusOSlabLeftTaperBordered = OffsetTaperComp(width1=PortWGWidth+2*SlabBorderExtensionOnTapers, width2=BusWGWidth+2*OSlabWidth+2*SlabBorderExtensionOnTapers, length=TaperLength+2*SlabBorderExtensionOnTapers, offset=BusRibCenterY-BusPortCenterY, layers=[layer_Si])
    cBusOSlabLeftTaperBordered.move(cBusOSlabLeftTaperBordered['in0'].center, cBusRibLeftTaper['in0'].center)
    cBusOSlabLeftTaperBordered.move((-SlabBorderExtensionOnTapers, 0))
    cBusOSlabLeftTaperBordered = BevelComponent(cBusOSlabLeftTaperBordered, MinAngle, BevelLength, layer=layer_Si)

    cBusOSlabRightTaperBordered = OffsetTaperComp(width1=BusWGWidth+2*OSlabWidth+2*SlabBorderExtensionOnTapers, width2=PortWGWidth+2*SlabBorderExtensionOnTapers, length=TaperLength+2*SlabBorderExtensionOnTapers, offset=BusPortCenterY-BusRibCenterY, layers=[layer_Si])
    cBusOSlabRightTaperBordered.move(cBusOSlabRightTaperBordered['in0'].center, cBusRibRightTaper['in0'].center)
    cBusOSlabRightTaperBordered.move((-SlabBorderExtensionOnTapers, 0))
    cBusOSlabRightTaperBordered = BevelComponent(cBusOSlabRightTaperBordered, MinAngle, BevelLength, layer=layer_Si)


    cBusCladdCenter = OffsetTaperComp(width1=BusWGWidth+2*OSlabWidth+2*CladdingWidth, width2=BusWGWidth+2*OSlabWidth+2*CladdingWidth, length=BusWGLength-2*TaperLength, offset=0, layers=[layer_Si])
    cBusCladdCenter.move(cBusCladdCenter['in0'].center, cBusRibCenter['in0'].center)
    cBusCladdLeftTaper = OffsetTaperComp(width1=PortWGWidth+2*CladdingWidth, width2=BusWGWidth+2*OSlabWidth+2*CladdingWidth, length=TaperLength, offset=BusRibCenterY-BusPortCenterY, layers=[layer_Si])
    cBusCladdLeftTaper.move(cBusCladdLeftTaper['in0'].center, cBusRibLeftTaper['in0'].center)
    cBusCladdRightTaper = OffsetTaperComp(width1=BusWGWidth+2*OSlabWidth+2*CladdingWidth, width2=PortWGWidth+2*CladdingWidth, length=TaperLength, offset=BusPortCenterY-BusRibCenterY, layers=[layer_Si])
    cBusCladdRightTaper.move(cBusCladdRightTaper['in0'].center, cBusRibRightTaper['in0'].center)

    cBusLoHiDopingRectangle = OffsetTaperComp(width1=BusWGWidth+2*HighDopingDistanceOutside, width2=BusWGWidth+2*HighDopingDistanceOutside, length=BusWGLength, offset=0, layers=[layer_Si])
    cBusLoHiDopingRectangle.move(cBusLoHiDopingRectangle['in0'].center, cBusRibLeftTaper['in0'].center)


    cBusRib = gf.boolean(A=cBusRibCenter, B=cBusRibLeftTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cBusRib = gf.boolean(A=cBusRib, B=cBusRibRightTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    #cBusOSlab = gf.boolean(A=cBusOSlabCenter, B=cBusOSlabLeftTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    #cBusOSlab = gf.boolean(A=cBusOSlab, B=cBusOSlabRightTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cBusOSlab = gf.boolean(A=cBusOSlabLeftTaper, B=cBusOSlabRightTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cBusOSlab = gf.boolean(A=cBusOSlab, B=cBusOSlabCenter, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cBusOSlabBordered = gf.boolean(A=cBusOSlab, B=cBusOSlabLeftTaperBordered, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cBusOSlabBordered = gf.boolean(A=cBusOSlabBordered, B=cBusOSlabRightTaperBordered, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cBusCladd = gf.boolean(A=cBusCladdCenter, B=cBusCladdLeftTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cBusCladd = gf.boolean(A=cBusCladd, B=cBusCladdRightTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    

    cDropRibCenter = OffsetTaperComp(width1=DropWGWidth, width2=DropWGWidth, length=DropWGLength-2*TaperLength, offset=0, layers=[layer_Si])
    cDropRibCenter.move((-DropWGLength/2.0+TaperLength, DropRibCenterY ))
    cDropRibLeftTaper = OffsetTaperComp(width1=PortWGWidth, width2=DropWGWidth, length=TaperLength, offset=DropRibCenterY-DropPortCenterY, layers=[layer_Si])
    cDropRibLeftTaper.move(cDropRibLeftTaper['out0'].center, cDropRibCenter['in0'].center)
    cDropRibRightTaper = OffsetTaperComp(width1=DropWGWidth, width2=PortWGWidth, length=TaperLength, offset=DropPortCenterY-DropRibCenterY, layers=[layer_Si])
    cDropRibRightTaper.move(cDropRibRightTaper['in0'].center, cDropRibCenter['out0'].center)
    cDropOSlabCenter = OffsetTaperComp(width1=DropWGWidth+2*OSlabWidth, width2=DropWGWidth+2*OSlabWidth, length=DropWGLength-2*TaperLength, offset=0, layers=[layer_Si])
    cDropOSlabCenter.move(cDropOSlabCenter['in0'].center, cDropRibCenter['in0'].center)
    cDropOSlabLeftTaper = OffsetTaperComp(width1=PortWGWidth, width2=DropWGWidth+2*OSlabWidth, length=TaperLength, offset=DropRibCenterY-DropPortCenterY, layers=[layer_Si])
    cDropOSlabLeftTaper.move(cDropOSlabLeftTaper['in0'].center, cDropRibLeftTaper['in0'].center)
    cDropOSlabRightTaper = OffsetTaperComp(width1=DropWGWidth+2*OSlabWidth, width2=PortWGWidth, length=TaperLength, offset=DropPortCenterY-DropRibCenterY, layers=[layer_Si])
    cDropOSlabRightTaper.move(cDropOSlabRightTaper['in0'].center, cDropRibRightTaper['in0'].center)

    cDropOSlabLeftTaperBordered = OffsetTaperComp(width1=PortWGWidth+2*SlabBorderExtensionOnTapers, width2=DropWGWidth+2*OSlabWidth+2*SlabBorderExtensionOnTapers, length=TaperLength+2*SlabBorderExtensionOnTapers, offset=DropRibCenterY-DropPortCenterY, layers=[layer_Si])
    cDropOSlabLeftTaperBordered.move(cDropOSlabLeftTaperBordered['in0'].center, cDropRibLeftTaper['in0'].center)
    cDropOSlabLeftTaperBordered.move((-SlabBorderExtensionOnTapers, 0))
    cDropOSlabLeftTaperBordered = BevelComponent(cDropOSlabLeftTaperBordered, MinAngle, BevelLength, layer=layer_Si)

    cDropOSlabRightTaperBordered = OffsetTaperComp(width1=DropWGWidth+2*OSlabWidth+2*SlabBorderExtensionOnTapers, width2=PortWGWidth+2*SlabBorderExtensionOnTapers, length=TaperLength+2*SlabBorderExtensionOnTapers, offset=DropPortCenterY-DropRibCenterY, layers=[layer_Si])
    cDropOSlabRightTaperBordered.move(cDropOSlabRightTaperBordered['in0'].center, cDropRibRightTaper['in0'].center)
    cDropOSlabRightTaperBordered.move((-SlabBorderExtensionOnTapers, 0))
    cDropOSlabRightTaperBordered = BevelComponent(cDropOSlabRightTaperBordered, MinAngle, BevelLength, layer=layer_Si)


    cDropCladdCenter = OffsetTaperComp(width1=DropWGWidth+2*OSlabWidth+2*CladdingWidth, width2=DropWGWidth+2*OSlabWidth+2*CladdingWidth, length=DropWGLength-2*TaperLength, offset=0, layers=[layer_Si])
    cDropCladdCenter.move(cDropCladdCenter['in0'].center, cDropRibCenter['in0'].center)
    cDropCladdLeftTaper = OffsetTaperComp(width1=PortWGWidth+2*CladdingWidth, width2=DropWGWidth+2*OSlabWidth+2*CladdingWidth, length=TaperLength, offset=DropRibCenterY-DropPortCenterY, layers=[layer_Si])
    cDropCladdLeftTaper.move(cDropCladdLeftTaper['in0'].center, cDropRibLeftTaper['in0'].center)
    cDropCladdRightTaper = OffsetTaperComp(width1=DropWGWidth+2*OSlabWidth+2*CladdingWidth, width2=PortWGWidth+2*CladdingWidth, length=TaperLength, offset=DropPortCenterY-DropRibCenterY, layers=[layer_Si])
    cDropCladdRightTaper.move(cDropCladdRightTaper['in0'].center, cDropRibRightTaper['in0'].center)

    cDropLoHiDopingRectangle = OffsetTaperComp(width1=DropWGWidth+2*HighDopingDistanceOutside, width2=DropWGWidth+2*HighDopingDistanceOutside, length=DropWGLength, offset=0, layers=[layer_Si])
    cDropLoHiDopingRectangle.move(cDropLoHiDopingRectangle['in0'].center, cDropRibLeftTaper['in0'].center)


    cDropRib = gf.boolean(A=cDropRibCenter, B=cDropRibLeftTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cDropRib = gf.boolean(A=cDropRib, B=cDropRibRightTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cDropOSlab = gf.boolean(A=cDropOSlabCenter, B=cDropOSlabLeftTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cDropOSlab = gf.boolean(A=cDropOSlab, B=cDropOSlabRightTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cDropOSlabBordered = gf.boolean(A=cDropOSlab, B=cDropOSlabLeftTaperBordered, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cDropOSlabBordered = gf.boolean(A=cDropOSlabBordered, B=cDropOSlabRightTaperBordered, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cDropCladd = gf.boolean(A=cDropCladdCenter, B=cDropCladdLeftTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cDropCladd = gf.boolean(A=cDropCladd, B=cDropCladdRightTaper, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)


    OuterESlabTopY = RingDiameter/2.0+RingWGWidth/2.0+DropGap+DropWGWidth+OSlabWidth+ESlabWidth
    OuterESlabBottomY = -(RingDiameter/2.0+RingWGWidth/2.0+BusGap+BusWGWidth+OSlabWidth+ESlabWidth)
    #cOuterESlabRectangle = gf.components.rectangle(size=(OuterESlabLength, OuterESlabTopY-OuterESlabBottomY) , layer=layer_Si)

    cOuterESlabRectangle = OffsetTaperComp(width1=OuterESlabTopY-OuterESlabBottomY, width2=OuterESlabTopY-OuterESlabBottomY, length=OuterESlabLength, offset=0, layers=[layer_Si])
    cOuterESlabRectangle.move((-OuterESlabLength/2.0, OuterESlabTopY-(OuterESlabTopY-OuterESlabBottomY)/2.0 ))
    cOuterESlabCladdRectangle = OffsetTaperComp(width1=OuterESlabTopY-OuterESlabBottomY+2*CladdingWidth, width2=OuterESlabTopY-OuterESlabBottomY+2*CladdingWidth, length=OuterESlabLength+2*CladdingWidth, offset=0, layers=[layer_Si])
    cOuterESlabCladdRectangle.move(cOuterESlabCladdRectangle['in0'].center, cOuterESlabRectangle['in0'].center)
    cOuterESlabCladdRectangle.move((-CladdingWidth, 0))

    cOuterESlabInnerRectangle = OffsetTaperComp(width1=OuterESlabTopY-OuterESlabBottomY-2*ESlabWidth, width2=OuterESlabTopY-OuterESlabBottomY-2*ESlabWidth, length=OuterESlabLength-2*ESlabWidth, offset=0, layers=[layer_Si])
    cOuterESlabInnerRectangle.move(cOuterESlabInnerRectangle['in0'].center, cOuterESlabRectangle['in0'].center)
    cOuterESlabInnerRectangle.move((ESlabWidth,0))

    cInnerM2Contact = OffsetTaperComp(width1=ElectricalLeadWidth, width2=ElectricalLeadWidth, length=ElectricalLeadProtrusion, offset=0, layers=[layer_Si])
    cInnerM2CutAwayTinyGap = OffsetTaperComp(width1=InnerElectricalMetal2RingMinimumGap, width2=InnerElectricalMetal2RingMinimumGap, length=ElectricalLeadProtrusion, offset=0, layers=[layer_Si])
    cInnerM2CutAwayTinyGap.move((-ElectricalLeadProtrusion, 0))
    cInnerM2CutAway = OffsetTaperComp(width1=RingDiameter, width2=RingDiameter, length=RingDiameter, offset=0, layers=[layer_Si])
    cInnerM2CutAway.move((-RingDiameter+InnerElectricalRingCutawayOffset, 0))
    PortPosRFInner = (ElectricalLeadProtrusion, 0)
    #plt.plot(PortPosRFInner[0], PortPosRFInner[1], 'bo')

    cOuterM2Contact = OffsetTaperComp(width1=ElectricalLeadWidth, width2=ElectricalLeadWidth, length=ElectricalLeadProtrusion, offset=0, layers=[layer_Si])
    cOuterM2Contact.move((-ElectricalLeadProtrusion, 0))
    cOuterM2CutAway = OffsetTaperComp(width1=ElectricalLeadWidth+2*ElectricalLeadGap, width2=ElectricalLeadWidth+2*ElectricalLeadGap, length=ElectricalLeadProtrusion, offset=0, layers=[layer_Si])
    PortPosRFOuter = (-ElectricalLeadProtrusion, 0)
    #plt.plot(PortPosRFOuter[0], PortPosRFOuter[1], 'bo')

    #cOuterESlabRectangle.move((-OuterESlabLength/2.0, OuterESlabTopY-(OuterESlabTopY-OuterESlabBottomY)/2.0   ))
    #cOuterESlabRectangle = gf.Component.add_polygon( [(-OuterESlabLength/2.0, OuterESlabTopY), (-OuterESlabLength/2.0, OuterESlabBottomY), (OuterESlabLength/2.0, OuterESlabBottomY), (OuterESlabLength/2.0, OuterESlabTopY)] , layer=layer_Si)
    #cOuterESlabCladdRectangle = gf.Component.add_polygon( [(-OuterESlabLength/2.0-CladdingWidth, OuterESlabTopY+CladdingWidth), (-OuterESlabLength/2.0-CladdingWidth, OuterESlabBottomY-CladdingWidth), (OuterESlabLength/2.0+CladdingWidth, OuterESlabBottomY-CladdingWidth), (OuterESlabLength/2.0+CladdingWidth, OuterESlabTopY+CladdingWidth)] , layer=layer_Si)


    cHeaterLegsRectWidth=DropRibCenterY+HeaterWidth/2.0
    #print(cHeaterLegsRectWidth)
    if round(cHeaterLegsRectWidth/scale)%2>0:
        cHeaterLegsRectWidth += GDS_DBU/GDS_UU
    #print(cHeaterLegsRectWidth)
    cHeaterLegsRectTop = DropRibCenterY+HeaterWidth/2.0
    cHeaterLegsRectYShift = -cHeaterLegsRectWidth/2+cHeaterLegsRectTop
    #print(cHeaterLegsRectYShift)
    cHeaterLegsRectWidth = round(cHeaterLegsRectWidth/scale)*scale
    #print("cHeaterLegsRectWidth=", cHeaterLegsRectWidth, "    cHeaterLegsRectYShift=", cHeaterLegsRectYShift)

    cRingInnerHeaterCircle = gf.components.circle(radius = RingDiameter/2.0-HeaterWidth/2.0, angle_resolution = circleAngleResolution, layer=layer_Si )
    cRingOuterHeaterCircle = gf.components.circle(radius = RingDiameter/2.0+HeaterWidth/2.0, angle_resolution = circleAngleResolution, layer=layer_Si )
    #print("1...")
    cHeaterLegsInnerRect = OffsetTaperComp(width1=cHeaterLegsRectWidth, width2=cHeaterLegsRectWidth, length=HeaterLegSpacing, offset=0, layers=[layer_Si])
    cHeaterLegsInnerRect.move((-HeaterLegSpacing/2.0, cHeaterLegsRectYShift))
    #print("2...")
    cHeaterLegsOuterRect = OffsetTaperComp(width1=cHeaterLegsRectWidth, width2=cHeaterLegsRectWidth, length=HeaterLegSpacing+2*HeaterWidth, offset=0, layers=[layer_Si])
    cHeaterLegsOuterRect.move((-HeaterLegSpacing/2.0-HeaterWidth, cHeaterLegsRectYShift))
    #print("3...")
    cHeaterContactRect = OffsetTaperComp(width1=HeaterWidth, width2=HeaterWidth, length=OuterESlabLength+2*HeaterContactExtension+2*HeaterMetal1Overlap, offset=0, layers=[layer_Si])
    cHeaterContactRect.move((-(OuterESlabLength+2*HeaterContactExtension+2*HeaterMetal1Overlap)/2.0, DropRibCenterY))

    HeaterMetal1Width = HeaterWidth+HeaterWidth/2+HeaterMetal1Protrusion
    if round(HeaterMetal1Width/scale)%2>0:
        HeaterMetal1Width += GDS_DBU/GDS_UU

    #print("4...")
    cHeaterMetal1Left = OffsetTaperComp(width1=HeaterMetal1Width, width2=HeaterMetal1Width, length=HeaterMetal1Overlap, offset=0, layers=[layer_Si])
    cHeaterMetal1Left.move((-(OuterESlabLength+2*HeaterContactExtension+2*HeaterMetal1Overlap)/2.0, DropRibCenterY+HeaterMetal1Width/2-HeaterWidth))
    #print("5...")
    cHeaterMetal1Right = OffsetTaperComp(width1=HeaterMetal1Width, width2=HeaterMetal1Width, length=HeaterMetal1Overlap, offset=0, layers=[layer_Si])
    cHeaterMetal1Right.move((OuterESlabLength/2.0+HeaterContactExtension, DropRibCenterY+HeaterMetal1Width/2-HeaterWidth))

    HeaterMetal2Width = HeaterWidth+HeaterWidth/2+HeaterMetal2Protrusion
    if round(HeaterMetal2Width/scale)%2>0:
        HeaterMetal2Width += GDS_DBU/GDS_UU

    #print("6...")
    cHeaterMetal2Left = OffsetTaperComp(width1=HeaterMetal2Width, width2=HeaterMetal2Width, length=HeaterMetal1Overlap, offset=0, layers=[layer_Si])
    cHeaterMetal2Left.move((-(OuterESlabLength+2*HeaterContactExtension+2*HeaterMetal1Overlap)/2.0, DropRibCenterY+HeaterMetal2Width/2-HeaterWidth))
    #print("7...")
    cHeaterMetal2Right = OffsetTaperComp(width1=HeaterMetal2Width, width2=HeaterMetal2Width, length=HeaterMetal1Overlap, offset=0, layers=[layer_Si])
    cHeaterMetal2Right.move((OuterESlabLength/2.0+HeaterContactExtension, DropRibCenterY+HeaterMetal2Width/2-HeaterWidth))

    PortPosHeaterLeft = (-OuterESlabLength/2.0-HeaterContactExtension-HeaterMetal1Overlap/2.0, DropRibCenterY+HeaterMetal2Width-HeaterWidth)
    PortPosHeaterRight = (OuterESlabLength/2.0+HeaterContactExtension+HeaterMetal1Overlap/2.0, DropRibCenterY+HeaterMetal2Width-HeaterWidth)
    #plt.plot(PortPosHeaterLeft[0], PortPosHeaterLeft[1], 'bo')
    #plt.plot(PortPosHeaterRight[0], PortPosHeaterRight[1], 'bo')

    #PlotComponent(component=cWGRingInnerCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=('black', 1.0), hatch='')
    if (global_vars['ShowConstructionLines'].get()==True):
        PlotComponent(component=cRingRibInnerCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cRingRibOuterCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cRingRibCenterCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.5, 1.0, 0.5), 0.5), hatch='')
        PlotComponent(component=cRingRibJunctionCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 1.0), 0.5), hatch='')
        PlotComponent(component=cRingInnerOSlabCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 1.0), 0.5), hatch='')
        PlotComponent(component=cRingOuterOSlabCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 1.0), 0.5), hatch='')
        PlotComponent(component=cRingInnerESlabCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cRingOuterESlabCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cRingInnerLoHiDopingCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 1.0), 0.5), hatch='')
        PlotComponent(component=cRingOuterLoHiDopingCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 1.0), 0.5), hatch='')
        PlotComponent(component=cRingInnerCladdCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.75, 0.75, 0.75), 0.5), hatch='')
        PlotComponent(component=cRingOuterCladdCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.75, 0.75, 0.75), 0.5), hatch='')
        #print("   CalculateGeometry function still invoked")
        #plt.show()

        PlotComponent(component=cBusRib, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cBusOSlab, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 1.0), 0.5), hatch='')
        PlotComponent(component=cBusCladd, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.75, 0.75, 0.75), 0.5), hatch='')
        PlotComponent(component=cBusLoHiDopingRectangle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 1.0), 0.5), hatch='')

        PlotComponent(component=cDropRib, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cDropOSlab, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 1.0), 0.5), hatch='')
        PlotComponent(component=cDropCladd, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.75, 0.75, 0.75), 0.5), hatch='')
        PlotComponent(component=cDropLoHiDopingRectangle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 1.0), 0.5), hatch='')

        PlotComponent(component=cOuterESlabRectangle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cOuterESlabInnerRectangle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cOuterESlabCladdRectangle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.75, 0.75, 0.75), 0.5), hatch='')

        PlotComponent(component=cInnerM2Contact, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.75, 0.75, 0.0), 0.5), hatch='')
        PlotComponent(component=cOuterM2Contact, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.85, 0.65, 0.0), 0.5), hatch='')

        PlotComponent(component=cRingInnerHeaterCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cRingOuterHeaterCircle, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cHeaterLegsInnerRect, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cHeaterLegsOuterRect, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cHeaterContactRect, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cHeaterMetal1Left, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cHeaterMetal1Right, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 0.0), 0.5), hatch='')
        PlotComponent(component=cHeaterMetal2Left, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 1.0), 0.5), hatch='')
        PlotComponent(component=cHeaterMetal2Right, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 1.0, 1.0), 0.5), hatch='')


    cAllCladd = gf.boolean(A=cRingOuterCladdCircle, B=cBusCladd, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cAllCladd = gf.boolean(A=cAllCladd, B=cDropCladd, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cAllCladd = gf.boolean(A=cAllCladd, B=cOuterESlabCladdRectangle, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cAllCladd = BevelComponent(cAllCladd, MinAngle, BevelLength, layer=layer_Si)
    cAllCladd = gf.boolean(A=cAllCladd, B=cRingInnerCladdCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    
    cAllSi = gf.boolean(A=cRingOuterESlabCircle, B=cBusOSlab, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cAllSi = gf.boolean(A=cAllSi, B=cDropOSlab, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cAllSi = gf.boolean(A=cAllSi, B=cOuterESlabRectangle, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cAllSi = BevelComponent(cAllSi, MinAngle, BevelLength, layer=layer_Si)
    cAllSi = gf.boolean(A=cAllSi, B=cRingInnerESlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    cOuterOSlab = gf.boolean(A=cRingOuterOSlabCircle, B=cBusOSlabBordered, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterOSlab = gf.boolean(A=cOuterOSlab, B=cDropOSlabBordered, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterOSlab = BevelComponent(cOuterOSlab, MinAngle, BevelLength, layer=layer_Si)
    cOuterOSlab = gf.boolean(A=cOuterOSlab, B=cRingRibOuterCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    cInnerOSlab = gf.boolean(A=cRingRibInnerCircle, B=cRingInnerOSlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    cWGRibs = gf.boolean(A=cBusRib, B=cDropRib, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    cInnerLowDoping = gf.boolean(A=cRingRibJunctionCircle, B=cRingInnerESlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    
    cOuterLowDoping = gf.boolean(A=cRingOuterESlabCircle, B=cOuterESlabRectangle, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterLowDoping = BevelComponent(cOuterLowDoping, MinAngle, BevelLength, layer=layer_Si)
    cOuterLowDoping = gf.boolean(A=cOuterLowDoping, B=cRingRibJunctionCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    cInnerHighDoping = gf.boolean(A=cRingInnerLoHiDopingCircle, B=cRingInnerESlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    cOuterHighDoping = gf.boolean(A=cRingOuterESlabCircle, B=cOuterESlabRectangle, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterHighDoping = gf.boolean(A=cOuterHighDoping, B=cRingOuterLoHiDopingCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterHighDoping = gf.boolean(A=cOuterHighDoping, B=cBusLoHiDopingRectangle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterHighDoping = gf.boolean(A=cOuterHighDoping, B=cDropLoHiDopingRectangle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterHighDoping = BevelComponent(cOuterHighDoping, MinAngle, BevelLength, layer=layer_Si)

    cInnerContactDoping = gf.boolean(A=cRingInnerOSlabCircle, B=cRingInnerESlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    cOuterContactDoping = gf.boolean(A=cRingOuterESlabCircle, B=cOuterESlabRectangle, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterContactDoping = gf.boolean(A=cOuterContactDoping, B=cRingOuterOSlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterContactDoping = gf.boolean(A=cOuterContactDoping, B=cBusOSlab, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterContactDoping = gf.boolean(A=cOuterContactDoping, B=cDropOSlab, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterContactDoping = BevelComponent(cOuterContactDoping, MinAngle, BevelLength, layer=layer_Si)

    # cInnerMetal1 is (currently) the same as cInnerContactDoping
    cInnerMetal1 = gf.boolean(A=cRingInnerOSlabCircle, B=cRingInnerESlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    #cInnerMetal1PerfBlock = cInnerMetal1

    # cOuterMetal1 is (currently) the same as cOuterContactDoping
    cOuterMetal1 = gf.boolean(A=cRingOuterESlabCircle, B=cOuterESlabRectangle, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterMetal1 = gf.boolean(A=cOuterMetal1, B=cRingOuterOSlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterMetal1 = gf.boolean(A=cOuterMetal1, B=cBusOSlab, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterMetal1 = gf.boolean(A=cOuterMetal1, B=cDropOSlab, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterMetal1 = BevelComponent(cOuterMetal1, MinAngle, BevelLength, layer=layer_Si)

    #cOuterMetal1PerfBlock = cOuterMetal1


    # make grid of possible contact plugs
    OpticsBBox=cAllCladd.bbox()
    #print(OpticsBBox)
    PlugPatternWidth=2*max([abs(OpticsBBox.left), abs(OpticsBBox.right)])
    PlugPatternHeight=2*max([abs(OpticsBBox.top), abs(OpticsBBox.bottom)])
    #print(PlugPatternWidth,'  ',PlugPatternHeight)
    cSquarePlugPattern = GenerateSquarePattern(ContactPlugSize, ContactPlugPitch, PlugPatternWidth, PlugPatternHeight, ContactPlugHex, ContactPlugPatternYOffset, layer_Si )
    #PlotComponent(component=cSquarePlugPattern, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 0.0), 0.6), hatch='')


    # calculate Contact Plugs
    rInnerContactDoping = cInnerContactDoping.get_region(layer=layer_Si)  # Regions are in DBU (1 nm in this case).
    rInnerContactPlugsRegion = rInnerContactDoping.sized(-ContactPlugRegionBorderWidth/scale)  # Regions are in DBU.
    cInnerContactPlugsRegion = gf.Component()
    cInnerContactPlugsRegion.add_polygon(rInnerContactPlugsRegion, layer=layer_Si)  # Add the region to the component.
    
    rOuterContactDoping = cOuterContactDoping.get_region(layer=layer_Si)  # Regions are in DBU (1 nm in this case).
    rOuterContactPlugsRegion = rOuterContactDoping.sized(-ContactPlugRegionBorderWidth/scale)  # Regions are in DBU.
    cOuterContactPlugsRegion = gf.Component()
    cOuterContactPlugsRegion.add_polygon(rOuterContactPlugsRegion, layer=layer_Si)  # Add the region to the component.

    cContactPlugsInnerBase = gf.boolean(A=cSquarePlugPattern, B=cInnerContactPlugsRegion, operation='and', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cContactPlugsOuterBase = gf.boolean(A=cSquarePlugPattern, B=cOuterContactPlugsRegion, operation='and', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cContactPlugsBase = gf.boolean(A=cContactPlugsInnerBase, B=cContactPlugsOuterBase, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cContactPlugs = RemoveDegenerateSquares(component=cContactPlugsBase, size=ContactPlugSize, layer=layer_Si)



    cInnerMetal2 = gf.boolean(A=cRingInnerOSlabCircle, B=cInnerM2Contact, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cInnerMetal2 = gf.boolean(A=cInnerMetal2, B=cRingInnerESlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cInnerMetal2 = gf.boolean(A=cInnerMetal2, B=cInnerM2CutAwayTinyGap, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cInnerMetal2 = gf.boolean(A=cInnerMetal2, B=cInnerM2CutAway, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cInnerMetal2 = BevelComponent(cInnerMetal2, MinAngle, BevelLength, layer=layer_Si)

    #cInnerMetal2PerfBlock = cInnerMetal2

    cOuterMetal2 = gf.boolean(A=cRingOuterESlabCircle, B=cOuterESlabRectangle, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterMetal2 = gf.boolean(A=cOuterMetal2, B=cOuterM2Contact, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterMetal2 = gf.boolean(A=cOuterMetal2, B=cRingOuterOSlabCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterMetal2 = gf.boolean(A=cOuterMetal2, B=cOuterESlabInnerRectangle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterMetal2 = gf.boolean(A=cOuterMetal2, B=cOuterM2CutAway, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cOuterMetal2 = BevelComponent(cOuterMetal2, MinAngle, BevelLength, layer=layer_Si)

    #cOuterMetal2PerfBlock = cOuterMetal2

    cHeaterWire = gf.boolean(A=cRingOuterHeaterCircle, B=cHeaterLegsOuterRect, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cHeaterWire = gf.boolean(A=cHeaterWire, B=cHeaterContactRect, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cHeaterWire = gf.boolean(A=cHeaterWire, B=cHeaterLegsInnerRect, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cHeaterWire = gf.boolean(A=cHeaterWire, B=cRingInnerHeaterCircle, operation='-', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    cHeaterMetal1 = gf.boolean(A=cHeaterMetal1Left, B=cHeaterMetal1Right, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cHeaterMetal2 = gf.boolean(A=cHeaterMetal2Left, B=cHeaterMetal2Right, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)

    cHeaterM1WireConnection = gf.boolean(A=cHeaterMetal1, B=cHeaterWire, operation='and', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)




    # make grid of possible Vias
    cM2BBox = gf.boolean(A=cInnerMetal2, B=cOuterMetal2, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cM2BBox = gf.boolean(A=cM2BBox, B=cHeaterMetal2, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    #ElectricsBBox=cOuterMetal2.bbox()
    ElectricsBBox=cM2BBox.bbox()
    #print(ElectricsBBox)
    ViaPatternWidth=2*max([abs(ElectricsBBox.left), abs(ElectricsBBox.right)])
    ViaPatternHeight=2*max([abs(ElectricsBBox.top), abs(ElectricsBBox.bottom)])
    #print(ViaPatternWidth,'  ',ViaPatternHeight)
    cSquareViaPattern = GenerateSquarePattern(ViaSize, ViaPitch, ViaPatternWidth, ViaPatternHeight, ViaHex, ViaPatternYOffset, layer_Si )
    #PlotComponent(component=cSquareViaPattern, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 0.0), 0.6), hatch='')


    # calculate Vias
    cInnerViaRegionCalc = gf.boolean(A=cInnerMetal1, B=cInnerMetal2, operation='and', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    rInnerViaRegion = cInnerViaRegionCalc.get_region(layer=layer_Si)  # Regions are in DBU (1 nm in this case).
    rInnerViaRegion = rInnerViaRegion.sized(-MetalViaBorderWidth/scale)  # Regions are in DBU.
    cInnerViaRegion = gf.Component()
    cInnerViaRegion.add_polygon(rInnerViaRegion, layer=layer_Si)  # Add the region to the component.

    cOuterViaRegionCalc = gf.boolean(A=cOuterMetal1, B=cOuterMetal2, operation='and', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    rOuterViaRegion = cOuterViaRegionCalc.get_region(layer=layer_Si)  # Regions are in DBU (1 nm in this case).
    rOuterViaRegion = rOuterViaRegion.sized(-MetalViaBorderWidth/scale)  # Regions are in DBU.
    cOuterViaRegion = gf.Component()
    cOuterViaRegion.add_polygon(rOuterViaRegion, layer=layer_Si)  # Add the region to the component.

    cViaRegion = gf.boolean(A=cInnerViaRegion, B=cOuterViaRegion, operation='or', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cViasBase = gf.boolean(A=cSquareViaPattern, B=cViaRegion, operation='and', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cVias = RemoveDegenerateSquares(component=cViasBase, size=ViaSize, layer=layer_Si)



    
    # calculate Heater Vias
    cHeaterViaRegionCalc = gf.boolean(A=cHeaterMetal1, B=cHeaterMetal2, operation='and', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    rHeaterViaRegion = cHeaterViaRegionCalc.get_region(layer=layer_Si)  # Regions are in DBU (1 nm in this case).
    rHeaterViaRegion = rHeaterViaRegion.sized(-MetalViaBorderWidth/scale)  # Regions are in DBU.
    cHeaterViaRegion = gf.Component()
    cHeaterViaRegion.add_polygon(rHeaterViaRegion, layer=layer_Si)  # Add the region to the component.

    cHeaterViasBase = gf.boolean(A=cSquareViaPattern, B=cHeaterViaRegion, operation='and', layer1=layer_Si, layer2=layer_Si, layer=layer_Si)
    cHeaterVias = RemoveDegenerateSquares(component=cHeaterViasBase, size=ViaSize, layer=layer_Si)



    """
    c2=gf.Component()
    test=cOuterContactDoping.get_polygons_points(by='tuple', layers=[layer_Si])
    for key in test:
        print(key)
        for polygon in test[key]:
        #if 1:
            #polygon = test[key][1]
            #print(polygon)
            #fastAxes.plot(polygon[:,0], polygon[:,1], 'r.')
            newPolygon = BevelContour(polygon, MinAngle, BevelLength)
            #fastAxes.plot(newPolygon[:,0], newPolygon[:,1], 'bo')
            #print(newPolygon)
            #print("-----")
            c2.add_polygon(newPolygon, layer=layer_Si)
    cOuterContactDoping=c2
    """



    #print(test[46][0])

    if (global_vars['ShowCladding'].get()==True):
        PlotComponent(component=cAllCladd, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.75, 0.75, 0.75), 0.5), hatch='......')

    if (global_vars['ShowSilicon'].get()==True):
        PlotComponent(component=cAllSi, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 0.2), 0.3), hatch='......')

    #cRingRib = gf.boolean(A=cRingRibOuterCircle, B=cRingRibInnerCircle, operation='-', layer=layer_Si)
    #PlotComponent(component=cRingRib, facecolor=((1.0, 0.0, 0.4), 0.2), edgecolor=('green', 0.5), hatch='')
    if (global_vars['ShowOSlabs'].get()==True):
        PlotComponent(component=cOuterOSlab, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.1, 0.0), 0.25), hatch='.....')
        PlotComponent(component=cInnerOSlab, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.5, 0.05, 0.0), 0.25), hatch='.....')
    if (global_vars['ShowWGRib'].get()==True):
        PlotComponent(component=cWGRibs, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.5, 0.0), 0.5), hatch='......')

    if (global_vars['ShowLowDoping'].get()==True):
        PlotComponent(component=cInnerLowDoping, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 0.0), 0.3), hatch='////')
        PlotComponent(component=cOuterLowDoping, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 1.0), 0.3), hatch='////')

    if (global_vars['ShowHighDoping'].get()==True):
        PlotComponent(component=cInnerHighDoping, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 0.0), 0.3), hatch='//////')
        PlotComponent(component=cOuterHighDoping, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 1.0), 0.3), hatch='//////')

    if (global_vars['ShowContactDoping'].get()==True):
        PlotComponent(component=cInnerContactDoping, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 0.0, 0.0), 0.3), hatch='\\\\\\\\\\\\')
        PlotComponent(component=cOuterContactDoping, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 1.0), 0.3), hatch='\\\\\\\\\\\\')
    
    if (global_vars['ShowMetal1'].get()==True):
        PlotComponent(component=cInnerMetal1, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.3, 1.0), 0.3), hatch='+++')
        PlotComponent(component=cOuterMetal1, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.3, 1.0), 0.3), hatch='+++')
        PlotComponent(component=cHeaterMetal1, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.3, 1.0), 0.3), hatch='+++')
    
    if (global_vars['ShowContactPlugRegion'].get()==True):
        PlotComponent(component=cInnerContactPlugsRegion, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 0.0), 0.6), hatch='')
        PlotComponent(component=cOuterContactPlugsRegion, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 0.0), 0.6), hatch='')
    
    if (global_vars['ShowContactPlugs'].get()==True):
        PlotComponent(component=cContactPlugs, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 0.0), 0.6), hatch='')

    if (global_vars['ShowMetal2'].get()==True):
        PlotComponent(component=cInnerMetal2, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.15, 0.5), 0.3), hatch='xxxxx')
        PlotComponent(component=cOuterMetal2, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.15, 0.5), 0.3), hatch='xxxxx')
        PlotComponent(component=cHeaterMetal2, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.15, 0.5), 0.3), hatch='xxxxx')

    if (global_vars['ShowMetalViaRegion'].get()==True):
        PlotComponent(component=cViaRegion, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 0.0), 0.6), hatch='')
        PlotComponent(component=cHeaterViaRegion, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 0.0), 0.6), hatch='')

    if (global_vars['ShowVias'].get()==True):
        PlotComponent(component=cVias, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 0.0), 0.6), hatch='')
        PlotComponent(component=cHeaterVias, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((0.0, 0.0, 0.0), 0.6), hatch='')

    if (global_vars['ShowHeater'].get()==True):
        PlotComponent(component=cHeaterWire, facecolor=((0.5, 0.1, 0.1), 0.7), edgecolor=((0.5, 0.1, 0.0), 0.7), hatch='')
        PlotComponent(component=cHeaterM1WireConnection, facecolor=((1.0, 1.0, 1.0), 0.0), edgecolor=((1.0, 1.0, 0.0), 0.5), hatch='xxx')


    plt.show(block=False)

    #print("Return from plt.show()...")

    


    #maincell = gf.Component("RM")
    maincell = gf.Component(maincellname)

    #m = maincell << .remap_layers({layer_Si: })
    #copyLayer(, layer_Si, )
    #cEmpty=gf.Component()
    # hier das eigentliche GDS2-Layout aufbauen
    mcAllCladd = maincell << copyLayer(cAllCladd, layer_Si, layer_Cladd)
    mcAllSi = maincell << copyLayer(cAllSi, layer_Si, layer_Si)
    mcOuterOSlab = maincell << copyLayer(cOuterOSlab, layer_Si, layer_OuterOSlabCladd)
    mcInnerOSlab = maincell << copyLayer(cInnerOSlab, layer_Si, layer_InnerOSlabCladd)
    mcWGRibs = maincell << copyLayer(cWGRibs, layer_Si, layer_OuterOSlabCore)
    mcInnerLowDoping = maincell << copyLayer(cInnerLowDoping, layer_Si, layer_InnerLowDope)
    mcOuterLowDoping= maincell << copyLayer(cOuterLowDoping, layer_Si, layer_OuterLowDope)
    mcInnerHighDoping = maincell << copyLayer(cInnerHighDoping, layer_Si, layer_InnerHiDope)
    mcOuterHighDoping = maincell << copyLayer(cOuterHighDoping, layer_Si, layer_OuterHiDope)
    mcInnerContactDoping = maincell << copyLayer(cInnerContactDoping, layer_Si, layer_InnerContactDope)
    mcOuterContactDoping = maincell << copyLayer(cOuterContactDoping, layer_Si, layer_OuterContactDope)
    mcInnerMetal1 = maincell << copyLayer(cInnerMetal1, layer_Si, layer_M1)
    mcInnerMetal1PerfBlock = maincell << copyLayer(cInnerMetal1, layer_Si, layer_M1PerfBlock)
    mcOuterMetal1 = maincell << copyLayer(cOuterMetal1, layer_Si, layer_M1)
    mcOuterMetal1PerfBlock = maincell << copyLayer(cOuterMetal1, layer_Si, layer_M1PerfBlock)
    mcHeaterMetal1 = maincell << copyLayer(cHeaterMetal1, layer_Si, layer_M1)
    mcHeaterMetal1PerfBlock = maincell << copyLayer(cHeaterMetal1, layer_Si, layer_M1PerfBlock)
    mcContactPlugs = maincell << copyLayer(cContactPlugs, layer_Si, layer_ContactPlugs)
    mcInnerMetal2 = maincell << copyLayer(cInnerMetal2, layer_Si, layer_M2)
    mcInnerMetal2PerfBlock = maincell << copyLayer(cInnerMetal2, layer_Si, layer_M2PerfBlock)
    mcOuterMetal2 = maincell << copyLayer(cOuterMetal2, layer_Si, layer_M2)
    mcOuterMetal2PerfBlock = maincell << copyLayer(cOuterMetal2, layer_Si, layer_M2PerfBlock)
    mcHeaterMetal2 = maincell << copyLayer(cHeaterMetal2, layer_Si, layer_M2)
    mcHeaterMetal2PerfBlock = maincell << copyLayer(cHeaterMetal2, layer_Si, layer_M2PerfBlock)
    mcVias = maincell << copyLayer(cVias, layer_Si, layer_Vias)
    mcHeaterVias = maincell << copyLayer(cHeaterVias, layer_Si, layer_Vias)
    mcHeaterWire = maincell << copyLayer(cHeaterWire, layer_Si, layer_Heater)
    mcHeaterM1WireConnection = maincell << copyLayer(cHeaterM1WireConnection, layer_Si, layer_HeaterContact)
    #m = maincell << copyLayer(, layer_Si, )
    #m = maincell << copyLayer(, layer_Si, )
    #m = maincell << copyLayer(, layer_Si, )

    # Ports...
    maincell.add_port(name="in0",  center=PortPosBusLeft, width=PortWGWidth, orientation=180, layer=layer_Si , port_type="optical")
    maincell.add_port(name="out0",  center=PortPosBusRight, width=PortWGWidth, orientation=0, layer=layer_Si , port_type="optical")
    maincell.add_port(name="in1",  center=PortPosDropRight, width=PortWGWidth, orientation=0, layer=layer_Si , port_type="optical")
    maincell.add_port(name="out1",  center=PortPosDropLeft, width=PortWGWidth, orientation=180, layer=layer_Si , port_type="optical")
   
    maincell.add_port(name="rf0",  center=PortPosRFOuter, width=ElectricalLeadWidth, orientation=180, layer=layer_M2 , port_type="electrical")
    maincell.add_port(name="rf1",  center=PortPosRFInner, width=ElectricalLeadWidth, orientation=0, layer=layer_M2 , port_type="electrical")

    maincell.add_port(name="dc0",  center=PortPosHeaterLeft, width=HeaterMetal1Overlap, orientation=90, layer=layer_M2 , port_type="electrical")
    maincell.add_port(name="dc1",  center=PortPosHeaterRight, width=HeaterMetal1Overlap, orientation=90, layer=layer_M2 , port_type="electrical")


    """
    #pointdict=cRing.get_polygons_points(merge=False)
    pointdict=cWGRingInnerCircle.get_polygons_points(merge=False)
    #print(WGRingInnerCircle)
    #print(pointdict)
    for layer, polygons in pointdict.items():
        for i, polygon in enumerate(polygons):
            plt.fill(polygon[:,0], polygon[:,1], facecolor=((1,0,0.4),0.2), edgecolor=('green', 1.0), hatch='O.')
    #plt.show()
    """

    return (maincell, ConfigString)

# --- End of CalculateGeometry ---



def MakePyProgram(c:gf.Component, filename):
    GDS_DBU = global_vars['GDS_DBU'].get()
    GDS_UU = global_vars['GDS_UU'].get()
    global scale
    scale = GDS_DBU/GDS_UU
    PortWGWidth = scale*global_vars['PortWGWidth'].get()
    ElectricalLeadWidth = scale*global_vars['ElectricalLeadWidth'].get()
    HeaterMetal1Overlap = scale*global_vars['HeaterMetal1Overlap'].get()


    PyProgramString = '# {filename}\n\n'
    PyProgramString += 'import gdsfactory as gf\n'
    PyProgramString += 'import conf\n'
    PyProgramString += 'import layers\n'
    PyProgramString += 'import xsections\n'
    PyProgramString += '\n'
    PyProgramString += '@gf.cell()\n'
    PyProgramString += 'def IPE_'+c.name+'():\n'
    PyProgramString += '    c = gf.Component()\n'
    PyProgramString += '    mygds = c.add_ref(gf.read.import_gds(gdspath=conf.Path_MyGDS+"\\'+filename+'"))\n'
    PyProgramString += '    c.add_port(name="in0" , center=('+str(c.ports['in0'].x)+', '+str(c.ports['in0'].y)+'), width='+str(PortWGWidth)+', orientation=180, layer=layers.layer_WG_COR, port_type="optical" )\n'
    PyProgramString += '    c.add_port(name="out0" , center=('+str(c.ports['out0'].x)+', '+str(c.ports['out0'].y)+'), width='+str(PortWGWidth)+', orientation=0, layer=layers.layer_WG_COR, port_type="optical" )\n'
    PyProgramString += '    c.add_port(name="in1" , center=('+str(c.ports['in1'].x)+', '+str(c.ports['in1'].y)+'), width='+str(PortWGWidth)+', orientation=0, layer=layers.layer_WG_COR, port_type="optical" )\n'
    PyProgramString += '    c.add_port(name="out1" , center=('+str(c.ports['out1'].x)+', '+str(c.ports['out1'].y)+'), width='+str(PortWGWidth)+', orientation=180, layer=layers.layer_WG_COR, port_type="optical" )\n'
    PyProgramString += '    c.add_port(name="rf0" , center=('+str(c.ports['rf0'].x)+', '+str(c.ports['rf0'].y)+'), width='+str(ElectricalLeadWidth)+', orientation=180, layer=layers.layer_M2_DRW, port_type="electrical" )\n'
    PyProgramString += '    c.add_port(name="rf1" , center=('+str(c.ports['rf1'].x)+', '+str(c.ports['rf1'].y)+'), width='+str(ElectricalLeadWidth)+', orientation=0, layer=layers.layer_M2_DRW, port_type="electrical" )\n'
    PyProgramString += '    c.add_port(name="dc0" , center=('+str(c.ports['dc0'].x)+', '+str(c.ports['dc0'].y)+'), width='+str(HeaterMetal1Overlap)+', orientation=90, layer=layers.layer_M2_DRW, port_type="electrical" )\n'
    PyProgramString += '    c.add_port(name="dc1" , center=('+str(c.ports['dc1'].x)+', '+str(c.ports['dc1'].y)+'), width='+str(HeaterMetal1Overlap)+', orientation=90, layer=layers.layer_M2_DRW, port_type="electrical" )\n'
    #PyProgramString += '\n'
    #PyProgramString += '\n'
    PyProgramString += '    return c\n'
    return PyProgramString







def on_closing():
    #print('Aaaaahhhhhh.....')
    plt.close('all')
    root.destroy()  # close the window
    root.quit() # stop the main loop



#def GraphWindowButtonPushed():
def OpenGraphWindow():
    global GraphWindowOpen
    global fastFig
    global fastAxes
    print('Opening graph window...')
    """
    if (GraphWindowOpen==False):
        # Create a new Toplevel window
        global Graph_window
        Graph_window = tk.Toplevel(root)
        Graph_window.title("Second Window")
        Graph_window.geometry("300x200")
        GraphWindowOpen=True
        
        # Add content to the second window
        label = tk.Label(Graph_window, text="This is the second window!")
        label.pack(pady=20)
    """
    """
    c = gf.Component()
    c.add_polygon([(-8, -6), (6, 8), (7, 17), (9, 5)], layer=(1, 0))
    c.add_ref(gf.components.analog.interdigital_capacitor(fingers=4, finger_length=20.0, finger_gap=2.0, thickness=5.0, layer='WG'))
    #c.plot()
    #c.show()

    pointdict=c.get_polygons_points(merge=False)
    #print(pointdict)


    x=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y=[random.randint(0, 10), random.randint(0, 10), random.randint(0, 10), random.randint(0, 10), random.randint(0, 10), random.randint(0, 10), random.randint(0, 10), random.randint(0, 10), random.randint(0, 10), random.randint(0, 10)]
    #print(x)
    #print(y)

    #plt.rcParams['figure.figsize'] = (12, 8)  # Set global default
    """

    existingFigures = plt.get_fignums()
    #print(existingFigures)
    if existingFigures==[]:
        fastFig, fastAxes = plt.subplots(nrows=1, ncols=1)
        fastFig.set_size_inches(7, 7)

    #111,figsize=(7,7))
    #fastFig.clear()
    fastAxes.clear()

    #fastAxes=fig.add_subplot(111)
    #...
    #ax.plot(...)
    #canvas.draw()

    #plt.cla()
    #plt.clf()
    #plt.xkcd()
    #plt.figure(figsize=(10, 6))  # Width=10 inches, Height=6 inches
    #plt.plot(x,y)

    
    """
    for layer, polygons in pointdict.items():
        #print(f"Layer: {layer}")
        for i, polygon in enumerate(polygons):
            #print(f"  Polygon {i}:")
            #print(f"    Vertices (shape {polygon.shape}):")
            #print(polygon)  # This is a numpy array of shape (N, 2)
            #plt.fill(polygon[:,0], polygon[:,1], facecolor=('red',0.2), edgecolor=('green', 1.0))
            plt.fill(polygon[:,0], polygon[:,1], facecolor=((1,0,0.4),0.2), edgecolor=('green', 1.0), hatch='O.')
    """
    """
    for layer, polys in pointdict.items():
        for poly in polys:
            for vertex in poly:
                x, y = vertex[0], vertex[1]
                # Do something with x, y
                print(x,y)
    """

    fastAxes.axis('equal')
    fastAxes.set_xlabel('x-dimension (µm)')
    fastAxes.set_ylabel('y-dimension (µm)')

    fastAxes.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    fastAxes.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    plt.minorticks_on()
    #plt.grid()

    
    #fastAxes.set_title('Ring Modulator')
    fastFig.suptitle('Ring Modulator')
    #fig = plt.gcf()  # Get current figure
    #fig.subplots_adjust(left=0.10, right=0.99, bottom=0.10, top=0.95)
    
    #fastFig.set_size_inches(7, 7)
    fastFig.tight_layout()
    #plt.ion()
    #plt.show()





def CalculateButtonPushed():
    print('Calculate button pushed:')
    #print('GDS database units = ', global_vars['GDS_DBU'].get())
    #print('GDS user units = ', global_vars['GDS_UU'].get())
    #print('')
    #for key,value in global_vars.items():
    #    print(f"{key}: ", global_vars[str(key)].get())

    #GraphWindowButtonPushed()
    OpenGraphWindow()
    CalculateGeometry(maincellname="RM")


def WriteGDS2ButtonPushed():
    print('Write GDS2 button pushed:')
    #print(ConfigString)

    #Dateinamen abfragen
    global savepath
    if (savepath==""):
        fullscriptnameandpath = os.path.abspath(__file__)
        fullscriptpath = os.path.dirname(fullscriptnameandpath) 
        savepath=fullscriptpath
    savefilename=tkf.asksaveasfilename(parent=root, initialdir=savepath, defaultextension="gds")
    if savefilename!="":
        print(savefilename)
        file=pathlib.Path(savefilename)
        gdsfilename = file.name
        gdsmaincellname = os.path.splitext(gdsfilename)[0]
        savepath = str(file.parent)
        print(savepath, "   ", gdsfilename)
        configfilename = gdsfilename+"_config.txt"
        programfilename = gdsfilename+"_component.py"

        OpenGraphWindow()
        maincell, ConfigString = CalculateGeometry(maincellname=gdsmaincellname)

        configfile = open(savepath+"/"+configfilename, "wt")
        configfile.write(ConfigString)
        configfile.close
        # ---------------------------------------------
        # hier noch die eigentliche GDS-Datei speichern
        # ---------------------------------------------
        #maincell.write_gds(gdsfilename, with_metadata=False)
        maincell.write_gds(gdsfilename, with_metadata=True)

        # start to write a small python-script-file for gdsfactory to load and show the RM
        PyProgramString = MakePyProgram(maincell, programfilename)
        programfile = open(savepath+"/"+programfilename, "wt")
        programfile.write(PyProgramString)
        programfile.close



    else:
        print("Not saved!")


    #OpenGraphWindow()
    #CalculateGeometry()
    #print(ConfigString)
    #options = kf.kcell.save_layout_options(write_context_info=False)
    #maincell.write_gds(r"dummy.gds", with_metadata=False,)




def LoadConfigButtonPushed():
    VariablesSection=False
    LayersSection=False

    print('Load Config button pushed:')
    #print(ConfigString)

    #Dateinamen abfragen
    global configpath
    if (configpath==""):
        fullscriptnameandpath = os.path.abspath(__file__)
        fullscriptpath = os.path.dirname(fullscriptnameandpath) 
        configpath=fullscriptpath
    configfilename=tkf.askopenfilename(parent=root, initialdir=configpath, defaultextension="txt")
    if configfilename!="":
        print(configfilename)
        configfile = open(configfilename, "rt")

        line = configfile.readline()        
        while line != '':
            #print(line, end='')

            # do something...
            if (line.find("[Variables]")>-1):
                VariablesSection=True
                LayersSection=False
            if (line.find("[Layers]")>-1):
                VariablesSection=False
                LayersSection=True
            if VariablesSection or LayersSection:
                EqualSignPos=line.find("=")
                if EqualSignPos>0:  #at least one character in front of equal sign
                    #read name
                    varname=line[0:EqualSignPos]
                    varvalue=line[EqualSignPos+1:-1]    #don't read the \n at the end of the line
                    print(varname, " = ", varvalue)
                    if varname in global_vars:
                        IntValueNames=["ContactPlugHex", "ViaHex", "ShowCladding", "ShowSilicon", "ShowWGRib",
                                    "ShowOSlabs", "ShowLowDoping", "ShowHighDoping", "ShowContactDoping", "ShowContactPlugRegion",
                                    "ShowContactPlugs", "ShowMetal1", "ShowMetal2", "ShowMetalViaRegion", "ShowVias", "ShowHeater", "ShowConstructionLines"]
                        if varname in IntValueNames:
                            global_vars[varname].set(int(varvalue))
                        else:
                            global_vars[varname].set(float(varvalue))
                    if varname in global_layers:
                        commapos=varvalue.find(",")
                        varvalue1 = varvalue[1:commapos]
                        varvalue2 = varvalue[commapos+1:-1]
                        global_layers[varname][0].set(int(varvalue1))
                        global_layers[varname][1].set(int(varvalue2))


            line = configfile.readline()
        configfile.close

    else:
        print("Nothing read!")







class MainWindow:
 


    def __init__(self, root):
        self.root = root
        self.root.title("RM-Layout Program")
        self.root.geometry("800x1000")
        self.root.minsize(300,300)

        KITLogoImg = ImageTk.PhotoImage(Image.open("KIT_Logo_XXXS.png"))
        IPELogoImg = ImageTk.PhotoImage(Image.open("IPE-Logo_XXXS.png"))
        
        # Configure grid weights for the root window
        self.root.grid_rowconfigure(1, weight=1)  # Middle frame takes remaining space
        self.root.grid_columnconfigure(0, weight=1)  # Full width
        
        # === Create top frame with fixed height
        self.top_logo_frame = tk.Frame(self.root, height=85)
        self.top_logo_frame.grid(row=0, column=0, sticky="ew")
        self.top_logo_frame.grid_propagate(False)  # Prevent frame from resizing based on content
        
        # Add label widgets to top logo frame
        program_name = tk.Label(self.top_logo_frame, text="Ring Modulator Layout Program 2", font=("Arial", 14, "bold"))
        program_name.place(x=10, y=10)
        
        program_author = tk.Label(self.top_logo_frame, text="by Dr.-Ing. Marc Schneider (KIT)", font=("Arial", 10))
        program_author.place(x=10, y=50)

        program_version = tk.Label(self.top_logo_frame, text=ProgramVersionString+" from "+ProgramVersionDateString, font=("Arial", 10))
        program_version.place(x=275, y=50)

        program_KITLogo = tk.Label(self.top_logo_frame, image=KITLogoImg)
        program_KITLogo.image = KITLogoImg  # don't know, why, but this anchors the image (https://stackoverflow.com/questions/13148975/tkinter-label-does-not-show-image)
        program_KITLogo.place(x=660, y=20)

        program_IPELogo = tk.Label(self.top_logo_frame, image=IPELogoImg)
        program_IPELogo.image = IPELogoImg
        program_IPELogo.place(x=460, y=20)
        
        # === Create (middle) main frame that spans remaining space
        self.main_frame = tk.Frame(self.root)
        self.main_frame.grid(row=1, column=0, sticky="nsew")
        
        # Configure bottom frame grid
        self.main_frame.grid_rowconfigure(0, weight=1)
        self.main_frame.grid_columnconfigure(0, weight=1)
        
        # Create notebook in bottom frame
        self.notebook = ttk.Notebook(self.main_frame)
        self.notebook.grid(row=0, column=0, sticky="nsew")
        
        # Add some sample tabs to the notebook
        self.create_notebook_tabs()
            
        
        # === Create bottom frame (fixed height)
        self.bottom_frame = tk.Frame(self.root, height=70)
        self.bottom_frame.grid(row=2, column=0, sticky="ew")
        self.bottom_frame.grid_propagate(False)  # Prevent resizing based on content
        
        # Add content to bottom frame
        calculate_button = ttk.Button(master=self.bottom_frame, text='     Calculate     ', bootstyle='success', padding=5, command=CalculateButtonPushed)
        calculate_button.pack(side='left', padx=10, pady=10)
        save_button = ttk.Button(master=self.bottom_frame, text='     Write GDS2     ', bootstyle='primary', padding=5, command=WriteGDS2ButtonPushed)
        save_button.pack(side='left', padx=10, pady=10)
        loadConfig_button = ttk.Button(master=self.bottom_frame, text='    Load Config.    ', bootstyle='info', padding=5, command=LoadConfigButtonPushed)
        loadConfig_button.pack(side='left', padx=10, pady=10)
        #graphwindow_button = ttk.Button(master=self.bottom_frame, text=' Open graph window ', bootstyle='secondary', padding=5, command=GraphWindowButtonPushed)
        #graphwindow_button.pack(side='left', padx=10, pady=10)
        #calculate_button.place(x=10, y=20)
        #status_label = tk.Label(self.bottom_frame, text="Status: Ready", font=("Arial", 10), bg="lightgreen")
        #status_label.pack(pady=10)
        
        


    def create_notebook_tabs(self):

        def create_checkbutton_element(tab, var_name, label_text, initial_value, row, column=0, command=None):
            """
            Creates a label and a checkbutton for a drawing parameter.
            
            :param tab: parent widget (e.g., tab_Drawing)
            :param var_name: string name of the variable (e.g., 'width', 'height')
            :param label_text: text for the label (e.g., 'Width (nm):')
            :param initial_value: initial value for IntVar
            :param row: row in grid
            :param column: optional column (default 0 for label, 1 for checkbutton)
            """
            # Create label
            label = tk.Label(tab, text=label_text)
            label.grid(column=column, row=row, sticky='e')
            
            # Create IntVar and store in dictionary
            global_vars[var_name] = ttk.IntVar(value=initial_value)
            
            # Create checkbutton
            check = ttk.Checkbutton(tab,variable=global_vars[var_name],style='Roundtoggle.Toolbutton', command=command)
            check.grid(column=column+1, row=row, padx=5, pady=10, sticky='w')

            return label, check  # optional: return widgets if you need to reference them later


        def create_entry_element(tab, var_name, label_text, initial_value, row, column=0):
            """
            Creates a label and a checkbutton for a drawing parameter.
            
            :param tab: parent widget (e.g., tab_Drawing)
            :param var_name: string name of the variable (e.g., 'width', 'height')
            :param label_text: text for the label (e.g., 'Width (nm):')
            :param initial_value: initial value for DoubleVar
            :param row: row in grid
            :param column: optional column (default 0 for label, 1 for checkbutton)
            """
            # Create label
            label = tk.Label(tab, text=label_text)
            label.grid(column=column, row=row, sticky='e')
            
            # Create DoubleVar and store in dictionary
            global_vars[var_name] = ttk.DoubleVar(value=initial_value)
            
            # Create entry
            entry = ttk.Entry(tab,textvariable=global_vars[var_name], width=10, justify='right')
            entry.grid(column=column+1, row=row, padx=5, pady=2)

            return label, entry  # optional: return widgets if you need to reference them later



        def create_doublelayerentry_element(tab, var_name, label_text, initial_value, row, column=0):
            """
            Creates a label and a checkbutton for a drawing parameter.
            
            :param tab: parent widget (e.g., tab_Drawing)
            :param var_name: string name of the variable (e.g., 'width', 'height')
            :param label_text: text for the label (e.g., 'Width (nm):')
            :param initial_value: tupel of initial value for IntVar
            :param row: row in grid
            :param column: optional column (default 0 for label, 1 for checkbutton)
            """
            # Create label
            label = tk.Label(tab, text=label_text)
            label.grid(column=column, row=row, sticky='e')
            
            # Create DoubleVar and store in dictionary
            global_layers[var_name] = (ttk.IntVar(value=initial_value[0]), ttk.IntVar(value=initial_value[1]))
            
            # Create entry
            entry1 = ttk.Entry(tab,textvariable=global_layers[var_name][0], width=10, justify='right')
            entry1.grid(column=column+1, row=row, padx=5, pady=2)
            entry2 = ttk.Entry(tab,textvariable=global_layers[var_name][1], width=10, justify='right')
            entry2.grid(column=column+2, row=row, padx=5, pady=2)

            return label, entry1, entry2  # optional: return widgets if you need to reference them later





        # Create tabs
        tab_Basics = tk.Frame(self.notebook, bg="white", padx=10, pady=10)
        tab_Doping = tk.Frame(self.notebook, bg="white", padx=10, pady=10)
        tab_Metallization = tk.Frame(self.notebook, bg="white", padx=10, pady=10)
        tab_Heater = tk.Frame(self.notebook, bg="white", padx=10, pady=10)
        tab_Drawing = tk.Frame(self.notebook, bg='white', padx=10, pady=10)
        tab_GDS = tk.Frame(self.notebook, bg="white", padx=10, pady=10)

        #tab_GDS.rowconfigure()
        
        # Add content to tabs

        # fill Basics tab
        #tk.Label(tab_Basics, text="This is Tab 1", font=("Arial", 14)).pack(pady=50)
        crow=0
        Basics_CladdingWidth_label,Basics_CladdingWidth_entry = create_entry_element(tab_Basics, var_name = 'CladdingWidth', label_text='Cladding width (nm):', initial_value = CladdingWidth_ini, row=crow)
        crow+=1
        Basics_RingDiscretization_label,Basics_RingDiscretization_entry = create_entry_element(tab_Basics, var_name = 'RingDiscretization', label_text='Ring discretization:', initial_value = RingDiscretization_ini, row=crow)
        crow+=1
        Basics_MinAngle_label,Basics_MinAngle_entry = create_entry_element(tab_Basics, var_name = 'MinAngle', label_text='Minimum angle (°):', initial_value = MinAngle_ini, row=crow)
        crow+=1
        Basics_BevelLength_label,Basics_BevelLength_entry = create_entry_element(tab_Basics, var_name = 'BevelLength', label_text='Bevel length (nm):', initial_value = BevelLength_ini, row=crow)

        crow+=1
        tk.Label(tab_Basics, text=' ').grid(column=0, row=crow, sticky='e')

        crow+=1
        Basics_RingDiameter_label,Basics_RingDiameter_entry = create_entry_element(tab_Basics, var_name = 'RingDiameter', label_text='Ring diameter (nm):', initial_value = RingDiameter_ini, row=crow)
        crow+=1
        Basics_RingWGWidth_label,Basics_RingWGWidth_entry = create_entry_element(tab_Basics, var_name = 'RingWGWidth', label_text='Ring waveguide width (nm):', initial_value = RingWGWidth_ini, row=crow)
        crow+=1
        Basics_BusGap_label,Basics_BusGap_entry = create_entry_element(tab_Basics, var_name = 'BusGap', label_text='Bus-ring gap (nm):', initial_value = BusGap_ini, row=crow)
        crow+=1
        Basics_DropGap_label,Basics_DropGap_entry = create_entry_element(tab_Basics, var_name = 'DropGap', label_text='Drop-ring gap (nm):', initial_value = DropGap_ini, row=crow)
        crow+=1
        Basics_BusWGWidth_label,Basics_BusWGWidth_entry = create_entry_element(tab_Basics, var_name = 'BusWGWidth', label_text='Bus waveguide width (nm):', initial_value = BusWGWidth_ini, row=crow)
        crow+=1
        Basics_DropWGWidth_label,Basics_DropWGWidth_entry = create_entry_element(tab_Basics, var_name = 'DropWGWidth', label_text='Drop waveguide width (nm):', initial_value = DropWGWidth_ini, row=crow)

        crow+=1
        tk.Label(tab_Basics, text=' ').grid(column=0, row=crow, sticky='e')

        crow+=1
        Basics_OSlabWidth_label,Basics_OSlabWidth_entry = create_entry_element(tab_Basics, var_name = 'OSlabWidth', label_text='Optics slab width (nm):', initial_value = OSlabWidth_ini, row=crow)
        crow+=1
        Basics_BusWGLength_label,Basics_BusWGLength_entry = create_entry_element(tab_Basics, var_name = 'BusWGLength', label_text='Bus waveguide length (nm):', initial_value = BusWGLength_ini, row=crow)
        crow+=1
        Basics_DropWGLength_label,Basics_DropWGLength_entry = create_entry_element(tab_Basics, var_name = 'DropWGLength', label_text='Drop waveguide length (nm):', initial_value = DropWGLength_ini, row=crow)
        crow+=1
        Basics_TaperLength_label,Basics_TaperLength_entry = create_entry_element(tab_Basics, var_name = 'TaperLength', label_text='Taper length (nm):', initial_value = TaperLength_ini, row=crow)
        crow+=1
        Basics_SlabBorderExtensionOnTapers_label,Basics_SlabBorderExtensionOnTapers_entry = create_entry_element(tab_Basics, var_name = 'SlabBorderExtensionOnTapers', label_text='OSlab border extension on Tapers (nm):', initial_value = SlabBorderExtensionOnTapers_ini, row=crow)
        crow+=1
        Basics_PortWGWidth_label,Basics_PortWGWidth_entry = create_entry_element(tab_Basics, var_name = 'PortWGWidth', label_text='Port waveguide width (nm):', initial_value = PortWGWidth_ini, row=crow)

        crow+=1
        tk.Label(tab_Basics, text=' ').grid(column=0, row=crow, sticky='e')

        crow+=1
        Basics_ESlabWidth_label,Basics_ESlabWidth_entry = create_entry_element(tab_Basics, var_name = 'ESlabWidth', label_text='Electrical contact slab width (nm):', initial_value = ESlabWidth_ini, row=crow)
        crow+=1
        Basics_OuterESlabLength_label,Basics_OuterESlabLength_entry = create_entry_element(tab_Basics, var_name = 'OuterESlabLength', label_text='Outer el. contact slab length (nm):', initial_value = OuterESlabLength_ini, row=crow)

        #crow+=1
        #Basics_xxx_label,Basics_xxx_entry = create_entry_element(tab_Basics, var_name = 'xxx', label_text='Blblub', initial_value = xxx_ini, row=crow)

        # fill Basics tab end


        # fill Doping tab
        #tk.Label(tab_Doping, text="This is Tab 2", font=("Arial", 14)).pack(pady=50)

        crow=0
        Doping_JunctionOffset_label,Doping_JunctionOffset_entry = create_entry_element(tab_Doping, var_name = 'JunctionOffset', label_text='Junction offset (+ larger radius) (nm):', initial_value = JunctionOffset_ini, row=crow)
        crow+=1
        Doping_HighDopingDistanceInside_label,Doping_HighDopingDistanceInside_entry = create_entry_element(tab_Doping, var_name = 'HighDopingDistanceInside', label_text='High doping distance to WG inside ring (nm):', initial_value = HighDopingDistanceInside_ini, row=crow)
        crow+=1
        Doping_HighDopingDistanceOutside_label,Doping_HighDopingDistanceOutside_entry = create_entry_element(tab_Doping, var_name = 'HighDopingDistanceOutside', label_text='High doping distance to WG outside ring (nm):', initial_value = HighDopingDistanceOutside_ini, row=crow)

        crow+=1
        tk.Label(tab_Doping, text=' ').grid(column=0, row=crow, sticky='e')

        crow+=1
        Doping_ContactPlugRegionBorderWidth_label,Doping_ContactPlugRegionBorderWidth_entry = create_entry_element(tab_Doping, var_name = 'ContactPlugRegionBorderWidth', label_text='Contact plug border width (nm):', initial_value = ContactPlugRegionBorderWidth_ini, row=crow)
        crow+=1
        Doping_ContactPlugSize_label,Doping_ContactPlugSize_entry = create_entry_element(tab_Doping, var_name = 'ContactPlugSize', label_text='Contact plug size (nm):', initial_value = ContactPlugSize_ini, row=crow)
        crow+=1
        Doping_ContactPlugPitch_label,Doping_ContactPlugPitch_entry = create_entry_element(tab_Doping, var_name = 'ContactPlugPitch', label_text='Contact plug pitch (nm):', initial_value = ContactPlugPitch_ini, row=crow)

        crow+=1
        Doping_ContactPlugHex_label,Doping_ContactPlugHex_check=create_checkbutton_element(tab_Doping, var_name='ContactPlugHex', label_text = 'Hexagonal contact plug pattern ', initial_value = ContactPlugHex_ini, row=crow)

        crow+=1
        Doping_ContactPlugPatternYOffset_label,Doping_ContactPlugPatternYOffset_entry = create_entry_element(tab_Doping, var_name = 'ContactPlugPatternYOffset', label_text='Contact plug pattern y-offset (nm):', initial_value = ContactPlugPatternYOffset_ini, row=crow)


        #crow+=1
        #Doping_xxx_label,Doping_xxx_entry = create_entry_element(tab_Doping, var_name = 'xxx', label_text='----', initial_value = xxx_ini, row=crow)

        # fill Doping tab end


        # fill Metallization tab
        #tk.Label(tab_Metallization, text="This is Tab 3", font=("Arial", 14)).pack(pady=50)


        crow=0
        Metallization_ElectricalLeadProtrusion_label,Metallization_ElectricalLeadProtrusion_entry = create_entry_element(tab_Metallization, var_name = 'ElectricalLeadProtrusion', label_text='Electrical lead protrusion (nm):', initial_value = ElectricalLeadProtrusion_ini, row=crow)
        crow+=1
        Metallization_ElectricalLeadWidth_label,Metallization_ElectricalLeadWidth_entry = create_entry_element(tab_Metallization, var_name = 'ElectricalLeadWidth', label_text='Electrical lead width (nm):', initial_value = ElectricalLeadWidth_ini, row=crow)
        crow+=1
        Metallization_ElectricalLeadGap_label,Metallization_ElectricalLeadGap_entry = create_entry_element(tab_Metallization, var_name = 'ElectricalLeadGap', label_text='Inner electrical lead - outer contact ring gap (nm):', initial_value = ElectricalLeadGap_ini, row=crow)
        crow+=1
        Metallization_InnerElectricalRingCutawayOffset_label,Metallization_InnerElectricalRingCutawayOffset_entry = create_entry_element(tab_Metallization, var_name = 'InnerElectricalRingCutawayOffset', label_text='Inner metal ring cut-away offset (nm):', initial_value = InnerElectricalRingCutawayOffset_ini, row=crow)
        crow+=1
        Metallization_MetalViaBorderWidth_label,Metallization_MetalViaBorderWidth_entry = create_entry_element(tab_Metallization, var_name = 'MetalViaBorderWidth', label_text='Metal via border width (nm):', initial_value = MetalViaBorderWidth_ini, row=crow)
        crow+=1
        Metallization_ViaSize_label,Metallization_ViaSize_entry = create_entry_element(tab_Metallization, var_name = 'ViaSize', label_text='Via size (nm):', initial_value = ViaSize_ini, row=crow)
        crow+=1
        Metallization_ViaPitch_label,Metallization_ViaPitch_entry = create_entry_element(tab_Metallization, var_name = 'ViaPitch', label_text='Via pitch (nm):', initial_value = ViaPitch_ini, row=crow)

        crow+=1
        Metallization_ViaHex_label,Metallization_ViaHex_check=create_checkbutton_element(tab_Metallization, var_name='ViaHex', label_text = 'Hexagonal via pattern ', initial_value = ViaHex_ini, row=crow)

        crow+=1
        Metallization_ViaPatternYOffset_label,Metallization_ViaPatternYOffset_entry = create_entry_element(tab_Metallization, var_name = 'ViaPatternYOffset', label_text='Via pattern y-offset (nm):', initial_value = ViaPatternYOffset_ini, row=crow)

        #crow+=1
        #Metallization_xxx_label,Metallization_xxx_entry = create_entry_element(tab_Metallization, var_name = 'xxx', label_text='----', initial_value = xxx_ini, row=crow)

        # fill Metallization tab end


        # fill Heater tab
        #tk.Label(tab_Heater, text="This is Tab 4", font=("Arial", 14)).pack(pady=50)


        crow=0
        Heater_HeaterWidth_label,Heater_HeaterWidth_entry = create_entry_element(tab_Heater, var_name = 'HeaterWidth', label_text='Heater width (nm):', initial_value = HeaterWidth_ini, row=crow)
        crow+=1
        Heater_HeaterLegSpacing_label,Heater_HeaterLegSpacing_entry = create_entry_element(tab_Heater, var_name = 'HeaterLegSpacing', label_text='Heater leg spacing (nm):', initial_value = HeaterLegSpacing_ini, row=crow)
        crow+=1
        Heater_HeaterContactExtension_label,Heater_HeaterContactExtension_entry = create_entry_element(tab_Heater, var_name = 'HeaterContactExtension', label_text='Heater contact extension (nm):', initial_value = HeaterContactExtension_ini, row=crow)
        crow+=1
        Heater_HeaterMetal1Overlap_label,Heater_HeaterMetal1Overlap_entry = create_entry_element(tab_Heater, var_name = 'HeaterMetal1Overlap', label_text='Metal 1 overlap (nm):', initial_value = HeaterMetal1Overlap_ini, row=crow)
        crow+=1
        Heater_HeaterMetal1Protrusion_label,Heater_HeaterMetal1Protrusion_entry = create_entry_element(tab_Heater, var_name = 'HeaterMetal1Protrusion', label_text='Metal 1 protrusion (nm):', initial_value = HeaterMetal1Protrusion_ini, row=crow)
        crow+=1
        Heater_HeaterMetal2Protrusion_label,Heater_HeaterMetal2Protrusion_entry = create_entry_element(tab_Heater, var_name = 'HeaterMetal2Protrusion', label_text='Metal 2 protrusion (nm):', initial_value = HeaterMetal2Protrusion_ini, row=crow)

        #crow+=1
        #Heater_xxx_label,Heater_xxx_entry = create_entry_element(tab_Heater, var_name = 'xxx', label_text='----', initial_value = xxx_ini, row=crow)

        """
        crow=crow+1
        Heater_xxx_label = tk.Label(tab_Heater, text='Blub (nm):')
        Heater_xxx_label.grid(column=0, row=crow, sticky='e')
        global xxx
        xxx = ttk.DoubleVar(value=xxx_ini)
        Heater_xxx_entry = ttk.Entry(tab_Heater, textvariable=xxx, width=10, justify='right')
        Heater_xxx_entry.grid(column=1, row=crow, padx=5, pady=2)
        """
        # fill Heater tab end


        # fill GDS tab
        #tk.Label(tab_GDS, text="This is the GDS tab", font=("Arial", 14)).pack(pady=50)

        crow=1
        GDS_DBU_label,GDS_DBU_entry = create_entry_element(tab_GDS, var_name = 'GDS_DBU', label_text='GDS database units:', initial_value = GDS_DBU_ini, row=crow, column=4)
        tk.Label(tab_GDS, text='               ').grid(column=3, row=crow, sticky='e')
        crow+=1
        GDS_UU_label,GDS_UU_entry = create_entry_element(tab_GDS, var_name = 'GDS_UU', label_text='GDS user units:', initial_value = GDS_UU_ini, row=crow, column=4)
        #crow+=1
        #tk.Label(tab_GDS, text=' ').grid(column=0, row=crow, sticky='e')
        crow=0
        tk.Label(tab_GDS, text=' ').grid(column=0, row=crow, sticky='e')
        tk.Label(tab_GDS, text='Layer').grid(column=1, row=crow, sticky='e')
        tk.Label(tab_GDS, text='Data type').grid(column=2, row=crow, sticky='e')
        crow+=1
        GDS_layer_Cladd_label, GDS_layer_Cladd_entry1, GDS_layer_Cladd_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_Cladd', label_text='SiO2 Cladding', initial_value=layer_Cladd_ini, row=crow)
        crow+=1
        GDS_layer_Si_label, GDS_layer_Si_entry1, GDS_layer_Si_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_Si', label_text='Any Silicon', initial_value=layer_Si_ini, row=crow)
        crow+=1
        GDS_layer_OuterOSlabCore_label, GDS_layer_OuterOSlabCore_entry1, GDS_layer_OuterOSlabCore_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_OuterOSlabCore', label_text='Outer Optical Slab Core', initial_value=layer_OuterOSlabCore_ini, row=crow)
        crow+=1
        GDS_layer_OuterOSlabCladd_label, GDS_layer_OuterOSlabCladd_entry1, GDS_layer_OuterOSlabCladd_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_OuterOSlabCladd', label_text='Outer Optical Slab Cladding', initial_value=layer_OuterOSlabCladd_ini, row=crow)
        crow+=1
        GDS_layer_InnerOSlabCladd_label, GDS_layer_InnerOSlabCladd_entry1, GDS_layer_InnerOSlabCladd_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_InnerOSlabCladd', label_text='Inner Optical Slab Cladding', initial_value=layer_InnerOSlabCladd_ini, row=crow)
        crow+=1
        GDS_layer_OuterLowDope_label, GDS_layer_OuterLowDope_entry1, GDS_layer_OuterLowDope_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_OuterLowDope', label_text='Outer Low Doping', initial_value=layer_OuterLowDope_ini, row=crow)
        crow+=1
        GDS_layer_InnerLowDope_label, GDS_layer_InnerLowDope_entry1, GDS_layer_InnerLowDope_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_InnerLowDope', label_text='Inner Low Doping', initial_value=layer_InnerLowDope_ini, row=crow)
        crow+=1
        GDS_layer_OuterHiDope_label, GDS_layer_OuterHiDope_entry1, GDS_layer_OuterHiDope_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_OuterHiDope', label_text='Outer Hi Doping', initial_value=layer_OuterHiDope_ini, row=crow)
        crow+=1
        GDS_layer_InnerHiDope_label, GDS_layer_InnerHiDope_entry1, GDS_layer_InnerHiDope_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_InnerHiDope', label_text='Inner Hi Doping', initial_value=layer_InnerHiDope_ini, row=crow)
        crow+=1
        GDS_layer_OuterContactDope_label, GDS_layer_OuterContactDope_entry1, GDS_layer_OuterContactDope_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_OuterContactDope', label_text='Outer Contact Doping', initial_value=layer_OuterContactDope_ini, row=crow)
        crow+=1
        GDS_layer_InnerContactDope_label, GDS_layer_InnerContactDope_entry1, GDS_layer_InnerContactDope_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_InnerContactDope', label_text='Inner Contact Doping', initial_value=layer_InnerContactDope_ini, row=crow)
        crow+=1
        GDS_layer_M1_label, GDS_layer_M1_entry1, GDS_layer_M1_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_M1', label_text='Metal 1', initial_value=layer_M1_ini, row=crow)
        crow+=1
        GDS_layer_M1PerfBlock_label, GDS_layer_M1PerfBlock_entry1, GDS_layer_M1PerfBlock_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_M1PerfBlock', label_text='Metal 1 Perforation Blocking', initial_value=layer_M1PerfBlock_ini, row=crow)
        crow+=1
        GDS_layer_M2_label, GDS_layer_M2_entry1, GDS_layer_M2_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_M2', label_text='Metal 2', initial_value=layer_M2_ini, row=crow)
        crow+=1
        GDS_layer_M2PerfBlock_label, GDS_layer_M2PerfBlock_entry1, GDS_layer_M2PerfBlock_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_M2PerfBlock', label_text='Metal 2 Perforation Blocking', initial_value=layer_M2PerfBlock_ini, row=crow)
        crow+=1
        GDS_layer_Heater_label, GDS_layer_Heater_entry1, GDS_layer_Heater_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_Heater', label_text='Heater', initial_value=layer_Heater_ini, row=crow)
        crow+=1
        GDS_layer_HeaterContact_label, GDS_layer_HeaterContact_entry1, GDS_layer_HeaterContact_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_HeaterContact', label_text='Heater Contact', initial_value=layer_HeaterContact_ini, row=crow)
        crow+=1
        GDS_layer_Vias_label, GDS_layer_Vias_entry1, GDS_layer_Vias_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_Vias', label_text='VIAs Metal 1 - Metal 2', initial_value=layer_Vias_ini, row=crow)
        crow+=1
        GDS_layer_ContactPlugs_label, GDS_layer_ContactPlugs_entry1, GDS_layer_ContactPlugs_entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_ContactPlugs', label_text='Silicon Contact Plugs', initial_value=layer_ContactPlugs_ini, row=crow)



        #crow+=1
        #GDS_layer__label, GDS_layer__entry1, GDS_layer__entry2 = create_doublelayerentry_element(tab_GDS, var_name='layer_', label_text='', initial_value=layer__ini, row=crow)

        
        #crow+=1
        #GDS_xxx_label,GDS_xxx_entry = create_entry_element(tab_GDS, var_name = 'xxx', label_text='----', initial_value = xxx_ini, row=crow)

        # fill GDS tab end


        # fill Drawing tab
        #tk.Label(tab_Drawing, text="This is the drawing configuration tab", font=("Arial", 14)).pack(pady=50)

        crow=0
        Drawing_ShowCladding_label,Drawing_ShowCladding_check=create_checkbutton_element(tab_Drawing, var_name='ShowCladding', label_text = 'Show cladding', initial_value = ShowCladding_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowSilicon_label,Drawing_ShowSilicon_check=create_checkbutton_element(tab_Drawing, var_name='ShowSilicon', label_text = 'Show silicon', initial_value = ShowSilicon_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowWGRib_label,Drawing_ShowWGRib_check=create_checkbutton_element(tab_Drawing, var_name='ShowWGRib', label_text = 'Show waveguide rib', initial_value = ShowWGRib_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowOSlabs_label,Drawing_ShowOSlabs_check=create_checkbutton_element(tab_Drawing, var_name='ShowOSlabs', label_text = 'Show optical slabs', initial_value = ShowOSlabs_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowLowDoping_label,Drawing_ShowLowDoping_check=create_checkbutton_element(tab_Drawing, var_name='ShowLowDoping', label_text = 'Show low doping', initial_value = ShowLowDoping_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowHighDoping_label,Drawing_ShowHighDoping_check=create_checkbutton_element(tab_Drawing, var_name='ShowHighDoping', label_text = 'Show high doping', initial_value = ShowHighDoping_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowContactDoping_label,Drawing_ShowContactDoping_check=create_checkbutton_element(tab_Drawing, var_name='ShowContactDoping', label_text = 'Show contact doping', initial_value = ShowContactDoping_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowContactPlugRegion_label,Drawing_ShowContactPlugRegion_check=create_checkbutton_element(tab_Drawing, var_name='ShowContactPlugRegion', label_text = 'Show contact plug region', initial_value = ShowContactPlugRegion_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowContactPlugs_label,Drawing_ShowContactPlugs_check=create_checkbutton_element(tab_Drawing, var_name='ShowContactPlugs', label_text = 'Show contact plugs', initial_value = ShowContactPlugs_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowMetal1_label,Drawing_ShowMetal1_check=create_checkbutton_element(tab_Drawing, var_name='ShowMetal1', label_text = 'Show metal 1', initial_value = ShowMetal1_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowMetal2_label,Drawing_ShowMetal2_check=create_checkbutton_element(tab_Drawing, var_name='ShowMetal2', label_text = 'Show metal 2', initial_value = ShowMetal2_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowMetalViaRegion_label,Drawing_ShowMetalViaRegion_check=create_checkbutton_element(tab_Drawing, var_name='ShowMetalViaRegion', label_text = 'Show metal via region', initial_value = ShowMetalViaRegion_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowVias_label,Drawing_ShowVias_check=create_checkbutton_element(tab_Drawing, var_name='ShowVias', label_text = 'Show vias', initial_value = ShowVias_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowHeater_label,Drawing_ShowHeater_check=create_checkbutton_element(tab_Drawing, var_name='ShowHeater', label_text = 'Show heater', initial_value = ShowHeater_ini, row=crow, command=CalculateButtonPushed)
        crow+=1
        Drawing_ShowConstructionLines_label,Drawing_ShowConstructionLines_check=create_checkbutton_element(tab_Drawing, var_name='ShowConstructionLines', label_text = 'Show Construction Lines', initial_value = ShowConstructionLines_ini, row=crow, command=CalculateButtonPushed)

        #crow+=1
        #Drawing_xxx_label,Drawing_xxx_check=create_checkbutton_element(tab_Drawing, var_name='xxx', label_text = 'Show cladding', initial_value = xxx_ini, row=crow)

        """
        crow=crow+1
        Drawing_xxx_label = tk.Label(tab_Drawing, text='Blub (nm):')
        Drawing_xxx_label.grid(column=0, row=crow, sticky='e')
        global xxx
        xxx = ttk.IntVar(value=xxx_ini)
        Drawing_xxx_check = ttk.Checkbutton(tab_Drawing, variable=xxx, style='Roundtoggle.Toolbutton')
        Drawing_xxx_check.grid(column=1, row=crow, padx=5, pady=10, sticky='w')
        """
        # fill Drawing tab end




        # Add tabs to notebook
        self.notebook.add(tab_Basics, text="Basics")
        self.notebook.add(tab_Doping, text="Doping")
        self.notebook.add(tab_Metallization, text="Metallization")
        self.notebook.add(tab_Heater, text="Heater")
        self.notebook.add(tab_GDS, text="GDS")
        self.notebook.add(tab_Drawing, text="Drawing")

# Create and run the application
if __name__ == "__main__":
    root = tk.Tk()
    style=ttk.Style()
    style.theme_use("yeti")
    app = MainWindow(root)
    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()
    


"""
# window
window = tk.Tk()
window.title('RM-Layout Program')
window.geometry('800x800')
window.minsize(300,300)

window.grid_columnconfigure(0,weight=1)

# main layout widgets
#content = ttk.Frame(window)
#content.grid(column=0, row=0)
top_logo_frame = ttk.Frame(window, relief='ridge', borderwidth=5, height=85)
#top_logo_frame = ttk.Frame(window, height=85)
top_logo_frame.grid(row=0, column=0, sticky='ew')

main_frame = tk.Frame(window, relief='ridge', borderwidth=5, height=700)
main_frame.grid(row=1, column=0, sticky='NSEW')


#ttk.Label(top_logo_frame, text='bla', background = 'red').pack(expand = True, fill = 'both')



# top program name and logos
program_name = ttk.Label(top_logo_frame, text='Ring Modulator Layout Program 2', font='Arial 18 bold')
program_name.place(x=10, y=10)
program_author = ttk.Label(top_logo_frame, text='by Dr.-Ing. Marc Schneider (KIT)', font='Arial 10')
program_author.place(x=10, y=50)
program_version = ttk.Label(top_logo_frame, text='v0.0 from 2025-12-11', font='Arial 10')
program_version.place(x=250, y=50)

ttk.Label(main_frame, text='main_frame', background = 'yellow').pack(expand=True, fill='both')




# run
window.mainloop()

"""