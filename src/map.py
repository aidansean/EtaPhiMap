##########################################################################################
# Written by Aidan Randle-Conde (aidan.randleconde@gmail.com)                            #
# Feel free to use, edit, and redistribute, but attribute me as the original creator     #
# 2014/11/03 21:54 UCT                                                                   #
##########################################################################################

import math

##########################################################################################
# Draw options                                                                           #
##########################################################################################
draw_crystals      = True # Draw outlines to crystals
fill_emptyCrystals = True # Colour in empty crystals
fill_dRCrystals    = True # Colour in crystals in the isolation cone
fill_5x5Crystals   = True # Colour in crystals in the 5x5

draw_dRCurves  = True  # Draw the isolation curves
draw_etaCurves = True  # Draw curves of constant eta
draw_phiLines  = True  # Draw lines of constant phi
draw_dRLabels  = True  # Draw labels showing 5x5/isolation crystals
draw_etaLabels = False # Draw eta labels (clutters up the plot)

draw_dRInfoLabel = True # Draw the label showing the size of the isolation curves
draw_zInfoLabel  = True # Show the z cooridinate

##########################################################################################
# Parameters                                                                             #
##########################################################################################
pi = 3.14159
dR       = 0.3   # Size in deltaR of "isolation" rings
DEta     = 0.3   # Separation of eta rings in the phi-eta grid
DPhi     = pi/10 # Separation of phi lines in the phi-eta grid
etaStart = 1.5   # Where to start making eta rings
etaStop  = 3.0   # And when to stop
nPoints  = 100   # How many points per ring
nCurves  =  20   # How many "isolation" rings to create

##########################################################################################
# Geometry                                                                               #
##########################################################################################
z0 = 3.14    # Distance in z from the centre of the CMS detector (Hey, it's pi m!)
eta0 = 1.479 # Start of endcap ecal in eta
eta1 = 3.000 # End of endcap ecal in eta (note this is the "inner" limit)
theta0 = 2*math.atan(math.exp(-eta0)) # Start of endcap ecal in theta
theta1 = 2*math.atan(math.exp(-eta1)) # End of endcap ecal in theta
inner_radius = z0*math.tan(theta1) # It's just trigonometry, dawg
outer_radius = z0*math.tan(theta0)

nCrystalsAcross = 53*2 # There are roughly 106 crystals across the endcap
nCrystalsAcross = 2*outer_radius/0.028
# Take 2*outer_radius and divide by crystal size, 2.28cm

##########################################################################################
# ROOT and style                                                                         #
# Standard boiler plate stuff by now                                                     #
##########################################################################################
import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetFillStyle(ROOT.kWhite)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetFrameBorderMode(ROOT.kWhite)
ROOT.gStyle.SetFrameFillColor(ROOT.kWhite)
ROOT.gStyle.SetCanvasBorderMode(ROOT.kWhite)
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetPadBorderMode(ROOT.kWhite)
ROOT.gStyle.SetPadColor(ROOT.kWhite)
ROOT.gStyle.SetStatColor(ROOT.kWhite)
ROOT.gStyle.SetErrorX(0)

lineColor = ROOT.kBlack # Colour of the eta-phi grid

##########################################################################################
# Make a pretty canvas                                                                   #
##########################################################################################
cw = 600
ch = 600
canvas = ROOT.TCanvas('canvas','',100,100,cw,ch)
canvas.SetGridx()
canvas.SetGridy()
canvas.SetFillColor(ROOT.kWhite)
canvas.SetBorderMode(0)

##########################################################################################
# Make the map from a large 2D histogram, with each cell being a crystal (or dead space) #
##########################################################################################
scale = 1.1 # Give the map some breathing space between it and the axes
size = scale*outer_radius
nBins = int(scale*nCrystalsAcross)
hMapBase = ROOT.TH2F('hMapBase', '', nBins, -size, size, nBins, -size, size)
hMapBase.GetXaxis().SetTitle('x [m]')
hMapBase.GetYaxis().SetTitle('y [m]')
hMapBase.GetXaxis().SetTitleOffset(1.25)
hMapBase.GetYaxis().SetTitleOffset(1.25)

# hMapBase will be square outlines, hMap will be colours
hMap = hMapBase.Clone('hMap')

# Fill hMapBase and hMap depending on whether they fall in the extent of the endcap ecal
for binX in range(1,1+hMapBase.GetNbinsX()):
    for binY in range(1,1+hMapBase.GetNbinsY()):
        xmax = hMapBase.GetXaxis().GetBinLowEdge(binX)+hMapBase.GetXaxis().GetBinWidth(binX)
        ymax = hMapBase.GetYaxis().GetBinLowEdge(binY)+hMapBase.GetYaxis().GetBinWidth(binY)
        r = math.sqrt(xmax*xmax+ymax*ymax)
        if r>inner_radius and r<outer_radius:
            hMapBase.SetBinContent(binX,binY,1)
            if fill_emptyCrystals:
                hMap.SetBinContent(binX,binY,1)

##########################################################################################
# Define some classes to hold useful infomrmation                                        #
# (Maybe put this in a separate file later)                                              #
##########################################################################################
class DR_curve_object:
    '''Class to hold information about the isolation curves'''
    def __init__(self, cx, cy, points, nIso, nSC):
        self.cx = cx
        self.cy = cy
        self.points = points
        self.nIso = nIso
        self.nSC = nSC
        
        self.r2max = -1e6
        self.r2min =  1e6
        self.xmax  = -1e6
        self.xmin  =  1e6
        self.ymax  = -1e6
        self.ymin  =  1e6
        for i in range(0,len(self.points)):
            x = self.points[i][0]
            x = self.points[i][1]
            if x*x+y*y > self.r2max:
                self.r2max = x*x + y*y
            if x*x+y*y < self.r2min:
                self.r2min = x*x + y*y
            if x>self.xmax:
                self.xmax = x
            if x<self.xmin:
                self.xmin = x
            if y>self.ymax:
                self.ymax = y
            if y<self.ymin:
                self.ymin = y
        
    def SC_label(self):
        label = ROOT.TLatex(self.cx,self.cy,'%d/%d'%(self.nSC,self.nIso))
        label.SetTextSize(0.02)
        label.SetTextAlign(22)
        label.SetTextColor(ROOT.kWhite)
        return label

class isoeta_curve:
    '''Class to hold information about the constant eta curves'''
    def __init__(self, eta, points, index):
        self.eta = eta
        self.points = points
        self.index = index
    def label(self):
        phi = -0.5*pi*self.index*DEta/(etaStop-etaStart)
        
        theta_tmp = 2*math.atan(math.exp(-self.eta))
        rho = z0*math.tan(theta_tmp)
        
        TV3 = ROOT.TVector3()
        TV3.SetPtEtaPhi(rho,self.eta,phi)
        x = TV3.X()
        y = TV3.Y()
        
        label = ROOT.TLatex(x,y,'#eta=%.1f'%self.eta)
        label.SetTextSize(0.02)
        label.SetTextAlign(22)
        label.SetTextColor(ROOT.kWhite)
        return label

class isophi_line:
    '''Class to hold information about the constant phi lines'''
    def __init__(self, phi, x1, y1, x2, y2):
        self.phi = phi
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.line = ROOT.TLine(self.x1,self.y1,self.x2,self.y2)
        self.line.SetLineWidth(1)
        self.line.SetLineColor(lineColor)

##########################################################################################
# Now make the isolation curves and update the histograms                                #
##########################################################################################
DR_curves = []
for j in range(0,nCurves):
    r = inner_radius + (j*1.0/nCurves)*(outer_radius-inner_radius)
    p = 2*2*pi*j/nCurves
    cx = r*math.cos(p)
    cy = r*math.sin(p)
    cp = ROOT.TVector3(cx,cy,z0)
    points = []
    for i in range(0,nPoints+1):
        psi  = 2*pi*(i+0)/nPoints
        dEta = dR*math.cos(psi)
        dPhi = dR*math.sin(psi)
        
        # This bit's tricky.  Get it wrong and you lose the eta variations
        eta_tmp   = cp.Eta()+dEta
        theta_tmp = 2*math.atan(math.exp(-eta_tmp))
        rho = z0*math.tan(theta_tmp)
        
        TV3 = ROOT.TVector3()
        TV3.SetPtEtaPhi(rho,cp.Eta()+dEta,cp.Phi()+dPhi)
        x = TV3.X()
        y = TV3.Y()
        points.append([x,y])
    
    # Count up how many crystals are in the 5x5 and how many are in the isolation cone
    # Update the histograms with this information
    nIso = 0
    nSC  = 0
    for binX in range(1,1+hMap.GetNbinsX()) :
        for binY in range(1,1+hMap.GetNbinsY()):
            crysx = hMap.GetXaxis().GetBinCenter(binX)
            crysy = hMap.GetYaxis().GetBinCenter(binY)
            TV3 = ROOT.TVector3(crysx,crysy,z0)
            if TV3.DeltaR(cp)<dR and TV3.Perp()>=inner_radius and TV3.Perp()<=outer_radius:
                if fill_dRCrystals:
                    hMap.SetBinContent(binX,binY,2)
                nIso += 1
    
    
    crysbinx = hMap.GetXaxis().FindBin(cx)
    crysbiny = hMap.GetYaxis().FindBin(cy)
    for binX in range(crysbinx-2,crysbinx+3):
        for binY in range(crysbiny-2,crysbiny+3):
            x = hMap.GetXaxis().GetBinCenter(binX)
            y = hMap.GetYaxis().GetBinCenter(binY)
            r = math.sqrt(x*x+y*y)
            if r>=inner_radius and r<=outer_radius:
                if fill_5x5Crystals:
                    hMap.SetBinContent(binX,binY,3)
                nSC += 1
    DCurve = DR_curve_object(cx,cy,points,nIso,nSC)
    DR_curves.append(DCurve)

##########################################################################################
# Now for the curves of constant eta and lines of constant phi                           #
##########################################################################################
eta_curves = []
eta = etaStart
index = 0
while eta < etaStop+1e-6:
    points = []
    for i in range(0,nPoints):
        phi = 2*pi*i/nPoints
        
        theta_tmp = 2*math.atan(math.exp(-eta))
        rho = z0*math.tan(theta_tmp)
        
        TV3 = ROOT.TVector3()
        TV3.SetPtEtaPhi(rho,eta,phi)
        x = TV3.X()
        y = TV3.Y()
        points.append([x,y])
    curve = isoeta_curve(eta, points, index)
    eta_curves.append(curve)
    eta += DEta
    index += 1

phi_lines = []
phi = 0
while phi < 2*pi:
    eta1 = etaStart
    eta2 = eta-DEta
    
    theta1 = 2*math.atan(math.exp(-eta1))
    theta2 = 2*math.atan(math.exp(-eta2))
    rho1 = z0*math.tan(theta1)
    rho2 = z0*math.tan(theta2)
    
    TV31 = ROOT.TVector3()
    TV32 = ROOT.TVector3()
    TV31.SetPtEtaPhi(rho1,eta1,phi)
    TV32.SetPtEtaPhi(rho2,eta2,phi)
    x1 = TV31.X()
    y1 = TV31.Y()
    x2 = TV32.X()
    y2 = TV32.Y()
    phi_lines.append(isophi_line(phi, x1, y1, x2, y2))
    
    phi += DPhi

##########################################################################################
# Put it all together and draw!                                                          #
##########################################################################################

# Make some arrays to contain objects so that python doesn't garbage collect them
lines   = []
markers = []
labels  = []

hMap.Draw('COL')
if draw_crystals:
    hMapBase.Draw('BOX:sames')

if draw_phiLines:
    for p in phi_lines:
        line = p.line
        lines.append(line)
        line.Draw()

if draw_etaCurves:
    for c in eta_curves:
        for i in range(0,len(c.points)):
            p0 = c.points[i+0]
            p1 = c.points[(i+1)%len(c.points)]
            line = ROOT.TLine(p0[0],p0[1],p1[0],p1[1])
            line.SetLineColor(lineColor)
            line.SetLineWidth(1)
            lines.append(line)
            line.Draw()
        label = c.label()
        labels.append(label)
        if draw_etaLabels:
            label.Draw()

if draw_dRCurves:
    for c in DR_curves:
        for i in range(0,len(c.points)-1):
            p0 = c.points[i+0]
            p1 = c.points[i+1]
            line = ROOT.TLine(p0[0],p0[1],p1[0],p1[1])
            line.SetLineColor(ROOT.kBlue)
            line.SetLineWidth(2)
            lines.append(line)
            line.Draw()
    
        label = c.SC_label()
        labels.append(label)
        if draw_dRLabels:
            label.Draw()

DR_label = ROOT.TLatex(1.25, 1.25, '#DeltaR=0.3')
DR_label.SetTextAlign(22)
if draw_dRInfoLabel:
    DR_label.Draw()

z_label = ROOT.TLatex(1.25, 1.5, 'z=3.14 m')
z_label.SetTextAlign(22)
if draw_zInfoLabel:
    z_label.Draw()

canvas.Print('baseMap.eps')

