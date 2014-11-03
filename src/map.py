import math

# Parameters
dR      = 0.3
nPoints = 100
nCurves =  20

##########################################################################################
# ROOT and style                                                                         #
##########################################################################################
import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
#ROOT.gROOT.ProcessLine('.L Loader.C+')

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

cw = 600
ch = 600
canvas = ROOT.TCanvas('canvas','',100,100,cw,ch)
canvas.SetGridx()
canvas.SetGridy()
canvas.SetFillColor(ROOT.kWhite)
canvas.SetBorderMode(0)

pi = 3.14159
z0 = 3.14
eta0 = 1.479
eta1 = 3.000
theta0 = 2*math.atan(math.exp(-eta0))
theta1 = 2*math.atan(math.exp(-eta1))
inner_radius = z0*math.tan(theta1)
outer_radius = z0*math.tan(theta0)

nCrystalsAcross = 53*2

size = 1.1*outer_radius
hMapBase = ROOT.TH2F('hMapBase', '', nCrystalsAcross, -size, size, nCrystalsAcross, -size, size)
hMapBase.GetXaxis().SetTitle('x [m]')
hMapBase.GetYaxis().SetTitle('y [m]')
hMapBase.GetXaxis().SetTitleOffset(1.25)
hMapBase.GetYaxis().SetTitleOffset(1.25)

hMap = hMapBase.Clone('hMap')

for binX in range(1,1+hMapBase.GetNbinsX()):
    for binY in range(1,1+hMapBase.GetNbinsY()):
        xmax = hMapBase.GetXaxis().GetBinLowEdge(binX)+hMapBase.GetXaxis().GetBinWidth(binX)
        ymax = hMapBase.GetYaxis().GetBinLowEdge(binY)+hMapBase.GetYaxis().GetBinWidth(binY)
        r = math.sqrt(xmax*xmax+ymax*ymax)
        if r>inner_radius and r<outer_radius:
            hMapBase.SetBinContent(binX,binY,1)
            hMap.SetBinContent(binX,binY,1)

class DR_curve_object:
    def __init__(self, cx, cy, points):
        self.cx = cx
        self.cy = cy
        self.points = points

class isoeta_curves:
    def __init__(self, eta, points):
        self.eta = eta
        self.points = points

DR_curves = []
for j in range(0,nCurves):
    r = inner_radius + (j*1.0/nCurves)*(outer_radius-inner_radius)
    p = 2*2*pi*j/nCurves
    cx = r*math.cos(p)
    cy = r*math.sin(p)
    cp = ROOT.TVector3(cx,cy,z0)
    r2min =  1e6
    r2max = -1e6
    xmin  =  1e6
    xmax  = -1e6
    ymin  =  1e6
    ymax  = -1e6
    points = []
    for i in range(0,nPoints+1):
        psi  = 2*pi*(i+0)/nPoints
        dEta = dR*math.cos(psi)
        dPhi = dR*math.sin(psi)
        
        eta_tmp = cp.Eta()+dEta
        theta_tmp = 2*math.atan(math.exp(-eta_tmp))
        rho = z0*math.tan(theta_tmp)
        
        TV3 = ROOT.TVector3()
        TV3.SetPtEtaPhi(rho,cp.Eta()+dEta,cp.Phi()+dPhi)
        x = TV3.X()
        y = TV3.Y()
        if x*x+y*y > r2max:
            r2max = x*x + y*y
        if x*x+y*y < r2min:
            r2min = x*x + y*y
        if x>xmax:
            xmax = x
        if x<xmin:
            xmin = x
        if y>ymax:
            ymax = y
        if y<ymin:
            ymin = y
        points.append([x,y])
    DCurve = DR_curve_object(cx,cy,points)
    
    nIso = 0
    nSC  = 0
    for binX in range(1,1+hMap.GetNbinsX()) :
        for binY in range(1,1+hMap.GetNbinsY()):
            crysx = hMap.GetXaxis().GetBinCenter(binX)
            crysy = hMap.GetYaxis().GetBinCenter(binY)
            TV3 = ROOT.TVector3(crysx,crysy,z0)
            if TV3.DeltaR(cp)<dR and TV3.Perp()>=inner_radius and TV3.Perp()<=outer_radius:
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
                hMap.SetBinContent(binX,binY,3)
                nSC += 1
    DCurve.nIso = nIso
    DCurve.nSC  = nSC
    DR_curves.append(DCurve)

hMap.Draw('COL')
hMapBase.Draw('BOX:sames')
lines   = []
markers = []
labels  = []
for c in DR_curves:
    for i in range(0,nPoints):
        p0 = c.points[i+0]
        p1 = c.points[i+1]
        line = ROOT.TLine(p0[0],p0[1],p1[0],p1[1])
        line.SetLineColor(ROOT.kBlue)
        line.SetLineWidth(2)
        lines.append(line)
        line.Draw()
        
        marker = ROOT.TMarker(p0[0],p0[1],22)
        marker.SetMarkerColor(50+i)
        markers.append(marker)
        #marker.Draw()
        
        #print '%10.7f  %10.7f  %10.7f  %10.7f'%(p0[0],p0[1],p1[0],p1[1])
        
    marker = ROOT.TMarker(c.cx,c.cy,20)
    markers.append(marker)
    #marker.Draw()
    
    label = ROOT.TLatex(c.cx,c.cy,'%d/%d'%(c.nSC,c.nIso))
    label.SetTextSize(0.02)
    label.SetTextAlign(22)
    label.SetTextColor(ROOT.kWhite)
    labels.append(label)
    label.Draw()

label = ROOT.TLatex(1.25, 1.25, '#DeltaR=0.3')
label.SetTextAlign(22)
label.Draw()

canvas.Print('plots/baseMap.eps')

