import math

# Parameters
pi = 3.14159
dR       = 0.3
DEta     = 0.3
DPhi     = pi/10
etaStart = 1.5
etaStop  = 3.0
nPoints  = 100
nCurves  =  20

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

lineColor = ROOT.kYellow+3
lineColor = ROOT.kBlack

##########################################################################################
# Geometry                                                                               #
##########################################################################################
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
    def __init__(self, phi, x1, y1, x2, y2):
        self.phi = phi
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.line = ROOT.TLine(self.x1,self.y1,self.x2,self.y2)
        self.line.SetLineWidth(1)
        self.line.SetLineColor(lineColor)

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
        
        eta_tmp   = cp.Eta()+dEta
        theta_tmp = 2*math.atan(math.exp(-eta_tmp))
        rho = z0*math.tan(theta_tmp)
        
        TV3 = ROOT.TVector3()
        TV3.SetPtEtaPhi(rho,cp.Eta()+dEta,cp.Phi()+dPhi)
        x = TV3.X()
        y = TV3.Y()
        points.append([x,y])
    
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
    DCurve = DR_curve_object(cx,cy,points,nIso,nSC)
    DR_curves.append(DCurve)

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

hMap.Draw('COL')
hMapBase.Draw('BOX:sames')
lines   = []
markers = []
labels  = []

for p in phi_lines:
    line = p.line
    lines.append(line)
    line.Draw()

for c in eta_curves:
    for i in range(0,len(c.points)):
        p0 = c.points[i+0]
        p1 = c.points[(i+1)%len(c.points)]
        line = ROOT.TLine(p0[0],p0[1],p1[0],p1[1])
        line.SetLineColor(lineColor)
        line.SetLineWidth(1)
        lines.append(line)
        line.Draw()
    #label = c.label()
    #labels.append(label)
    #label.Draw()
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
    label.Draw()


DR_label = ROOT.TLatex(1.25, 1.25, '#DeltaR=0.3')
DR_label.SetTextAlign(22)
DR_label.Draw()

z_label = ROOT.TLatex(1.25, 1.5, 'z=3.14 m')
z_label.SetTextAlign(22)
z_label.Draw()

canvas.Print('plots/baseMap.eps')

