# -*- coding: utf-8 -*-
"""
"""
__author__          =   "Stephan Schlamminger and Alireza Panna"
__email__           =   "apanna1236@gmail.com"
__status__          =   "Stable"
__date__            =   "06/26/13"
__version__         =   "0.1"

import xlrd
import numpy as np

class hyreader:
    def __init__(self, bpath, rpath, fname):

        self.bpath = bpath
        self.rpath = rpath
        self.fname = fname
        self.filename = bpath + '/' + rpath + '/' + fname
        self.wb = xlrd.open_workbook(self.filename, "r")
        self.sh1 = self.wb.sheet_by_index(0)
        self.num_rows = self.sh1.nrows
        H = []
        B = []
        tc = []
        mu = []
        M = []
        temp = []
        Bint = []
        Hint = []
        self.Htot = []
        self.Btot = []
        self.tctot = []
        self.temptot = []
        self.BR = []
        self.HC = []
        self.Mu = []
        self.MuH = []
        self.Hmax = []
        self.Bmax = []
        self.mumean = []
        self.muhmean = []
        self.BHE = []
        self.B2 = []
        self.hd = []
        self.bd = []
        self.BRratio = []
        self.Bm = []
        self.Hm = []
        self.hmean = []
        self.bmean = []
        self.Mtot = []
        self.Mudiff = []
        self.Hdiff = []
        co = 0
        mu_0 = 4*np.pi*1e-7

        for lines in range(self.num_rows):
            sr = self.sh1.row_values(lines)
            if (sr != ['', u' ', ''] and sr != [u' ', '', '' ] \
                and sr != [u' ', u' ', u' '] and sr != ['', '', u' ']):
                tc.append(co)
                H.append((sr[0]))
                B.append((sr[1]))
                M.append(((sr[1]/mu_0) - sr[0]))
                temp.append((sr[2]))
                co += 1
            elif(H != [] and B != []):
                self.Htot.append(H)
                self.Btot.append(B)
                self.tctot.append(tc)
                self.temptot.append(temp)
                self.Mtot.append(M)
                H = []
                B = []
                M = []
                temp = []
                co = 0
                tc = []

        self.N = len(self.Htot)
        self.tt = []

        for tem in self.temptot:
            for j in tem:
                if(j != 0):
                    self.tt.append(j)

        for x, i, j in zip(self.tctot, self.Htot, self.Btot):
            Hsum = 0
            Hco = 0
            Bsum = 0
            Bco = 0
            x1= np.linspace(min(x), max(x), num=1000)
            Bint = np.interp(x1, x, j)
            Hint = np.interp(x1, x, i)
            for Bi1, Bi2, xi1, xi2 in \
            zip(Bint[0:-1], Bint[1:], x1[0:-1], x1[1:]):
                if Bi1*Bi2 <= 0:
                    x0 = (0.5)*(xi1 + xi2)
                    Hint1 = np.interp(x0, x, i)
                    if Hint1 > 0:
                        Hsum += Hint1
                        Hco += 1
            if Hco > 0:
                self.HC.append(Hsum/Hco)
            for Hi1, Hi2, xi3, xi4 in \
            zip(Hint[0:-1], Hint[1:], x1[0:-1], x1[1:]):
                if Hi1*Hi2 <= 0:
                    x0 = (0.5)*(xi3 + xi4)
                    Bint1 = np.interp(x0, x, j)
                    if Bint1 > 0:
                        Bsum += Bint1
                        Bco += 1
            if Bco > 0:
                self.BR.append(Bsum/Bco)

        for Hone, Bone in zip(self.Htot, self.Btot):
            mu = []
            muH = []
            for H, B in zip(Hone, Bone):
                if H > 100:
                    mu.append(B/mu_0/H)
                    muH.append(H)
            self.Mu.append(mu)
            self.MuH.append(muH)

        for HH, BB in zip(self.Htot, self.Btot):
                self.Hma = max(HH)
                self.Bma = max(BB)
                self.Bm.append(self.Bma)
                self.Hm.append(self.Hma)
                self.mumean.append(self.Bma/mu_0/self.Hma)
                self.muhmean.append(self.Hma)
        for h1, h2, b1, b2 in \
        zip(self.Hm[0:-1], self.Hm[1:], self.Bm[0:-1], self.Bm[1:]):
                if((h2 - h1) != 0):
                    self.Mudiff.append(((b2 - b1)/(h2 - h1)/mu_0))
                    self.Hdiff.append((h1 + h2)/2)
        BHE = []
        hd = []
        bd = []
        B2 = []
        for HH, BB in zip(self.Htot, self.Btot):
               for hh, bb in zip(HH, BB):
                   if(hh < 0 and bb > 0):
                       BHE.append(hh*bb)
                       hd.append(hh)
                       bd.append(bb)
                       B2.append(bb)
               self.BHE.append(BHE)
               self.hd.append(hd)
               self.bd.append(bd)
               self.B2.append(B2)
               BHE = []
               hd = []
               bd = []
               B2 = []

        for br, bm in zip(self.BR, self.Bm):
            self.BRratio.append(br/bm)
    # Get individual BH curves
    def getHyN(self, n):
        if n >= 0:
            n = n
        else:
            n = len(self.Htot) + n
        return self.Htot[n], self.Btot[n]
    def getHyM(self, n):
         if n >= 0:
            n = n
         else:
            n = len(self.Htot) + n
         return self.Htot[n], self.Mtot[n]
    # Get all BH curves
    def getallHB(self):
        return self.Htot, self.Btot

    def getallHM(self):
        return self.Htot, self.Mtot
    # Get all Coercivity
    def getallCoercivity(self):
        return self.HC
    # Get all Retentivity
    def getallRetentivity(self):
        return self.BR
    # Max of each curve
    def getMaxHB(self, n):
        return max(self.Htot[n]), max(self.Btot[n])
    # Mean line of BH curves
    def getallMax(self):
        for i in range(self.N):
            Ht,Bt = self.getMaxHB(i)
            self.Hmax.append(Ht)
            self.Bmax.append(Bt)
        return self.Hmax, self.Bmax
    # Get temperature data
    def getallTemp(self):
        return self.tt
    # Get all mu
    def getallMu(self):
        return self.Mu, self.MuH
    def getMumean(self):
         return self.mumean, self.muhmean
    def getEproduct(self):
         return self.BHE, self.B2
    def getBrRatio(self):
         return self.BRratio, self.Bm
    def getMeanCurve(self):
        return self.hmean, self.bmean
    def getDemagCurve(self):
        return self.hd, self.bd
    def getallMudiffMean(self):
        return self.Hdiff, self.Mudiff
