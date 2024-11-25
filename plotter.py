from event import Event
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

class Plotter:

    def __init__(self, event):
        print("plotter: initilization")
        self.event = event
        self.Jump(0)

    def Collect(self):
        self.xx = np.array([])
        self.yy = np.array([])
        self.zz = np.array([])
        self.cc = np.array([])
        self.tt = np.array([])
        self.ee = np.array([])
        self.ll = np.array([])
        self.li = np.array([])
        self.pdg = np.array([])
        self.ans = np.array([])

        self.USER_COLORS = ['black', 'red', 'blue', 'magenta']

        nColor = len(self.USER_COLORS)
        mm2m = 0.001
        mm2cm = 0.1
        for i, track in enumerate(self.event.tracks):
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            pdg = track.GetPDGCode()
            for di in depoList:
                depo = self.event.depos[di]
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                t = (depo.GetStart().T() + depo.GetStop().T()) / 2 # ns
                e = depo.GetEnergyDeposit()
                l = depo.GetTrackLength() *mm2cm # most are 0.5 mm
                self.xx = np.append(self.xx, x)
                self.yy = np.append(self.yy, y)
                self.zz = np.append(self.zz, z)
                self.tt = np.append(self.tt, t)
                self.ee = np.append(self.ee, e)
                self.ll = np.append(self.ll, l)
                self.pdg = np.append(self.pdg, pdg)
                self.ans = np.append(self.ans, ancestor)
                self.cc = np.append(self.cc, self.USER_COLORS[ancestor % nColor])

    def Collect_detected(self,threshold=0):
        self.xx = np.array([])
        self.yy = np.array([])
        self.zz = np.array([])
        self.cc = np.array([])
        self.tt = np.array([])
        self.ee = np.array([])
        self.qq = np.array([])
        self.ll = np.array([])
        self.li = np.array([])
        self.pdg = np.array([])
        self.ans = np.array([])

        self.USER_COLORS = ['black', 'red', 'blue', 'magenta']

        nColor = len(self.USER_COLORS)
        mm2m = 0.001
        mm2cm = 0.1
        for i, track in enumerate(self.event.tracks):
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            pdg = track.GetPDGCode()
            # print("********************")
            for di in depoList:
                # print(di)
                depo = self.event.depos[di]
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                t = (depo.GetStart().T() + depo.GetStop().T()) / 2 # ns
                e = depo.GetEnergyDeposit()
                l = depo.GetTrackLength() *mm2cm # most are 0.5 mm
                dedx=0
                q = 0
                if(l>=0.01):
                    dedx=e/l
                    q=self.event.GetdQdx(dedx)*l
                else:
                    q=e*0.7
                light = e-q*0.83
                self.xx = np.append(self.xx, x)
                self.yy = np.append(self.yy, y)
                self.zz = np.append(self.zz, z)
                self.tt = np.append(self.tt, t)
                self.ee = np.append(self.ee, e)
                self.qq = np.append(self.qq, q)
                self.ll = np.append(self.ll, l)
                self.li = np.append(self.li, light)
                self.pdg = np.append(self.pdg, pdg)
                self.ans = np.append(self.ans, ancestor)
                self.cc = np.append(self.cc, self.USER_COLORS[ancestor % nColor])
        self.index = np.where(self.qq>threshold)
        #self.index = np.where(self.ll<threshold)
        # self.index = np.where(self.qq>(1.3/0.6*self.ll+0.2))
        # if threshold==1:
        #     self.index = np.where((self.ll<0.1) & (self.qq>0.075))
        #     # print(self.index)
        # if threshold==0:
        #     self.index = np.where(self.qq>0.075)

        # if threshold==2:
        #     self.index = np.where((self.ee>5/0.6*self.ll+1)&(self.ee<8/0.6*self.ll+2))


    def Jump(self, entryNo,threshold=0):
        self.event.Jump(entryNo)
        self.Collect_detected(threshold)

    def Next(self,threshold=0):
        self.event.Next()
        self.Collect_detected(threshold)

    def Prev(self,threshold=0):
        self.event.Prev()
        self.Collect_detected(threshold)

    def Draw(self, axis='yz', value='time', markerSize=0.2, cmap='jet', vmax=2000):
        # particle, timing, dE/dx

        mapping = {'x': self.xx[self.index], 'y': self.yy[self.index], 'z': self.zz[self.index]}

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5*2, 4), dpi=100)
        # fig.suptitle(self.event.vertex.GetReaction())
        cb_ax = fig.add_axes([.94,.124,.03,.754])

        # particle plot
        ax1.scatter(mapping[axis[1]], mapping[axis[0]], c=self.cc[self.index], s=markerSize)
        if value == 'time':
            # timing plot
            plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.tt[self.index], cmap=cmap, vmin=1, vmax=vmax, norm=matplotlib.colors.LogNorm(), s=markerSize)
            cb_ax.set_xlabel('ns')

        elif value == 'charge':
            # charge plot
            plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.ee[self.index], cmap=cmap, vmax=4, s=markerSize)
            cb_ax.set_xlabel('MeV/cm')

        elif value == 'length':
            # length plot
            plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.ll[self.index], cmap=cmap, vmax=0.2, s=markerSize)
            cb_ax.set_xlabel('cm')

        elif value == 'light':
            # length plot
            plot_12 = ax2.scatter(mapping[axis[1]], mapping[axis[0]], c=self.li[self.index], cmap=cmap, vmax=0.2, s=markerSize)
            cb_ax.set_xlabel('MeV')

        fig.colorbar(plot_12, orientation='vertical', cax=cb_ax)
        # fig.tight_layout()

        for ax in (ax1, ax2):
            ax.set_ylim(-4, 4)
            ax.set_xlim(-2, 8)
            ax.tick_params(axis='y', direction='in', length=2)
            ax.tick_params(axis='x', direction='in', length=2)
            ax.set_xlabel(f'{axis[1]} [m]')
            if ax == ax1:
                ax.set_ylabel(f'{axis[0]} [m]')


        xpos = -1.8
        ypos = 3.6
        nColor = len(self.USER_COLORS)
        countnegId = 0
        for i, particle in enumerate(self.event.vertex.Particles):
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            trkId = particle.GetTrackId()
            if trkId < 0:
                # Because we skipped the trk id, need to subtract the index for the proper coloring of deposits
                countnegId += 1
                continue
            name = particle.GetName()
            color = self.USER_COLORS[(i-countnegId) % nColor]
            # pdg = particle.GetPDGCode()
            # name = particle.GetName()
            # trkId = particle.GetTrackId()
            mom = particle.GetMomentum()
            KE = mom.E() - mom.M()
            name = '%s: %.1f MeV' % (name, KE)

            ax1.text(xpos, ypos, name, color=color)
            ypos -= 0.35

        #ax2.plot()
        #fig.savefig(self.event.plotpath + '/particle_%s_evt_%d.pdf' % (value, self.event.currentEntry) )
        fig.savefig(self.event.plotpath + '/particle_%s_evt_%d.pdf' % (value, self.event.currentEntry) )

        # extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        # fig.savefig(self.event.plotpath + '/particle_%s_evt_%d_ax1.pdf' % (value, self.event.currentEntry) , bbox_inches=extent.expanded(1.2, 1.5))
#         plt.clf() # important to clear figure
#         plt.close()
        plt.show()

    def Draw_sp(self, axis='yz', markerSize=0.2, cmap='jet', vmax=2000):
        # particle, timing, dE/dx
        mapping = {'x': self.xx[self.index], 'y': self.yy[self.index], 'z': self.zz[self.index]}

        fig, ax1 = plt.subplots(1, 1, figsize=(5, 4), dpi=100)
        # fig.suptitle(self.event.vertex.GetReaction())
        # particle plot
        ax1.scatter(mapping[axis[1]], mapping[axis[0]], c=self.cc[self.index], s=markerSize)

        # fig.tight_layout()

     
        ax1.set_ylim(-4, 4)
        ax1.set_xlim(-2, 8)
        ax1.tick_params(axis='y', direction='in', length=2)
        ax1.tick_params(axis='x', direction='in', length=2)
        ax1.set_xlabel(f'{axis[1]} [m]')
        ax1.set_ylabel(f'{axis[0]} [m]')


        xpos = -1.8
        ypos = 3.6
        nColor = len(self.USER_COLORS)
        countnegId = 0
        for i, particle in enumerate(self.event.vertex.Particles):
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            trkId = particle.GetTrackId()
            if trkId < 0:
                # Because we skipped the trk id, need to subtract the index for the proper coloring of deposits
                countnegId += 1
                continue
            name = particle.GetName()
            color = self.USER_COLORS[(i-countnegId) % nColor]
            # pdg = particle.GetPDGCode()
            # name = particle.GetName()
            # trkId = particle.GetTrackId()
            mom = particle.GetMomentum()
            KE = mom.E() - mom.M()
            name = '%s: %.1f MeV' % (name, KE)

            ax1.text(xpos, ypos, name, color=color)
            ypos -= 0.35

        #ax2.plot()
        #fig.savefig(self.event.plotpath + '/particle_%s_evt_%d.pdf' % (value, self.event.currentEntry) )
        fig.savefig(self.event.plotpath + '/particle_evt_%d.pdf' % (self.event.currentEntry) )

        # extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        # fig.savefig(self.event.plotpath + '/particle_%s_evt_%d_ax1.pdf' % (value, self.event.currentEntry) , bbox_inches=extent.expanded(1.2, 1.5))
#         plt.clf() # important to clear figure
#         plt.close()
        plt.show()

    def scatter_dedl(self):
        index = np.where((self.pdg==2212))
        plt.scatter(self.ll[index],self.ee[index],s=0.5)
        plt.draw()

    def hist_dedx(self):
        plt.hist(self.ee, range=(0,6), bins=100)
        plt.draw()
        # plt.savefig(self.event.plotpath + '/dE_dx_evt_%d.pdf' % self.event.currentEntry)
#         plt.clf() # important to clear figure
#         plt.close()
        plt.show()

    def hist_dqdx(self):
        plt.hist(self.qq, range=(0,6), bins=100)
        plt.draw()
        # plt.savefig(self.event.plotpath + '/dQ_dx_evt_%d.pdf' % self.event.currentEntry)
#         plt.clf() # important to clear figure
#         plt.close()
        plt.show()

    def hist_dqdx_thre(self):
        plt.hist(self.qq[self.index], range=(0,6), bins=100)
        plt.draw()
        # plt.savefig(self.event.plotpath + '/dQ_dx_thre_evt_%d.pdf' % self.event.currentEntry)
#         plt.clf() # important to clear figure
#         plt.close()
        plt.show()

    def hist_tracklength(self):
        plt.hist(self.ll, range=(0,0.6), bins=100)
        plt.draw()
        # plt.savefig(self.event.plotpath + '/tracklength_evt_%d.pdf' % self.event.currentEntry)
#         plt.clf() # important to clear figure
#         plt.close()
        plt.show()

    def hist_tracklength_thre(self):
        plt.hist(self.ll[self.index], range=(0,0.6), bins=100)
        plt.draw()
        # plt.savefig(self.event.plotpath + '/tracklength_evt_%d.pdf' % self.event.currentEntry)
#         plt.clf() # important to clear figure
#         plt.close()
        plt.show()
    #---------------------------------------------
    def DrawROOT(self, dim2d='yz', markerSize=0.2):
        import ROOT
        from ROOT import TH2F, TMarker, TCanvas, TLatex
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetMarkerStyle(24)
        ROOT.gStyle.SetMarkerSize(markerSize)
        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta]
        nColor = len(colors)

        c1 = TCanvas("c1", self.event.vertex.GetReaction(), 800, 800)

        # canvas of 5m x 5m
        dummy = TH2F("dummy", "", 100, -5, 5, 100, -5, 5)
        dummy.GetXaxis().SetTitle('[m]')
        dummy.GetYaxis().SetTitle('[m]')
        dummy.Draw()

        mm2m = 0.001
        txts = []
        markers = []
        txtX = 0.15
        txtY = 0.85
        txtSize = 0.03
        # nDepo = 0
        for i, track in enumerate(self.event.tracks):
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            color = colors[ancestor % nColor]


            if ancestor == i:
                txt = TLatex()
                txt = txt.DrawLatexNDC(txtX, txtY, track.GetName())
                txt.SetTextColor(color)
                txt.SetTextSize(txtSize)
                txts.append(txt)
                txtY -= txtSize

            for j, di in enumerate(depoList):
                depo = self.event.depos[di]
                x = (depo.GetStart().X() + depo.GetStop().X()) / 2 * mm2m
                y = (depo.GetStart().Y() + depo.GetStop().Y()) / 2 * mm2m
                z = (depo.GetStart().Z() + depo.GetStop().Z()) / 2 * mm2m
                mapping = {'x': x, 'y': y, 'z': z}

                m = TMarker(mapping[dim2d[1]], mapping[dim2d[0]], 24)
                m.SetMarkerColor(color)

                m.Draw()
                markers.append(m)
                # nDepo += 1

        # print('depo points drawn: ', nDepo, '| total depo: ', self.event.depos.size)

        ROOT.gPad.Update()
        return c1

if __name__ == "__main__":
    event = Event(sys.argv[1])
    p = Plotter(event)
    p.Next()
    p.Draw('yz')
    # c1 = p.DrawROOT('xz', 0.2)
    # input('press a key to continue ...')
