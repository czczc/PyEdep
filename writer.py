from event import Event
import sys
import numpy as np
from array import array
from ROOT import TFile, TTree, TH1D,TH2D

class Writer:

    def __init__(self, event, outfile='output.root'):
        self.event = event
        self.h_dedx_e_total = TH1D("h_dedx_e_total","dEdx_e",1000,0,20)
        self.h_dedx_p_total = TH1D("h_dedx_p_total","dEdx_e",1000,0,20)
        self.h_dedx_n_total = TH1D("h_dedx_n_total","dEdx_e",1000,0,20)
        self.h_dedx_pi_total = TH1D("h_dedx_pi_total","dEdx_e",1000,0,20)
        self.h_dedx_e_2d_total = TH2D("h_dedx_e_2d_total","dEdx_e",1000,0,20,1000,0,20)
        self.h_dedx_p_2d_total = TH2D("h_dedx_p_2d_total","dEdx_e",1000,0,20,1000,0,20)
        self.h_dedx_n_2d_total = TH2D("h_dedx_n_2d_total","dEdx_e",1000,0,20,1000,0,20)
        self.h_dedx_pi_2d_total = TH2D("h_dedx_pi_2d_total","dEdx_e",1000,0,20,1000,0,20)
        self.f_out = TFile(outfile, 'RECREATE')
        self.initOutputTree2()

    def initOutputTree2(self):
        self.T_out = TTree('Sim', 'Sim') # output tree

        self.nu_pdg = array('i', [0]) #
        self.T_out.Branch('nu_pdg', self.nu_pdg, 'nu_pdg/I')

        self.nu_xs = array('f', [0])
        self.T_out.Branch('nu_xs', self.nu_xs, 'nu_xs/F')

        self.nu_proc = array('i', [0]) # genie process
        self.T_out.Branch('nu_proc', self.nu_proc, 'nu_proc/I')

        self.nu_nucl = array('i', [0]) # which nucleon
        self.T_out.Branch('nu_nucl', self.nu_nucl, 'nu_nucl/I')

        self.E_nu = array('f', [0]) # true neutrino energy
        self.T_out.Branch('E_nu', self.E_nu, 'E_nu/F')

        self.E_depoTotal_recoil = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('E_depoTotal_recoil', self.E_depoTotal_recoil, 'E_depoTotal_recoil/F')

    def initOutputTree(self):
        self.T_out = TTree('Sim', 'Sim') # output tree

        self.nu_pdg = array('i', [0]) #
        self.T_out.Branch('nu_pdg', self.nu_pdg, 'nu_pdg/I')

        self.nu_xs = array('f', [0])
        self.T_out.Branch('nu_xs', self.nu_xs, 'nu_xs/F')

        self.nu_proc = array('i', [0]) # genie process
        self.T_out.Branch('nu_proc', self.nu_proc, 'nu_proc/I')

        self.nu_nucl = array('i', [0]) # which nucleon
        self.T_out.Branch('nu_nucl', self.nu_nucl, 'nu_nucl/I')

        self.E_nu = array('f', [0]) # true neutrino energy
        self.T_out.Branch('E_nu', self.E_nu, 'E_nu/F')

        self.E_avail = array('f', [0]) # energy availabe from vertex interaction (excluding energy lost inside nuclei)
        self.T_out.Branch('E_avail', self.E_avail, 'E_avail/F')

        # self.E_availList = np.zeros((6,), dtype=np.float32) # E avail for: lepton, proton, neutron, pi+-, pi0, others.
        # self.T_out.Branch('E_availList', self.E_availList, 'E_availList[6]/F')

        self.E_availList = np.zeros((7,), dtype=np.float32) # E avail for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_availList', self.E_availList, 'E_availList[7]/F')

        self.E_depoTotal = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('E_depoTotal', self.E_depoTotal, 'E_depoTotal/F')

        self.E_depoTotal_recoil = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('E_depoTotal_recoil', self.E_depoTotal_recoil, 'E_depoTotal_recoil/F')

        self.Q_depoTotal = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('Q_depoTotal', self.Q_depoTotal, 'Q_depoTotal/F')

        self.Q_depoTotal_thre = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('Q_depoTotal_thre', self.Q_depoTotal_thre, 'Q_depoTotal_thre/F')

        self.E_depoTotal_re = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('E_depoTotal_re', self.E_depoTotal_re, 'E_depoTotal_re/F')

        self.E_depoTotal_l = array('f', [0]) # total energy deposit from all (charged) tracks
        self.T_out.Branch('E_depoTotal_l', self.E_depoTotal_l, 'E_depoTotal_l/F')

        # self.E_depoList = np.zeros((6,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        # self.T_out.Branch('E_depoList', self.E_depoList, 'E_depoList[6]/F')

        self.E_depoList = np.zeros((7,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_depoList', self.E_depoList, 'E_depoList[7]/F')

        self.Q_depoList = np.zeros((7,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('Q_depoList', self.Q_depoList, 'Q_depoList[7]/F')

        self.Q_depoList_thre = np.zeros((7,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('Q_depoList_thre', self.Q_depoList_thre, 'Q_depoList_thre[7]/F')

        self.E_depoList_re = np.zeros((7,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_depoList_re', self.E_depoList_re, 'E_depoList_re[7]/F')

        self.E_depoList_re_track = np.zeros((7,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_depoList_re_track', self.E_depoList_re_track, 'E_depoList_re_track[7]/F')

        self.E_depoList_re_lep = np.zeros((7,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_depoList_re_lep', self.E_depoList_re_lep, 'E_depoList_re_lep[7]/F')

        self.E_depoList_re_had = np.zeros((7,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_depoList_re_had', self.E_depoList_re_had, 'E_depoList_re_had[7]/F')

        self.E_depoList_l = np.zeros((7,), dtype=np.float32) # depo for: lepton, proton, neutron, pi+-, pi0, others.
        self.T_out.Branch('E_depoList_l', self.E_depoList_l, 'E_depoList_l[7]/F')


    def Write(self):
        # self.stat = {
        # }
        self.f_out.cd()

        for i in range(self.event.nEntry):
        # for i in range(2):
            self.event.Jump(i)

            # proc = str(self.event.info['nu_proc']) + '-' + str(self.event.info['nu_nucl'])
            # v = self.stat.setdefault(proc, 0)
            # self.stat[proc] = v + 1
            self.nu_pdg[0] = self.event.info['nu_pdg']
            self.nu_xs[0] = self.event.info['nu_xs']
            self.nu_proc[0] = self.event.info['nu_proc']
            self.nu_nucl[0] = self.event.info['nu_nucl']
            self.E_nu[0] = self.event.info['E_nu']

            self.E_depoTotal[0] = self.event.info['E_depoTotal']
            self.Q_depoTotal[0] = self.event.info['Q_depoTotal']
            self.Q_depoTotal_thre[0] = self.event.info['Q_depoTotal_thre']
            self.E_depoTotal_re[0] = self.event.info['E_depoTotal_re']
            self.E_depoTotal_l[0] = self.event.info['E_depoTotal_l']
            self.E_avail[0] = self.event.info['E_avail']
            self.E_availList[:] = self.event.info['E_availList']
            self.E_depoList[:] = self.event.info['E_depoList']
            self.Q_depoList[:] = self.event.info['Q_depoList']
            self.Q_depoList_thre[:] = self.event.info['Q_depoList_thre']
            self.E_depoList_re[:] = self.event.info['E_depoList_re']
            self.E_depoList_re_track[:] = self.event.info['E_depoList_re_track']
            self.E_depoList_re_had[:] = self.event.info['E_depoList_re_had']
            self.E_depoList_re_lep[:] = self.event.info['E_depoList_re_lep']
            self.E_depoList_l[:] = self.event.info['E_depoList_l']
            self.E_depoTotal_recoil[0] = self.event.info['E_depoTotal_recoil']


            self.T_out.Fill()

        self.T_out.Write()
        # self.event.CheckTH2D()
        # print(self.stat)

    def Write_dedx(self):
        # self.stat = {
        # }
        self.f_out.cd()

        for i in range(self.event.nEntry):
        # for i in range(100):
            self.event.Jump(i)
            self.h_dedx_e_2d_total.Add(self.event.h_dedx_2d_e)
            self.h_dedx_e_total.Fill(self.event.h_dedx_e.GetMean())
            self.h_dedx_p_2d_total.Add(self.event.h_dedx_2d_p)
            self.h_dedx_p_total.Fill(self.event.h_dedx_p.GetMean())
            self.h_dedx_n_2d_total.Add(self.event.h_dedx_2d_n)
            self.h_dedx_n_total.Fill(self.event.h_dedx_n.GetMean())
            self.h_dedx_pi_2d_total.Add(self.event.h_dedx_2d_pi)
            self.h_dedx_pi_total.Fill(self.event.h_dedx_pi.GetMean())


            self.nu_pdg[0] = self.event.info['nu_pdg']
            self.nu_xs[0] = self.event.info['nu_xs']
            self.nu_proc[0] = self.event.info['nu_proc']
            self.nu_nucl[0] = self.event.info['nu_nucl']
            self.E_nu[0] = self.event.info['E_nu']
            self.E_depoTotal_recoil[0] = self.event.info['E_depoTotal_recoil']

            self.T_out.Fill()

        self.h_dedx_e_total.Write()
        self.h_dedx_e_2d_total.Write()
        self.h_dedx_p_total.Write()
        self.h_dedx_p_2d_total.Write()
        self.h_dedx_n_total.Write()
        self.h_dedx_n_2d_total.Write()
        self.h_dedx_pi_total.Write()
        self.h_dedx_pi_2d_total.Write()
        self.T_out.Write()

if __name__ == "__main__":
    if (len(sys.argv)>2):
        outfile = sys.argv[2]
    else:
        outfile = '/home/xning/output/output_test.root'
    E_field=float(sys.argv[3])
    event = Event(sys.argv[1],E_field)
    w = Writer(event, outfile)
    # w.Write()
    w.Write_dedx()
