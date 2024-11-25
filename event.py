import ROOT
from ROOT import TG4Event, TFile, TChain, TF1, TH2D, TH1D

import sys, os
import numpy as np

class Event:

    def __init__(self, fileName, eField=0.5):
        print("event: initilization")
        self.fileName = fileName
        self.ReadTree()
        # self.h = TH2D("h","length vs energy",500,0,0.6,500,0,20)
        # self.h2 = TH2D("h2","ratio vs energy",500,0,1,500,0,20)
        # self.h4 = TH2D("h4","ratio vs time",500,0,1,500,0,100)
        # self.h3 = TH2D("h3","dEdx vs residual",500,0,10,500,0,50)
        self.h_dedx_2d_e = TH2D("h_dedx_2d_e","dEdx_e",1000,0,20,1000,0,20)
        self.h_dedx_2d_p = TH2D("h_dedx_2d_p","dEdx_p",1000,0,20,1000,0,20)
        self.h_dedx_2d_n = TH2D("h_dedx_2d_n","dEdx_n",1000,0,20,1000,0,20)
        self.h_dedx_2d_pi = TH2D("h_dedx_2d_pi","dEdx_pi",1000,0,20,1000,0,20)
        self.h_dedx_e = TH1D("h_dedx_e","dEdx_e",1000,0,20)
        self.h_dedx_p = TH1D("h_dedx_p","dEdx_p",1000,0,20)
        self.h_dedx_n = TH1D("h_dedx_n","dEdx_n",1000,0,20)
        self.h_dedx_pi = TH1D("h_dedx_pi","dEdx_pi",1000,0,20)

        self.ly_array=np.asarray([156.599, 182.889, 187.899, 181.286, 159.422, 195.915, 180.879, 224.883, 215.799, 175.474, 179.844, 152.537, 158.442, 173.718, 164.869, 170.57, 206.106, 190.386, 189.51, 175.023, 172.214, 214.587, 168.861, 212.289, 192.466, 192.867, 176.751, 200.706, 206.586, 179.376, 200.951, 201.63, 200.69, 210.991, 170.035, 181.236, 192.317, 158.458, 193.091, 181.738, 179.475, 146.137, 186.001, 162.501, 180.432, 195.794, 184.454, 176.903, 176.154, 165.126])
        self.xbin=np.linspace(-6,6,num=(self.ly_array.ndim+1))

        self.currentEntry = 0
        self.A=0
        self.B=0
        self.E_field=eField #kV/cm
        # create a folder to store plots
        self.plotpath = "./plots"
        if not os.path.exists( self.plotpath):
            os.makedirs( self.plotpath)
            print("plotpath '" + self.plotpath + "' did not exist. It has been created!")
        # self.Jump(self.currentEntry)

    # ------------------------
    def ReadTree(self):
        # self.rootFile = TFile(self.fileName)
        self.simTree = TChain("EDepSimEvents")
        self.simTree.Add(self.fileName)
        self.nEntry = self.simTree.GetEntries()
        # Only Genie has gRooTracker, Marley doesn't
        try:
            edep_file = TFile(self.fileName, "OPEN")
            test = edep_file.Get("DetSimPassThru/gRooTracker")
            if not test:
                raise Exception("Sorry, no gRooTracker in file")
        except:
            # MARLEY
            print("Skip looking for gRooTracker directory.")
        else:
            # Genie
            self.genieTree = TChain("DetSimPassThru/gRooTracker")
            self.genieTree.Add(self.fileName)
            # self.genieTree = self.rootFile.Get("DetSimPassThru/gRooTracker")
            self.nGenieEntry = self.genieTree.GetEntries()
            if self.nEntry != self.nGenieEntry:
                print("Edep-sim tree and GENIE tree number of entries do not match!")
                sys.exit()

        self.event = TG4Event()
        self.simTree.SetBranchAddress("Event", self.event)

    # ------------------------
    def Jump(self, entryNo):
        # print(f'reading event {entryNo}/{self.nEntry}')

        self.currentEntry = entryNo
        self.simTree.GetEntry(entryNo)

        self.ReadVertex()
        self.ReadTracks()
        self.ReadEnergyDepo('SimEnergyDeposit')
        # Find tracks, correct energy along tracks.
        # self.ReadEnergyDepoByTrack()
        self.GetdEdxByTrack()

        self.info = {}
        self.info['E_nu'] = 0
        self.info['E_avail'] = 0
        #self.info['E_availList'] = np.zeros(6) # lepton, proton, neutron, pi+-, pi0, others.
        self.info['E_availList'] = np.zeros(7) # electron, proton, neutron, pi+-, pi0, others, muon.
        self.info['E_depoTotal'] = 0
        self.info['Q_depoTotal'] = 0 #XN
        self.info['Q_depoTotal_thre'] = 0
        self.info['E_depoTotal_re'] = 0
        self.info['E_depoTotal_l'] = 0
        self.info['E_depoTotal_recoil'] = 0
        # self.info['E_depoList'] = np.zeros(6) # lepton, proton, neutron, pi+-, pi0, others.
        # self.info['E_depoList_thre'] = np.zeros(6) #XN lepton, proton, neutron, pi+-, pi0, others.
        self.info['E_depoList'] = np.zeros(7) # electron, proton, neutron, pi+-, pi0, others, muon.
        self.info['Q_depoList'] = np.zeros(7) # electron, proton, neutron, pi+-, pi0, others, muon.
        self.info['Q_depoList_thre'] = np.zeros(7) #XN electron, proton, neutron, pi+-, pi0, others, muon.
        self.info['E_depoList_re'] = np.zeros(7) #XN electron, proton, neutron, pi+-, pi0, others, muon.
        self.info['E_depoList_re_track'] = np.zeros(7) #XN electron, proton, neutron, pi+-, pi0, others, muon.
        self.info['E_depoList_re_lep'] = np.zeros(7) #XN electron, proton, neutron, pi+-, pi0, others, muon.
        self.info['E_depoList_re_had'] = np.zeros(7) #XN electron, proton, neutron, pi+-, pi0, others, muon.
        self.info['E_depoList_l'] = np.zeros(7) #XN electron, proton, neutron, pi+-, pi0, others, muon.
        self.info['nu_pdg'] = 0
        self.info['nu_xs'] = self.vertex.GetCrossSection()
        self.info['nu_proc'], self.info['nu_nucl'] = self.GetReaction()
        self.FillEnergyInfo()
        try:
            self.ReadGenie()
        except:
            # MARLEY
            print("Jump: Marley events, assert info from file name")
            self.ReadMarley()

    # ------------------------
    def ReadGenie(self):
        # gRooTracker info
        # https://github.com/GENIE-MC/Generator/blob/master/src/Apps/gNtpConv.cxx#L1837
        # StdHepStatus: 0: initial state; 1: final state particles; others: intermediate transport
        # following assumes the first particle is always the neutrino.
        self.genieTree.GetEntry(self.currentEntry)
        self.info['nu_pdg'] = self.genieTree.StdHepPdg[0]
        self.info['E_nu'] = self.genieTree.StdHepP4[3]*1000

    # ------------------------
    def ReadMarley(self):
        # we are mostly looking at nue anyway
        self.info['nu_pdg'] = 12
        if self.GetnuPDGFromFileName() == 'nue': self.info['nu_pdg'] = 12
        if self.GetnuPDGFromFileName() == 'numu': self.info['nu_pdg'] = 14
        if self.GetnuPDGFromFileName() == 'anue': self.info['nu_pdg'] = -12
        if self.GetnuPDGFromFileName() == 'anumu': self.info['nu_pdg'] = -14
        self.info['E_nu'] = self.GetEnuFromFileName()

    # ------------------------
    def ReadVertex(self):
        primaries = np.array(self.event.Primaries)
        if (primaries.size != 1):
            print("Number of primaries not equal to 1 (not neutrino vertex)!")
            return

        self.vertex = primaries[0]

    #--------------------------
    def GetReaction(self):
        txt_list = self.vertex.GetReaction().split(';')
        proc = ''
        nucl = 0
        for x in txt_list[2:]:
            if 'proc:' in x:
                proc = x.replace('proc:', '').replace('Weak[','').replace('],', '')
            elif 'N:' in x:
                nucl = int(x.replace('N:', ''))
        # print(proc, nucl)
        proc_num = 0
        if (proc.startswith('CC')):
            proc_num = 10
            proc = proc.replace('CC', '')
        elif (proc.startswith('NC')):
            proc_num = 20
            proc = proc.replace('NC', '')
        else:
            proc_num = 30
        proc_dict = {
            'QES' : 1,
            'RES' : 2,
            'DIS' : 3,
            'COH' : 4,
            'MEC' : 5,
        }
        proc_num += proc_dict.get(proc, 0)
        return proc_num, nucl

    # ------------------------
    def ReadTracks(self):
        self.tracks = np.array(self.event.Trajectories)
        # print(self.tracks[100])
        self.tracks_tag = np.zeros(self.tracks.size)
        for i in range(self.tracks.size):
            self.tracks[i].energy = {}
            self.tracks[i].association = {}
            self.tracks[i].energy['depoTotal'] = 0
            self.tracks[i].energy['depoTotal_q'] = 0 #XN
            self.tracks[i].energy['depoTotal_q_thre'] = 0 #XN
            self.tracks[i].energy['depoTotal_e_re'] = 0 #XN
            self.tracks[i].energy['depoTotal_e_re_track'] = 0 #XN
            self.tracks[i].energy['depoTotal_e_re_had'] = 0 #XN
            self.tracks[i].energy['depoTotal_e_re_lep'] = 0 #XN
            self.tracks[i].energy['depoTotal_l'] = 0 #XN
            self.tracks[i].energy['depoTotal_recoil'] = 0 #XN

            self.tracks[i].association['depoList'] = []
            self.tracks[i].association['children'] = []
            self.tracks[i].association['ancestor'] = i
            self.tracks[i].association['steplength'] = []
            self.tracks[i].association['tracklength'] = 0

        for i in range(self.tracks.size):
            track = self.tracks[i]
            parId = track.GetParentId()
            if parId == -1: continue
            self.tracks[parId].association['children'].append(i)
            while parId != -1:
                track.association['ancestor'] = parId
                parId = self.tracks[parId].GetParentId()

    # ------------------------
    def ReadEnergyDepo(self, detName):
        self.depos = np.array(self.event.SegmentDetectors[detName])
        # E1=0
        # E2=0

        # add depo info to tracks
        # depoList = self.FindDepoListFromTrack(1)
        # depoEnergy = np.sum([depo.GetEnergyDeposit() for depo in self.depos[depoList]])
        # print('debug: muon deposit energy:', depoEnergy)

        for i, depo in enumerate(self.depos):
            # trkId_list = depo.Contrib
            trkId = depo.Contrib[0]
            edep = depo.GetEnergyDeposit()
            mm2cm=0.1
            step_length=depo.GetTrackLength()*mm2cm # tracklength in cm
            # dEdx = 0
            # if(step_length>0.01):
            #     dEdx=edep/step_length
            #     qdep=self.GetdQdx(dEdx)*step_length
            # else:
            #     qdep=edep*0.7 # using R(18[MeV/cm])=0.35
            
            # qdep_thre = qdep
            # thresthold=0.075 #MeV
            # if qdep < thresthold:
            #     qdep_thre=0
            
            # # reconstruct dE after threshold
            # edep_re=0
            # if step_length>=0.48 or ((qdep_thre<(1.3/0.6*step_length+0.2)) & (step_length>0.01)):
            #     E1+=qdep_thre
            #     dqdx_re=qdep_thre/step_length
            #     edep_re=self.GetdEdx(dqdx_re)*step_length
            #     # if(edep_re>0):
            #         # print(qdep_thre/edep_re)
            # else:
            #     E2+=qdep_thre
            #     edep_re=qdep_thre/0.7 # using R(MIP)=0.7
            # edep_l=0
            # #self.h.Fill(step_length,qdep_thre)
            # #self.h.Fill(step_length,edep/step_length)
                
            # edep_l=edep-qdep*0.83

            # edepSecond = depo.GetSecondaryDeposit()
            # trkLength = depo.GetTrackLength()
            track = self.tracks[trkId]
            track.association['depoList'].append(i)
            track.association['steplength'].append(depo.GetTrackLength()*mm2cm)
            track.association['tracklength'] += step_length 
            track.energy['depoTotal'] += edep
            # track.energy['depoTotal_q'] += qdep #XN
            # track.energy['depoTotal_q_thre'] += qdep_thre #XN
            # track.energy['depoTotal_e_re'] += edep_re #XN
            # track.energy['depoTotal_l'] += edep_l #XN
            # for signal_track in trkId_list:
            #         track = self.tracks[signal_track]
            #         track.association['depoList'].append(i)
        
        #print("legnth of depo: ",len(self.depos))
        # E_tot = np.sum([depo.GetEnergyDeposit() for depo in self.depos])
        # print('total deposit energy: ', E_tot)
        # print(E1,"vs",E2) 

    def ReadEnergyDepoByTrack(self):
        #find particle
        #trkId_proton=999999999999999
        trkId_lepton=[]
        for particle in self.vertex.Particles:
            if particle.GetPDGCode() in [11,-11,13,-13]:
                trkId_lepton.append(particle.GetTrackId())        

        threshold=0.075 # threshold is 75keV
        mm2cm = 0.1
        for i,track in enumerate(self.tracks):
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            step = np.asarray(track.association['steplength'])
            length = np.asarray(track.association['tracklength'])
            trackpdg = track.GetPDGCode()
            if len(depoList)>0:
                for j,di in enumerate(depoList):
                    depo = self.depos[di]
                    e = depo.GetEnergyDeposit()
                    steplength = depo.GetTrackLength()*mm2cm
                    residual = np.sum(step[j:])
                    # convert dE to dQ
                    if steplength>0.01:
                        dEdx = e/steplength
                        dQdx =self.GetdQdx(dEdx)
                        dQ = dQdx*steplength
                    else: 
                        dQ = 0.7*e

                    # Apply threshold
                    dQ_thre = dQ
                    if dQ<threshold:
                        dQ_thre=0
                    # a=0.66 # for lepton
                    # b=0.41 # for hadron
                    if length>2: #cm  identify tracks
                    # if steplength>0.01: #cm
                            dE_re = self.GetdEdx(dQ_thre/steplength)*steplength
                            track.energy['depoTotal_e_re_track'] += dE_re
                    elif ancestor in trkId_lepton:
                            dE_re =dQ_thre
                            track.energy['depoTotal_e_re_lep'] += dE_re
                    else:
                            dE_re = dQ_thre
                            track.energy['depoTotal_e_re_had'] += dE_re

                    # if trackpdg == 2212:
                    #     dedx = e/steplength
                    #     self.h3.Fill(residual,dedx)
                    # if ancestor== trkId_proton:
                    #     ratio = dQ/e
                    #     self.h.Fill(steplength,e)
                    #     self.h2.Fill(dQ/e,e)
                    #     if (e>5/0.6*steplength+1) & (e<8/0.6*steplength+2):
                    #         self.h4.Fill(dQ/e,depo.GetStart().T())
                    #     if dQ/e<0.6:
                    #         self.A +=e
                    #     else:
                    #         self.B +=e 

                    dL = e-dQ*0.83 #energy goes to light
                    track.energy['depoTotal_q'] += dQ #XN charge
                    track.energy['depoTotal_q_thre'] += dQ_thre #XN charge with threshold
                    track.energy['depoTotal_e_re'] += dE_re #XN Total reconstructed energy
                    track.energy['depoTotal_l'] += dL #XNlength  Light
                    # print(dQ,dL,e)
                    # print("proton: ",track.energy['depoTotal_q'],track.energy['depoTotal'],track.energy['depoTotal_q']/track.energy['depoTotal'])

    def GetdEdxByTrack(self):
        #find particle
        trkId_proton=[]
        trkId_neutron=[]
        trkId_cpion=[]
        trkId_lepton=[]
        for particle in self.vertex.Particles:
            # print("all", particle.GetTrackId())        
            if particle.GetPDGCode() in [11,-11,13,-13,111]:
                trkId_lepton.append(particle.GetTrackId())
                # print(particle.GetTrackId())    
            if particle.GetPDGCode() == 2212:
                trkId_proton.append(particle.GetTrackId()) 
                # print("proton", particle.GetTrackId())    
            if particle.GetPDGCode() == 2112:
                trkId_neutron.append(particle.GetTrackId()) 
            if particle.GetPDGCode() in [211, -211]:
                trkId_cpion.append(particle.GetTrackId()) 

        mm2cm = 0.1
        for i,track in enumerate(self.tracks):
            depoList = track.association['depoList']
            ancestor = track.association['ancestor']
            # print(ancestor)
            # step = np.asarray(track.association['steplength'])
            # trackpdg = track.GetPDGCode()
            if track.GetPDGCode()==2112:
                selfDepo = track.energy['depoTotal']
                track.energy['depoTotal_recoil'] += selfDepo
                if selfDepo!=0:
                    continue
            elif track.GetPDGCode()>1e9:
                mom = track.GetInitialMomentum()
                mass = mom.M()
                KE = mom.E() - mass 
                track.energy['depoTotal_recoil'] += KE



            if len(depoList)>0:
                for j,di in enumerate(depoList):
                    depo = self.depos[di]
                    e = depo.GetEnergyDeposit()
                    steplength = depo.GetTrackLength()*mm2cm
                    dEdx= e/steplength
                    #if e>0.075:
                    if steplength>0.01:
                        if ancestor in trkId_lepton:
                            self.h_dedx_e.Fill(dEdx)
                            self.h_dedx_2d_e.Fill(dEdx,e)
                            # print(dEdx,e,steplength) 
                        elif ancestor in trkId_proton :
                            self.h_dedx_p.Fill(dEdx)
                            self.h_dedx_2d_p.Fill(dEdx,e)
                        elif ancestor in trkId_neutron :
                            self.h_dedx_n.Fill(dEdx)
                            self.h_dedx_2d_n.Fill(dEdx,e)
                        elif ancestor in trkId_cpion :
                            self.h_dedx_pi.Fill(dEdx)
                            self.h_dedx_2d_pi.Fill(dEdx,e)
                            # print(dEdx,e,steplength) 
                        # elif ancestor in trkId_proton:
                        #     self.h_dedx_h.Fill(dEdx)
                        #     print(dEdx,e,steplength) 
      
    def CheckTH2D(self):
        outfile = ROOT.TFile.Open("proton_resi_dedx.root","RECREATE")

        # self.h.Write()
        # self.h2.Write()
        # self.h4.Write()
        # self.h3.Write()
        self.h_dedx_e.Write()
        self.h_dedx_h.Write()
        outfile.Close()
        # print(self.A, self.B)

    # ------------------------
    def FindDepoListFromTrack(self, trkId):
        x = [i for (i, depo)
             in enumerate(self.depos)
             if (trkId in depo.Contrib)]
            # if (trkId == depo.GetPrimaryId())]
        return x

    # ------------------------
    def PrintVertex(self):
        try:
            self.info['nu_pdg'] = self.genieTree.StdHepPdg[0]
            self.info['E_nu'] = self.genieTree.StdHepP4[3]*1000
            print(f"neutrino {self.info['nu_pdg']}: {self.info['E_nu']} MeV")
        except:
            print("PrintVertex: Marley events, assert info from file name")
            self.ReadMarley()

        posx = self.vertex.GetPosition().X()
        posy = self.vertex.GetPosition().Y()
        posz = self.vertex.GetPosition().Z()
        print(f"vertex @: ({posx}, {posy}, {posz}) [mm]")
        print(f"reaction: {self.vertex.GetReaction()}")
        # print(f"interaction #: {self.vertex.GetInteractionNumber()}")

        # print(f"{self.vertex.Particles.size()} particles at the vertex", )
        print(f'{"pdg":>8}{"name":>8}{"trkId":>6}{"mass":>10}{"KE":>10}')
        print(f'{"":>8}{"":>8}{"":>6}{"[MeV]":>10}{"[MeV]":>10}')
        print('-'*(8+8+6+10+10))
        for particle in self.vertex.Particles:
            trkId = particle.GetTrackId()
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            if trkId < 0: continue
            pdg = particle.GetPDGCode()
            name = particle.GetName()
            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            print(f'{pdg:>8d}{name:>8s}{trkId:>6d}{mass:>10.2f}{KE:>10.2f}')
        print('-'*(8+8+6+10+10))

        print(f'{self.info}')

    # ------------------------
    def FillEnergyInfo(self):

        for particle in self.vertex.Particles:
            trkId = particle.GetTrackId()
            # Skip negative trk id: in the case of Marley events,
            # this usually is the final nucleus before deexcitation that G4 doesn't track
            # the kinematics are not correct either
            if trkId < 0: continue
            pdg = particle.GetPDGCode()
            depoE = self.GetEnergyDepoWithDesendents(trkId)
            depoQ = self.GetEnergyDepoWithDesendents_Q(trkId)
            depoQ_thre = self.GetEnergyDepoWithDesendents_Q_thre(trkId) #XN
            depoE_re = self.GetEnergyDepoWithDesendents_e_re(trkId) #XN
            depoE_re_track = self.GetEnergyDepoWithDesendents_key(trkId,"depoTotal_e_re_track") #XN
            depoE_re_lep = self.GetEnergyDepoWithDesendents_key(trkId,"depoTotal_e_re_lep") #XN
            depoE_re_had = self.GetEnergyDepoWithDesendents_key(trkId,"depoTotal_e_re_had") #XN
            depoE_recoil = self.GetEnergyDepoWithDesendents_key(trkId,"depoTotal_recoil") #XN
            depoE_l = self.GetEnergyDepoWithDesendents_l(trkId) #XN
            mom = particle.GetMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            self.info['E_depoTotal'] += depoE
            self.info['Q_depoTotal'] += depoQ
            self.info['Q_depoTotal_thre'] += depoQ_thre
            self.info['E_depoTotal_re'] += depoE_re
            self.info['E_depoTotal_l'] += depoE_l
            self.info['E_depoTotal_recoil'] += depoE_recoil
            # fill E_availList: lepton, proton, neutron, pi+-, pi0, others.
            # if (pdg in [13, -13, 11, -11]):
            #     self.info['E_avail'] += (KE + mass)
            #     self.info['E_availList'][0] += (KE + mass)
            #     self.info['E_depoList'][0] += depoE
            # fill E_availList: electron, muon, proton, neutron, pi+-, pi0, others.
            if (pdg in [11, -11]):        #XN
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][0] += (KE + mass)
                self.info['E_depoList'][0] += depoE
                self.info['Q_depoList'][0] += depoQ
                self.info['Q_depoList_thre'][0] += depoQ_thre
                self.info['E_depoList_re'][0] += depoE_re
                self.info['E_depoList_re_track'][0] += depoE_re_track
                self.info['E_depoList_re_lep'][0] += depoE_re_lep
                self.info['E_depoList_re_had'][0] += depoE_re_had
                self.info['E_depoList_l'][0] += depoE_l
                # print("electron",depoQ_thre,depoE_re,KE,depoE_re/KE)
            elif (pdg in [13,-13]):                
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][6] += (KE + mass)
                self.info['E_depoList'][6] += depoE
                self.info['Q_depoList'][6] += depoQ
                self.info['Q_depoList_thre'][6] += depoQ_thre
                self.info['E_depoList_re'][6] += depoE_re
                self.info['E_depoList_re_track'][6] += depoE_re_track
                self.info['E_depoList_re_lep'][6] += depoE_re_lep
                self.info['E_depoList_re_had'][6] += depoE_re_had
                self.info['E_depoList_l'][6] += depoE_l
            elif (pdg == 2212):
                self.info['E_avail'] += KE
                self.info['E_availList'][1] += KE
                self.info['E_depoList'][1] += depoE
                self.info['Q_depoList'][1] += depoQ
                self.info['Q_depoList_thre'][1] += depoQ_thre
                self.info['E_depoList_re'][1] += depoE_re
                self.info['E_depoList_re_track'][1] += depoE_re_track
                self.info['E_depoList_re_lep'][1] += depoE_re_lep
                self.info['E_depoList_re_had'][1] += depoE_re_had
                self.info['E_depoList_l'][1] += depoE_l
                # if depoE!=0:
                    # print("in seperate particle:",depoQ/depoE,depoE,depoQ,self.info['E_depoList'][1],self.info['Q_depoList'][1])
            elif (pdg == 2112):
                self.info['E_avail'] += KE
                self.info['E_availList'][2] += KE
                self.info['E_depoList'][2] += depoE
                self.info['Q_depoList'][2] += depoQ
                self.info['Q_depoList_thre'][2] += depoQ_thre
                self.info['E_depoList_re'][2] += depoE_re
                self.info['E_depoList_re_track'][2] += depoE_re_track
                self.info['E_depoList_re_lep'][2] += depoE_re_lep
                self.info['E_depoList_re_had'][2] += depoE_re_had
                self.info['E_depoList_l'][2] += depoE_l
            elif (pdg in [211, -211]):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][3] += (KE + mass)
                self.info['E_depoList'][3] += depoE
                self.info['Q_depoList'][3] += depoQ
                self.info['Q_depoList_thre'][3] += depoQ_thre
                self.info['E_depoList_re'][3] += depoE_re
                self.info['E_depoList_re_track'][3] += depoE_re_track
                self.info['E_depoList_re_lep'][3] += depoE_re_lep
                self.info['E_depoList_re_had'][3] += depoE_re_had
                self.info['E_depoList_l'][3] += depoE_l
            elif (pdg == 111):
                self.info['E_avail'] += (KE + mass)
                self.info['E_availList'][4] += (KE + mass)
                self.info['E_depoList'][4] += depoE
                self.info['Q_depoList'][4] += depoQ
                self.info['Q_depoList_thre'][4] += depoQ_thre
                self.info['E_depoList_re'][4] += depoE_re
                self.info['E_depoList_re_track'][4] += depoE_re_track
                self.info['E_depoList_re_lep'][4] += depoE_re_lep
                self.info['E_depoList_re_had'][4] += depoE_re_had
                self.info['E_depoList_l'][4] += depoE_l
            else:
                self.info['E_avail'] += KE
                self.info['E_availList'][5] += KE
                self.info['E_depoList'][5] += depoE
                self.info['Q_depoList'][5] += depoQ
                self.info['Q_depoList_thre'][5] += depoQ_thre
                self.info['E_depoList_re'][5] += depoE_re
                self.info['E_depoList_re_track'][5] += depoE_re_track
                self.info['E_depoList_re_lep'][5] += depoE_re_lep
                self.info['E_depoList_re_had'][5] += depoE_re_had
                self.info['E_depoList_l'][5] += depoE_l

        # for track in self.tracks:
        #     pdg = track.GetPDGCode()
        #     depoE = track.energy['depoTotal']
        #     if (pdg in [13, -13, 11, -11]):
        #         self.info['E_depoList'][0] += depoE
        #     elif (pdg == 2212):
        #         self.info['E_depoList'][1] += depoE
        #     elif (pdg == 2112):
        #         self.info['E_depoList'][2] += depoE
        #     elif (pdg in [211, -211]):
        #         self.info['E_depoList'][3] += depoE
        #     elif (pdg == 111):
        #         self.info['E_depoList'][4] += depoE
        #     else:
        #         self.info['E_depoList'][5] += depoE


    # ------------------------
    def PrintDepo(self, i):
        print(f"{self.depos.size} depo points stored in total")
        depo = self.depos[i]
        contrib = depo.GetContributors()
        primaryId = depo.GetPrimaryId()
        edep = depo.GetEnergyDeposit()
        # edepSecond = depo.GetSecondaryDeposit()
        trkLength = depo.GetTrackLength()
        print(contrib, "%6d %.2f %.2f" %(primaryId, edep, trkLength))

    # ------------------------
    def PrintTrack(self, trkId):
        self.PrintTracks(trkId, trkId+1)
        track = self.tracks[trkId]
        children = track.association['children']
        print(f'children: {children}')
        for childId in children:
            child = self.tracks[childId]
            name = child.GetName()
            mom = child.GetInitialMomentum()
            KE = mom.E() - mom.M()
            print(f"{childId} {name}: {KE:.2f} MeV, {child.energy['depoTotal']:.2f} MeV,{child.energy['depoTotal_q']:.2f} MeV ,{child.association['children']}")

        print(f"{track.Points.size()} points stored in track {trkId}")
        for point in track.Points:
            x = point.GetPosition().X()
            y = point.GetPosition().Y()
            z = point.GetPosition().Z()
            t = point.GetPosition().T()
            print(f"{point.GetProcess()}, {point.GetSubprocess()}, {x}, {y}, {z}, {t}")


        depoList = track.association['depoList']     
        mm2m = 0.001
        mm2cm = 0.1
        print("depolist:",depoList) 
        # if len(depoList)>0:
        #         depo0=self.depos[depoList[0]]
        #         depoN=self.depos[depoList[-1]]
        #         X0=depo0.GetStart().X()
        #         Xn=depoN.GetStop().X()
        #         Y0=depo0.GetStart().Y()
        #         Yn=depoN.GetStop().Y()
        #         Z0=depo0.GetStart().Z()
        #         Zn=depoN.GetStop().Z()
        #         length_sq=(X0-Xn)*(X0-Xn)+(Y0-Yn)*(Y0-Yn)+(Z0-Zn)*(Z0-Zn)
        #         length = np.sqrt(length_sq)*mm2cm #cm
        #         print("tracklength = ",length,"cm") 
        for di in depoList:
            depo = self.depos[di]
            # x = (depo.GetStart().X() ) 
            # y = (depo.GetStart().Y() ) 
            # z = (depo.GetStart().Z() ) 
            # x_e= (depo.GetStop().X() ) 
            # y_e= (depo.GetStop().Y() ) 
            # z_e= (depo.GetStop().Z() ) 
            # t = (depo.GetStart().T() ) # ns
            # e = depo.GetEnergyDeposit()
            l = depo.GetTrackLength() *mm2cm # most are 0.5 mm
            print("tracklength = ",l,"cm") 
        #     print(e,l)  

    # ------------------------
    def PrintTracks(self, start=0, stop=-1):
        # print(f"{self.tracks.size} trajectories stored", )
        print(f"{'pdg':>8}{'name':>8}{'trkId':>6}{'parId':>6}{'acId':>6}{'KE':>10}{'selfDepo':>10}{'allDepo':>10}")
        print(f"{'':>8}{'':>8}{'':>6}{'':>6}{'':>6}{'[MeV]':>10}{'[MeV]':>10}{'[MeV]':>10}")
        print('-'*(8+8+6+6+6+10+10+10))

        for track in self.tracks[start:stop]:
            pdg = track.GetPDGCode()
            name = track.GetName()
            trkId = track.GetTrackId()
            parId = track.GetParentId()
            mom = track.GetInitialMomentum()
            mass = mom.M()
            KE = mom.E() - mass
            ancestor = track.association['ancestor']
            selfDepo = track.energy['depoTotal']
            allDepo = self.GetEnergyDepoWithDesendents(trkId)
            print(f"{pdg:>8d}{name:>8s}{trkId:>6d}{parId:>6d}{ancestor:>6d}{KE:>10.2f}{selfDepo:>10.2f}{allDepo:>10.2f}")

        print('-'*(8+8+6+6+6+10+10+10))


    # ------------------------
    def Next(self):
        if self.currentEntry != self.nEntry -1:
            self.currentEntry += 1
        else:
            self.currentEntry = 0
        self.Jump(self.currentEntry)

    # ------------------------
    def Prev(self):
        if self.currentEntry != 0:
            self.currentEntry -= 1
        else:
            self.currentEntry = self.nEntry -1
        self.Jump(self.currentEntry)

    #-------------------------
    def GetEnergyDepoWithDesendents(self, trkId):
        track = self.tracks[trkId]
        energy = track.energy['depoTotal']
        children = track.association['children']
        for childId in children:
            energy += self.GetEnergyDepoWithDesendents(childId)
        return energy

    #XN-----------------------------
    def GetEnergyDepoWithDesendents_Q_thre(self, trkId):
        track = self.tracks[trkId]
        energy = track.energy['depoTotal_q_thre']
        children = track.association['children']
        for childId in children:
            energy += self.GetEnergyDepoWithDesendents_Q_thre(childId)
        return energy

    def GetEnergyDepoWithDesendents_e_re(self, trkId):
        track = self.tracks[trkId]
        energy = track.energy['depoTotal_e_re']
        children = track.association['children']
        for childId in children:
            energy += self.GetEnergyDepoWithDesendents_e_re(childId)
        return energy
    def GetEnergyDepoWithDesendents_key(self, trkId,key):
        track = self.tracks[trkId]
        energy = track.energy[key]
        children = track.association['children']
        for childId in children:
            energy += self.GetEnergyDepoWithDesendents_key(childId,key)
        return energy

    def GetEnergyDepoWithDesendents_l(self, trkId):
        track = self.tracks[trkId]
        energy = track.energy['depoTotal_l']
        children = track.association['children']
        for childId in children:
            energy += self.GetEnergyDepoWithDesendents_l(childId)
        return energy


    def GetEnergyDepoWithDesendents_Q(self, trkId):
        track = self.tracks[trkId]
        energy = track.energy['depoTotal_q']
        children = track.association['children']
        for childId in children:
            energy += self.GetEnergyDepoWithDesendents_Q(childId)
        return energy
    
    
    def BirksModel(self, dEdx):
        A_3t=0.8
        k_3t=0.0486 # [(g/MeVcm2)(kV/cm)]
        # E=0.5 # [kV/cm]
        E=self.E_field
        # print("e field = ",E)

        rho=1.38 # [g/cm3]
        ips=E*rho
        R=A_3t/(1+k_3t/ips*dEdx)
        return R
    
    def BirksModel_inverse(self, dQdx):
        A_3t=0.8
        k_3t=0.0486 # [(g/MeVcm2)(kV/cm)]
        # E=0.5 # [kV/cm]
        E=self.E_field
        rho=1.38 # [g/cm3]
        ips=E*rho
        R=1/(A_3t-k_3t/ips*dQdx)
        return R
    
    def GetdQdx(self,dEdx):
        return dEdx*self.BirksModel(dEdx)
        # return dEdx*0.7

    def GetdEdx(self,dQdx):
        return dQdx*self.BirksModel_inverse(dQdx)
        # return dQdx/0.7

    #-------------------------
    def GetEnergyDepoWithAncestor(self, acId):
        energy = 0
        for track in self.tracks:
            ancestor = track.association['ancestor']
            # print(ancestor)
            if ancestor == acId:
                energy += track.energy['depoTotal']
        return energy

    #-----------------------
    def GetEnuFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-2]
        filename = filename.replace('MeV', '')
        return float(filename)  # Marley file names already in MeV

    #-----------------------
    def GetnuPDGFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-3]
        return filename

    #-----------------------
    def GetFileName(self):
        return self.simTree.GetFile().GetName()
# ------------------------
    #-------------------------
    # def GetEnergyDepoWithAncestor(self, acId):
    #     energy = 0
    #     for track in self.tracks:
    #         ancestor = track.association['ancestor']
    #         if ancestor == acId:
    #             energy += track.energy['depoTotal']
    #     return energy

    #-----------------------
    def GetEnuFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-2]
        filename = filename.replace('MeV', '')
        return float(filename)  # Marley file names already in MeV

    #-----------------------
    def GetnuPDGFromFileName(self):
        # This is used for Marley low energy files
        # No need for Genie, can read from gRooTracker tree directly
        filename = self.GetFileName().split('/')[-1]
        filename = filename.split('_')[-3]
        return filename

    #-----------------------
    def GetFileName(self):
        return self.simTree.GetFile().GetName()
# ------------------------
if __name__ == "__main__":
    event = Event(sys.argv[1])
    event.Jump(0)
    event.PrintVertex()
    event.Next()
    event.PrintVertex()
    event.PrintTracks(0,8)
    # event.PrintTrack(1)
