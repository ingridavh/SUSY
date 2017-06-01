import AnalysisHelpers as AH
import ROOT
import Analysis
import math

#======================================================================
        
class susyAnalysis(Analysis.Analysis):
  """Semileptonic susyAnalysis loosely based on the ATLAS analyses of top pair events where 
  one W boson decays to leptons and one decays to hadrons.
  """
  def __init__(self, store):
    super(susyAnalysis, self).__init__(store)
  
  def initialize(self):
      self.hist_WtMass      =  self.addStandardHistogram("WtMass")

      self.hist_leptn       =  self.addStandardHistogram("lep_n")
      self.hist_leptpt      =  self.addStandardHistogram("lep_pt")
      self.hist_lepteta     =  self.addStandardHistogram("lep_eta")
      self.hist_leptE       =  self.addStandardHistogram("lep_E")
      self.hist_leptphi     =  self.addStandardHistogram("lep_phi")
      self.hist_leptch      =  self.addStandardHistogram("lep_charge")
      self.hist_leptID      =  self.addStandardHistogram("lep_type")
      self.hist_leptptc     =  self.addStandardHistogram("lep_ptconerel30")
      self.hist_leptetc     =  self.addStandardHistogram("lep_etconerel20")
      self.hist_lepz0       =  self.addStandardHistogram("lep_z0")
      self.hist_lepd0       =  self.addStandardHistogram("lep_d0")

      self.hist_njets       =  self.addStandardHistogram("n_jets")       
      self.hist_jetspt      =  self.addStandardHistogram("jet_pt")       
      self.hist_jetm        =  self.addStandardHistogram("jet_m")        
      self.hist_jetJVF      =  self.addStandardHistogram("jet_jvf")      
      self.hist_jeteta      =  self.addStandardHistogram("jet_eta")      
      self.hist_jetmv1      =  self.addStandardHistogram("jet_MV1")
      self.hist_jet1pt      =  self.addStandardHistogram("jet1_pt")
      self.hist_jet2pt      =  self.addStandardHistogram("jet2_pt")
      self.hist_jet3pt      =  self.addStandardHistogram("jet3_pt")

      self.hist_etmiss      = self.addStandardHistogram("etmiss")
      self.hist_vxp_z       = self.addStandardHistogram("vxp_z")
      self.hist_pvxp_n      = self.addStandardHistogram("pvxp_n")

      #Self-added
      self.hist_mt          = self.addStandardHistogram("mt")
      self.hist_meff        = self.addStandardHistogram("meff")
      self.hist_meratio     = self.addStandardHistogram("meratio")
  
  
  def analyze(self):
      # retrieving objects
      eventinfo = self.Store.getEventInfo()
      weight = eventinfo.scalefactor()*eventinfo.eventWeight() if not self.getIsData() else 1
      self.countEvent("all", weight)

      # apply standard event based selection
      if not AH.StandardEventCuts(eventinfo): return False
      self.countEvent("EventCuts", weight)

      # neutrinos are expected, so cut on missing transverse momentum
		#-----------------------------------Change cut according to signal 
      etmiss = self.Store.getEtMiss()

      #self.countEvent("MET", weight)
      
      # one good lepton from one of the W boson decays is expected, so require exactly one good lepton
      goodLeptons = AH.selectAndSortContainer(self.Store.getLeptons(), AH.isGoodLepton, lambda p: p.pt())
      if not (len(goodLeptons) == 1): return False
      self.countEvent("1 Lepton", weight)

      leadlepton = goodLeptons[0]
      
      # two jets from one of the W boson decays as well as two b-jets from the top pair decays are expected
      #-------------------------------------------------------
      #HER kan du bestemme antall jets
      goodJets = AH.selectAndSortContainer(self.Store.getJets(), AH.isGoodJet, lambda p: p.pt())

      #Require at least 3 jets
      if not (len(goodJets) >= 3) : return False
      self.countEvent("Jets", weight)

      #PT requirement on leading lepton
      
      if not (leadlepton.pt() > 25) : return False
      
      #lepton veto
      
      if len(goodLeptons) > 1 :
        if not (goodLeptons[-1].pt() > 10): return False
      

      #----------------------------------------------------
      #PT requirement on leading jets

      
      if not goodJets[0].pt() > 80: return False
      if not goodJets[1].pt() > 80: return False
      if not goodJets[2].pt() > 30: return False
      
      #Jet veto
      
      if len(goodJets) > 4:
        if goodJets[4].pt() > 40: return False
      

      #Cut on missing energy
      #change to 400GeV because of plot
      if not (etmiss.et() > 500.0): return False
      
      
      # apply the b-tagging requirement using the MV1 algorithm at 80% efficiency
		#-------------------------Not looking for b-jets      
      btags = sum([1 for jet in goodJets if jet.mv1() > 0.7892])
      #if not (btags >= 2): return False
      #self.countEvent("btags", weight)

      # apply a cut on the transverse mass of the W boson decaying to leptons
		#-----------------Change transverse mass mass according to sample
      #if not (AH.WTransverseMass(leadlepton, etmiss) > 30.0): return False

      #trasverse mass
      mt = math.sqrt( 2*leadlepton.pt() * etmiss.et() * (1-math.cos(leadlepton.phi()) ) )

      #inclusive effective mass and exclusive effective mass
      meffincl = etmiss.et()
      meffexcl = etmiss.et()
      
      for i in range(len(goodLeptons)) :
        meffincl += goodLeptons[i].pt()
        meffexcl +=  goodLeptons[i].pt()
        
      for j in range(len(goodJets)) : meffincl += goodJets[j].pt()
      for k in range(3) : meffexcl += goodJets[k].pt()

      #Ratio
      emratio = 0
      if not meffexcl == 0 :
        emratio = etmiss.et()/float(meffexcl)

      # cut on transverse mass
      
      if not (mt > 150) : return False
      

      #cut on ratio
      
      if not (emratio > 0.3) : return False
      

      #cut on effective inclusive mass
      
      if not (meffincl > 1400) : return False
      
      

      # Histograms detailing event information
      self.hist_vxp_z.Fill(eventinfo.primaryVertexPosition(), weight)
      self.hist_pvxp_n.Fill(eventinfo.numberOfVertices(), weight)

      # histograms for the W boson properties
      self.hist_WtMass.Fill(AH.WTransverseMass(leadlepton, etmiss), weight)

      # histograms for missing et
      self.hist_etmiss.Fill(etmiss.et(),weight)  

      # histograms detailing lepton information
      self.hist_leptn.Fill(len(goodLeptons), weight)
      self.hist_leptpt.Fill(leadlepton.pt(), weight)
      self.hist_lepteta.Fill(leadlepton.eta(), weight)
      self.hist_leptE.Fill(leadlepton.e(), weight)
      self.hist_leptphi.Fill(leadlepton.phi(), weight)
      self.hist_leptch.Fill(leadlepton.charge(), weight)
      self.hist_leptID.Fill(leadlepton.pdgId(), weight)
      self.hist_lepz0.Fill(leadlepton.z0(), weight)
      self.hist_lepd0.Fill(leadlepton.d0(), weight)      
      self.hist_leptptc.Fill(leadlepton.isoptconerel30(), weight)
      self.hist_leptetc.Fill(leadlepton.isoetconerel20(), weight)
      
      # histograms detailing jet information
      self.hist_njets.Fill(len(goodJets), weight)
      [self.hist_jetm.Fill(jet.m(), weight) for jet in goodJets]
      [self.hist_jetspt.Fill(jet.pt(), weight) for jet in goodJets]
      [self.hist_jetJVF.Fill(jet.jvf(), weight) for jet in goodJets]
      [self.hist_jeteta.Fill(jet.eta(), weight) for jet in goodJets]
      [self.hist_jetmv1.Fill(jet.mv1(), weight) for jet in goodJets]
      self.hist_jet1pt.Fill(goodJets[0].pt(), weight)
      self.hist_jet2pt.Fill(goodJets[1].pt(), weight)
      self.hist_jet3pt.Fill(goodJets[2].pt(), weight)

      #Histograms detailing self-added
      self.hist_mt.Fill(mt, weight) 
      self.hist_meff.Fill(meffincl, weight)
      self.hist_meratio.Fill(emratio, weight)   
      
      return True

  def finalize(self):
      pass
