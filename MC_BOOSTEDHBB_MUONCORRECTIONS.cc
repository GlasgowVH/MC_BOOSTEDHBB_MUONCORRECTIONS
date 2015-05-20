// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


    class MC_BOOSTEDHBB_MUONCORRECTIONS : public Analysis {
        public:

            /// Constructor
            MC_BOOSTEDHBB_MUONCORRECTIONS()
                : Analysis("MC_BOOSTEDHBB_MUONCORRECTIONS")
            {    }


            /// Book histograms and initialise projections before the run
            void init() {
                // what goes into complete jets:
                // hadrons, charged leptons, neutrinos, photons
                // with pt > 100 MeV
                FinalState allPartFS(Cuts::abseta < 3.0);

                // all muons
                IdentifiedFinalState allMuonFS(allPartFS);
                allMuonFS.acceptIdPair(13);
                addProjection(allMuonFS, "AllMuons");

                // prompt muons
                PromptFinalState promptMuonFS(allMuonFS);
                promptMuonFS.acceptTauDecays(true);
                addProjection(promptMuonFS, "PromptMuons");

                // nonprompt muons
                VetoedFinalState nonpromptMuonFS(allMuonFS);
                nonpromptMuonFS.addVetoOnThisFinalState(promptMuonFS);
                addProjection(nonpromptMuonFS, "NonpromptMuons");

                // all neutrinos
                IdentifiedFinalState allNeutrinoFS(allPartFS);
                allNeutrinoFS.acceptNeutrinos();
                addProjection(allNeutrinoFS, "AllNeutrinos");

                // prompt neutrinos
                PromptFinalState promptNeutrinoFS(allNeutrinoFS);
                promptNeutrinoFS.acceptTauDecays(true);
                addProjection(promptNeutrinoFS, "PromptNeutrinos");

                // nonprompt neutrinos
                VetoedFinalState nonpromptNeutrinoFS(allNeutrinoFS);
                nonpromptNeutrinoFS.addVetoOnThisFinalState(promptNeutrinoFS);
                addProjection(nonpromptNeutrinoFS, "NonpromptNeutrinos");


                // allJetPartFS: exclude only prompt neutrinos and muons
                VetoedFinalState allJetPartFS(allPartFS);
                allJetPartFS.addVetoOnThisFinalState(promptMuonFS);
                allJetPartFS.addVetoOnThisFinalState(promptNeutrinoFS);

                // caloJetPartFS: also exclude non-prompt neutrinos
                // and muons
                VetoedFinalState caloJetPartFS(allPartFS);
                caloJetPartFS.addVetoOnThisFinalState(promptMuonFS);
                caloJetPartFS.addVetoOnThisFinalState(promptNeutrinoFS);
                caloJetPartFS.addVetoOnThisFinalState(nonpromptMuonFS);
                caloJetPartFS.addVetoOnThisFinalState(nonpromptNeutrinoFS);

                // all charged particles with pt > 100 MeV
                ChargedFinalState trackPartFS(-2.5, 2.5, 100*MeV);

                // TODO
                // there's a problem with the second projection
                // veto of veto FS causing it?
                // various large-R jet collections
                addProjection(FastJets(allJetPartFS, FastJets::ANTIKT, 1.0), "Akt10All");
                addProjection(FastJets(caloJetPartFS, FastJets::ANTIKT, 1.0), "Akt10Calo");

                // TODO
                // there's a problem with the second projection
                // veto of veto FS causing it?
                // various small-R jet collections
                addProjection(FastJets(allJetPartFS, FastJets::ANTIKT, 0.4), "Akt04All");
                addProjection(FastJets(caloJetPartFS, FastJets::ANTIKT, 0.4), "Akt04Calo");

                // small-R track jets
                addProjection(FastJets(trackPartFS, FastJets::ANTIKT, 0.2), "Akt02Track");


                // book histograms
                akt10AllJetPt =
                    bookHisto1D("akt10AllJetPt", 50, 0, 1000*GeV,
                            "anti-$k_t$ $R=1.0$ jets",
                            "jet $p_T$ / GeV", "$\\frac{d\\sigma}{dp_T} / \\frac{\\rm pb}{\\rm GeV}$");

                akt10AllJetMass =
                    bookHisto1D("akt10AllJetMass", 50, 0, 500*GeV,
                            "anti-$k_t$ $R=1.0$ jets",
                            "jet mass / GeV", "$\\frac{d\\sigma}{dm} / \\frac{\\rm pb}{\\rm GeV}$");


                uncorrAkt10CaloJetPt =
                    bookHisto1D("uncorrAkt10CaloJetPt", 50, 0, 1000*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "jet $p_T$ / GeV", "$\\frac{d\\sigma}{dp_T} / \\frac{\\rm pb}{\\rm GeV}$");

                uncorrAkt10CaloJetMass =
                    bookHisto1D("uncorrAkt10CaloJetMass", 50, 0, 500*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "jet mass / GeV", "$\\frac{d\\sigma}{dm} / \\frac{\\rm pb}{\\rm GeV}$");

                corrAkt10CaloJetPt =
                    bookHisto1D("corrAkt10CaloJetPt", 50, 0, 1000*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "jet $p_T$ / GeV", "$\\frac{d\\sigma}{dp_T} / \\frac{\\rm pb}{\\rm GeV}$");

                corrAkt10CaloJetMass =
                    bookHisto1D("corrAkt10CaloJetMass", 50, 0, 500*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "jet mass / GeV", "$\\frac{d\\sigma}{dm} / \\frac{\\rm pb}{\\rm GeV}$");


                // book profiles
                uncorrAkt10CaloJetMeanPtResponseVsPt =
                    bookProfile1D("uncorrAkt10CaloJetMeanPtResponseVsPt", 50, 0, 1000*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< p_T\\ {\\rm response} \\right>$");

                uncorrAkt10CaloJetMeanMassResponseVsPt =
                    bookProfile1D("uncorrAkt10CaloJetMeanMassResponseVsPt", 50, 0, 1000*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< {\\rm mass\\ response} \\right>$");

                corrAkt10CaloJetMeanPtResponseVsPt =
                    bookProfile1D("corrAkt10CaloJetMeanPtResponseVsPt", 50, 0, 1000*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< p_T\\ {\\rm response} \\right>$");

                corrAkt10CaloJetMeanMassResponseVsPt =
                    bookProfile1D("corrAkt10CaloJetMeanMassResponseVsPt", 50, 0, 1000*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< {\\rm mass\\ response} \\right>$");


                uncorrAkt10CaloJetMeanPtResponseVsMass =
                    bookProfile1D("uncorrAkt10CaloJetMeanPtResponseVsMass", 50, 0, 500*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "full jet mass / GeV", "$\\left< p_T\\ {\\rm response} \\right>$");

                uncorrAkt10CaloJetMeanMassResponseVsMass =
                    bookProfile1D("uncorrAkt10CaloJetMeanMassResponseVsMass", 50, 0, 500*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "full jet mass / GeV", "$\\left< {\\rm mass\\ response} \\right>$");

                corrAkt10CaloJetMeanPtResponseVsMass =
                    bookProfile1D("corrAkt10CaloJetMeanPtResponseVsMass", 50, 0, 500*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "full jet mass / GeV", "$\\left< p_T\\ {\\rm response} \\right>$");

                corrAkt10CaloJetMeanMassResponseVsMass =
                    bookProfile1D("corrAkt10CaloJetMeanMassResponseVsMass", 50, 0, 500*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "full jet mass / GeV", "$\\left< {\\rm mass\\ response} \\right>$");


            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {

                const double weight = event.weight();

                const Particles& allMuons =
                    applyProjection<IdentifiedFinalState>(event, "AllMuons").particlesByPt();

                const Particles& promptMuons =
                    applyProjection<PromptFinalState>(event, "PromptMuons").particlesByPt();

                const Particles& nonpromptMuons =
                    applyProjection<VetoedFinalState>(event, "NonpromptMuons").particlesByPt();

                const Particles& allNeutrinos =
                    applyProjection<IdentifiedFinalState>(event, "AllNeutrinos").particlesByPt();

                const Particles& promptNeutrinos =
                    applyProjection<PromptFinalState>(event, "PromptNeutrinos").particlesByPt();

                const Particles& nonpromptNeutrinos =
                    applyProjection<VetoedFinalState>(event, "NonpromptNeutrinos").particlesByPt();


                // large-R jets
                const Jets& akt10AllJets =
                    applyProjection<FastJets>(event, "Akt10All").jetsByPt(250*GeV);

                const Jets& akt10CaloJets =
                    applyProjection<FastJets>(event, "Akt10Calo").jetsByPt(250*GeV);


                // small-R jets
                const Jets& akt04AllJets =
                    applyProjection<FastJets>(event, "Akt04All").jetsByPt(25*GeV);

                const Jets& akt04CaloJets =
                    applyProjection<FastJets>(event, "Akt04Calo").jetsByPt(25*GeV);


                // track jets
                const Jets& akt02TrackJets =
                    applyProjection<FastJets>(event, "Akt02Track").jetsByPt(7*GeV);

                // TODO!!
                // - match the Calo jets to the All jets
                // - store jet kinematics
                // - store muon kinematics
                // - store neutrino kinematics
                // - (all of the above relative to the full jet)
                // - find non-prompt muons in track jet using Rivet::Jet.containsParticle()
                // - add non-prompt muon momentum to calo jet momentum
                // - try adding a multiple of muon momentum to calo jet momentum
                // - etc.


                foreach (const Jet& uncorrJet, akt10CaloJets) {

                    // find the corresponding all jet
                    Jet allJet;
                    foreach (const Jet& jet, akt10AllJets) {
                        // TODO
                        // not really the best way to match
                        if (deltaR(uncorrJet, jet) < 0.2) {
                            allJet = jet;
                            break;
                        }
                    }

                    // TODO
                    // fine muons through trackjets

                    // find all corresponding trackjets
                    Jets assocTrackjets;
                    foreach (const Jet& jet, akt02TrackJets) {
                        if (deltaR(uncorrJet, jet) < 1.0)
                            assocTrackjets.push_back(jet);
                    }


                    // find all corresponding muons
                    Particles assocMuons;
                    foreach (const Particle& muon, nonpromptMuons) {
                        // TODO
                        // not really the best way to match
                        if (deltaR(uncorrJet, muon) < 1.0)
                            assocMuons.push_back(muon);
                    }

                    // find all corresponding neutrinos
                    Particles assocNeutrinos;
                    foreach (const Particle& neutrino, nonpromptNeutrinos) {
                        // TODO
                        // not really the best way to match
                        if (deltaR(uncorrJet, neutrino) < 1.0)
                            assocNeutrinos.push_back(neutrino);
                    }

                    // TODO
                    // how many of these muons and neutrinos aren't
                    // jet constituents?

                    // TODO
                    // for now just add the muons back and see what
                    // happens.

                    FourMomentum corrJetP4 = uncorrJet.mom();
                    foreach (const Particle& muon, assocMuons)
                        corrJetP4 += muon;


                    const Jet corrJet = Jet(corrJetP4);

                    double uncorrPt = uncorrJet.pt();
                    double corrPt = corrJet.pt();
                    double uncorrMass = uncorrJet.mass();
                    double corrMass = corrJet.mass();

                    // fill uncorrected and corrected histograms
                    uncorrAkt10CaloJetPt->fill(uncorrPt, weight);
                    corrAkt10CaloJetPt->fill(corrPt, weight);
                    uncorrAkt10CaloJetMass->fill(uncorrMass, weight);
                    corrAkt10CaloJetMass->fill(corrMass, weight);

                    // we didn't find a matched jet...
                    if (!allJet.pt())
                        continue;

                    double allPt = allJet.pt();
                    double allMass = allJet.mass();

                    double uncorrPtResponse = (uncorrPt - allPt) / allPt;
                    double corrPtResponse = (corrPt - allPt) / allPt;
                    double uncorrMassResponse = (uncorrMass - allMass) / allMass;
                    double corrMassResponse = (corrMass - allMass) / allMass;

                    uncorrAkt10CaloJetMeanPtResponseVsPt->fill(allPt, uncorrPtResponse, weight);
                    uncorrAkt10CaloJetMeanMassResponseVsPt->fill(allPt, uncorrMassResponse, weight);
                    corrAkt10CaloJetMeanPtResponseVsPt->fill(allPt, corrPtResponse, weight);
                    corrAkt10CaloJetMeanMassResponseVsPt->fill(allPt, corrMassResponse, weight);

                    uncorrAkt10CaloJetMeanPtResponseVsMass->fill(allMass, uncorrPtResponse, weight);
                    uncorrAkt10CaloJetMeanMassResponseVsMass->fill(allMass, uncorrMassResponse, weight);
                    corrAkt10CaloJetMeanPtResponseVsMass->fill(allMass, corrPtResponse, weight);
                    corrAkt10CaloJetMeanMassResponseVsMass->fill(allMass, corrMassResponse, weight);
                }

                foreach (const Jet& jet, akt10AllJets) {
                    // fill full jet histograms
                    akt10AllJetPt->fill(jet.pt(), weight);
                    akt10AllJetMass->fill(jet.mass(), weight);
                }

            }


            /// Normalise histograms etc., after the run
            void finalize() {

                // normalize histograms to cross section
                scale(akt10AllJetPt, crossSection()/sumOfWeights());
                scale(akt10AllJetMass, crossSection()/sumOfWeights());
                scale(uncorrAkt10CaloJetPt, crossSection()/sumOfWeights());
                scale(uncorrAkt10CaloJetMass, crossSection()/sumOfWeights());
                scale(corrAkt10CaloJetPt, crossSection()/sumOfWeights());
                scale(corrAkt10CaloJetMass, crossSection()/sumOfWeights());

            }

            //@}


    private:

            // Data members like post-cuts event weight counters go here

            Histo1DPtr akt10AllJetPt;
            Histo1DPtr akt10AllJetMass;
            Histo1DPtr uncorrAkt10CaloJetPt; 
            Histo1DPtr uncorrAkt10CaloJetMass; 
            Histo1DPtr corrAkt10CaloJetPt; 
            Histo1DPtr corrAkt10CaloJetMass; 

            Profile1DPtr uncorrAkt10CaloJetMeanPtResponseVsPt;
            Profile1DPtr uncorrAkt10CaloJetMeanMassResponseVsPt;
            Profile1DPtr corrAkt10CaloJetMeanPtResponseVsPt;
            Profile1DPtr corrAkt10CaloJetMeanMassResponseVsPt;
            Profile1DPtr uncorrAkt10CaloJetMeanPtResponseVsMass;
            Profile1DPtr uncorrAkt10CaloJetMeanMassResponseVsMass;
            Profile1DPtr corrAkt10CaloJetMeanPtResponseVsMass;
            Profile1DPtr corrAkt10CaloJetMeanMassResponseVsMass;

};



// The hook for the plugin system
DECLARE_RIVET_PLUGIN(MC_BOOSTEDHBB_MUONCORRECTIONS);


}
