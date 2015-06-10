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

                // various large-R jet collections
                FastJets akt10InclusiveJetProj (allJetPartFS, FastJets::ANTIKT, 1.0);
                akt10InclusiveJetProj.useInvisibles(true);
                addProjection(akt10InclusiveJetProj, "Akt10Inclusive");

                FastJets akt10CaloJetProj (caloJetPartFS, FastJets::ANTIKT, 1.0);
                akt10CaloJetProj.useInvisibles(true);
                addProjection(akt10CaloJetProj, "Akt10Calo");


                // various small-R jet collections
                FastJets akt04InclusiveJetProj (allJetPartFS, FastJets::ANTIKT, 0.4);
                akt04InclusiveJetProj.useInvisibles(true);
                addProjection(akt04InclusiveJetProj, "Akt04Inclusive");

                FastJets akt04CaloJetProj (caloJetPartFS, FastJets::ANTIKT, 0.4);
                akt04CaloJetProj.useInvisibles(true);
                addProjection(akt04CaloJetProj, "Akt04Calo");

                // small-R track jets
                addProjection(FastJets(trackPartFS, FastJets::ANTIKT, 0.2), "Akt02Track");

                // book histograms
                akt10InclusiveJetPt =
                    bookHisto1D("akt10InclusiveJetPt", 50, 0, 1000*GeV,
                            "anti-$k_t$ $R=1.0$ jets",
                            "jet $p_T$ / GeV", "$\\frac{d\\sigma}{dp_T} / \\frac{\\rm pb}{\\rm GeV}$");

                akt10InclusiveJetMass =
                    bookHisto1D("akt10InclusiveJetMass", 50, 0, 500*GeV,
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
                const Jets& akt10InclusiveJets =
                    applyProjection<FastJets>(event, "Akt10Inclusive").jetsByPt(250*GeV);

                const Jets& akt10CaloJets =
                    applyProjection<FastJets>(event, "Akt10Calo").jetsByPt(250*GeV);


                // small-R jets
                const Jets& akt04InclusiveJets =
                    applyProjection<FastJets>(event, "Akt04Inclusive").jetsByPt(25*GeV);

                const Jets& akt04CaloJets =
                    applyProjection<FastJets>(event, "Akt04Calo").jetsByPt(25*GeV);


                // track jets
                const Jets& akt02TrackJets =
                    applyProjection<FastJets>(event, "Akt02Track").jetsByPt(7*GeV);

                // TODO!!
                // - match the Calo jets to the Inclusive jets
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
                    Jet inclusiveJet;
                    foreach (const Jet& jet, akt10InclusiveJets) {
                        // TODO
                        // not really the best way to match
                        if (deltaR(uncorrJet, jet) < 0.2) {
                            inclusiveJet = jet;
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
                        // only identify muons with |eta| < 2.5
                        if (muon.abseta() > 2.5)
                            continue;

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
                    if (!inclusiveJet.pt() || !inclusiveJet.mass())
                        continue;

                    double inclusivePt = inclusiveJet.pt();
                    double inclusiveMass = inclusiveJet.mass();

                    double uncorrPtResponse = (uncorrPt - inclusivePt) / inclusivePt;
                    double corrPtResponse = (corrPt - inclusivePt) / inclusivePt;
                    double uncorrMassResponse = (uncorrMass - inclusiveMass) / inclusiveMass;
                    double corrMassResponse = (corrMass - inclusiveMass) / inclusiveMass;

                    uncorrAkt10CaloJetMeanPtResponseVsPt->fill(inclusivePt, uncorrPtResponse, weight);
                    uncorrAkt10CaloJetMeanMassResponseVsPt->fill(inclusivePt, uncorrMassResponse, weight);
                    corrAkt10CaloJetMeanPtResponseVsPt->fill(inclusivePt, corrPtResponse, weight);
                    corrAkt10CaloJetMeanMassResponseVsPt->fill(inclusivePt, corrMassResponse, weight);

                    uncorrAkt10CaloJetMeanPtResponseVsMass->fill(inclusiveMass, uncorrPtResponse, weight);
                    uncorrAkt10CaloJetMeanMassResponseVsMass->fill(inclusiveMass, uncorrMassResponse, weight);
                    corrAkt10CaloJetMeanPtResponseVsMass->fill(inclusiveMass, corrPtResponse, weight);
                    corrAkt10CaloJetMeanMassResponseVsMass->fill(inclusiveMass, corrMassResponse, weight);
                }

                foreach (const Jet& jet, akt10InclusiveJets) {
                    // fill full jet histograms
                    akt10InclusiveJetPt->fill(jet.pt(), weight);
                    akt10InclusiveJetMass->fill(jet.mass(), weight);
                }

            }


            /// Normalise histograms etc., after the run
            void finalize() {

                // normalize histograms to cross section
                scale(akt10InclusiveJetPt, crossSection()/sumOfWeights());
                scale(akt10InclusiveJetMass, crossSection()/sumOfWeights());
                scale(uncorrAkt10CaloJetPt, crossSection()/sumOfWeights());
                scale(uncorrAkt10CaloJetMass, crossSection()/sumOfWeights());
                scale(corrAkt10CaloJetPt, crossSection()/sumOfWeights());
                scale(corrAkt10CaloJetMass, crossSection()/sumOfWeights());

            }

            //@}


    private:

            // Data members like post-cuts event weight counters go here

            Histo1DPtr akt10InclusiveJetPt;
            Histo1DPtr akt10InclusiveJetMass;
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
