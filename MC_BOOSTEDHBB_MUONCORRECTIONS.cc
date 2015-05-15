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

                // only prompt particles (not coming from hadron
                // decays)
                // particles from prompt taus are considered prompt.
                PromptFinalState promptFS(allPartFS);
                promptFS.acceptTauDecays(true);

                // only non-prompt particles (coming from hadron
                // decays)
                VetoedFinalState nonpromptFS(allPartFS);
                nonpromptFS.addVetoOnThisFinalState(promptFS);

                // prompt muons
                IdentifiedFinalState promptMuonFS(promptFS);
                promptMuonFS.acceptIdPair(13);
                addProjection(promptMuonFS, "PromptMuons");

                // prompt neutrinos
                IdentifiedFinalState promptNeutrinoFS(promptFS);
                promptNeutrinoFS.acceptNeutrinos();
                addProjection(promptNeutrinoFS, "PromptNeutrinos");

                // nonprompt muons
                IdentifiedFinalState nonpromptMuonFS(nonpromptFS);
                nonpromptMuonFS.acceptIdPair(13);
                addProjection(nonpromptMuonFS, "NonpromptMuons");

                // nonprompt neutrinos
                IdentifiedFinalState nonpromptNeutrinoFS(nonpromptFS);
                nonpromptNeutrinoFS.acceptNeutrinos();
                addProjection(nonpromptNeutrinoFS, "NonpromptNeutrinos");

                MergedFinalState allMuonFS(promptMuonFS, nonpromptMuonFS);
                addProjection(allMuonFS, "AllMuons");

                MergedFinalState allNeutrinoFS(promptNeutrinoFS, nonpromptNeutrinoFS);
                addProjection(allNeutrinoFS, "AllNeutrinos");


                // allJetPartFS: exclude only prompt neutrinos and muons
                VetoedFinalState allJetPartFS(allPartFS);
                allJetPartFS.addVetoOnThisFinalState(promptMuonFS);
                allJetPartFS.addVetoOnThisFinalState(promptNeutrinoFS);

                // caloJetPartFS: also exclude non-prompt neutrinos
                // and muons
                VetoedFinalState caloJetPartFS(allJetPartFS);
                caloJetPartFS.addVetoOnThisFinalState(nonpromptMuonFS);
                caloJetPartFS.addVetoOnThisFinalState(nonpromptNeutrinoFS);

                // all charged particles with pt > 100 MeV
                ChargedFinalState trackPartFS(-2.5, 2.5, 100*MeV);

                // TODO
                // there's a problem with the second projection
                // veto of veto FS causing it?
                // various large-R jet collections
                addProjection(FastJets(allJetPartFS, FastJets::ANTIKT, 1.0), "AktAll10");
                addProjection(FastJets(caloJetPartFS, FastJets::ANTIKT, 1.0), "AktCalo10");

                // TODO
                // there's a problem with the second projection
                // veto of veto FS causing it?
                // various small-R jet collections
                addProjection(FastJets(allJetPartFS, FastJets::ANTIKT, 0.4), "AktAll04");
                addProjection(FastJets(caloJetPartFS, FastJets::ANTIKT, 0.4), "AktCalo04");

                // small-R track jets
                addProjection(FastJets(trackPartFS, FastJets::ANTIKT, 0.2), "AktTrack02");


                // book histograms
                uncorrAkt10CaloJetMeanPtResponseVsPt =
                    bookProfile1D("uncorrAkt10CaloJetMeanPtResponseVsPt", 50, 0, 500*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< p_{T,{\\rm jet}} \\right>$ / GeV");

                uncorrAkt10CaloJetMeanMassResponseVsPt =
                    bookProfile1D("uncorrAkt10CaloJetMeanMassResponseVsPt", 50, 0, 500*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< m_{\\rm jet} \\right>$ / GeV");

                uncorrAkt10CaloJetMeanPtResponseVsMass =
                    bookProfile1D("uncorrAkt10CaloJetMeanPtResponseVsMass", 50, 0, 500*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< p_{T,{\\rm jet}} \\right>$ / GeV");

                uncorrAkt10CaloJetMeanMassResponseVsMass =
                    bookProfile1D("uncorrAkt10CaloJetMeanMassResponseVsMass", 50, 0, 500*GeV,
                            "uncorrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< m_{\\rm jet} \\right>$ / GeV");

                corrAkt10CaloJetMeanPtResponseVsPt =
                    bookProfile1D("corrAkt10CaloJetMeanPtResponseVsPt", 50, 0, 500*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< p_{T,{\\rm jet}} \\right>$ / GeV");

                corrAkt10CaloJetMeanMassResponseVsPt =
                    bookProfile1D("corrAkt10CaloJetMeanMassResponseVsPt", 50, 0, 500*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< m_{\\rm jet} \\right>$ / GeV");

                corrAkt10CaloJetMeanPtResponseVsMass =
                    bookProfile1D("corrAkt10CaloJetMeanPtResponseVsMass", 50, 0, 500*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< p_{T,{\\rm jet}} \\right>$ / GeV");

                corrAkt10CaloJetMeanMassResponseVsMass =
                    bookProfile1D("corrAkt10CaloJetMeanMassResponseVsMass", 50, 0, 500*GeV,
                            "corrected anti-$k_t$ $R=1.0$ jets",
                            "full jet $p_T$ / GeV", "$\\left< m_{\\rm jet} \\right>$ / GeV");


            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {

                const double weight = event.weight();

                const Particles& allMuons =
                    applyProjection<IdentifiedFinalState>(event, "AllMuons").particlesByPt();

                const Particles& promptMuons =
                    applyProjection<IdentifiedFinalState>(event, "PromptMuons").particlesByPt();

                const Particles& nonpromptMuons =
                    applyProjection<IdentifiedFinalState>(event, "NonpromptMuons").particlesByPt();

                const Particles& allNeutrinos =
                    applyProjection<IdentifiedFinalState>(event, "AllNeutrinos").particlesByPt();

                const Particles& promptNeutrinos =
                    applyProjection<IdentifiedFinalState>(event, "PromptNeutrinos").particlesByPt();

                const Particles& nonpromptNeutrinos =
                    applyProjection<IdentifiedFinalState>(event, "NonpromptNeutrinos").particlesByPt();


                // large-R jets
                const Jets& aktAll10Jets =
                    applyProjection<FastJets>(event, "AktAll10").jetsByPt(250*GeV);

                const Jets& aktCalo10Jets =
                    applyProjection<FastJets>(event, "AktCalo10").jetsByPt(250*GeV);


                // small-R jets
                const Jets& aktAll04Jets =
                    applyProjection<FastJets>(event, "AktAll04").jetsByPt(25*GeV);

                const Jets& aktCalo04Jets =
                    applyProjection<FastJets>(event, "AktCalo04").jetsByPt(25*GeV);


                // track jets
                const Jets& aktTrack02Jets =
                    applyProjection<FastJets>(event, "AktTrack02").jetsByPt(7*GeV);

                // TODO!!
                // - match the Calo jets to the All jets
                // - find non-prompt muons in track jet using Rivet::Jet.containsParticle()
                // - add non-prompt muon momentum to calo jet momentum
                // - try adding a multiple of muon momentum to calo jet momentum
                // - etc.


                foreach (const Jet& uncorrJet, aktCalo10Jets) {

                    // find the corresponding full jet
                    Jet fullJet;
                    foreach (const Jet& jet, aktAll10Jets) {
                        // TODO
                        // not really the best way to match
                        if (deltaR(uncorrJet, jet) < 0.2) {
                            fullJet = jet;
                            break;
                        }
                    }

                    // TODO
                    // fine muons through trackjets

                    // find all corresponding trackjets
                    Jets assocTrackjets;
                    foreach (const Jet& jet, aktTrack02Jets) {
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

                    // find all corresponding muons
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
                    double fullPt = fullJet.pt();
                    double uncorrMass = uncorrJet.mass();
                    double corrMass = corrJet.mass();
                    double fullMass = fullJet.mass();

                    uncorrAkt10CaloJetMeanPtResponseVsPt->fill(fullPt, uncorrPt);
                    uncorrAkt10CaloJetMeanMassResponseVsPt->fill(fullPt, uncorrMass);
                    corrAkt10CaloJetMeanPtResponseVsPt->fill(fullPt, corrPt);
                    corrAkt10CaloJetMeanMassResponseVsPt->fill(fullPt, corrMass);

                    uncorrAkt10CaloJetMeanPtResponseVsMass->fill(fullMass, uncorrPt);
                    uncorrAkt10CaloJetMeanMassResponseVsMass->fill(fullMass, uncorrMass);
                    corrAkt10CaloJetMeanPtResponseVsMass->fill(fullMass, corrPt);
                    corrAkt10CaloJetMeanMassResponseVsMass->fill(fullMass, corrMass);
                }

            }


            /// Normalise histograms etc., after the run
            void finalize() {

                /// @todo Normalise, scale and otherwise manipulate histograms here

                // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
                // normalize(_h_YYYY); // normalize to unity

            }

            //@}


    private:

            // Data members like post-cuts event weight counters go here

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
