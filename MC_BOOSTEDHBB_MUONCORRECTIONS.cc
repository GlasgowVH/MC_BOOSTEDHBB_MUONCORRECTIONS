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

                // various large-R jet collections
                addProjection(FastJets(allJetPartFS, FastJets::ANTIKT, 1.0), "AktAll10");
                addProjection(FastJets(caloJetPartFS, FastJets::ANTIKT, 1.0), "AktCalo10");

                // various small-R jet collections
                addProjection(FastJets(allJetPartFS, FastJets::ANTIKT, 0.4), "AktAll04");
                addProjection(FastJets(caloJetPartFS, FastJets::ANTIKT, 0.4), "AktCalo04");

                // small-R track jets
                addProjection(FastJets(trackPartFS, FastJets::ANTIKT, 0.2), "AktTrack02");
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


            /// @name Histograms
            //@{
            // some examples
            Profile1DPtr UncorrAkt10CaloJet_Mean_Pt_Response_Vs_Akt10AllJet_Pt;
            Profile1DPtr UncorrAkt10CaloJet_Mean_Mass_Response_Vs_Akt10AllJet_Pt;
            Profile1DPtr UncorrAkt04CaloJet_Mean_Pt_Response_Vs_Akt10AllJet_Pt;
            Profile1DPtr CorrAkt10CaloJet_Mean_Pt_Response_Vs_Akt10AllJet_Pt;
            Profile1DPtr CorrAkt10CaloJet_Mean_Mass_Response_Vs_Akt10AllJet_Pt;
            Profile1DPtr CorrAkt04CaloJet_Mean_Pt_Response_Vs_Akt04AllJet_Pt;
            //@}


};



// The hook for the plugin system
DECLARE_RIVET_PLUGIN(MC_BOOSTEDHBB_MUONCORRECTIONS);


}
