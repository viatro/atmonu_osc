#include "MyEarthAtm.h"
#include "nuSQuIDS/nuSQuIDS.h"

#include "Math/Functor.h"
#include "Math/Integrator.h"

#include "TFile.h"
#include "TNtupleD.h"
#include "TString.h"
#include "TStopwatch.h"

#include <array>
#include <map>
#include <iostream>
#include <string>
#include <vector>

#include "/home/victor/include/ezOptionParser-0.2.2/ezOptionParser.hpp"

using namespace nusquids;

squids::Const units;

namespace std
{
    template<typename T, size_t N>
    struct hash<array<T, N> >
    {
        typedef array<T, N> argument_type;
        typedef size_t result_type;

        result_type operator()(const argument_type& a) const
        {
            hash<T> hasher;
            result_type h = 0;
            for (result_type i = 0; i < N; ++i)
            {
                h = h * 31 + hasher(a[i]);
            }
            return h;
        }
    };
}

class MyEarthAtmIntegrator {
  
  public:
    
    MyEarthAtmIntegrator(NeutrinoType NT = neutrino) {
        nus = std::make_shared<nuSQUIDS>(3, NT);
        //Declaration of the body
        body = std::make_shared<MyEarthAtm>();
        //Definition of the track, in encodes the trajectory inside the body
        track = std::make_shared<MyEarthAtm::Track>(0., 0., 0.);
        //We set this in the nusSQuID object.
        nus->Set_Body(body);
        nus->Set_Track(track);
        
        // mixing angles
        nus->Set_MixingAngle(0, 1, 0.583995872);
        nus->Set_MixingAngle(0, 2, 0.148532030);
        nus->Set_MixingAngle(1, 2, 0.799399993);
        // square mass differences
        nus->Set_SquareMassDifference(1, 7.53e-5);
        nus->Set_SquareMassDifference(2, 2.42e-3);
        // CP phase
        nus->Set_CPPhase(0, 2, 0.0);
        
        //Set the maximum size for the integration step, important for fast or sharp variations of the density.
        //nus->Set_h_max( 1000.*units.km );
        //Set the GSL step function
        //nus->Set_GSL_step(gsl_odeiv2_step_rk2); // 45.27 s
        //nus->Set_GSL_step(gsl_odeiv2_step_rk4); // 138.19 s
        //nus->Set_GSL_step(gsl_odeiv2_step_rkf45); // 77.25 s
        //nus->Set_GSL_step(gsl_odeiv2_step_rkck); // 81.56 s
        //nus->Set_GSL_step(gsl_odeiv2_step_rk8pd); // 166.01 s
        //nus->Set_GSL_step(gsl_odeiv2_step_rk2imp); // *** Break *** segmentation violation
        //nus->Set_GSL_step(gsl_odeiv2_step_rk4imp); // *** Break *** segmentation violation
        //nus->Set_GSL_step(gsl_odeiv2_step_bsimp); // *** Break *** segmentation violation
        //nus->Set_GSL_step(gsl_odeiv2_step_rk1imp); // *** Break *** segmentation violation
        nus->Set_GSL_step(gsl_odeiv2_step_msadams); // 40.41 s
        //nus->Set_GSL_step(gsl_odeiv2_step_msbdf); // *** Break *** segmentation violation
        //Set the numerical precision of gsl integrator.
        nus->Set_rel_error(1.0e-10);
        nus->Set_abs_error(1.0e-10);
    }
    
    ~MyEarthAtmIntegrator() {}
    
    double EvalFlavorSingle(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.},
        double zenith_angle = 0., double production_height = 15.*units.km, double detector_depth = 0.)
    {
        static double storage_energy;
        static std::array<double, 3> storage_inistate;
        static std::array<double, 3> storage_key;
        static std::map<std::array<double, 3>, std::array<double, 3> > storage;
        
        storage_key = {zenith_angle, production_height, detector_depth};
        
        if (energy == storage_energy) {
            if (ini_state == storage_inistate) {
                if (storage.count(storage_key)) {
                    return storage[storage_key][out_flv];
                }
            } else {
                storage_inistate = ini_state;
                storage.clear();
            }
        } else {
            storage_energy = energy;
            storage_inistate = ini_state;
            storage.clear();
        }
        
        //Set energy
        nus->Set_E(energy);
        //Set initial state
        for (unsigned int i = 0; i < 3; ++i) {
            inistate[i] = ini_state[i];
        }
        //Set the initial state in nuSQuIDS object
        nus->Set_initial_state(inistate, flavor);
        //Set track parameters
        track->SetParams(zenith_angle, production_height, detector_depth);
        nus->Set_Track(track);
        //Propagate the neutrinos in the earth for the path defined in path
        nus->EvolveState();
        
        if (storage.size() < 50000000) {
            storage.emplace(storage_key, std::array<double, 3>{nus->EvalFlavor(0), nus->EvalFlavor(1), nus->EvalFlavor(2)});
        }
        //Return probability for out_flv neutrino flavor
        return nus->EvalFlavor(out_flv);
    }
    
    double ProductionHeightKmPDF(double height) {
        if (height < 0. || height > 86.) return 0.;
        return 0.0017024776552744395*(
                    86.49523161945041*exp(
                        0.15324051641391565*(14.484440464510044 - height)
                        - exp(0.15324051641391565*(14.484440464510044 - height))
                    ) + 1.4868453125378192*exp(-0.06437711874083203*height)
                );
    }
    
    double EvalFlavorIntegratedOverProductionHeight(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.}, double zenith_angle = 0., double detector_depth = 0.) {
        ROOT::Math::Functor1D wf([&](double production_height_in_km) {
            double hprob = ProductionHeightKmPDF(production_height_in_km);
            double prob_single = EvalFlavorSingle(out_flv, energy, ini_state, zenith_angle, production_height_in_km*units.km, detector_depth);
            double result = hprob*prob_single;
            return result;
        });
        ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR, 0, 0, 1e6);
        ig.SetFunction(wf);
        return ig.Integral(0.,86.);
    }
    
    double EvalFlavorIntegratedOverZenithAngle(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.}, double production_height = 15.*units.km, double detector_depth = 0.) {
        ROOT::Math::Functor1D wf([&](double zenith_angle) {
            double prob_single = EvalFlavorSingle(out_flv, energy, ini_state, zenith_angle, production_height, detector_depth);
            double result = prob_single/units.pi;
            return result;
        });
        ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR, 0, 0, 1e6);
        ig.SetFunction(wf);
        return ig.Integral(0.,units.pi);
    }
    
    double EvalFlavorIntegrated(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.}, double detector_depth = 0.)
    {
        ROOT::Math::Functor1D wf([&](double production_height_in_km) {
            double hprob = ProductionHeightKmPDF(production_height_in_km);
            return hprob*EvalFlavorIntegratedOverZenithAngle(out_flv, energy, ini_state, production_height_in_km*units.km, detector_depth);
        });
        ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR, 0, 0, 1e6);
        ig.SetFunction(wf);
        return ig.Integral(0.,86.);
    }
    
    double EvalFlavorIntegrated2(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.}, double detector_depth = 0.)
    {
        ROOT::Math::Functor1D wf([&](double zenith_angle) {
            return EvalFlavorIntegratedOverProductionHeight(out_flv, energy, ini_state, zenith_angle, detector_depth);
        });
        ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR, 0, 0, 1e6);
        ig.SetFunction(wf);
        return ig.Integral(0.,units.pi)/units.pi;
    }
    
    
  private:
    std::shared_ptr<nuSQUIDS> nus;
    std::shared_ptr<MyEarthAtm> body;
    std::shared_ptr<MyEarthAtm::Track> track;
    marray<double, 1> inistate = marray<double, 1>({3});
    
};

ez::ezOptionParser opt;

void CmdLineOptions(int argc, const char** argv) {
    opt.overview = ""; // General description in human language on what the user's tool does. (1st section)
    opt.syntax = "atmonu_osc [OPTIONS]"; // A synopsis of command and options usage to show expected order of input arguments. (2nd section)
    opt.example = "EXAMPLE"; // Example. (3rd section)
    opt.footer = "FOOTER"; // Final section printed in usage message. For contact, copyrights, version info.
    
    /*
    opt.add(
		Default,
		Required?,
		Number of args expected,
		Delimiter if expecting multiple args,
		Help description,
		Flag token,
		Flag token,
		Flag token,
		Flag token
	);
    */
    
    opt.add("", 0, 0,   0, "Display usage instructions.", "-h", "-help", "--help", "-?" );
    opt.add("", 1, 3, ',', "Set initial probabilities state.", "-i", "--inistate", "--initial-state" );
    opt.add("", 1, 3, ',', "Set energy for-loop parameters [in MeV].", "-e", "--energy", "--energy-loop" );
    opt.add("", 0, 2, ',', "Output ROOT file & TTree names.", "-o", "--output", "--output-file" );
    
    opt.parse(argc, argv);
    
    if (opt.isSet("-h")) {
		std::string usage;
        opt.getUsage(usage);
        std::cout << usage << std::endl;
		exit(0);
	}
    
    std::vector<std::string> badOptions;
	if(!opt.gotRequired(badOptions)) {
		for(size_t i = 0; i < badOptions.size(); ++i)
			std::cerr << "ERROR: Missing required option " << badOptions[i] << ".\n\n";
		std::string usage;
        opt.getUsage(usage);
        std::cout << usage << std::endl;
		exit(1);
	}
    
    if(!opt.gotExpected(badOptions)) {
		for(size_t i = 0; i < badOptions.size(); ++i)
			std::cerr << "ERROR: Got unexpected number of arguments for option " << badOptions[i] << ".\n\n";
		std::string usage;
        opt.getUsage(usage);
        std::cout << usage << std::endl;
		exit(1);
	}
}

int main(int argc, const char** argv) {
    
    CmdLineOptions(argc, argv);
    
    std::vector<double> OPT_inistate;
    opt.get("-i")->getDoubles(OPT_inistate);
    
    std::vector<double> OPT_energy_loop;
    opt.get("-e")->getDoubles(OPT_energy_loop);
    
    TString filename = TString::Format("atmonu_osc_integrated_inistate-%g,%g,%g_energies-%g-%g-%g_MeV.root",
                                        OPT_inistate[0], OPT_inistate[1], OPT_inistate[2],
                                        OPT_energy_loop[0], OPT_energy_loop[1], OPT_energy_loop[2]
                                      );
    TString treename = "atmonu_osc";
    
    if (opt.isSet("-o")) {
        std::vector<std::string> OPT_file_and_tree;
        opt.get("-o")->getStrings(OPT_file_and_tree);
        filename = OPT_file_and_tree[0];
        treename = OPT_file_and_tree[1];
    }
    
    TFile *file = TFile::Open(filename, "recreate");
    TNtupleD *tree = new TNtupleD(treename, "", "energy:ini_prob_nu_e:ini_prob_nu_mu:ini_prob_nu_tau:final_prob_nu_e:final_prob_nu_mu:final_prob_nu_tau:cpu_time");
    tree->SetAutoSave(1);
    
    MyEarthAtmIntegrator meai(neutrino);
    TStopwatch timer;
    std::array<double, 3> final_state;
    
    for (double EE = OPT_energy_loop[1]; EE >= OPT_energy_loop[0]; EE -= OPT_energy_loop[2]) {
        timer.Start();
        for(unsigned int i = 0; i < 3; i++){
            final_state[i] = meai.EvalFlavorIntegrated(i, EE*units.MeV, {OPT_inistate[0], OPT_inistate[1], OPT_inistate[2]});
        }
        timer.Stop();
        std::cout << TString::Format("%8.3f MeV:    %12.10f  %12.10f  %12.10f    =>  cputime %12.5f", EE, final_state[0], final_state[1], final_state[2], timer.CpuTime()) << std::endl;
        tree->Fill(EE, OPT_inistate[0], OPT_inistate[1], OPT_inistate[2], final_state[0], final_state[1], final_state[2], timer.CpuTime());
    }
    
    tree->Write(0,TObject::kOverwrite);
    file->Close();

    return 0;
}