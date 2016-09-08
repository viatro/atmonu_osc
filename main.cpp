#include "MyEarthAtm.h"
#include "nuSQuIDS/nuSQuIDS.h"

#include "Math/Functor.h"
#include "Math/Integrator.h"

#include "TStopwatch.h"

#include <array>
#include <tuple>
#include <map>
#include <iostream>

using namespace nusquids;

squids::Const units;

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
        static std::tuple<double, double, double> storage_key;
        static std::map<std::tuple<double, double, double>, std::array<double, 3> > storage;
        
        storage_key = std::make_tuple(zenith_angle, production_height, detector_depth);
        
        if (energy == storage_energy) {
            if (ini_state == storage_inistate) {
                if (storage.count(storage_key)) {
//std::cout << "." << std::flush;
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
        //Return probability for out_flv neutrino flavor
        
        for (unsigned int i = 0; i < 3; ++i) {
            storage[storage_key][i] = nus->EvalFlavor(i);
        }
//std::cout << " " << std::flush;
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

int main() {
    
    MyEarthAtmIntegrator meai(neutrino);
    TStopwatch timer;
    
    for (double EE : linspace(100.*units.MeV, 1000.*units.MeV, 9)) {
        std::cout << EE/units.MeV << " MeV:\t\t" << std::flush;
        timer.Start();
        for(unsigned int i = 0; i < 3; i++){
            std::cout << meai.EvalFlavorIntegrated(i, EE, {0., 1., 0.}) << " " << std::flush;
        }
        timer.Stop();
        std::cout << "\t\t=>\tcpu_time: " << timer.CpuTime() << "\treal_time: " << timer.RealTime() << std::endl;
    }
    

    return 0;
}
