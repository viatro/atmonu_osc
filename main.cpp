#include "MyEarthAtm.h"
#include "nuSQuIDS/nuSQuIDS.h"

#include "Math/Functor.h"
#include "Math/IntegratorMultiDim.h"

#include "Math/Integrator.h"

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
        nus->Set_h_max( 1000.*units.km );
        //Set the GSL step function
        //nus.Set_GSL_step(gsl_odeiv2_step_rk4);
        //Set the numerical precision of gsl integrator.
        nus->Set_rel_error(1.0e-10);
        nus->Set_abs_error(1.0e-10);
    }
    
    ~MyEarthAtmIntegrator() {}
    
    double EvalFlavorSingle(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.},
        double zenith_angle = units.pi, double production_height = 22.*units.km, double detector_depth = 0.)
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
    
    double EvalFlavorIntegrated(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.}, double detector_depth = 0.)
    {
        double limits_zenith_angle[2] = {0.,units.pi};
        double limits_production_height[2] = {0., 86.};
        ROOT::Math::Functor wf([&](const double *x) {
            double hprob = ProductionHeightKmPDF(x[1]/units.km);
            double zang_range = limits_zenith_angle[1]-limits_zenith_angle[0];
            double prob_single = EvalFlavorSingle(out_flv, energy, ini_state, x[0], x[1], detector_depth);
            double result = hprob/zang_range*prob_single;
std::cout << "lambda: " << "zenith_angle=" << x[0]*180./units.pi << "\t\tproduction_height=" << x[1]/units.km << "\t\thprob=" << hprob << "\tprob_single=" << prob_single << "\tresult=" << result << std::endl;
            return result;
        }, 2);
        ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
        ig.SetFunction(wf);
        return ig.Integral(limits_zenith_angle, limits_production_height);
    }
    
    double EvalFlavorIntegratedOverProductionHeight(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.}, double zenith_angle = units.pi, double detector_depth = 0.) {
        ROOT::Math::Functor1D wf([&](double production_height_in_km) {
            double hprob = ProductionHeightKmPDF(production_height_in_km);
            double prob_single = EvalFlavorSingle(out_flv, energy, ini_state, zenith_angle, production_height_in_km*units.km, detector_depth);
            double result = hprob*prob_single;
            return result;
        });
        ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
        ig.SetFunction(wf);
        return ig.Integral(0.,100.);
    }
    
    double EvalFlavorIntegratedOverZenithAngle(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.}, double production_height = 22.*units.km, double detector_depth = 0.) {
        ROOT::Math::Functor1D wf([&](double zenith_angle) {
            double prob_single = EvalFlavorSingle(out_flv, energy, ini_state, zenith_angle, production_height, detector_depth);
            double result = prob_single/units.pi;
            return result;
        });
        ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR, 0, 0, 1e4);
        ig.SetFunction(wf);
        return ig.Integral(0.,units.pi);
    }
    
    double EvalFlavorIntegrated2(unsigned int out_flv, double energy, const std::array<double, 3>& ini_state = {0., 1., 0.}, double detector_depth = 0.)
    {
        ROOT::Math::Functor1D wf([&](double zenith_angle) {
            double prob_single = EvalFlavorIntegratedOverProductionHeight(out_flv, energy, ini_state, zenith_angle, detector_depth);
            double result = prob_single;
            return result;
        });
        ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
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
    
    nuSQUIDS nus(3, neutrino);
    // mixing angles
    nus.Set_MixingAngle(0, 1, 0.583995872);
    nus.Set_MixingAngle(0, 2, 0.148532030);
    nus.Set_MixingAngle(1, 2, 0.799399993);
    // square mass differences
    nus.Set_SquareMassDifference(1, 7.53e-5);
    nus.Set_SquareMassDifference(2, 2.42e-3);
    // CP phase
    nus.Set_CPPhase(0, 2, 0.0);

    //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
    nus.Set_h_max( 100.*units.km );

    //We set the GSL step function
    //nus.Set_GSL_step(gsl_odeiv2_step_rk4);

    //Setting the numerical precision of gsl integrator.
    nus.Set_rel_error(1.0e-10);
    nus.Set_abs_error(1.0e-10);
    
    
    //Declaration of the body, EarthAtm is one of the predefined bodies
    auto earth_atm = std::make_shared<MyEarthAtm>();
    //Definition of the track, in encodes the trajectory inside the body, here is declared with the zenith angle.
    auto track_atm = std::make_shared<MyEarthAtm::Track>(units.pi*0.5, 22.*units.km, 0.);
    //We set this in the nusSQuID object.
    nus.Set_Body(earth_atm);
    nus.Set_Track(track_atm);
    
    //Set true the progress bar during the evolution.
    //nus.Set_ProgressBar(true);
    
    nus.Set_E(1000.*units.MeV);
    
    //Construct the initial state
    marray<double, 1> inistate({3}, {1000., 2000., 0.});
    
    //Set the initial state in nuSQuIDS object
    nus.Set_initial_state(inistate, flavor);
    
    //Propagate the neutrinos in the earth for the path defined in path
    nus.EvolveState();
    
    std::cout << "Out state" << std::endl;
    for (double EE : nus.GetERange()){
        std::cout << EE/units.MeV << " MeV: ";
        for(int i = 0; i < 3; i++){
            std::cout << nus.EvalFlavor(i) << " ";
        }
        std::cout << std::endl;
    }
    
    MyEarthAtmIntegrator meai;
    double EE = 1000.*units.MeV;
    std::cout << "Out state from integrator: single" << std::endl;
    std::cout << EE/units.MeV << " MeV: ";
    for(unsigned int i = 0; i < 3; i++){
        std::cout << meai.EvalFlavorSingle(i, EE, {1000., 2000., 0.}, units.pi*0.5, 22.*units.km) << " ";
    }
    std::cout << std::endl;
    
    std::cout << "Out state from integrator: integrated" << std::endl;
    std::cout << EE/units.MeV << " MeV: ";
    for(unsigned int i = 0; i < 3; i++){
        //std::cout << meai.EvalFlavorIntegratedOverProductionHeight(i, EE, {0., 1., 0.}, units.pi*0.5) << " " << std::flush;
        //std::cout << meai.EvalFlavorIntegratedOverZenithAngle(i, EE, {0., 1., 0.}, 22.*units.km) << " " << std::flush;
        std::cout << meai.EvalFlavorIntegrated2(i, EE, {1000., 2000., 0.}) << " " << std::flush;
    }
    std::cout << std::endl;

    return 0;
}
