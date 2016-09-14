#include "MyEarthAtm.h"

namespace nusquids{
    static squids::Const units;
    
    /*
    ----------------------------------------------------------------------
        MyEarthAtm CLASS DEFINITIONS
    ----------------------------------------------------------------------
    */
    
    double MyEarthAtm::earth_radius = 6371.*units.km;
    
    // track constructor
    //MyEarthAtm::Track::Track() {};
    
    MyEarthAtm::Track::Track(double zenith_angle, double production_height, double detector_depth)
    : Body::Track(0.,0.)
    , zenith_angle(zenith_angle)
    , production_height(production_height)
    , detector_depth(detector_depth)
    , L()
    {
        SetParams(zenith_angle, production_height, detector_depth);
    }
    
    void MyEarthAtm::Track::SetParams(double zenith_angle, double production_height, double detector_depth) {
        this->zenith_angle = zenith_angle;
        this->production_height = production_height;
        this->detector_depth = detector_depth;
        
        if (zenith_angle == 0.0) {
            L = production_height + detector_depth;
        } else if (zenith_angle == units.pi) {
            L = 2*MyEarthAtm::GetEarthRadius() + production_height - detector_depth;
        } else {
            double singamma = sin(units.pi - zenith_angle);
            double sinbeta = (MyEarthAtm::GetEarthRadius() - detector_depth) / (MyEarthAtm::GetEarthRadius() + production_height) * singamma;
            double sinalpha = sin(zenith_angle - asin(sinbeta)); // pi - (pi - zenith_angle) - beta
            L = (MyEarthAtm::GetEarthRadius() + production_height) * sinalpha / singamma;
        }
        
        x = xini = 0.0;
        xend = L;
    }
    
    double MyEarthAtm::Track::GetCurrentRadius() const {
        double c = MyEarthAtm::GetEarthRadius() + production_height;
        double b = MyEarthAtm::GetEarthRadius() - detector_depth;
        
        double cosbeta = (L*L + c*c - b*b) / (2.*L*c);
        double r = sqrt(x*x + c*c -2.*x*c*cosbeta);
        return r;
    }
    
    void MyEarthAtm::Track::FillDerivedParams(std::vector<double>& TrackParams) const{
        TrackParams.push_back(zenith_angle);
        TrackParams.push_back(production_height);
        TrackParams.push_back(detector_depth);
    }
    
    // MyEarthAtm constructor
    
    MyEarthAtm::MyEarthAtm()
    : Body(7,"MyEarthAtm")
    {
        ye_map[1221.5*units.km] = 0.4656;
        ye_map[3480.0*units.km] = 0.4656;
        ye_map[   earth_radius] = 0.4957;
        ye_map[2* earth_radius] = 0.5   ;
        
        // r in meters, a[0] + a[1]*r + a[2]*r^2 density in kg/m^3
        earth_density_map[1.2215e6] = {1.3088e4,  1.9110e-8, -2.1773e-10};
        earth_density_map[3.4800e6] = {1.2346e4,  1.3976e-4, -2.4123e-10};
        earth_density_map[3.6300e6] = {7.3067e3, -5.0007e-4,  0.        };
        earth_density_map[5.7010e6] = {6.7823e3, -2.4441e-4, -3.0922e-11};
        earth_density_map[5.7710e6] = {5.3197e3, -2.3286e-4,  0.        };
        earth_density_map[5.9710e6] = {1.1249e4, -1.2603e-3,  0.        };
        earth_density_map[6.1510e6] = {7.1083e3, -5.9706e-4,  0.        };
        earth_density_map[6.3466e6] = {2.6910e3,  1.0869e-4,  0.        };
        earth_density_map[6.3560e6] = {2.9000e3,  0.       ,  0.        };
        earth_density_map[6.3680e6] = {2.6000e3,  0.       ,  0.        };
        earth_density_map[6.3710e6] = {1.0200e3,  0.       ,  0.        };
        
        // geopotential h in meters, T in Kelvin, T lapse in K/m, density in g/cm^3
        atmo_temperature_map[    0.           ] = { 15.0  + 273.15, -0.0065, 1.2250e-3           };
        atmo_temperature_map[11000.           ] = {-56.5  + 273.15,  0.0   , 3.639271988833481e-4};
        atmo_temperature_map[20000.           ] = {-56.5  + 273.15,  0.001 , 8.80391839515545e-5 };
        atmo_temperature_map[32000.           ] = {-44.5  + 273.15,  0.0028, 1.3226067235598e-5  };
        atmo_temperature_map[47000.           ] = {- 2.5  + 273.15,  0.0   , 1.427313410448156e-6};
        atmo_temperature_map[51000.           ] = {- 2.5  + 273.15, -0.0028, 8.61479984576618e-7 };
        atmo_temperature_map[71000.           ] = {-58.5  + 273.15, -0.002 , 6.420472983251111e-8};
        atmo_temperature_map[84854.57642868205] = {-86.28 + 273.15,  0.0   , 6.954393506100748e-9};
        
        
    }
    
    MyEarthAtm::~MyEarthAtm() {}
    
    void MyEarthAtm::Serialize(hid_t group) const {}
    std::shared_ptr<EarthAtm> EarthAtm::Deserialize(hid_t group) { return nullptr; }
    
    double MyEarthAtm::density(const GenericTrack& track_input) const {
        const MyEarthAtm::Track& track_MyEarthAtm = static_cast<const MyEarthAtm::Track&>(track_input);
        
        double r = track_MyEarthAtm.GetCurrentRadius();
        double r_in_m = r/units.meter;

        if (r <= earth_radius) {
            std::array<double,3> pars = earth_density_map.lower_bound(r_in_m)->second;
            return (pars[0] + pars[1]*r_in_m + pars[2]*r_in_m*r_in_m)/1000.;
        } else if (r <= earth_radius + 86.*units.km) {
            double h = r_in_m - earth_radius/units.meter;
            h = h*earth_radius/(h + earth_radius); // convertion from geometric to geopotential
            
            auto iter = atmo_temperature_map.lower_bound(h);
            --iter;
            
            double g0 = 9.80665; // m/s^2
            double R = 8.3144598; // N·m/(mol·K)
            double M = 0.0289644; // kg/mol
            
            double hb   = iter->first;
            double Tb   = iter->second[0];
            double Lb   = iter->second[1];
            double rhob = iter->second[2];
            
            if (Lb == 0.0) {
                return rhob * exp(-g0*M*(h - hb)/(R*Tb));
            } else {
                return rhob * pow(Tb/(Tb + Lb*(h - hb)), 1 + g0*M/(R*Lb));
            }
        } else return 0;
    }

    double MyEarthAtm::ye(const GenericTrack& track_input) const {
        const MyEarthAtm::Track& track_MyEarthAtm = static_cast<const MyEarthAtm::Track&>(track_input);
        double r = track_MyEarthAtm.GetCurrentRadius();
        return ye_map.lower_bound(r)->second;
    }
    
}
