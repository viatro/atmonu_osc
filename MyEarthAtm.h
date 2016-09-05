#ifndef MyEarthAtm_H
#define MyEarthAtm_H

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>

#include <vector>
#include <map>



namespace nusquids{
    
    /// \class MyEarthAtm
    /// \brief A model of the Earth with atmospheric neutrinos geometry.
    class MyEarthAtm: public Body{
      protected:
        /// \brief Radius of the Earth.
        static double earth_radius;
      public:
        static double GetEarthRadius() { return earth_radius; }
        static void SetEarthRadius(double r) { earth_radius = r; }
    
      protected:
        /// \brief Earth density parameters map
        std::map<double, std::array<double, 3> > earth_density_map;
        /// \brief Atmosphere temperature parameters map
        std::map<double, std::array<double, 3> > atmo_temperature_map;
        /// \brief Earth electron fraction map
        std::map<double, double> ye_map; // !!! Use map::lower_bound !!!
      
      public:
        /// \brief Default constructor using supplied PREM.
        MyEarthAtm();
        
        ~MyEarthAtm();
    
        /// \brief Serialization function
        void Serialize(hid_t group) const;
        /// \brief Deserialization function
        static std::shared_ptr<MyEarthAtm> Deserialize(hid_t group);
    
        /// \class Track
        /// \brief MyEarthAtm trajectory
        class Track: public Body::Track {
          friend class MyEarthAtm;
          private:
            /// \brief Zenith angle.
            double zenith_angle;
            /// \brief Neutrino production heihgt.
            double production_height;
            /// \brief Detector depth.
            double detector_depth;
            /// \brief Baseline.
            double L;
          public :
            /// \brief Construct trajectory.
            /// @param zenith_angle Zenith angle in radians.
            Track(double zenith_angle = 0., double production_height = 0., double detector_depth = 0.);
            /// 
            void SetParams(double zenith_angle, double production_height, double detector_depth);
            /// \brief Returns the neutrino baseline in natural units.
            double GetBaseline() const {return L;}
            /// \brief Returns current distance from center of Earth
            double GetCurrentRadius() const;
            virtual void FillDerivedParams(std::vector<double>& TrackParams) const;
        };
    
        /// \brief Returns the density in g/cm^3
        double density(const GenericTrack&) const;
        /// \brief Returns the electron fraction
        double ye(const GenericTrack&) const;
        /// \brief Returns the radius of the Earth in natural units.
    };
    
}

#endif
