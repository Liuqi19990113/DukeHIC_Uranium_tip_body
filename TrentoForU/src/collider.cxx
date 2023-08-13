// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "collider.h"

#include <cmath>
#include <string>
#include <vector>

#include <boost/program_options/variables_map.hpp>

#include "fwd_decl.h"
#include "nucleus.h"

#include <iostream>


namespace trento {

namespace {

// Helper functions for Collider ctor.

// Create one nucleus from the configuration.
NucleusPtr create_nucleus(const VarMap& var_map, std::size_t index) {
  const auto& species = var_map["projectile"]
                        .as<std::vector<std::string>>().at(index);
  const auto& nucleon_dmin = var_map["nucleon-min-dist"].as<double>();
  return Nucleus::create(species, nucleon_dmin);
}

// Determine the maximum impact parameter.  If the configuration contains a
// non-negative value for bmax, use it; otherwise, fall back to the minimum-bias
// default.
double determine_bmax(const VarMap& var_map,
    const Nucleus& A, const Nucleus& B, const NucleonCommon& nc) {
  auto bmax = var_map["b-max"].as<double>();
  if (bmax < 0.)
    bmax = A.radius() + B.radius() + nc.max_impact();
  return bmax;
}

// Determine the asymmetry parameter (Collider::asymmetry_) for a pair of
// nuclei.  It's just rA/(rA+rB), falling back to 1/2 if both radii are zero
// (i.e. for proton-proton).
double determine_asym(const Nucleus& A, const Nucleus& B) {
  double rA = A.radius();
  double rB = B.radius();
  double sum = rA + rB;
  if (sum < 0.1)
    return 0.5;
  else
    return rA/sum;
}

}  // unnamed namespace

// Lots of members to initialize...
// Several helper functions are defined above.
Collider::Collider(const VarMap& var_map)
    : nucleusA_(create_nucleus(var_map, 0)),
      nucleusB_(create_nucleus(var_map, 1)),
      nucleon_common_(var_map),
      nevents_(var_map["number-events"].as<int>()),
      calc_ncoll_(var_map["ncoll"].as<bool>()),
      bmin_(var_map["b-min"].as<double>()),
      bmax_(determine_bmax(var_map, *nucleusA_, *nucleusB_, nucleon_common_)),
      asymmetry_(determine_asym(*nucleusA_, *nucleusB_)),
      event_(var_map),
      output_(var_map), 
      ///LiuQi change it
      projectileA_name_(var_map["projectile"].as<std::vector<std::string>>().at(0)),
      projectileB_name_(var_map["projectile"].as<std::vector<std::string>>().at(1))
      {
  // Constructor body begins here.
  // Set random seed if requested.
  auto seed = var_map["random-seed"].as<int64_t>();
  if (seed > 0)
    random::engine.seed(static_cast<random::Engine::result_type>(seed));
}

// See header for explanation.
Collider::~Collider() = default;

///LiuQi change it
/// LiuQi change it
void Collider::run_events() {
  // The main event loop.
  std::string special_pro = "U2LiuQi";
  if(projectileA_name_==special_pro and projectileB_name_==special_pro)
  {
    for (int n = 0; n < nevents_; ++n) {
      // Sampling the impact parameter also implicitly prepares the nuclei for
      // event computation, i.e. by sampling nucleon positions and participants.
      std::tuple<double, int, double, double, double, double> collision_attr = sample_collision("U2LiuQi");
      double b = std::get<0>(collision_attr);
      int ncoll = std::get<1>(collision_attr);
      double spin_a = std::get<2>(collision_attr);
      double tilt_a = std::get<3>(collision_attr);
      double spin_b = std::get<4>(collision_attr);
      double tilt_b = std::get<5>(collision_attr);
      // Pass the prepared nuclei to the Event.  It computes the entropy profile
      // (thickness grid) and other event observables.
      event_.compute(*nucleusA_, *nucleusB_, nucleon_common_);

      // Write event data.
      output_(n, b, ncoll, event_, spin_a, tilt_a, spin_b, tilt_b);
    }
  }
  else
  {
    for (int n = 0; n < nevents_; ++n) {
      // Sampling the impact parameter also implicitly prepares the nuclei for
      // event computation, i.e. by sampling nucleon positions and participants.
      std::tuple<double, int> collision_attr = sample_collision();
      double b = std::get<0>(collision_attr);
      int ncoll = std::get<1>(collision_attr);
      double spin_a = -1;
      double tilt_a = -1;
      double spin_b = -1;
      double tilt_b = -1;
      // Pass the prepared nuclei to the Event.  It computes the entropy profile
      // (thickness grid) and other event observables.
      event_.compute(*nucleusA_, *nucleusB_, nucleon_common_);
      // Write event data.
      output_(n, b, ncoll, event_, spin_a, tilt_a, spin_b, tilt_b);
    }
  }
}

//LiuQi reload this function to special collision of U+U
std::tuple<double, int, double, double, double, double> Collider::sample_collision(std::string name) {
  // Sample impact parameters until at least one nucleon-nucleon pair
  // participates.  The bool 'collision' keeps track -- it is effectively a
  // logical OR over all possible participant pairs.
  // Returns the sampled impact parameter b, and binary collision number ncoll.
  double b;
  int ncoll = 0;
  bool collision = false;
  double spin_a = 0;
  double tilt_a = 0;
  double spin_b = 0;
  double tilt_b = 0;

  do {
    // Sample b from P(b)db = 2*pi*b.
    b = std::sqrt(bmin_ * bmin_ + (bmax_ * bmax_ - bmin_ * bmin_) * random::canonical<double>());

    // Offset each nucleus depending on the asymmetry parameter (see header).
    nucleusA_->sample_nucleons(asymmetry_ * b);
    nucleusB_->sample_nucleons((asymmetry_ - 1.) * b);
    
    // Check each nucleon-nucleon pair.
    for (auto&& A : *nucleusA_) {
      for (auto&& B : *nucleusB_) {
        auto new_collision = nucleon_common_.participate(A, B);
        if (new_collision && calc_ncoll_) ++ncoll;
        collision = new_collision || collision;
      }
    }
  } while (!collision);
  auto nucleus_A = dynamic_cast<DeformedNucleusSpecialConfig*>(nucleusA_.get());
  auto nucleus_B = dynamic_cast<DeformedNucleusSpecialConfig*>(nucleusB_.get());
  spin_a = nucleus_A->spin_angle_;
  tilt_a = nucleus_A->tilt_angle_;
  spin_b = nucleus_B->spin_angle_;
  tilt_b = nucleus_B->tilt_angle_;

  return std::make_tuple(b, ncoll, spin_a, tilt_a, spin_b, tilt_b);
}

std::tuple<double, int> Collider::sample_collision() {
  /// Sample impact parameters until at least one nucleon-nucleon pair
  /// participates.  The bool 'collision' keeps track -- it is effectively a
  /// logical OR over all possible participant pairs.
  /// Returns the sampled impact parameter b, and binary collision number ncoll.
    double b;
    int ncoll = 0;
    bool collision = false;
    do {
  // Sample b from P(b)db = 2*pi*b.
        b = std::sqrt(bmin_ * bmin_ + (bmax_ * bmax_ - bmin_ * bmin_) * random::canonical<double>());
  //
  //    Offset each nucleus depending on the asymmetry parameter (see header).
        nucleusA_->sample_nucleons(asymmetry_ * b);
        nucleusB_->sample_nucleons((asymmetry_ - 1.) * b);
  //                                      
  // Check each nucleon-nucleon pair.
        for (auto&& A : *nucleusA_) {
           for (auto&& B : *nucleusB_) {
              auto new_collision = nucleon_common_.participate(A, B);
              if (new_collision && calc_ncoll_) ++ncoll;
              collision = new_collision || collision;
                                               }
                                                 }
       } while (!collision);
     return std::make_tuple(b, ncoll);
                               }

}  // namespace trento
