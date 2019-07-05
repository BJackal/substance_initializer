// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#ifndef INTEGRATION_SUBSTANCE_INITIALIZATION_H_
#define INTEGRATION_SUBSTANCE_INITIALIZATION_H_

#include <vector>

#include "biodynamo.h"
#include "core/substance_initializers.h"
//#include "diffusion_biology_modules.h"


// setup simulation parameters
// set numer of simulation steps
const int simulation_steps = 500; // 130000 Time between two simulation steps equals: 0.01hours (default)

// size of 2d layer for precursor cells
const double x_range = 150, y_range = 150, z_range = 4500; // set dims of simulation space
const double z_pos_precursor = 10; // (x/y range + z postion)
const double simulation_cube_dim = std::max(x_range, y_range);

// number of precursor cells
const size_t num_precursor_cells = 15;  // number of precursor cells (S1) in the simulation
const int default_cell_diameter = 6;


namespace bdm {

/// An initializer that creates a linear concentration along one axis given two points of the line
/// i.e. (startpos_, startvalue_) and (endpos_, endvalue_)
/// The linear concentration results in a constant gradient along the axis
///
struct LinearConcentration {
    double startvalue_;
    double endvalue_;
    double startpos_;
    double endpos_;
    double slope_;
    double intercept_;
    uint8_t axis_;

    /// @brief      The constructor ... computes slope and intercept of the linear function
    ///
    /// @param[in]  startvalue value at start postion
    /// @param[in]  endvalue   value at end postion
    /// @param[in]  startpos   start postion on axis_
    /// @param[in]  endpos     end postion on axis_
    /// @param[in]  axis       The axis along which you want linear concentration to be oriented to
    ///
    LinearConcentration(double startvalue, double endvalue, double startpos, double endpos, uint8_t axis) {
        startvalue_ = startvalue;
        endvalue_ = endvalue;
        startpos_ = startpos;
        endpos_ = endpos;
        axis_ = axis;

        // compute slope of the linear function
        slope_ = (endvalue_ - startvalue_) / (endpos_ - startpos_);
        // and its intercept
        intercept_ = startvalue_ - (slope_ * startpos_);
    }

    /// @brief      The model that we want to apply for substance initialization.
    ///             The operator is called for the entire space
    ///
    /// @param[in]  x     The x coordinate
    /// @param[in]  y     The y coordinate
    /// @param[in]  z     The z coordinate
    ///
    double operator()(double x, double y, double z) {
        switch(axis_) {
            case Axis::kXAxis: return (slope_ * x) + intercept_;
            case Axis::kYAxis: return (slope_ * y) + intercept_;
            case Axis::kZAxis: return (slope_ * z) + intercept_;
            default: throw std::logic_error("You have chosen an non-existing axis!");
        }
    }
};

// -----------------------------------------------------------------------------
// In this integration test we should how to make use of the 'substance
// initializers', in order to initialize the concentration of a particular
// substance. We create a gaussian distribution along each axis.
// -----------------------------------------------------------------------------

// Create list of substances
enum Substances { kSubstance };

// Define displacement behavior:

const double concentrationThreshold = 50;  // stop ovement a this concentraion

// Cells move along the diffusion gradient (from low concentration to high)
struct MoveAlongGradient : public BaseBiologyModule {
        BDM_STATELESS_BM_HEADER(MoveAlongGradient, BaseBiologyModule, 1);

    public:
        MoveAlongGradient() : BaseBiologyModule(gAllEventIds) {}

        void Run(SimObject* so) override {
            auto* sim = Simulation::GetActive();
            auto* rm = sim->GetResourceManager();
            static auto* kDg = rm->GetDiffusionGrid(kSubstance);

            kDg->SetConcentrationThreshold(1e15);

            if (auto* cell = dynamic_cast<Cell*>(so)) {

                const auto& current_position = so->GetPosition();  // get current cell postion

                double current_concentration;  // get concentraion at current cell position
                current_concentration = kDg->GetConcentration(current_position);

                if (current_concentration < concentrationThreshold) {
                    Double3 gradient;
                    kDg->GetGradient(current_position, &gradient);

                    gradient[0] *= 0.5;
                    gradient[1] *= 0.5;
                    gradient[2] *= 0.5;

                    cell->UpdatePosition(gradient);
                } else {
                    std::cout <<  "Concentraion " <<  current_concentration << " to high. Stopped migration!\n";
                }
            }
        }
    };



inline int Simulate(int argc, const char** argv) {

    // set space parameters of the simulation
    auto set_param = [](auto* param) {
        param->bound_space_ = true;
        param->min_bound_ = -(simulation_cube_dim/2);
        param->max_bound_ = (simulation_cube_dim/2);  // cube of 4500*4500*4500
        // param->run_mechanical_interactions_ = true; tru by default
    };

    Simulation simulation(argc, argv, set_param);
    auto* rm = simulation.GetResourceManager();  // get pointer to resource manager
    auto* random = simulation.GetRandom();  // get thread of local random number generator.

    double x_coord, y_coord, z_coord;

    // 2D plate for precursor cells (150x150)
    double x_min = 0 - (x_range/2);  // set position of the plate with (0,0) at the center of the simulation space
    double x_max = 0 + (x_range/2);
    double y_min = 0 - (y_range/2);
    double y_max = 0 + (y_range/2);

    // create 2d Layer of cells
    for (size_t i = 0; i < num_precursor_cells; ++i) {
        // create coordinates for cells in 2D plate
        x_coord = random->Uniform(x_min, x_max);
        y_coord = random->Uniform(y_min, y_max);
        z_coord = z_pos_precursor;

        // creating the cell at position x, y, z
        Cell* cell = new Cell({x_coord, y_coord, z_coord});
        // set cell parameters
        cell->SetDiameter(default_cell_diameter);
        cell->AddBiologyModule(new MoveAlongGradient());
        rm->push_back(cell);  // put the created cell in our cells structure
    }

    // 3. Define the substances in our simulation
    // Order: substance id, substance_name, diffusion_coefficient, decay_constant,
    // resolution
    ModelInitializer::DefineSubstance(kSubstance, "Substance", 0, 0, 20);

    // Init substance with linear concentration distribution
    //  LinearGradiend(double startvalue, double endvalue, double startpos, double endpos, uint8_t axis)
    ModelInitializer::InitializeSubstance(kSubstance, "Substance",
                                          LinearConcentration(0, 100, 0, simulation_cube_dim/2, Axis::kZAxis));


    // 4. Run simulation for N timesteps
    simulation.GetScheduler()->Simulate(simulation_steps);

    std::cout << "Simulation completed successfully!\n";
    return 0;
    }

}  // namespace bdm

#endif  // INTEGRATION_SUBSTANCE_INITIALIZATION_H_