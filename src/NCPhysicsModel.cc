
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCPhysicsModel.hh"

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include <iostream>
#include <iomanip>
#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/sabscatter/NCSABScatter.hh"
#include "NCrystal/internal/sab/NCSABFactory.hh"
#include "NCrystal/internal/sab/NCSABExtender.hh"
#include "NCrystal/internal/sab/NCSABIntegrator.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/freegas/NCFreeGas.hh"

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/core/NCDefs.hh"

#include <cmath>
#include <fstream>

using namespace std;

// Forward declarations
double egelstaffDiffusion(double alpha, double beta, double c, double wt_diffusion);
struct ConvolutionParams;
NC::SABData convoluteSABData(const NC::SABData& sabdata, const ConvolutionParams& params);
NC::SABData applyYoungKoppel(const NC::SABData& sabdata, const std::string& yk_model, const NC::AtomData& atomData, const NCP::LiquidInfo& liq);
std::vector<double> getColumn(const std::vector<std::vector<double>>& matrix, size_t j);
// Applies the Skold or Vineyard correction for a given alpha point.
//   Skold:    S(Q)·S_s(α',β) − S_s(α,β),  α' = α/S(Q)
//   Vineyard: (S(Q)−1)·S_s(α,β)
std::vector<double> coherentCorrection(const NC::SABData& s_s, const NCP::LiquidInfo& liq,
                                       int ia, double alpha, double Q);

// Output (alpha, beta) grids on which the final convolved S(alpha,beta) is stored.
// Note: the actual convolution integral is performed on the finer beta_prime grid inside processConvolutionIteration.
struct ConvolutionGrids {
    std::vector<double> alpha_out;
    std::vector<double> beta_out;
    std::vector<double> posbeta;  // positive beta values, mirrored to construct the full output beta grid
    double beta_max_phys;
    double alpha_max_phys;
};

// Constructs logarithmically-spaced output (alpha, beta) grids on which the convolved S(alpha,beta) is returned.
// grid_points = 400 by default; 
ConvolutionGrids createConvolutionGrids(const NC::VectD& betaGrid, size_t grid_points = 400) {
    ConvolutionGrids grids;
    
    // Calculate physical limits for grid construction
    // The maximum alpha value is computed as 4*beta, see discussion on pg. 674 of the NJOY manual
    grids.beta_max_phys = 0.0;
    for (const auto& beta : betaGrid) {
        grids.beta_max_phys = std::max(grids.beta_max_phys, std::abs(beta));
    }
    grids.alpha_max_phys = 4.0 * grids.beta_max_phys; // Maximum alpha
    
    // Construct the positive side of the beta grid: start from 1e-5 to beta_max_phys
    grids.posbeta.resize(grid_points);
    double log_start = std::log10(1e-5);
    double log_end = std::log10(grids.beta_max_phys);
    double log_step = (log_end - log_start) / (grid_points - 1);
    for (size_t i = 0; i < grid_points; ++i) {
        grids.posbeta[i] = std::pow(10.0, log_start + i * log_step);
    }
    
    // Construct output alpha grid: logarithmically spaced from 1e-4 to alpha_max_phys
    grids.alpha_out.resize(grid_points);
    log_start = std::log10(1e-4);
    log_end = std::log10(grids.alpha_max_phys);
    log_step = (log_end - log_start) / (grid_points - 1);
    for (size_t i = 0; i < grid_points; ++i) {
        grids.alpha_out[i] = std::pow(10.0, log_start + i * log_step);
    }

    // Construct beta_out by mirroring positive side to negative
    std::vector<double> negbeta(grid_points);
    for (size_t i = 0; i < grid_points; ++i) {
        negbeta[i] = -grids.posbeta[grid_points - 1 - i]; // Negative values, reversed order
    }
    
    // Combine: negbeta + 0 + posbeta
    grids.beta_out.reserve(2 * grid_points + 1);
    grids.beta_out.insert(grids.beta_out.end(), negbeta.begin(), negbeta.end());
    grids.beta_out.push_back(0.0);
    grids.beta_out.insert(grids.beta_out.end(), grids.posbeta.begin(), grids.posbeta.end());
    
    return grids;
}

// Pre-computed quantities for the Egelstaff convolution, cached per alpha_out point
struct ConvolutionPrecomputed {
    std::vector<double> aws_values;  // ws * alpha_out[i]: solid-scaled alpha for VDOS interpolation
    std::vector<int>    i1_values;   // lower VDOS alpha grid index for each alpha_out point
    std::vector<int>    i2_values;   // upper VDOS alpha grid index for each alpha_out point
    std::vector<double> weights;     // linear interpolation weight between i1 and i2
    std::vector<double> exp_terms;   // exp(-alpha * dwi): zeroth-phonon prefactor
    std::vector<std::vector<double>> Sab; // convolved S(alpha,beta) on the (alpha_out, beta_out) grid
};

// Struct to hold interpolation indices and weights
struct InterpolationData {
    int i1, i2;
    double weight;
    bool valid;
};

// Find interpolation indices and weight for a value in a grid
InterpolationData findInterpolationData(const std::vector<double>& grid, double value) {
    InterpolationData data;
    
    // Find insertion point (equivalent to searchsorted)
    auto it = std::lower_bound(grid.begin(), grid.end(), value);
    data.i2 = std::distance(grid.begin(), it);
    data.i1 = data.i2 - 1; // Can be negative
    
    // Check if interpolation is valid
    data.valid = (data.i1 >= 0 && data.i2 < static_cast<int>(grid.size()));
    
    if (data.valid) {
        double grid_diff = grid[data.i2] - grid[data.i1];
        if (grid_diff > 0) {
            data.weight = (value - grid[data.i1]) / grid_diff;
        } else {
            data.weight = 0.0;
        }
    } else {
        data.weight = 0.0;
    }
    
    return data;
}

// Helper function for linear interpolation using pre-computed indices
double interpolateWithIndices(const std::vector<double>& values, const InterpolationData& data) {
    if (!data.valid) {
        // Use boundary values
        if (data.i1 < 0) return values.front();
        if (data.i2 >= static_cast<int>(values.size())) return values.back();
    }
    
    return (1.0 - data.weight) * values[data.i1] + data.weight * values[data.i2];
}

// Struct to hold convolution physics parameters
struct ConvolutionParams {
    double dwi;              // used for the Debye-Waller factor, exp(-alpha*dwi)
    double delta_beta_vdos; // beta grid spacing from VDOS, sets integration step in convolution
    double awr;              // atomic weight ratio
    NCP::LiquidInfo liquid_info; // coherent model parameters and S(Q) data
};

// For each point on the convolution alpha grid, finds where ws*alpha falls on the
// original VDOS alpha grid (storing the two bracketing indices and interpolation weight)
// and computes the Debye-Waller prefactor exp(-alpha*dwi).
ConvolutionPrecomputed precomputeConvolutionData(const ConvolutionGrids& grids, 
                                                  const NC::VectD& alphaGrid,
                                                  const ConvolutionParams& params) {
    ConvolutionPrecomputed precomp;
    size_t n_alpha_out = grids.alpha_out.size();
    size_t n_beta_out = grids.beta_out.size();
    
    // Pre-allocate output array
    precomp.Sab.resize(n_alpha_out, std::vector<double>(n_beta_out, 0.0));
    
    // Pre-compute aws values and interpolation indices
    precomp.aws_values.resize(n_alpha_out);
    precomp.i2_values.resize(n_alpha_out);
    precomp.i1_values.resize(n_alpha_out);
    
    precomp.weights.resize(n_alpha_out, 0.0);
    for (size_t i = 0; i < n_alpha_out; ++i) {
        precomp.aws_values[i] = params.liquid_info.ws * grids.alpha_out[i];
        
        // Use unified interpolation logic
        InterpolationData interpData = findInterpolationData(alphaGrid, precomp.aws_values[i]);
        precomp.i1_values[i] = interpData.i1;
        precomp.i2_values[i] = interpData.i2;
        precomp.weights[i] = interpData.weight;
    }
    
    // Pre-compute the zeroth-phonon term
    precomp.exp_terms.resize(n_alpha_out);
    for (size_t i = 0; i < n_alpha_out; ++i) {
        precomp.exp_terms[i] = std::exp(-grids.alpha_out[i] * params.dwi);
    }
    
    return precomp;
}

// Helper function for 1D convolution using trapezoidal integration
double convolve1d(const std::vector<double>& y, const std::vector<double>& f2_values, double delta_beta) {
    if (y.size() != f2_values.size()) return 0.0;
    
    // Direct trapezoidal integration without temporary array allocation
    double result = 0.0;
    for (size_t k = 0; k < y.size() - 1; ++k) {
        double integrand_left = y[k] * f2_values[k];
        double integrand_right = y[k + 1] * f2_values[k + 1];
        result += (integrand_left + integrand_right) * delta_beta * 0.5;
    }
    return result;
}

// Linear interpolation
double interpolate_solid(const std::vector<double>& grid, const std::vector<double>& values, double x) {
    if (x <= grid.front()) return values.front();
    if (x >= grid.back()) return values.back();
    
    InterpolationData data = findInterpolationData(grid, x);
    return interpolateWithIndices(values, data);
}

// Struct to hold parameters for single convolution iteration
struct ConvolutionIterationParams {
    size_t alpha_idx;
    double alpha_val;
    const std::vector<double>& alpha_out;
    const std::vector<double>& beta_out;
    const NC::VectD& alphaGrid;
    const NC::VectD& betaGrid; 
    const NC::VectD& sabValues;
    const ConvolutionPrecomputed& precomp;
    const ConvolutionParams& convParams;
    // Working arrays (passed by reference to avoid reallocation)
    std::vector<double>& f1_beta_out;
    std::vector<double>& beta_prime;
    std::vector<double>& f1_beta_prime;
};

// Set to true to print convolution progress every 30 alpha steps
static constexpr bool verbose_convolution = false;

// Print progress information for convolution iteration
void printConvolutionProgress(size_t iteration_idx, double alpha_val,
                              double beta_min, double beta_max,
                              size_t n_beta_points, double f1_integral,
                              const std::vector<double>& beta_out,
                              const std::vector<std::vector<double>>& Sab) {
    if (!verbose_convolution) return;
    if (iteration_idx % 30 != 0) return; // Only print every 30 iterations

    size_t n_beta_out = beta_out.size();
    size_t i = iteration_idx;

    // Calculate integral of Sab[i] using trapezoidal rule
    double sab_integral = 0.0;
    for (size_t j = 0; j < n_beta_out - 1; ++j) {
        double dx = beta_out[j + 1] - beta_out[j];
        double integrand = (Sab[i][j] + Sab[i][j + 1]) * dx * 0.5;
        if (std::isfinite(integrand)) {
            sab_integral += integrand;
        }
    }

    std::cout << i << ", α=" << std::fixed << std::setprecision(3) << alpha_val
              << ", β_range=[" << std::setprecision(1) << beta_min << ", " << beta_max << "]"
              << ", n_pts=" << n_beta_points
              << ", ∫f1=" << std::setprecision(6) << f1_integral
              << ", ∫Sab=" << std::setprecision(6) << sab_integral << std::endl;
}

// Computes the Egelstaff convolution S(alpha,beta) = f1(beta)*exp(-alpha*dwi)
// + integral[ f1(beta') * S_solid(alpha, beta-beta') dbeta' ] for a single alpha point.
void processConvolutionIteration(const ConvolutionIterationParams& params) {
    size_t i = params.alpha_idx;
    double a = params.alpha_val;
    const auto& alpha_out = params.alpha_out;
    const auto& beta_out = params.beta_out;
    const auto& alphaGrid = params.alphaGrid;
    const auto& betaGrid = params.betaGrid;
    const auto& sabValues = params.sabValues;
    const auto& precomp = params.precomp;
    auto& Sab = const_cast<std::vector<std::vector<double>>&>(precomp.Sab);
    auto& f1_beta_out = params.f1_beta_out;
    auto& beta_prime = params.beta_prime;
    auto& f1_beta_prime = params.f1_beta_prime;
    
    size_t n_beta_out = beta_out.size();
    
    // Evaluate diffusion on the existing beta_out grid (evaluated at alpha/awr for the NCrystal convention)
    for (size_t j = 0; j < n_beta_out; ++j) {
        double diff_val = egelstaffDiffusion(a/params.convParams.awr, beta_out[j], params.convParams.liquid_info.c, params.convParams.liquid_info.wt);
        if (!std::isfinite(diff_val)) {
            diff_val = 0.0;
        }
        f1_beta_out[j] = diff_val;
    }
    
    // Find where it drops below threshold
    double threshold = 1e-3;
    std::vector<size_t> indices_above;
    for (size_t j = 0; j < n_beta_out; ++j) {
        if (f1_beta_out[j] > threshold) {
            indices_above.push_back(j);
        }
    }
    
    double beta_min, beta_max;
    if (!indices_above.empty()) {
        beta_min = beta_out[indices_above.front()];
        beta_max = beta_out[indices_above.back()];
    } else {
        beta_min = beta_out[0];
        beta_max = beta_out.back();
    }
    
    // The following sections build a fine uniform integration grid (beta_prime) over the range 
    // where f1 is significant, which is used for the convolution integral.
    double delta_beta_vdos = std::min(params.convParams.delta_beta_vdos, a);
    
    // Safety check: prevent huge arrays by setting minimum delta_beta_vdos
    double beta_range = beta_max - beta_min;
    double min_delta = beta_range / 1000000.0; // Max 1M points
    delta_beta_vdos = std::max(delta_beta_vdos, min_delta);
    
    beta_prime.clear();
    for (double b = beta_min; b <= beta_max + delta_beta_vdos; b += delta_beta_vdos) {
        beta_prime.push_back(b);
        // Safety break to prevent infinite loops
        if (beta_prime.size() > 5000000) {
            std::cout << "Warning: beta_prime array too large, truncating at " << beta_prime.size() << " points" << std::endl;
            break;
        }
    }
    
    // Evaluate Egelstaff diffusion on beta_prime (evaluated at alpha/awr for the NCrystal convention)
    f1_beta_prime.resize(beta_prime.size());
    for (size_t j = 0; j < beta_prime.size(); ++j) {
        double diff_val = egelstaffDiffusion(a/params.convParams.awr, beta_prime[j], params.convParams.liquid_info.c, params.convParams.liquid_info.wt);
        // Check for NaN/inf values
        if (!std::isfinite(diff_val)) {
            std::cout << "Warning: NaN/inf in egelstaff at α=" << a << ", β=" << beta_prime[j] << std::endl;
            diff_val = 0.0;
        }
        f1_beta_prime[j] = diff_val;
    }
    
    // Check if the f1 integral is negligible; if so, skip the convolution.
    double f1_integral = 0.0;
    for (size_t j = 0; j < beta_prime.size() - 1; ++j) {
        double dx = beta_prime[j + 1] - beta_prime[j];
        double integrand = (f1_beta_prime[j] + f1_beta_prime[j + 1]) * dx * 0.5;
        if (std::isfinite(integrand)) {
            f1_integral += integrand;
        }
    }
    
    if (f1_integral < 1e-3 ) {
        // Skip convolution, just use direct term
        for (size_t j = 0; j < n_beta_out; ++j) {
            Sab[i][j] = f1_beta_out[j] * precomp.exp_terms[i];
        }
    } else {
        // Interpolate solid part
        int i1 = precomp.i1_values[i];
        int i2 = precomp.i2_values[i];
        std::vector<double> ss(betaGrid.size());
        
        if (i1 < 0) {
            // i1 < 0: use first column
            for (size_t j = 0; j < betaGrid.size(); ++j) {
                ss[j] = sabValues[j * alphaGrid.size() + 0];
            }
        } else if (i2 >= static_cast<int>(alphaGrid.size())) {
            // i2 >= len(alpha): use last column
            for (size_t j = 0; j < betaGrid.size(); ++j) {
                ss[j] = sabValues[j * alphaGrid.size() + (alphaGrid.size() - 1)];
            }
        } else {
            // Interpolate
            for (size_t j = 0; j < betaGrid.size(); ++j) {
                double s1 = sabValues[j * alphaGrid.size() + i1];
                double s2 = sabValues[j * alphaGrid.size() + i2];
                ss[j] = s1 + precomp.weights[i] * (s2 - s1);
            }
        }
        
        // Vectorized computation of convolution
        std::vector<double> f2_values(beta_prime.size());
        for (size_t j = 0; j < n_beta_out; ++j) {
            double b = beta_out[j];
            
            // Reuse f2_values vector for this beta
            for (size_t k = 0; k < beta_prime.size(); ++k) {
                double beta_diff = b - beta_prime[k];
                f2_values[k] = interpolate_solid(betaGrid, ss, beta_diff);
            }
            
            double conv_result = convolve1d(f1_beta_prime, f2_values, delta_beta_vdos);
            Sab[i][j] = f1_beta_out[j] * precomp.exp_terms[i] + conv_result;
        }
    }
    
    // Progress printing
    printConvolutionProgress(i, a, beta_min, beta_max, beta_prime.size(), 
                             f1_integral, beta_out, Sab);
}

bool NCP::PhysicsModel::isApplicable( const NC::Info& info )
{
  //Accept if input is NCMAT data with @CUSTOM_<pluginname> section:
  return info.countCustomSections(pluginNameUpperCase()) > 0;
}

NCP::PhysicsModel NCP::PhysicsModel::createFromInfo( const NC::Info& info )
{
  //Parse the content of our custom section. In case of syntax errors, we should
  //raise BadInput exceptions, to make sure users gets understandable error
  //messages. We should try to avoid other types of exceptions.

  //Get the relevant custom section data (and verify that there are not multiple
  //such sections in the input data):
  if ( info.countCustomSections( pluginNameUpperCase() ) != 1 )
    NCRYSTAL_THROW2(BadInput,"Multiple @CUSTOM_"<<pluginNameUpperCase()
                    <<" sections are not allowed");
  auto data = info.getCustomSection( pluginNameUpperCase() );

  // data is here a vector of lines, and each line is a vector of words. In our
  // case, we want to accept sections of the form (units are barn and angstrom
  // as is usual in NCrystal):
  //
  // @CUSTOM_<ourpluginname>
  //    <c_diff> <wt_diff> <wt_solid>
  //

  //Parse and validate values:
  std::unordered_map<std::string, LiquidInfo> liquid_data;
  std::string yk_model = "";

   for (int i=0;i!=data.size();i++)
    {
      auto line = data.at(i);
      if (line.empty()) continue;

      // Handle standalone yk_model keyword line (INCOHERENT + YK case)
      if (line.at(0) == "yk_model") {
        if (line.size() < 2)
          NCRYSTAL_THROW2(BadInput,"yk_model keyword requires a value (ortho or para) in @CUSTOM_"
                          <<pluginNameUpperCase()<<" section");
        std::string ykm = line.at(1);
        std::transform(ykm.begin(), ykm.end(), ykm.begin(), ::toupper);
        if (ykm != "ORTHO" && ykm != "PARA")
          NCRYSTAL_THROW2(BadInput,"Invalid yk_model value '"<<line.at(1)<<"' in @CUSTOM_"
                          <<pluginNameUpperCase()<<" section: must be 'ortho' or 'para'");
        yk_model = ykm;
        continue;
      }

      auto lbl = line.at(0);
      LiquidInfo elem;
      std::string model = line.at(1);
      std::transform(model.begin(), model.end(), model.begin(), ::toupper);
      if (model != "INCOHERENT" && model != "VINEYARD" && model != "SKOLD")
          NCRYSTAL_THROW2( BadInput,"Invalid coherent model in @CUSTOM_"
                           <<pluginNameUpperCase()<<" section: "<< model );
      elem.coherent_model = model;
      if ( ! NC::safe_str2dbl( line.at(2), elem.c )
        || ! NC::safe_str2dbl( line.at(3), elem.wt )
        || ! NC::safe_str2dbl( line.at(4), elem.ws )
        || ! (elem.c>0.0) || !(elem.wt>=0.0) || !(elem.ws>=0.0) )
            NCRYSTAL_THROW2( BadInput,"Invalid values specified in @CUSTOM_"
                     <<pluginNameUpperCase()<<" section: c, wt and ws must be"
                                            " positive floating point values)" );
      // FIXME: Add checks for values
      if (model == "VINEYARD" || model == "SKOLD")
      {
        // Optional yk_model line before S(Q) table (SKOLD or VINEYARD + YK case)
        if (i+1 != (int)data.size() && data.at(i+1).size() >= 2
            && data.at(i+1).at(0) == "yk_model") {
          std::string ykm = data.at(i+1).at(1);
          std::transform(ykm.begin(), ykm.end(), ykm.begin(), ::toupper);
          if (ykm != "ORTHO" && ykm != "PARA")
            NCRYSTAL_THROW2(BadInput,"Invalid yk_model value in @CUSTOM_"
                            <<pluginNameUpperCase()<<" section: must be 'ortho' or 'para'");
          yk_model = ykm;
          i++;
        }
        double old_x = -1;
        while (i+1!=(int)data.size()){
           if (data.at(i+1).size() != 2) break;
           double temp_x, temp_y;
           if ( ! NC::safe_str2dbl( data.at(i+1).at(0), temp_x )
             || ! NC::safe_str2dbl( data.at(i+1).at(1), temp_y )
             || ! (temp_x>=0.0) || !(temp_y>=0.0) )
            NCRYSTAL_THROW2( BadInput,"Invalid values specified in @CUSTOM_"
                     <<pluginNameUpperCase()<<" section: q and S(Q) must be"
                                            " positive floating point values)" );
            if ( !(temp_x > old_x) )
            NCRYSTAL_THROW2( BadInput,"Invalid values specified in @CUSTOM_"
                     <<pluginNameUpperCase()<<" section: q must be incremental" );
           old_x = temp_x;
           elem.qval.push_back(temp_x);
           elem.Sval.push_back(temp_y);
           i++;
         }
      }
      liquid_data[lbl] = elem;
    }
  // Assign yk_model to all parsed elements
  if (!yk_model.empty())
    for (auto& kv : liquid_data)
      kv.second.yk_model = yk_model;

  //Parsing done! Create and return our model:
  return PhysicsModel(liquid_data,info);
}

NCP::PhysicsModel::PhysicsModel( std::unordered_map<std::string, LiquidInfo> liquid_data, const NC::Info& info )
  : m_liquid_data(liquid_data)
{
  NC::ProcImpl::ProcComposition::ComponentList components;

   for ( auto& di : info.getDynamicInfoList() )
   {
     auto di_vdos = dynamic_cast<const NC::DI_VDOS*>(di.get());
     auto lbl = info.displayLabel(di->atom().index);
     if (liquid_data.find(lbl) != liquid_data.end())
     {
       const double mass_neutron_amu = NC::const_neutron_mass_amu;  // neutron mass in AMU
       const double mass_neutron = (mass_neutron_amu * NC::constant_dalton2eVc2) /
                                   ((NC::constant_c * 1e-12) * (NC::constant_c * 1e-12));  // eV·ps²/Å²
       const double hbar = NC::constant_planck / NC::k2Pi * 1e12;  // eV·ps
       const double bk = NC::constant_boltzmann;  // eV/K
       const double temperature = info.getTemperature().dbl(); // K
       const double awr = (double) di->atomData().averageMassAMU().relativeToNeutronMass();
       const double sigma_coh = di->atomData().scatteringXS().dbl() - di->atomData().incoherentXS().dbl();
       liquid_data[lbl].cfrac = sigma_coh / di->atomData().scatteringXS().dbl();  // sigma_coh/(sigma_coh+sigma_inc)
       
       // Extract VDOS data
       const auto& vdosData = di_vdos->vdosData();
       double emax = vdosData.vdos_egrid().second;
       const auto& rho = vdosData.vdos_density();

       // Calculate delta_beta from VDOS data
       double delta_beta_vdos = (emax / rho.size()) / (bk * std::abs(temperature));  // dimensionless beta spacing

       // Calculate MSD from VDOS data
       NC::VDOSEval vdosEval(di_vdos->vdosData());
       double msd = vdosEval.getMSD();
       
       // Calculate dwi, used for the Debye-Waller factor exp(-alpha*dwi)
       double dwi = msd * 2.0 * mass_neutron / (hbar * hbar) * bk * temperature * liquid_data[lbl].ws;  // unitless
       
       NC::ScatKnlData skd = NC::createScatteringKernel(di_vdos->vdosData(), 3,  0, NC::VDOSGn::TruncAndThinningChoices::Default);
       NC::SABData sabdata = NC::SABUtils::transformKernelToStdFormat(std::move(skd));

       ConvolutionParams convParams = {dwi, delta_beta_vdos, awr, liquid_data[lbl]};

       // Convolve with the diffusion model
       NC::SABData s_s = convoluteSABData(sabdata, convParams);

       // Young-Koppel (+ coherent correction if Skold or Vineyard)
       if (!liquid_data[lbl].yk_model.empty()) {
         NC::SABData yk_data = applyYoungKoppel(s_s, liquid_data[lbl].yk_model, di->atomData(), liquid_data[lbl]);
         components.push_back({di->fraction(), NC::makeSO<NC::SABScatter>(std::move(yk_data))});

       // Coherent correction only (no Young-Koppel)
       } else if (liquid_data[lbl].coherent_model == "SKOLD" || liquid_data[lbl].coherent_model == "VINEYARD") {
         const auto& alpha_grid = s_s.alphaGrid();
         const auto& beta_grid  = s_s.betaGrid();
         const int   n_alpha    = (int)alpha_grid.size();
         const int   n_beta     = (int)beta_grid.size();
         NC::VectD sab_out(s_s.sab().begin(), s_s.sab().end());
         for (int ia = 0; ia < n_alpha; ++ia) {
           const double alpha   = alpha_grid[ia];
           const double Q       = std::sqrt(2.0*mass_neutron*s_s.temperature().kT()*alpha)/hbar;
           const auto   bracket = coherentCorrection(s_s, liquid_data[lbl], ia, alpha, Q);
           for (int ib = 0; ib < n_beta; ++ib) {
             sab_out[ib * n_alpha + ia] += liquid_data[lbl].cfrac * bracket[ib];
             sab_out[ib * n_alpha + ia]  = std::max(0.0, sab_out[ib * n_alpha + ia]);
           }
         }
         NC::SABData corrected(NC::VectD(alpha_grid), NC::VectD(beta_grid), std::move(sab_out),
                               s_s.temperature(), s_s.boundXS(), s_s.elementMassAMU(), s_s.suggestedEmax());
         components.push_back({di->fraction(), NC::makeSO<NC::SABScatter>(std::move(corrected))});

       // Translational kernel only
       } else {
         components.push_back({di->fraction(), NC::makeSO<NC::SABScatter>(std::move(s_s))});
       }

     }
     else {
       components.push_back({di->fraction(), NC::makeSO<NC::FreeGas>(di->temperature(), di->atomData())});
    }
   }

  NC::ProcImpl::ProcPtr procptr = NC::ProcImpl::ProcComposition::consumeAndCombine( std::move(components), NC::ProcessType::Scatter );

  auto rngproducer = NC::getDefaultRNGProducer();
  auto rng = rngproducer->produce();

  m_proc = std::make_shared<NC::Scatter>( std::move(rngproducer), std::move(rng), std::move(procptr));
}

double NCP::PhysicsModel::calcCrossSection( double neutron_ekin ) const
{
  return m_proc->crossSectionIsotropic(NC::NeutronEnergy{neutron_ekin}).dbl();
}

NCP::PhysicsModel::ScatEvent
NCP::PhysicsModel::sampleScatteringEvent( NC::RNG& rng,
                                          double neutron_ekin ) const
{
  ScatEvent result;
  auto res = m_proc->sampleScatterIsotropic(NC::NeutronEnergy{neutron_ekin});
  result.ekin_final = res.ekin.dbl();
  result.mu = res.mu.dbl();
  return result;
}

// Numerical evaluaion of the modified Bessel function I1
// Adapted from polynimial expansion in Abramowitz and Stegun, vol. 55, sec. 9.8
double bessel_i1( double x ){
    double t = x / 3.75;
    double t2 = t*t;
    double res;
    if (t < 1.0){
        res = (0.5
               + t2*(0.87890594
               + t2*(0.51498869
               + t2*(0.15084934
               + t2*(0.02658733
               + t2*(0.00301532
               +  t2*0.00032411))))));
        return res * x;
    } else {
        double u = 1.0/t;
        res = (     0.39894228
                - u*(0.03988024
                + u*(0.00362018
                - u*(0.00163801
                - u*(0.01031555
                - u*(0.02282967
                - u*(0.02895312
                - u*(0.01787654
                - u*0.00420059))))))));
        return res*std::exp(x)/std::sqrt(x);
    }
}

// Numerical evaluaion of the modified Bessel function K1
// Adapted from polynimial expansion in Abramowitz and Stegun, vol. 55, sec. 9.8
double bessel_k1( double x ){
    double t = 0.5*x;
    double t2 = t*t;
    double res;
    if (t < 1.0){
        res = ( x * std::log(t) * bessel_i1(x)
                + 1.0
                + t2*(0.15443144
                - t2*(0.67278579
                + t2*(0.18156897
                + t2*(0.01919402
                + t2*(0.00110404
                +  t2*0.00004686))))));
        return res/x;
    } else {
        double u = 1.0 / t;
        res = (      1.25331414
                + u*(0.23498619
                - u*(0.03655620
                - u*(0.01504268
                - u*(0.00780353
                - u*(0.00325614
                -  u*0.00068245))))));
        return res * std::exp(-x) / std::sqrt(x);
    }
}

double spherical_jn(int n, double x) {
    /*
    Implementation of the spherical Bessel function j_n(x).
    */
    if (n < 0 || x < 0.0)
        return 0.0;
    if (x == 0.0)
        return n == 0 ? 1.0 : 0.0;

    // Small-x power series from j_n(x)=sqrt(pi/(2x))*J_{n+1/2}(x):
    //   https://dlmf.nist.gov/10.47.E3
    // with the Bessel-J power series:
    //   https://dlmf.nist.gov/10.2.E2
    // Near x=0, evaluate the power series directly instead of using recurrence,
    // which is numerically unstable for higher orders.
    if (x < 0.1) {
        // Leading small-x term j_n(x) ~ x^n/(2n+1)!!:
        //   https://dlmf.nist.gov/10.52.E1
        double leading = 1.0;
        for (int k = 1; k <= n; ++k)
            leading *= x / double(2 * k + 1);

        // Higher-order correction terms from the J_{n+1/2} power series:
        //   https://dlmf.nist.gov/10.2.E2
        double term = 1.0;
        double series = 1.0;
        double y = -0.25 * x * x;
        for (int k = 1; k < 32; ++k) {
            term *= y / (double(k) * (double(n + k) + 0.5));
            series += term;
            if (std::abs(term) <= 1.0e-16 * std::max(1.0, std::abs(series)))
                break;
        }
        return leading * series;
    }

    // Exact base values j_0(x)=sin(x)/x and
    // j_1(x)=sin(x)/x^2-cos(x)/x.
    double j0 = std::sin(x) / x;
    if (n == 0)
        return j0;
    double j1 = std::sin(x) / (x * x) - std::cos(x) / x;
    if (n == 1)
        return j1;

    // Upward use of the three-term recurrence:
    //   j_{n+1}(x)=(2n+1)j_n(x)/x-j_{n-1}(x)
    //   https://dlmf.nist.gov/10.51.E1
    // This is stable when x is larger than n.
    if (x > n) {
        double jm1 = j0;
        double j = j1;
        for (int ell = 1; ell < n; ++ell) {
            double jp1 = (2 * ell + 1) * j / x - jm1;
            jm1 = j;
            j = jp1;
        }
        return j;
    }

    // Downward use of the same three-term recurrence:
    //   https://dlmf.nist.gov/10.51.E1
    // This is stable when x is less than or comparable to n. The result is
    // normalized to the exact j_0 value computed above.
    int m = std::max(n + 32, int(x) + 32);
    double jp1 = 0.0;
    double j = 1.0;
    double target = 0.0;
    for (int ell = m; ell > 0; --ell) {
        double jm1 = (2 * ell + 1) * j / x - jp1;
        if (ell - 1 == n)
            target = jm1;
        jp1 = j;
        j = jm1;
        if (std::abs(j) > 1.0e100) {
            jp1 *= 1.0e-100;
            j *= 1.0e-100;
            target *= 1.0e-100;
        }
    }
    return j0 * target / j;
}


// Helper: stable evaluation of exp(var) * K1(x)
static inline double stable_exp_times_K1( double var, double x ) {
  // Platform-derived safe exponent limits
  static const double EXP_MAX = std::log(std::numeric_limits<double>::max()) - 5.0; // safety margin
  static const double EXP_MIN = std::log( (std::numeric_limits<double>::denorm_min()>0.0?
                                           std::numeric_limits<double>::denorm_min()
                                           : std::numeric_limits<double>::min()) ) + 5.0;
  auto clampExp = [&](double v){ return std::max(EXP_MIN, std::min(EXP_MAX, v)); };

  // 1) Direct evaluation
  double direct = 0.0;
  double k1v = bessel_k1(x);
  if ( std::isfinite(k1v) && k1v > 0.0 ) {
    double val = std::exp( clampExp(var) ) * k1v;
    if ( std::isfinite(val) && val > 0.0 )
      direct = val;
  }

  // 2) Asymptotic evaluation (no thresholds): K1(x) ~ sqrt(pi/(2x)) * exp(-x)
  double asymp = 0.0;
  if ( x > 0.0 ) {
    double val = std::exp( clampExp(var - x) ) * std::sqrt( NC::kPi / ( 2.0 * x ) );
    if ( std::isfinite(val) && val > 0.0 )
      asymp = val;
  }

  // Choose the stable, positive result with larger magnitude
  if (direct>0.0 && asymp>0.0)
    return std::max(direct,asymp);
  return (direct>0.0) ? direct : asymp;
}

// Egelstaff diffusion model (NJOY16 manual, Eq. 536, p.660)
double egelstaffDiffusion(double alpha, double beta, double c, double wt_diffusion) {
  const double pi = NC::kPi;
  
  // Early return for extreme alpha values to prevent overflow 
  if (alpha > 10000.0) {
    return 0.0;  // For very large alpha, diffusion becomes negligible
  }
  
  // Check exponential argument to prevent overflow 
  double exp_arg = 2.0 * c * c * wt_diffusion * alpha - beta / 2.0;
  
  // Calculate sqrt(c² + 0.25)
  double sqrt_c_term = std::sqrt(c * c + 0.25);
  
  // Calculate sqrt(β² + 4c²wt_diffusion²α²)
  double sqrt_beta_term = std::sqrt(beta * beta + 4.0 * c * c * wt_diffusion * wt_diffusion * alpha * alpha);
  
  // Calculate the argument for K₁: sqrt(c² + 0.25) * sqrt(β² + 4c²wt_diffusion²α²)
  double bessel_arg = sqrt_c_term * sqrt_beta_term;
  
  // Evaluate exp(exp_arg) * K1(bessel_arg) in a stable way 
  double expk1 = stable_exp_times_K1( exp_arg, bessel_arg );
  
  // Combine all terms: (2c*wt_diffusion*alpha/pi) * expk1 * (sqrt_c_term/sqrt_beta_term)
  double result = (2.0 * c * wt_diffusion * alpha / pi) * expk1 * (sqrt_c_term / sqrt_beta_term);
  
  // Final safety check
  if (!std::isfinite(result)) return 0.0;
  return result;
}

std::vector<double> getColumn(const std::vector<std::vector<double>>& matrix, size_t j) {
    std::vector<double> result;
    result.reserve(matrix.size());
    
    for (const auto& row : matrix) {
        if (row.size() <= j) {
            throw std::out_of_range("Column index out of bounds");
        }
        result.push_back(row[j]);
    }
    return result;
}

// Convolution of SAB data with Egelstaff diffusion model
NC::SABData convoluteSABData(const NC::SABData& sabdata, const ConvolutionParams& params) {
  const auto& alphaGrid = sabdata.alphaGrid();
  const auto& betaGrid = sabdata.betaGrid();
  const auto& sabValues = sabdata.sab();
  
  // Create convolution grids
  ConvolutionGrids grids = createConvolutionGrids(betaGrid);
  const auto& alpha_out = grids.alpha_out;
  const auto& beta_out = grids.beta_out;
  
  // Pre-compute convolution data
  ConvolutionPrecomputed precomp = precomputeConvolutionData(grids, alphaGrid, params);
  
  // Extract commonly used values
  size_t n_alpha_out = alpha_out.size();
  size_t n_beta_out = beta_out.size();
  auto& Sab = precomp.Sab;
  
  // delta_beta_vdos is passed as parameter from VDOS calculation

  // Useful constants
  const double mass_neutron_amu = NC::const_neutron_mass_amu;  // neutron mass in AMU
  const double mass_neutron = (mass_neutron_amu * NC::constant_dalton2eVc2) / 
                               ((NC::constant_c * 1e-12) * (NC::constant_c * 1e-12));  // eV·ps²/Å²
  const double hbar = NC::constant_planck / NC::k2Pi * 1e12;  // eV·ps
  const double bk = NC::constant_boltzmann;  // eV/K
  
  // Pre-allocate reusable vectors outside main loop
  std::vector<double> f1_beta_out(n_beta_out);
  std::vector<double> beta_prime;
  std::vector<double> f1_beta_prime;
  beta_prime.reserve(1000); // Reserve space to avoid reallocations
  f1_beta_prime.reserve(1000);
  
  // Main convolution loop
  for (size_t i = 0; i < n_alpha_out; ++i) {
    ConvolutionIterationParams iterParams = {
      i, alpha_out[i], alpha_out, beta_out, alphaGrid, betaGrid, sabValues, 
      precomp, params,
      f1_beta_out, beta_prime, f1_beta_prime
    };
    processConvolutionIteration(iterParams);
  }
  


  
  // Convert 2D Sab array to 1D format for SABData: S(a1,b1), S(a2,b1), S(a3,b1)... S(aA,b1), S(a1,b2), S(a2,b2)...
  NC::VectD convolvedSAB(n_alpha_out * n_beta_out);
  for (size_t j = 0; j < n_beta_out; ++j) {
    for (size_t i = 0; i < n_alpha_out; ++i) {
      convolvedSAB[j * n_alpha_out + i] = Sab[i][j];  // Beta varies slowest, alpha varies fastest
    }
  }
  
  // Create new SABData with convolved values using new grids
  return NC::SABData(NC::VectD(alpha_out), NC::VectD(beta_out), std::move(convolvedSAB),
                     sabdata.temperature(), sabdata.boundXS(), sabdata.elementMassAMU(),
                     sabdata.suggestedEmax());
}

// Squared Clebsch-Gordan coefficient C²(j1,j2,l;0,0,0) = |<j1,0;j2,0|l,0>|²
// Computed via the log-gamma form of DLMF §34.3, Eq. 34.3.5,
// together with the standard Clebsch-Gordan / 3j-symbol relation.
// See https://dlmf.nist.gov/34.3
// Returns 0 if triangle rule or parity rule (j1+j2+l odd) is violated.
static double cg00sq(int j1, int j2, int l)
{
  if (l < std::abs(j1-j2) || l > j1+j2) return 0.0;
  if ((j1+j2+l) % 2 != 0)               return 0.0;
  const int s = (j1+j2+l) / 2;
  const double log_C2 =
      std::log(2*l+1)
    + std::lgamma(j1+j2-l+1) + std::lgamma(j1-j2+l+1) + std::lgamma(-j1+j2+l+1)
    + 2.0*std::lgamma(s+1)
    - std::lgamma(j1+j2+l+2)
    - 2.0*std::lgamma(s-j1+1) - 2.0*std::lgamma(s-j2+1) - 2.0*std::lgamma(s-l+1);
  return std::exp(log_C2);
}

// Applies the Skold or Vineyard correction for a given alpha point.
//   Skold:    S(Q)·S_s(α',β) − S_s(α,β),  α' = α/S(Q)
//   Vineyard: (S(Q)−1)·S_s(α,β)
std::vector<double> coherentCorrection(const NC::SABData& s_s, const NCP::LiquidInfo& liq,
                                       int ia, double alpha, double Q)
{
  const auto& alpha_grid = s_s.alphaGrid();
  const auto& beta_grid  = s_s.betaGrid();
  const auto& sab_ss     = s_s.sab();
  const int   n_alpha    = (int)alpha_grid.size();
  const int   n_beta     = (int)beta_grid.size();

  // Interpolate S(Q) from the input table
  const double SQ = interpolate_solid(liq.qval, liq.Sval, Q);
  const double SQ_eff = std::max(SQ, 1e-8);  // clamp to avoid divergence of alpha' when S(Q)->0

  // alpha' = alpha/S(Q) for Skold; alpha' = alpha for Vineyard (no shift)
  const bool   is_skold = (liq.coherent_model == "SKOLD");
  const double alpha_p  = is_skold ? alpha / SQ_eff : alpha;

  // Interpolate S_s at shifted alpha for fixed beta row ib
  auto interp_ss = [&](double alpha_val, int ib) -> double {
    if (alpha_val <= alpha_grid.front()) return sab_ss[ib * n_alpha];
    if (alpha_val >= alpha_grid.back())  return sab_ss[ib * n_alpha + (n_alpha - 1)];
    int lo = 0, hi = n_alpha - 1;
    while (hi - lo > 1) { int mid = (lo+hi)/2; if (alpha_grid[mid] <= alpha_val) lo=mid; else hi=mid; }
    const double t = (alpha_val - alpha_grid[lo]) / (alpha_grid[hi] - alpha_grid[lo]);
    return (1.0-t)*sab_ss[ib*n_alpha+lo] + t*sab_ss[ib*n_alpha+hi];
  };

  std::vector<double> bracket(n_beta);
  for (int ib = 0; ib < n_beta; ++ib) {
    const double Ss_alpha  = sab_ss[ib * n_alpha + ia];
    const double Ss_alphap = interp_ss(alpha_p, ib);
    bracket[ib] = SQ * Ss_alphap - Ss_alpha;
  }
  return bracket;
}

// Applies the Young-Koppel rotational model (Young & Koppel, Phys. Rev. 135, A603, 1964)
// to compute S_rot(alpha,beta) for liquid ortho/para H2 or D2, convolved with
// the translational S(alpha,beta). The coherent intermolecular correction (Skold or
// Vineyard) is applied separately via coherentCorrection().
NC::SABData applyYoungKoppel(const NC::SABData& sabdata,
                              const std::string& yk_model,
                              const NC::AtomData& atomData,
                              const NCP::LiquidInfo& liq)
{
  // Determine H2 or D2 from the atomic mass (H ~1 amu, D ~2 amu)
  const double mass_atom_amu = atomData.averageMassAMU().dbl();
  if (mass_atom_amu > 2.5)
    NCRYSTAL_THROW2(BadInput,"yk_model is only supported for H2 or D2, got atom mass "
                    <<mass_atom_amu<<" amu");
  const bool   is_H2   = (mass_atom_amu < 1.5);
  const bool   is_ortho = (yk_model == "ORTHO");
  const std::string species = is_H2 ? "H2" : "D2";
  const std::string spin    = is_ortho ? "ortho" : "para";


  // Physical parameters

  const double T_K   = sabdata.temperature().dbl();
  const double kT    = NC::constant_boltzmann * T_K;   // eV

  const double B_rot = is_H2 ? 7.35e-3 : 3.71e-3;     // eV,  E_J = J(J+1)*B_rot
  const double M_amu = 2.0 * mass_atom_amu;            // molecular mass [amu]
  const double a_cm  = 0.3707;                         // atom-to-CM distance [Å] (half the H-H bond length)
  const int    J_max = 10;                             // truncation limit for rotational states
  const double hbar         = NC::constant_planck / NC::k2Pi * 1e12;              // eV·ps
  const double mass_neutron = (NC::const_neutron_mass_amu * NC::constant_dalton2eVc2)
                              / ((NC::constant_c*1e-12) * (NC::constant_c*1e-12)); // eV·ps²/Å²

  // Allowed J values: ortho-H2: odd, para-H2: even, ortho-D2: even, para-D2: odd
  const int J_start = ((is_H2 && is_ortho) || (!is_H2 && !is_ortho)) ? 1 : 0;

  // P_J: statistical weight of rotational state J (Young & Koppel, Eq. 3),
  // summed over allowed J values only (ortho or para selection rules).
  std::vector<double> P_J(J_max + 1, 0.0);
  double Z = 0.0;
  for (int J = J_start; J <= J_max; J += 2) {
    P_J[J] = (2*J + 1) * std::exp(-J*(J+1)*B_rot / kT);
    Z += P_J[J];
  }
  for (int J = J_start; J <= J_max; J += 2)
    P_J[J] /= Z;

  // Scattering lengths squared from NCrystal atom data [barn]
  const double sigma_coh = atomData.scatteringXS().dbl() - atomData.incoherentXS().dbl();
  const double sigma_inc = atomData.incoherentXS().dbl();
  const double bc2 = sigma_coh / (2.0 * NC::k2Pi);  // b_c² [barn]
  const double bi2 = sigma_inc / (2.0 * NC::k2Pi);  // b_i² [barn]

  // Spin-statistics weights (NJOY LEAPR table, p.667)
  // A: coefficient for even-J' transitions
  // B: coefficient for odd-J'  transitions
  double A_coeff, B_coeff;
  if      ( is_H2 && !is_ortho) { A_coeff = bc2;                  B_coeff = bi2;                 }  // para-H2
  else if ( is_H2 &&  is_ortho) { A_coeff = bi2/3.0;              B_coeff = bc2+(2.0/3.0)*bi2;   }  // ortho-H2
  else if (!is_H2 && !is_ortho) { A_coeff = (3.0/4.0)*bi2;        B_coeff = bc2+(1.0/4.0)*bi2;   }  // para-D2
  else                          { A_coeff = bc2+(5.0/8.0)*bi2;    B_coeff = (3.0/8.0)*bi2;        }  // ortho-D2

  // Pre-compute C²(J,J',l;0,0,0) for all J, J', l into a lookup table.
  // cg_sq[J][Jp][l],  J,J' in [0,J_max],  l in [0, 2*J_max]
  const int l_max = 2 * J_max;
  std::vector<std::vector<std::vector<double>>> cg_sq(
    J_max+1, std::vector<std::vector<double>>(
      J_max+1, std::vector<double>(l_max+1, 0.0)));

  for (int J = 0; J <= J_max; ++J)
    for (int Jp = 0; Jp <= J_max; ++Jp)
      for (int l = std::abs(J-Jp); l <= J+Jp; l += 2)
        cg_sq[J][Jp][l] = cg00sq(J, Jp, l);

  const auto& alpha_grid = sabdata.alphaGrid();
  const int   n_alpha    = (int)alpha_grid.size();

  // Momentum-transfer dependent factor for the J→J' rotational transition.
  // G_fac[J][Jp][ia] = sum_l 4*j_l²(y)*C²(J,J',l;0,0,0) (NJOY16, Eqs. 572-573, p.667)
  // where y = a_cm*sqrt(2*m_n*kT*alpha)/hbar  
  std::vector<std::vector<std::vector<double>>> G_fac(
    J_max+1, std::vector<std::vector<double>>(J_max+1, std::vector<double>(n_alpha, 0.0)));

  for (int ia = 0; ia < n_alpha; ++ia) {
    const double y = a_cm * std::sqrt(2.0 * mass_neutron * kT * alpha_grid[ia]) / hbar;
    for (int J = 0; J <= J_max; ++J)
      for (int Jp = 0; Jp <= J_max; ++Jp)
        for (int l = std::abs(J-Jp); l <= J+Jp; l += 2) {
          double jl = spherical_jn(l, y);
          G_fac[J][Jp][ia] += 4.0 * jl*jl * cg_sq[J][Jp][l];
        }
  }

  // beta_shift[J][Jp]: dimensionless rotational energy transfer for J→J' transition,
  // used to shift the translational S(alpha,beta) in the Young-Koppel sum.
  std::vector<std::vector<double>> beta_shift(J_max+1, std::vector<double>(J_max+1, 0.0));
  for (int J  = 0; J  <= J_max; ++J)
    for (int Jp = 0; Jp <= J_max; ++Jp)
      beta_shift[J][Jp] = (Jp*(Jp+1) - J*(J+1)) * B_rot / kT;

  const auto& beta_grid = sabdata.betaGrid();
  const auto& sab_vals  = sabdata.sab();
  const int   n_beta    = (int)beta_grid.size();

  // Returns S_conv(alpha[ia], beta_target) by linear interpolation.
  // Returns 0 if beta_target is outside the grid.
  auto interp_sconv = [&](int ia, double beta_target) -> double {
    if (beta_target <= beta_grid.front() || beta_target >= beta_grid.back())
      return 0.0;
    int lo = 0, hi = n_beta - 1;
    while (hi - lo > 1) {
      int mid = (lo + hi) / 2;
      if (beta_grid[mid] <= beta_target) lo = mid;
      else                               hi = mid;
    }
    const double t = (beta_target - beta_grid[lo]) / (beta_grid[hi] - beta_grid[lo]);
    return (1.0-t) * sab_vals[lo*n_alpha + ia]
           +    t  * sab_vals[hi*n_alpha + ia];
  };

  // Summation over rotational transitions J→J'
  // S_YK(alpha,beta) = sum_J P_J * sum_J' W_JJ' * S_conv(alpha, beta+beta_JJ') * G(J,J',alpha)
  // Divide by b_sq_sum so the result is the dimensionless S(alpha,beta)

  // b_sq_sum = bc2 + bi2: normalization factor for the YK sum
  const double b_sq_sum = bc2 + bi2;

  NC::VectD sab_yk(n_alpha * n_beta, 0.0);

  for (int ia = 0; ia < n_alpha; ++ia) {
    for (int ib = 0; ib < n_beta; ++ib) {
      double S_yk = 0.0;
      for (int J = J_start; J <= J_max; J += 2) {
        for (int Jp = 0; Jp <= J_max; ++Jp) {
          const double coeff = (Jp % 2 == 0) ? A_coeff : B_coeff;
          const double W     = coeff * (2*Jp + 1);
          const double S     = interp_sconv(ia, beta_grid[ib] + beta_shift[J][Jp]);
          S_yk += 0.5 * P_J[J] * W * S * G_fac[J][Jp][ia];  // 0.5: per-atom (2 atoms per molecule)
        }
      }
      sab_yk[ib * n_alpha + ia] = S_yk / b_sq_sum;
    }
  }

  // Skold/Vineyard coherent correction: adds the intermolecular (S^out) contribution to the
  // intramolecular (S^in) term already computed above, giving S_total = S^in + S^out.
  // (IAEA International Evaluation Co-operation Vol. 42,
  // "Thermal Scattering Law S(alpha,beta): Measurement, Evaluation and Application", eq. 7, p.20)
  // S^out(α,β) = coh_prefactor · j₀²(Q·a_cm) · { S(Q)·S_s(α',β) − S_s(α,β) }
  // coh_prefactor = 2 b_c²/b_sq_sum,  α' = α/S(Q),  Q = √(2 m_n kT α)/ħ
  if (!liq.qval.empty()) {
    const double coh_prefactor = 2.0 * bc2 / b_sq_sum;  // = 2*cfrac
    for (int ia = 0; ia < n_alpha; ++ia) {
      const double alpha  = alpha_grid[ia];
      const double Q      = std::sqrt(2.0 * mass_neutron * kT * alpha) / hbar;
      const double Qb     = Q * a_cm;
      const double j0_Qb  = (Qb > 1.0e-10) ? std::sin(Qb) / Qb : 1.0;
      const auto   bracket = coherentCorrection(sabdata, liq, ia, alpha, Q);
      for (int ib = 0; ib < n_beta; ++ib) {
        sab_yk[ib * n_alpha + ia] += coh_prefactor * j0_Qb * j0_Qb * bracket[ib];
        sab_yk[ib * n_alpha + ia]  = std::max(0.0, sab_yk[ib * n_alpha + ia]);
      }
    }
  }

  return NC::SABData( NC::VectD(alpha_grid),
                      NC::VectD(beta_grid),
                      std::move(sab_yk),
                      sabdata.temperature(),
                      sabdata.boundXS(),
                      sabdata.elementMassAMU(),
                      sabdata.suggestedEmax() );
}
