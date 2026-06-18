
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

#include "NCTestPlugin.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCMath.hh"

#include <cmath>
#include <string>

void NCP::customPluginTest()
{
  //This function is called by NCrystal after the plugin is loaded, but only if
  //the NCRYSTAL_PLUGIN_RUNTESTS environment variable is set to "1". In case of
  //errors or anything about the plugin not working, simply throw an exception
  //(or use nc_assert_always).

  //Note, emit all messages here and elsewhere in plugin code with NCPLUGIN_MSG
  //(or NCPLUGIN_WARN for warnings), never raw usage of std::cout or printf!

  // Cross-section check against v1.0.0 reference values.
  NCPLUGIN_MSG("Checking cross sections against v1.0.0 reference values");
  const double xs_reldiff_tol = 1.0e-4;
  const double energies[] = {
    1.0000000000000001e-05,
    3.3274211919892401e-05,
    0.00011071731788899118,
    0.00036840314986403866,
    0.0012258324480531914,
    0.0040788608654802379,
    0.013572088082974531,
    0.045160053506834122,
    0.15026651907000801,
    0.49999999999999994
  };

  struct XSTestCase {
    const char * material;
    double xs[10];
  };

  const XSTestCase testcases[] = {
    { "H2O_liquid.ncmat",
      { 349.63368106043191, 216.92701831950976, 143.304952570854,
        101.98729111570167, 77.066382309953624, 59.161536267710062,
        44.099775070284387, 28.819135478212306, 20.93788548109589,
        16.385722077650069 } },
    { "D2O_liquid.ncmat",
      { 35.000032425630152, 20.05794238740739, 11.919137569856556,
        7.6495006994483621, 6.2891032378549472, 7.6019422475262326,
        5.2975859579285869, 4.3941121829418268, 3.8602930574185108,
        3.5766275832827348 } },
    { "H2_liquid.ncmat",
      { 332.7293730565309, 196.71724840733688, 122.7345232454622,
        82.734649956490003, 60.658419580704305, 46.762143693582999,
        39.171424853950114, 37.120001510641792, 36.483155578544689,
        36.042487151863 } },
    { "tsl-para-H_20K.ncmat",
      { 1.5056216799131767, 0.83696387717080101, 0.48575807940673849,
        0.33404121416828531, 0.47363009943197687, 1.2621217491809509,
        1.7390717677008061, 21.492674645196203, 22.856932302802399,
        21.919755773610184 } },
    { "tsl-ortho-D_20K.ncmat",
      { 9.6598934667225702, 6.2190956219250095, 4.312220787620225,
        3.2770772531909249, 3.0627644386289195, 6.9229719660621667,
        5.4377855858195758, 3.6501749847573697, 3.856990864046784,
        3.4305998025225555 } }
  };

  for ( const auto& testcase : testcases ) {
    auto scatter = NC::createScatter( std::string("plugins::Liquids/")
                                      + testcase.material );
    for ( unsigned i = 0; i < 10; ++i ) {
      const double xs = scatter.crossSectionIsotropic(
        NC::NeutronEnergy{ energies[i] } ).dbl();
      const double ref = testcase.xs[i];
      const double reldiff = std::abs(xs/ref - 1.0);
      if ( !std::isfinite(xs) || !std::isfinite(ref) || ref <= 0.0
           || reldiff > xs_reldiff_tol ) {
        NCPLUGIN_MSG( "Cross-section check failed for "
                      << testcase.material
                      << " at E=" << energies[i] << " eV: xs=" << xs
                      << ", ref=" << ref
                      << ", reldiff=" << reldiff
                      << ", tolerance=" << xs_reldiff_tol );
        nc_assert_always( false );
      }
    }
  }

  NCPLUGIN_MSG("Cross-section regression checks passed");

  // Sampling sanity check.
  NCPLUGIN_MSG("Checking sampled scattering events");
  const unsigned nsamples = 20;

  for ( const auto& testcase : testcases ) {
    auto scatter = NC::createScatter( std::string("plugins::Liquids/")
                                      + testcase.material );
    for ( const auto& energy : energies ) {
      for ( unsigned i = 0; i < nsamples; ++i ) {
        auto event = scatter.sampleScatterIsotropic(
          NC::NeutronEnergy{ energy } );
        const double eout = event.ekin.dbl();
        const double mu = event.mu.dbl();
        if ( !std::isfinite(eout) || eout < 0.0
             || !std::isfinite(mu) || mu < -1.0 || mu > 1.0 ) {
          NCPLUGIN_MSG( "Sampling check failed for "
                        << testcase.material
                        << " at E=" << energy << " eV: eout=" << eout
                        << ", mu=" << mu );
          nc_assert_always( false );
        }
      }
    }
  }

  NCPLUGIN_MSG("Sampling sanity checks passed");
  NCPLUGIN_MSG("All tests of plugin were successful!");
}
