
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

double spherical_jn(int n, double x);
double bessel_k1(double x);

void NCP::customPluginTest()
{
  //This function is called by NCrystal after the plugin is loaded, but only if
  //the NCRYSTAL_PLUGIN_RUNTESTS environment variable is set to "1". In case of
  //errors or anything about the plugin not working, simply throw an exception
  //or use nc_assert_always.

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

  NCPLUGIN_MSG("Cross-section reference checks passed");

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

  const double specialfunc_abs_tol = 1.0e-12;
  const double specialfunc_rel_tol = 1.0e-6;

  // Spherical Bessel check against reference values.
  // Reference values generated with Python scipy.special.spherical_jn.
  NCPLUGIN_MSG("Checking spherical Bessel values against reference values");

  struct SphericalJnTestCase {
    int n;
    double x;
    double ref;
  };

  const SphericalJnTestCase spherical_jn_testcases[] = {
    { 0, 0.0, 1.0 },
    { 0, 9.9999999999999995e-07, 0.99999999999983336 },
    { 0, 0.0001, 0.99999999833333342 },
    { 0, 0.01, 0.99998333341666645 },
    { 0, 0.10000000000000001, 0.99833416646828155 },
    { 0, 1.0, 0.8414709848078965 },
    { 0, 5.0, -0.1917848549326277 },
    { 0, 10.0, -0.054402111088936979 },
    { 0, 20.0, 0.045647262536381385 },
    { 1, 0.0, 0.0 },
    { 1, 9.9999999999999995e-07, 3.3333333333330117e-07 },
    { 1, 0.0001, 3.3333333300000061e-05 },
    { 1, 0.01, 0.0033333000001190523 },
    { 1, 0.10000000000000001, 0.033300011902557616 },
    { 1, 1.0, 0.30116867893975707 },
    { 1, 5.0, -0.095089408079170795 },
    { 1, 10.0, 0.078466941798751549 },
    { 1, 20.0, -0.018121739963850528 },
    { 2, 0.0, 0.0 },
    { 2, 0.0001, 6.666666661904783e-10 },
    { 2, 0.01, 6.6666190477513343e-06 },
    { 2, 0.10000000000000001, 0.00066619060844556942 },
    { 2, 1.0, 0.062035052011373916 },
    { 2, 5.0, 0.13473121008512523 },
    { 2, 10.0, 0.077942193628562445 },
    { 2, 20.0, -0.048365523530958965 },
    { 5, 0.0, 0.0 },
    { 5, 0.10000000000000001, 9.6163102329164501e-10 },
    { 5, 1.0, 9.2561158611258252e-05 },
    { 5, 5.0, 0.1068111614565046 },
    { 5, 10.0, -0.055534511621452162 },
    { 5, 20.0, 0.016683908063095696 },
    { 10, 0.0, 0.0 },
    { 10, 1.0, 7.1165526400473407e-11 },
    { 10, 5.0, 0.00040734424424946199 },
    { 10, 10.0, 0.064605154492563974 },
    { 10, 20.0, 0.039686698644626373 },
    { 20, 0.0, 0.0 },
    { 20, 5.0, 5.4277267607932615e-12 },
    { 20, 10.0, 2.3083719613194551e-06 },
    { 20, 20.0, 0.038324851639805507 }
  };

  for ( const auto& testcase : spherical_jn_testcases ) {
    const double val = spherical_jn(testcase.n, testcase.x);
    const double ref = testcase.ref;
    const double diff = std::abs(val-ref);
    if ( !std::isfinite(val)
         || ( diff > specialfunc_abs_tol
              && diff > specialfunc_rel_tol * std::abs(ref) ) ) {
      NCPLUGIN_MSG( "Spherical Bessel check failed for n=" << testcase.n
                    << ", x=" << testcase.x
                    << ": value=" << val
                    << ", ref=" << ref
                    << ", absdiff=" << diff
                    << ", abs_tol=" << specialfunc_abs_tol
                    << ", rel_tol=" << specialfunc_rel_tol );
      nc_assert_always( false );
    }
  }

  NCPLUGIN_MSG("Spherical Bessel checks passed");

  // Modified Bessel K1 check against reference values.
  // Reference values generated with Python scipy.special.k1.
  NCPLUGIN_MSG("Checking K1 values against reference values");

  struct K1TestCase {
    double x;
    double ref;
  };

  const K1TestCase k1_testcases[] = {
    { 9.9999999999999995e-07, 999999.99999278435 },
    { 0.0001, 9999.9995086864055 },
    { 0.01, 99.973894118296229 },
    { 0.10000000000000001, 9.853844780870606 },
    { 1.0, 0.60190723019723458 },
    { 2.0, 0.13986588181652246 },
    { 5.0, 0.0040446134454521629 },
    { 10.0, 1.8648773453825585e-05 },
    { 20.0, 5.8830579695570384e-10 }
  };

  for ( const auto& testcase : k1_testcases ) {
    const double val = bessel_k1(testcase.x);
    const double ref = testcase.ref;
    const double diff = std::abs(val-ref);
    if ( !std::isfinite(val)
         || ( diff > specialfunc_abs_tol
              && diff > specialfunc_rel_tol * std::abs(ref) ) ) {
      NCPLUGIN_MSG( "K1 check failed for x=" << testcase.x
                    << ": value=" << val
                    << ", ref=" << ref
                    << ", absdiff=" << diff
                    << ", abs_tol=" << specialfunc_abs_tol
                    << ", rel_tol=" << specialfunc_rel_tol );
      nc_assert_always( false );
    }
  }

  NCPLUGIN_MSG("K1 checks passed");
  NCPLUGIN_MSG("All tests of plugin were successful!");
}
