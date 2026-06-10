
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

void NCP::customPluginTest()
{
  //This function is called by NCrystal after the plugin is loaded, but only if
  //the NCRYSTAL_PLUGIN_RUNTESTS environment variable is set to "1". In case of
  //errors or anything about the plugin not working, simply throw an exception
  //(which is what the nc_assert_always function does below, but feel free to
  //simply throw an exception directly).

  //Note, emit all messages here and elsewhere in plugin code with NCPLUGIN_MSG
  //(or NCPLUGIN_WARN for warnings), never raw usage of std::cout or printf!

  NCPLUGIN_MSG("Testing plugin");

  //Create some test NCMAT data. For simplicity we will here base it on an
  //existing file, but add our @CUSTOMPLUGIN section (see NCPhysicsModel.cc for
  //a description of the format).:

  auto scat_somefile
    = NC::createScatter( "plugins::Liquids/D2O_liquid.ncmat");

  for ( auto& wlval : NC::linspace(0.1,3.1,16) ) {
    auto wl = NC::NeutronWavelength( wlval );
    auto xs = scat_somefile.crossSectionIsotropic( wl );
    NCPLUGIN_MSG( "D2O_liquid.ncmat xs @ "<<wl<<" : "<<xs);
  }

  //Note: We do not test the scattering here for this simple dummy plugin.

  NCPLUGIN_MSG("All tests of plugin were successful!");
}
