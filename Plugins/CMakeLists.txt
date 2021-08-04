# plugins are manually ordered such that dependent plugins appear after their
# dependencies in the component list.
# independent plugins
add_component_if(Autodiff PluginAutodiff ACTS_BUILD_PLUGIN_AUTODIFF)
add_component_if(Cuda PluginCuda ACTS_BUILD_PLUGIN_CUDA)
add_component_if(Digitization PluginDigitization ACTS_BUILD_PLUGIN_DIGITIZATION)
add_component_if(Identification PluginIdentification ACTS_BUILD_PLUGIN_IDENTIFICATION)
add_component_if(Json PluginJson ACTS_BUILD_PLUGIN_JSON)
add_component_if(Legacy PluginLegacy ACTS_BUILD_PLUGIN_LEGACY)
add_component_if(Onnx PluginOnnx ACTS_BUILD_PLUGIN_ONNX)
add_component_if(Python PluginPython ACTS_BUILD_PLUGIN_PYTHON)
add_component_if(Sycl PluginSycl ACTS_BUILD_PLUGIN_SYCL)
# dependent plugins. depend either on a independent plugins or on one another
add_component_if(TGeo PluginTGeo ACTS_BUILD_PLUGIN_TGEO)
add_component_if(DD4hep PluginDD4hep ACTS_BUILD_PLUGIN_DD4HEP)

propagate_components_to_parent()
