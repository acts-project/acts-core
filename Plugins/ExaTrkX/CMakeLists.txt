add_library(ActsPluginExaTrkX SHARED src/buildEdges.cpp src/ExaTrkXPipeline.cpp)

if(ACTS_EXATRKX_ENABLE_CUDA)
    target_compile_definitions(ActsPluginExaTrkX PUBLIC ACTS_EXATRKX_WITH_CUDA)
    target_sources(ActsPluginExaTrkX PRIVATE src/CudaTrackBuilding.cu)
endif()

if(ACTS_EXATRKX_ENABLE_ONNX)
    target_sources(
        ActsPluginExaTrkX
        PRIVATE src/OnnxEdgeClassifier.cpp src/OnnxMetricLearning.cpp
    )
endif()

if(ACTS_EXATRKX_ENABLE_TORCH)
    target_sources(
        ActsPluginExaTrkX
        PRIVATE
            src/TorchEdgeClassifier.cpp
            src/TorchMetricLearning.cpp
            src/BoostTrackBuilding.cpp
            src/TorchTruthGraphMetricsHook.cpp
            src/TorchGraphStoreHook.cpp
    )
endif()

if(ACTS_EXATRKX_ENABLE_TENSORRT)
    find_package(TensorRT REQUIRED)
    message(STATUS "Found TensorRT ${TensorRT_VERSION}")
    target_link_libraries(
        ActsPluginExaTrkX
        PUBLIC trt::nvinfer trt::nvinfer_plugin
    )
    target_sources(ActsPluginExaTrkX PRIVATE src/TensorRTEdgeClassifier.cpp)
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC ACTS_EXATRKX_WITH_TENSORRT
    )
endif()

target_include_directories(
    ActsPluginExaTrkX
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
    ActsPluginExaTrkX
    PUBLIC
        ActsCore
        Boost::boost
        ${TORCH_LIBRARIES} # TODO try to make this private again (torch::Device is leaking)
        std::filesystem
)

if(ACTS_EXATRKX_ENABLE_CUDA)
    target_link_libraries(ActsPluginExaTrkX PRIVATE frnn)
    target_compile_features(ActsPluginExaTrkX PUBLIC cuda_std_20)
    set_target_properties(
        ActsPluginExaTrkX
        PROPERTIES CUDA_STANDARD_REQUIRED ON CUDA_SEPARABLE_COMPILATION ON
    )
    target_compile_options(
        ActsPluginExaTrkX
        PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:-g
            --generate-line-info
            --extended-lambda>
    )
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC CUDA_API_PER_THREAD_DEFAULT_STREAM
    )
else()
    target_compile_definitions(ActsPluginExaTrkX PUBLIC ACTS_EXATRKX_CPUONLY)
endif()

if(ACTS_EXATRKX_ENABLE_ONNX)
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC ACTS_EXATRKX_ONNX_BACKEND
    )

    target_link_libraries(ActsPluginExaTrkX PRIVATE OnnxRuntime)
endif()

if(ACTS_EXATRKX_ENABLE_TORCH)
    target_compile_definitions(
        ActsPluginExaTrkX
        PUBLIC ACTS_EXATRKX_TORCH_BACKEND
    )

    target_link_libraries(ActsPluginExaTrkX PRIVATE TorchScatter::TorchScatter)

    # Should not discard TorchScatter even if its not needed at this point
    # since we need the scatter_max operation in the torch script later
    target_link_options(ActsPluginExaTrkX PUBLIC "-Wl,-no-as-needed")
endif()

install(
    TARGETS ActsPluginExaTrkX
    EXPORT ActsPluginExaTrkXTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
