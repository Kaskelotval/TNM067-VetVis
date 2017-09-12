#--------------------------------------------------------------------
# Inviwo TNM067Lab1 Module
ivw_module(TNM067Lab1)

#--------------------------------------------------------------------
# Add header files
set(HEADER_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/tnm067lab1processor.h
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/imagemappingcpu.h
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/imagetoheightfield.h
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/imageupsampler.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/interpolationmethods.h
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/scalartocolormapping.h
)
ivw_group("Header Files" ${HEADER_FILES})

#--------------------------------------------------------------------
# Add source files
set(SOURCE_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/tnm067lab1processor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/imagemappingcpu.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/imagetoheightfield.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/imageupsampler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/scalartocolormapping.cpp
)
ivw_group("Source Files" ${SOURCE_FILES})


#--------------------------------------------------------------------
# Add shaders
set(SHADER_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/glsl/tnm067lab1processor.frag
)
ivw_group("Shader Files" ${SHADER_FILES})


#--------------------------------------------------------------------
# Add Unittests
set(TEST_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/unittests/imageupsampler-test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/unittests/interploation-test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/unittests/tnm067lab1-unittest-main.cpp
)
ivw_add_unittest(${TEST_FILES})


#--------------------------------------------------------------------
# Create module
ivw_create_module(${SOURCE_FILES} ${HEADER_FILES} ${SHADER_FILES})

#--------------------------------------------------------------------
# Add shader directory to pack
# ivw_add_to_module_pack(${CMAKE_CURRENT_SOURCE_DIR}/glsl)

ivw_folder(inviwo-module-tnm067lab1 TNM067)