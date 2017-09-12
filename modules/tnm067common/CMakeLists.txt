#--------------------------------------------------------------------
# Inviwo TNM067Common Module
ivw_module(TNM067Common)

#--------------------------------------------------------------------
# Add header files
set(HEADER_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/tnm067commonprocessor.h
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/gauss2dfunction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/test2by2image.h
)
ivw_group("Header Files" ${HEADER_FILES})

#--------------------------------------------------------------------
# Add source files
set(SOURCE_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/tnm067commonprocessor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/gauss2dfunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/processors/test2by2image.cpp
)
ivw_group("Source Files" ${SOURCE_FILES})


#--------------------------------------------------------------------
# Add shaders
set(SHADER_FILES
    #${CMAKE_CURRENT_SOURCE_DIR}/glsl/tnm067commonprocessor.frag
)
ivw_group("Shader Files" ${SHADER_FILES})


#--------------------------------------------------------------------
# Add Unittests
set(TEST_FILES
	#${CMAKE_CURRENT_SOURCE_DIR}/tests/tnm067common-test.cpp
)
ivw_add_unittest(${TEST_FILES})

#--------------------------------------------------------------------
# Create module
ivw_create_module(${SOURCE_FILES} ${HEADER_FILES} ${SHADER_FILES})

#--------------------------------------------------------------------
# Add shader directory to pack
# ivw_add_to_module_pack(${CMAKE_CURRENT_SOURCE_DIR}/glsl)

ivw_folder(inviwo-module-tnm067common TNM067)