set (lib_sources "")
list (APPEND lib_sources
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/zfunct.f>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/zfun0.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fzeta.f>

)
add_library (math_functions ${lib_sources})

#********************************************

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../tests")

add_executable (test_zfun)

target_sources (test_zfun
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_zfun.f90>
)

target_link_libraries (test_zfun
                       PUBLIC
                       math_functions
)


