SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../tests")
add_executable (compare_analyt_2_interp)

target_sources (compare_analyt_2_interp
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/compare_analyt_2_interp.f90>
 )

target_include_directories (compare_analyt_2_interp
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
)
target_link_libraries (compare_analyt_2_interp
                       PUBLIC
                       RAYS_lib
)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin")
add_executable (solovev_2_eqdsk)

target_sources (solovev_2_eqdsk
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/solovev_2_eqdsk.f90>
 )

target_include_directories (solovev_2_eqdsk
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
)
target_link_libraries (solovev_2_eqdsk
                       PUBLIC
                       RAYS_lib
)
