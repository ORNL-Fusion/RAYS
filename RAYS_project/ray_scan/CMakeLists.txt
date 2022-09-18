add_executable (ray_scan)


target_sources (ray_scan
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ray_scan.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ray_scanner_m.f90>
)
target_include_directories (ray_scan
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
)
target_link_libraries (ray_scan
                       PUBLIC
                       RAYS_lib
)
