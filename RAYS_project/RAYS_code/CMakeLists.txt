add_executable (RAYS)


target_sources (RAYS
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/RAYS.f90>
)
target_include_directories (RAYS
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
)
target_link_libraries (RAYS
                       PUBLIC
                       RAYS_lib
)
