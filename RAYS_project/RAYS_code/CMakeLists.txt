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

add_executable (RAYS_P)


target_sources (RAYS_P
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/RAYS_P.f90>
)
target_include_directories (RAYS_P
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
                            $<TARGET_PROPERTY:post_process_lib,BINARY_DIR>
)
target_link_libraries (RAYS_P
                       PUBLIC
                       RAYS_lib
                       post_process_lib
)
