add_executable (post_process_RAYS)

target_sources (post_process_RAYS
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/post_process_RAYS.f90>
)

target_include_directories (post_process_RAYS
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
                            $<TARGET_PROPERTY:post_process_lib,BINARY_DIR>
                            $<TARGET_PROPERTY:math_functions,BINARY_DIR>
)
target_link_libraries (post_process_RAYS
                      PUBLIC
                      RAYS_lib math_functions post_process_lib)
