add_executable (post_process_RAYS)

target_sources (post_process_RAYS
                PRIVATE
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/post_process_RAYS.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/post_processing_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/intialize_post_process_RAYS.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/slab_processor_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/solovev_processor_m.f90>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/finalize_post_process_RAYS.f90>
)

target_include_directories (post_process_RAYS
                            PUBLIC
                            $<TARGET_PROPERTY:RAYS_lib,BINARY_DIR>
)
target_link_libraries (post_process_RAYS
                       PUBLIC
                       RAYS_lib
)
