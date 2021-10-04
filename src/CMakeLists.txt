add_library(mc2err
            mc2err_analyze.c
            mc2err_append.c
            mc2err_begin.c
            mc2err_end.c
            mc2err_expand.c
            mc2err_input.c
            mc2err_likelihood.c
            mc2err_load.c
            mc2err_output.c
            mc2err_save.c)

target_include_directories(mc2err PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(mc2err PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
