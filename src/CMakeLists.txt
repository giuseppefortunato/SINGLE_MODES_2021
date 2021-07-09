#CMake file for Adding executables

function( ADD_RBM_TESTCASE arg1)

    add_executable(${arg1} ${arg1}.cpp)
    if (BUILD_STATIC)
        target_link_libraries(${arg1} MODES_lib_static ${MANDATORY_LIBRARIES})
    else ()
        target_link_libraries(${arg1} MODES ${MANDATORY_LIBRARIES})
    endif()
endfunction()


ADD_RBM_TESTCASE(SingleMODES)
ADD_RBM_TESTCASE(SingleMODES_leave_one_out)
ADD_RBM_TESTCASE(OptimalISO)
