cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_COMPILER "g++-9")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -O -Wall -Wextra -lstdc++fs -g -lm")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -O -Wall -Wextra -g -lm")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -O -Wall -Wextra -g -lm")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -lm")
#link_libraries(stdc++fs)

find_package(Boost 1.73.0 REQUIRED)
if (Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")

    enable_language(CXX)

    # enabling test
    enable_testing()

    # main project
    project(master)

    # master 
    # add_executable(master_exe 
    #             src/Main.cpp
    # )

    # ringleb 
    add_executable(master_exe 
                src/ringleb.cpp
    )

    target_link_libraries(master_exe
                        master_lib   
    )

    add_library(master_lib
                src/Element.cpp
                src/Ghost.cpp
                src/Mesh.cpp
                src/Node.cpp
                #src/Pipeline.cpp
                #src/Posproc.cpp
                #src/Preproc.cpp
                src/Poly.cpp
                src/SD.cpp
                #src/Solver.cpp
                include/Dict.h
                include/DVector.h
                include/Dummies.h
                include/Edge.h
                include/Element.h
                include/Ghost.h
                include/Helpers.h
                include/Mesh.h
                include/Node.h
                include/Property.h
                #include/Pipeline.h
                #include/Posproc.h
                #include/Preproc.h
                include/Poly.h
                include/SD.h
                #include/Solver.h

    )

    target_include_directories(master_lib PUBLIC include)
    target_include_directories(master_lib PUBLIC src)
    target_include_directories(master_lib PUBLIC ${Boost_INCLUDE_DIRS})

    target_link_libraries(master_exe ${Boost_LIBRARY_DIRS})


    # tests
    add_executable(master_test 
                tests/main.cpp
                tests/tests.cpp
                tests/test_dict.cpp
                tests/test_dvector.cpp
                tests/test_element.cpp
                tests/test_helpers.cpp
                tests/test_sd.cpp
                tests/test_posproc.cpp
                #tests/test_pipeline.cpp
                #tests/test_sd.cpp
                #tests/test_ublas.cpp

    )

    target_link_libraries(master_test master_lib)
    target_link_libraries(master_test catch_lib)
    target_link_libraries(master_test ${Boost_LIBRARY_DIRS})

    add_library(catch_lib
                thirdparty/catch/catch.hpp
    )

    set_target_properties(catch_lib PROPERTIES LINKER_LANGUAGE CXX)

    target_include_directories(catch_lib PUBLIC tests)
    target_include_directories(catch_lib PUBLIC thirdparty)


    # add test
    add_test(master_test master_test)
else()
    message(STATUS "Boost Library not found")
endif()