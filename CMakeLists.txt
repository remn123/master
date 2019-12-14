cmake_minimum_required(VERSION 3.10.1)

set(CMAKE_CXX_COMPILER "g++-7")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -O -Wall -Wextra -lstdc++fs -g -lm")
link_libraries(stdc++fs)

enable_language(CXX)

# enabling test
enable_testing()

# main project
project(master)

# master 
add_executable(master_exe 
               src/Main.cpp
)

target_link_libraries(master_exe
                      master_lib   
)

add_library(master_lib
            src/Mesh.cpp
            src/Element.cpp
            src/Pipeline.cpp
            src/Poly.cpp
            include/Dict.h
            include/DVector.h
            include/Helpers.h
            include/Pipeline.h
            include/Poly.h
            include/Mesh.h

)

target_include_directories(master_lib PUBLIC include)
target_include_directories(master_lib PUBLIC src)

# tests
add_executable(master_test 
               tests/main.cpp
               tests/tests.cpp
               tests/test_helpers.cpp
               tests/test_dict.cpp
               tests/test_dvector.cpp

)

target_link_libraries(master_test master_lib)
target_link_libraries(master_test catch_lib)

add_library(catch_lib
            thirdparty/catch/catch.hpp
)

set_target_properties(catch_lib PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(catch_lib PUBLIC tests)
target_include_directories(catch_lib PUBLIC thirdparty)


# add test
add_test(master_test master_test)
