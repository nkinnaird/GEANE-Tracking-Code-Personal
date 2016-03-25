#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(toymodelstrawerrprop)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include/eigen3")

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
#find_package(Geant4 REQUIRED)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#

find_package(CLHEP REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR}
                    ${CLHEP_INCLUDE_DIR}
                    /usr/local/include/eigen3)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(basic basic.cc ${sources} ${headers})
target_link_libraries(basic ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES})

add_executable(errprop errprop.cc ${sources} ${headers})
target_link_libraries(errprop ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build toymodelbasicerrprop. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(toymodelstrawerrprop_SCRIPTS
   
  )

foreach(_script ${toymodelstrawerrprop_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(toymodelstrawerrprop DEPENDS basic errprop)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS basic DESTINATION bin)
install(TARGETS errprop DESTINATION bin)

