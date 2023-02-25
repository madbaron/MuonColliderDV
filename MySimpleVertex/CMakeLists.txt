########################################################
# cmake file for building Marlin example Package
# @author Federico Meloni
CMAKE_MINIMUM_REQUIRED(VERSION 3.15 FATAL_ERROR)
########################################################

# project name
PROJECT( MySimpleVertex )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
#INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

find_package( ROOT REQUIRED )
set( ROOT_COMPONENT_LIBRARIES Geom Reflex)
INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

FIND_PACKAGE( Acts REQUIRED COMPONENTS Core PluginJson PluginTGeo )
INCLUDE_DIRECTORIES( ${Acts_INCLUDE_DIRS} )

# optional package
#FIND_PACKAGE( AIDA )
#IF( AIDA_FOUND )
#    INCLUDE_DIRECTORIES( ${AIDA_INCLUDE_DIRS} )
#    LINK_LIBRARIES( ${AIDA_LIBRARIES} )
#    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
#    MESSAGE( STATUS "AIDA -- found" )
#ELSE()
#    MESSAGE( STATUS "AIDA -- not found" )
#ENDIF()

INCLUDE(GNUInstallDirs)

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

ADD_LIBRARY(${PROJECT_NAME} SHARED
    ./src/MyVertexFinder.cc  
    ./src/VertexFitterActs.cc	
    ./src/VertexFitterSimple.cc  
    ./src/VertexingUtils.cc
    ./src/DVFinder.cc
    ./src/Sel2TrkVertices.cc
    ./src/Helpers.cc
    ./src/IterativeVertexFinder.cc
)

TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME}
  PRIVATE
    ./include
    ${CMAKE_BINARY_DIR}
    ${Marlin_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS})

TARGET_LINK_LIBRARIES(${PROJECT_NAME}
  ${ROOT_LIBRARIES} DD4hep::DDCore
  ActsCore ActsPluginTGeo ActsPluginJson)

INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib )


