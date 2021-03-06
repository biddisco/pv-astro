PROJECT("pv-astro")

cmake_minimum_required(VERSION 2.8)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
endif(COMMAND cmake_policy)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

IF (VTK_USE_MPI)
  INCLUDE_DIRECTORIES(${VTK_MPI_INCLUDE_DIR})
  ADD_DEFINITIONS(-DUSE_MPI)
ENDIF (VTK_USE_MPI)

###############
# pv_astro plugin
# INCLUDES
#   o Tipsy binary reader
#   o Ramses binary reader
#   o Grafic binary reader
#   o color bars equivalent to those available in tipsy
#   o smoothing filter (recommended to use pv-meshless SPH Smooth utility)
#   o center of mass filter 
#   o moments of inertia filter
#   o profile filter
#   o add additional attribute filter
#   o friends-of-friends halo finder filter
# Author: Christine Corbett Moran, contributions by Rafael Kueng and John Biddiscombe
###

###############
# StrangeAttractors plugin
# o strange attractors random/binary visualization
# Author: Christine Corbett Moran
###



INCLUDE_DIRECTORIES(AstroVizHelpersLib)
INCLUDE_DIRECTORIES(fio)

# for special gui
IF(PARAVIEW_BUILD_QT_GUI)
  QT4_WRAP_CPP(MOC_SRCS pqTrackFilterPanel.h)

  ADD_PARAVIEW_OBJECT_PANEL(
    IFACES IFACE_SRCS
    CLASS_NAME pqTrackFilterPanel
    XML_NAME TrackFilter
    XML_GROUP filters
  )
  SET(PV_ASTRO_GUI_FILES pqTrackFilterPanel.cxx)

ENDIF(PARAVIEW_BUILD_QT_GUI)

#Defines the pv_astro plugin
IF (NOT WIN32)
  SET(GRAPHIC_SRCS
    vtkGraficReader.cxx    
  )
ENDIF (NOT WIN32)

  #-----------------------------------------------------------------------------
  # OpenMP
  #-----------------------------------------------------------------------------
  OPTION(ASTRO_USE_OPENMP "Compile with OpenMP support" ON)
  IF (ASTRO_USE_OPENMP)
    ADD_DEFINITIONS(-DHAVE_OPENMP)
    IF (WIN32 AND MSVC) 
      SET(PV_ASTRO_OPENMP_CXX_FLAGS "/openmp")
    ENDIF (WIN32 AND MSVC) 
    if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
      SET(PV_ASTRO_OPENMP_CXX_FLAGS "-fopenmp" )
    endif (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
  ENDIF (ASTRO_USE_OPENMP)

#--------------------------------------------------
# pv_astro Plugin
#--------------------------------------------------
ADD_PARAVIEW_PLUGIN(pv_astro "1.1"
  SERVER_MANAGER_SOURCES 
    vtkTipsyReader.cxx
    vtkNSmoothFilter.cxx
    vtkCenterOfMassFilter.cxx 
    vtkProfileFilter.cxx
    vtkMomentsOfInertiaFilter.cxx
    vtkVirialRadiusFilter.cxx
    vtkAddAdditionalAttribute.cxx 
    vtkFriendsOfFriendsHaloFinder.cxx 
    vtkRamsesReader.cxx
		vtkGadgetReader.cxx
    vtkSQLiteReader.cxx
    vtkTrackFilter.cxx
    vtkSimpleBin.cxx
    ${GRAPHIC_SRCS}

  SERVER_MANAGER_XML
    TipsyReaderSM.xml
    NSmoothFilterSM.xml
    CenterOfMassFilter.xml
    ProfileFilter.xml
    MomentsOfInertiaFilter.xml 
    VirialRadiusFilter.xml
    AddAdditionalAttribute.xml 
    FriendsOfFriendsHaloFinder.xml  
    GraficReaderSM.xml
    SQLiteReaderSM.xml
    TrackFilterSM.xml
    SimpleBinSM.xml
    RamsesReaderSM.xml
    #GadgetReaderSM.xml

  GUI_RESOURCES
    TipsyReader.qrc
    SQLiteReader.qrc
    RamsesReader.qrc
    #GadgetReader.qrc
    GraficReader.qrc    
  GUI_RESOURCE_FILES
    AstroVizFilterMenu.xml
  GUI_INTERFACES 
    ${IFACES} 
  GUI_SOURCES  
    ${PV_ASTRO_GUI_FILES}
    ${MOC_SRCS} 
    ${IFACE_SRCS}
  )
#PARAVIEW_INCLUDE_WRAPPED_SOURCES(fio/fio.c fio/romberg.c)
#Defines the StrangeAttractors plugin

ADD_PARAVIEW_PLUGIN(StrangeAttractorsPlugin "1.0"
  SERVER_MANAGER_SOURCES 
    vtkStrangeAttractors.cxx
  SERVER_MANAGER_XML
    StrangeAttractorsSM.xml
  GUI_RESOURCES
    StrangeAttractors.qrc 
  GUI_INTERFACES 
  SOURCES 
  )
  
# For helper functions often used, will later include these in a single
# VTK class.
ADD_LIBRARY(AstroVizHelpers STATIC AstroVizHelpersLib/AstroVizHelpers.cxx)

SET_TARGET_PROPERTIES(AstroVizHelpers PROPERTIES COMPILE_FLAGS "-fPIC")
TARGET_LINK_LIBRARIES(AstroVizHelpers vtkParallelMPI vtkFiltersSources)

TARGET_LINK_LIBRARIES(pv_astro AstroVizHelpers ) 
SET_TARGET_PROPERTIES(pv_astro PROPERTIES COMPILE_FLAGS "${PV_ASTRO_OPENMP_CXX_FLAGS}")

# For the Tipsy reader plugin I use some of Doug Potter's Tipsy lib.
ADD_LIBRARY(
  TipsyHelpers
  STATIC
    tipsylib/adapter.cpp
    tipsylib/binner.cpp 
    tipsylib/ftipsy.cpp tipsylib/native.cpp 
    tipsylib/standard.cpp
    tipsylib/vtipsy.cpp
  )
  
SET_TARGET_PROPERTIES(TipsyHelpers PROPERTIES COMPILE_FLAGS "-fPIC")  

# add the winsock2 library for net lookup names
IF (WIN32)
  TARGET_LINK_LIBRARIES(TipsyHelpers ws2_32)  
ENDIF (WIN32)


# For SQLiteReader Plugin, SQLite lib
ADD_LIBRARY(
  SQLitelib
  STATIC
    sqlitelib/sqlite3.c
  )

SET_TARGET_PROPERTIES(SQLitelib PROPERTIES COMPILE_FLAGS "-fPIC")
#For the Grafic reader plugin (and later) integrated tipsy HDF, tipsy binary, and grafic reader I use some of Doug Potter's/Joachim Stadel's fio
# For the Tipsy reader plugin I use some of Doug Potter's Tipsy lib.
# TODO: this is not working, hardcoding library path for now

IF (NOT WIN32)
  ADD_LIBRARY(
    fio
    STATIC
    fio/fio.c
    fio/romberg.c
  )
  SET_TARGET_PROPERTIES(fio PROPERTIES LINKER_LANGUAGE C)
  SET_TARGET_PROPERTIES(fio PROPERTIES COMPILE_FLAGS "-fPIC")
  TARGET_LINK_LIBRARIES(pv_astro fio)
ENDIF (NOT WIN32)


TARGET_LINK_LIBRARIES(
  pv_astro
  TipsyHelpers
  SQLitelib
  #gmp
  )
#TARGET_LINK_LIBRARIES(pv_astro #/Users/corbett/Documents/Projects/pvaddons/ParaViz/ParaViz_src/fio/libFio.so)

IF (NOT WIN32)
  TARGET_LINK_LIBRARIES(pv_astro
    fio
  )
ENDIF (NOT WIN32)
