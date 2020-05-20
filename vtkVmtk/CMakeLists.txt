project(VTK_VMTK)

include(${VMTK_SOURCE_DIR}/CMake/vmtkFunctions.cmake)

set(dirs
  Utilities
  Common
  ComputationalGeometry
  DifferentialGeometry
  IO
  Misc
  Segmentation
  )

if (VMTK_USE_RENDERING)
  list(APPEND dirs Rendering)
endif (VMTK_USE_RENDERING)

option(VTK_VMTK_CONTRIB "Build and install classes in the vtkVmtk/Contrib directory." OFF)

if (VTK_VMTK_CONTRIB)
  list(APPEND dirs Contrib)
endif (VTK_VMTK_CONTRIB)

include (${VMTK_SOURCE_DIR}/CMake/CMakeOptions.cmake)

if(APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-dylib_file,/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib")
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,-dylib_file,/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib")
endif(APPLE)

option(VTK_VMTK_BUILD_TETGEN "Build TetGen and TetGen wrapper. Check TetGen license before you activate this." ON)

option(VTK_VMTK_BUILD_STREAMTRACER "Build static temporal stream tracer." ON)

if (VTK_USE_COCOA)
  option(VTK_VMTK_USE_COCOA "Build the Cocoa vmtk classes." ON)
endif (VTK_USE_COCOA)

if(NOT DEFINED VTK_VMTK_INSTALL_BIN_DIR)
  set(VTK_VMTK_INSTALL_BIN_DIR bin)
endif(NOT DEFINED VTK_VMTK_INSTALL_BIN_DIR)

if(NOT DEFINED VTK_VMTK_INSTALL_INCLUDE_DIR)
  set(VTK_VMTK_INSTALL_INCLUDE_DIR include/vmtk)
endif(NOT DEFINED VTK_VMTK_INSTALL_INCLUDE_DIR)

if(NOT DEFINED VTK_VMTK_INSTALL_LIB_DIR)
  set(VTK_VMTK_INSTALL_LIB_DIR lib)
endif(NOT DEFINED VTK_VMTK_INSTALL_LIB_DIR)

if(NOT DEFINED VTK_VMTK_MODULE_INSTALL_LIB_DIR)
  set(VTK_VMTK_MODULE_INSTALL_LIB_DIR ${VTK_VMTK_INSTALL_LIB_DIR}/${VMTK_PYTHON_VERSION}/site-packages/vmtk)
endif(NOT DEFINED VTK_VMTK_MODULE_INSTALL_LIB_DIR)

include_directories (${VTK_VMTK_SOURCE_DIR})
include_directories (${VTK_VMTK_SOURCE_DIR}/Common)
include_directories (${VTK_VMTK_SOURCE_DIR}/ComputationalGeometry)
include_directories (${VTK_VMTK_SOURCE_DIR}/DifferentialGeometry)
include_directories (${VTK_VMTK_SOURCE_DIR}/IO)
include_directories (${VTK_VMTK_SOURCE_DIR}/Misc)
if (VMTK_USE_RENDERING)
  include_directories (${VTK_VMTK_SOURCE_DIR}/Rendering)
endif (VMTK_USE_RENDERING)
include_directories (${VTK_VMTK_SOURCE_DIR}/Segmentation)
if (VTK_VMTK_CONTRIB)
  include_directories (${VTK_VMTK_SOURCE_DIR}/Contrib)
endif (VTK_VMTK_CONTRIB)
include_directories (${VTK_VMTK_SOURCE_DIR}/Utilities/vtkvmtkITK)
include_directories (${VTK_VMTK_BINARY_DIR})

set(VTK_DOXYGEN_HOME ${VTK_VMTK_SOURCE_DIR}/Utilities/Doxygen)

if (VTK_WRAP_PYTHON AND VTK_VMTK_WRAP_PYTHON)
  include_directories (${VTK_PYTHON_INCLUDE_DIR})
  configure_file(${VTK_VMTK_SOURCE_DIR}/vtkvmtk.py ${VTK_VMTK_BINARY_DIR}/vtkvmtk.py COPYONLY)
  install(FILES ${VTK_VMTK_BINARY_DIR}/vtkvmtk.py
    DESTINATION ${VTK_VMTK_MODULE_INSTALL_LIB_DIR}
    COMPONENT PythonRuntimeLibraries)
endif (VTK_WRAP_PYTHON AND VTK_VMTK_WRAP_PYTHON)

set(TETGEN_SOURCE_DIR "${VTK_VMTK_SOURCE_DIR}/Utilities/tetgen1.4.3")

foreach(dir ${dirs})
  add_subdirectory(${dir})
endforeach()


if(NOT VTK_VERSION VERSION_LESS 8.90)
    # When building as a stand-alone module for VTK 9
    message(STATUS "vtkvmtkCommontest: Building as an External VTK Module")

    # Clear VTK_LIBS, linking will be handled by modular build
    set(VTK_LIBS)
    set(VTK_Commontest_LINK_TARGET VMTK::Commontest)

    vtk_module_scan(
      MODULE_FILES        "${CMAKE_CURRENT_SOURCE_DIR}/Commontest/Common.module"
      REQUEST_MODULES     "VMTK::Commontest"
      PROVIDES_MODULES    vmtk_modules
      ENABLE_TESTS        "${BUILD_TESTING}")

    vtk_module_build(
      MODULES             ${vmtk_modules}
      INSTALL_EXPORT      Commontest
      ARCHIVE_DESTINATION "lib"
      HEADERS_DESTINATION "include/Commontest${DICOM_SHORT_VERSION}"
      CMAKE_DESTINATION   "${CMAKE_INSTALL_LIBDIR}/cmake/vtk"
      LICENSE_DESTINATION "${CMAKE_INSTALL_LIBDIR}/vtk"
      HIERARCHY_DESTINATION "${CMAKE_INSTALL_LIBDIR}/vtk/hierarchy/${CMAKE_PROJECT_NAME}"
      LIBRARY_NAME_SUFFIX "${Commontest_CUSTOM_LIBRARY_SUFFIX}"
      VERSION             "${Commontest_VERSION}"
      SOVERSION           "1"
      TEST_DATA_TARGET    CommontestData
      USE_EXTERNAL        "${Commontest_USE_EXTERNAL}"
      TEST_INPUT_DATA_DIRECTORY   "${Commontest_test_data_directory_input}"
      TEST_OUTPUT_DATA_DIRECTORY  "${Commontest_test_data_directory_output}")

    if(BUILD_PYTHON_WRAPPERS)
      find_package(PythonInterp ${VTK_PYTHON_VERSION} QUIET)

      vtk_module_wrap_python(
        MODULES         ${vmtk_modules}
        INSTALL_EXPORT  CommonPython
        PYTHON_PACKAGE  "vmtk.modules"
        #MODULE_DESTINATION  ""
        CMAKE_DESTINATION   "${CMAKE_INSTALL_LIBDIR}/cmake/vtk"
        LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        BUILD_STATIC    OFF)
    endif()
  endif() #(NOT VTK_VERSION VERSION_LESS 8.90)