#
# Source files
#
# Here is where you can add the name of your local common classes.
#
# set(classes
#   vtkvmtkMath
# )

#
# Create the vtkmyCommon C++ library.
#
# vtk_module_add_module(VMTK::Common
#   CLASSES ${classes})



set(VTK_VMTK_COMMON_TARGET_LINK_LIBRARIES )

set ( VTK_VMTK_COMMON_COMPONENTS
vtkCommonCore
)
if (VTK_WRAP_PYTHON AND VTK_VMTK_WRAP_PYTHON)
list(APPEND VTK_VMTK_COMMON_COMPONENTS
  vtkWrappingPythonCore
  )
endif(VTK_WRAP_PYTHON AND VTK_VMTK_WRAP_PYTHON)
find_package( VTK ${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR} COMPONENTS ${VTK_VMTK_COMMON_COMPONENTS} REQUIRED )
include( ${VTK_USE_FILE} )
list( APPEND VTK_VMTK_COMMON_TARGET_LINK_LIBRARIES ${VTK_LIBRARIES} )

set (VTK_VMTK_COMMON_SRCS
  vtkvmtkMath.cxx
  )

vmtk_build_library(
  NAME Common
  SRCS ${VTK_VMTK_COMMON_SRCS}
  TARGET_LINK_LIBRARIES ${VTK_VMTK_COMMON_TARGET_LINK_LIBRARIES}
  )


# vtkvmtkConfigure.h

configure_file(
  ${VTK_VMTK_SOURCE_DIR}/vtkvmtkConfigure.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/vtkvmtkConfigure.h
  )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/vtkvmtkConfigure.h
  DESTINATION ${VTK_VMTK_INSTALL_INCLUDE_DIR}
  COMPONENT Development
  )
target_include_directories(vtkvmtkCommon PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # For vtkvmtkConfigure.h
  )
