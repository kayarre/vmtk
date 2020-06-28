#
# Build the documentation
#
#include (${CMAKE_ROOT}/Modules/Documentation.cmake OPTIONAL)

# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
vmtkDocumentation
-------------

This module provides support for the VTK documentation framework.  It
relies on several tools (Doxygen, Perl, etc).
#]=======================================================================]

#
# Build the documentation ?
#
option(BUILD_DOCUMENTATION "Build the documentation (Doxygen)." OFF)

if (BUILD_DOCUMENTATION)

  #
  # Check for the tools
  #
  find_package(UnixCommands)
  find_package(Doxygen)
  find_package(Gnuplot)
  find_package(HTMLHelp)
  find_package(Perl)
  find_package(Wget)

  option(DOCUMENTATION_HTML_HELP
    "Build the HTML Help file (CHM)." OFF)

  option(DOCUMENTATION_HTML_TARZ
    "Build a compressed tar archive of the HTML doc." OFF)

  mark_as_advanced(
    DOCUMENTATION_HTML_HELP
    DOCUMENTATION_HTML_TARZ
    )

  #
  # The documentation process is controlled by a batch file.
  # We will probably need bash to create the custom target
  #

endif ()

if (BUILD_DOCUMENTATION)
                                                                                                               
  #
  # Configure the script and the doxyfile, then add target
  #
  configure_file(
    ${VTK_VMTK_SOURCE_DIR}/Utilities/Doxygen/doxyfile.in
    ${VTK_VMTK_BINARY_DIR}/Utilities/Doxygen/doxyfile)
                                                                                                               
  configure_file(
    ${VTK_VMTK_SOURCE_DIR}/Utilities/Doxygen/doc_makeall.sh.in
    ${VTK_VMTK_BINARY_DIR}/Utilities/Doxygen/doc_makeall.sh)

  add_custom_target(DoxygenDoc
    ${BASH}
    ${VTK_VMTK_BINARY_DIR}/Utilities/Doxygen/doc_makeall.sh)
                                                                                                               
  include (${VTK_VMTK_SOURCE_DIR}/Utilities/Doxygen/LocalUserOptions.cmake OPTIONAL)
  
  install_files(${VTK_VMTK_INSTALL_LIB_DIR}/doxygen "\\.css$")
  install_files(${VTK_VMTK_INSTALL_LIB_DIR}/doxygen "\\.gif$")
  install_files(${VTK_VMTK_INSTALL_LIB_DIR}/doxygen "\\.html$")
  install_files(${VTK_VMTK_INSTALL_LIB_DIR}/doxygen "\\.png$")
  #install_files(${VTK_VMTK_INSTALL_LIB_DIR}/doxygen "\\.pl$")
  #install_files(${VTK_VMTK_INSTALL_LIB_DIR}/doxygen "\\.stop$")
  #install_files(${VTK_VMTK_INSTALL_LIB_DIR}/doxygen .txt )
  #install_files(${VTK_VMTK_INSTALL_LIB_DIR}/doxygen .txt authors doc_readme)

endif (BUILD_DOCUMENTATION)
                                                                                                               
