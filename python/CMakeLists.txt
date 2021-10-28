add_subdirectory(dune)

configure_file(setup.py.in setup.py)

dune_python_install_package(
  PATH "."
  CMAKE_METADATA_FILE dune/foamgrid/metadata.cmake
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)

if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()
