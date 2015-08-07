add_executable(elementparametrization elementparametrization.cc)
target_link_libraries(elementparametrization ${DUNE_LIBS})
dune_symlink_to_source_files(FILES "heat.ini")
