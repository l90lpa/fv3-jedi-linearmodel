
# Grab Python, 3.10 or newer
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module NumPy)
set(pyver python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR})

# Grab the variables from a local Python installation F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Print out the discovered paths
include(CMakePrintHelpers)
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)


#### Arguments for generating wrappers and python extension module

set(package_name "pyfv3jedilm")
set(src_files
    "${CMAKE_CURRENT_SOURCE_DIR}/runtime_mod.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/fv3jedi_lm_mod.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/utils/fv3jedi_lm_kinds_mod.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/utils/fv3jedi_lm_utils_mod.F90"
)

##### Generation target and associated custom command
set(f90wrap_module_py "${package_name}.py")
set(f2py_module_c "_${package_name}module.c")

set(generated_python_files
    "${CMAKE_CURRENT_BINARY_DIR}/${f90wrap_module_py}"
)
set(generated_src_files 
    "${CMAKE_CURRENT_BINARY_DIR}/f90wrap_runtime_mod.f90"
    "${CMAKE_CURRENT_BINARY_DIR}/f90wrap_fv3jedi_lm_mod.f90"
    "${CMAKE_CURRENT_BINARY_DIR}/f90wrap_fv3jedi_lm_kinds_mod.f90"
    "${CMAKE_CURRENT_BINARY_DIR}/f90wrap_fv3jedi_lm_utils_mod.f90"
    "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
)

add_custom_target(
  "_${package_name}_generator"
  DEPENDS ${generated_src_files} ${generated_python_files}
)
add_custom_command(
  OUTPUT ${generated_src_files} ${generated_python_files}
  COMMAND f90wrap -m ${package_name} --f90-mod-name "${INSTALL_LIB_DIR}._${package_name}" -k "${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap_default_real_is_double" ${src_files}
  COMMAND f2py-f90wrap -m "_${package_name}" "${CMAKE_CURRENT_BINARY_DIR}/f90wrap_runtime_mod.f90" "${CMAKE_CURRENT_BINARY_DIR}/f90wrap_fv3jedi_lm_mod.f90" "${CMAKE_CURRENT_BINARY_DIR}/f90wrap_fv3jedi_lm_kinds_mod.f90" "${CMAKE_CURRENT_BINARY_DIR}/f90wrap_fv3jedi_lm_utils_mod.f90"
  DEPENDS ${src_files}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Define the extension module
Python_add_library("_${package_name}" MODULE WITH_SOABI
  "${CMAKE_CURRENT_SOURCE_DIR}/runtime_mod.F90"
  ${generated_src_files}                  # <-- Generated source
  "${F2PY_INCLUDE_DIR}/fortranobject.c"   # <-- From NumPy
)
target_link_libraries("_${package_name}" PRIVATE fv3jedilm)
if(FV3_FORECAST_MODEL MATCHES "UFS")
  target_link_libraries(fv3jedilm PUBLIC FMS::fms_r8)
elseif(NOT FV3_FORECAST_MODEL MATCHES "GEOS")
  target_link_libraries("_${package_name}" PRIVATE fms)
endif()
target_link_libraries("_${package_name}" PRIVATE Python::NumPy)
target_include_directories("_${package_name}" PRIVATE "${F2PY_INCLUDE_DIR}")
target_include_directories("_${package_name}" PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
set_target_properties("_${package_name}" PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${INSTALL_LIB_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${INSTALL_LIB_DIR}"
)

# Add dependency from extension module to generation target
add_dependencies("_${package_name}" "_${package_name}_generator")

# Copy Python wrapper interface to library directory once the extension module has been built
set(files_to_copy ${generated_python_files})
list(APPEND files_to_copy "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")
add_custom_command(TARGET "_${package_name}" POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy ${files_to_copy} "${CMAKE_BINARY_DIR}/${INSTALL_LIB_DIR}")

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${CMAKE_CURRENT_BINARY_DIR}/${package_name}.py "${CMAKE_BINARY_DIR}/${INSTALL_LIB_DIR}/_${package_name}.${Python_SOABI}${CMAKE_SHARED_LIBRARY_SUFFIX}"
	      DESTINATION ${INSTALL_LIB_DIR}
        COMPONENT Python)