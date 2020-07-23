include(functions.cmake)

set(CMAKE_CXX_STANDARD 11)

# Set a predefined build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release'.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

option(TTK_ENABLE_64BIT_IDS "Enable processing on large datasets" OFF)
mark_as_advanced(TTK_ENABLE_64BIT_IDS)

option(TTK_ENABLE_KAMIKAZE "Enable Kamikaze compilation mode" ON)
mark_as_advanced(TTK_ENABLE_KAMIKAZE)

option(TTK_ENABLE_CPU_OPTIMIZATION "Enable native CPU optimizations" ON)
mark_as_advanced(TTK_ENABLE_CPU_OPTIMIZATION)

option(TTK_BUILD_DOCUMENTATION "Build doxygen developer documentation" OFF)
if(TTK_BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ttk.doxygen
      ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen)
    add_custom_target(doc
      ALL
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/ttk.doxygen
      WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT
        "Generating API documentation with Doxygen"
      VERBATIM
      )
    install(
      DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/doc/html
      DESTINATION
        ${CMAKE_INSTALL_PREFIX}/share/doc/ttk
        )
    install(
      DIRECTORY
        ${CMAKE_SOURCE_DIR}/doc/img
      DESTINATION
        ${CMAKE_INSTALL_PREFIX}/share/doc/ttk
        )
  endif()
endif()

find_package(Boost COMPONENTS system)
if(NOT Boost_FOUND)
  find_package(Boost REQUIRED)
  if(Boost_FOUND)
    message(STATUS "BOOST_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
  endif()
endif()

if (NOT APPLE)
  if (MSVC)
    option(TTK_ENABLE_OPENMP "Enable OpenMP support" FALSE)
  else()
    option(TTK_ENABLE_OPENMP "Enable OpenMP support" TRUE)
  endif()
endif()

if (TTK_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OPENMP_FOUND)
    option(TTK_ENABLE_OMP_PRIORITY "Gives tasks priority, high perf improvement" OFF)

    if (OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 4 AND OpenMP_CXX_VERSION_MINOR GREATER_EQUAL 5)
      set(TTK_ENABLE_OMP_PRIORITY OFF CACHE BOOL "Enable priorities on opnemp tasks" FORCE)
    endif()

    mark_as_advanced(TTK_ENABLE_OMP_PRIORITY)

  endif()
else()
  if(TTK_ENABLE_OMP_PRIORITY)
    # priorities are only meaningful when openmp is on
    set(TTK_ENABLE_OMP_PRIORITY OFF CACHE BOOL "Enable priorities on opnemp tasks" FORCE)
  endif()
endif()


add_subdirectory(base)
if (TTK_BUILD_VTK_WRAPPERS)
  add_subdirectory(vtk)
endif()
