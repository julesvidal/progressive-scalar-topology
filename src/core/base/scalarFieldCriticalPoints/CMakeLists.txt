ttk_add_base_library(scalarFieldCriticalPoints
  SOURCES
    ScalarFieldCriticalPoints.cpp
  HEADERS
    ScalarFieldCriticalPoints.h
  LINK
    triangulation
    unionFind
    multiresTriangulation
    dynamicTree
    )

option(TTK_ENABLE_DYNAMIC_TREES "Enabled dynamic trees" ON)
if (TTK_ENABLE_DYNAMIC_TREES)
  target_compile_definitions(scalarFieldCriticalPoints PUBLIC TTK_ENABLE_DYNAMIC_TREES)
endif()

option(TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC "Enable implicit link" ON)
if (TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC)
  target_compile_definitions(scalarFieldCriticalPoints PUBLIC TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC)
endif()

option(TTK_USE_PRECALCULATED_LINKS_FOR_UF "Use the precalculated for the union finds" ON)
if (TTK_USE_PRECALCULATED_LINKS_FOR_UF)
  target_compile_definitions(scalarFieldCriticalPoints PUBLIC TTK_USE_PRECALCULATED_LINKS_FOR_UF)
endif()
