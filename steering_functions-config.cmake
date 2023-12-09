include(CMakeFindDependencyMacro)

# Same syntax as find_package
find_dependency(Eigen3 REQUIRED)

# Add the targets file
include("${CMAKE_CURRENT_LIST_DIR}/steering_functions-targets.cmake")
