cmake_minimum_required(VERSION 3.0.2)
project(steering_functions)

set(PROJECT_VERSION 0.3.0)

## Build with ROS 1 by default, if it is available
if($ENV{ROS_VERSION} EQUAL 1)
    option(BUILD_WITH_ROS "Build with ROS 1 support" ON)
else()
    option(BUILD_WITH_ROS "Build with ROS 1 support" OFF)
endif()

## System dependencies are found with CMake's conventions
find_package(Eigen3 REQUIRED)

if(BUILD_WITH_ROS)
    find_package(catkin REQUIRED COMPONENTS
        roscpp
        geometry_msgs
        nav_msgs
        visualization_msgs
        costmap_2d
        tf
    )

    catkin_package(
        INCLUDE_DIRS include
        LIBRARIES steering_functions
        DEPENDS EIGEN3
    )
endif()

## Get all source files
file(GLOB SRC_HC_CC_SS src/hc_cc_state_space/*.cpp)

add_library(${PROJECT_NAME}
    src/dubins_state_space/dubins_state_space.cpp
    src/reeds_shepp_state_space/reeds_shepp_state_space.cpp
    src/utilities/utilities.cpp
    src/filter/ekf.cpp
    ${SRC_HC_CC_SS}
)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
)
set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${PROJECT_NAME}-config-version.cmake
    COMPATIBILITY AnyNewerVersion
)
install(DIRECTORY include/${PROJECT_NAME} DESTINATION include)
install(EXPORT ${PROJECT_NAME}-targets
    FILE "${PROJECT_NAME}-targets.cmake"
    NAMESPACE steering_functions::
    DESTINATION "lib/cmake/${PROJECT_NAME}"
)
install(
    FILES
        "${PROJECT_NAME}-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
    DESTINATION "lib/cmake/${PROJECT_NAME}"
)

if(BUILD_WITH_ROS)
    add_executable(${PROJECT_NAME}_node src/${PROJECT_NAME}_node.cpp)

    target_link_libraries(${PROJECT_NAME}_node
        PRIVATE
            steering_functions
            Eigen3::Eigen
            ${catkin_LIBRARIES}
    )

    target_include_directories(${PROJECT_NAME}_node
        PRIVATE
            ${catkin_INCLUDE_DIRS}
    )

    install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_node
        ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
    )

    install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
    )

    install(DIRECTORY
        config
        launch
        DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
    )

    if(CATKIN_ENABLE_TESTING)
        ## Add gtest based cpp test target and link libraries
        catkin_add_gtest(steering_functions_test test/steering_functions_test.cpp)
        if(TARGET steering_functions_test)
            target_link_libraries(steering_functions_test ${catkin_LIBRARIES} ${PROJECT_NAME})
            target_include_directories(steering_functions_test PRIVATE ${catkin_INCLUDE_DIRS})
        endif()

        catkin_add_gtest(fresnel_test test/fresnel_test.cpp)
        if(TARGET fresnel_test)
            target_link_libraries(fresnel_test ${catkin_LIBRARIES} ${PROJECT_NAME})
            target_include_directories(fresnel_test PRIVATE ${catkin_INCLUDE_DIRS})
        endif()

        catkin_add_gtest(timing_test test/timing_test.cpp)
        if(TARGET timing_test)
            target_link_libraries(timing_test ${catkin_LIBRARIES} ${PROJECT_NAME})
            target_include_directories(timing_test PRIVATE ${catkin_INCLUDE_DIRS})
        endif()

        catkin_add_gtest(jacobian_test test/jacobian_test.cpp)
        if(TARGET jacobian_test)
            target_link_libraries(jacobian_test ${catkin_LIBRARIES} ${PROJECT_NAME})
            target_include_directories(jacobian_test PRIVATE ${catkin_INCLUDE_DIRS})
        endif()

        catkin_add_gtest(hc_cc_circle_test test/hc_cc_circle_test.cpp)
        if(TARGET hc_cc_circle_test)
            target_link_libraries(hc_cc_circle_test ${catkin_LIBRARIES} ${PROJECT_NAME})
            target_include_directories(hc_cc_circle_test PRIVATE ${catkin_INCLUDE_DIRS})
        endif()
    endif()
endif()
