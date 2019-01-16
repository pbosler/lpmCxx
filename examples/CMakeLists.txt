include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src" "${PROJECT_BINARY_DIR}" ${VTK_INCLUDE_DIRS} ${KokkosInclude})

ADD_EXECUTABLE(swe1d.exe SWEOneD.cpp)
TARGET_LINK_LIBRARIES(swe1d.exe lpm)