foreach(component IN ITEMS core common arrays interpolate integrate io linalg
        optimize polynomial stats)
    add_subdirectory(${component})
    if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${component}/include")
        target_include_directories(numfort PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/${component}/include")
    endif()
endforeach()
