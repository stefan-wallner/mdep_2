cmake_minimum_required (VERSION 2.8.8)

project (Mdep)
set(CMAKE_AUTOMOC ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Werror -Wno-reorder -fopenmp")

set(use_auto_diff ON) 
set(use_YAML ON)
#set(use_multiNest ON)
#set(use_anchor_t ON)

set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR} )

find_package(ROOT 5.34.20 COMPONENTS "Minuit2" REQUIRED)

set( INCLUDE_DIRS
    ${PROJECT_BINARY_DIR} 
    ${PROJECT_SOURCE_DIR} 
    ${ROOT_INCLUDE_DIR}
)

 
#### Switch on single modules

add_definitions( -DSTORE_ACTIVE )

#if ( use_anchor_t )
#  add_definitions( -DUSE_ANCHOR_T )
#endif()

if( use_YAML )
  add_definitions( -DUSE_YAML )
  set( INCLUDE_DIRS
       ${INCLUDE_DIRS}
       /nfs/hicran/project/compass/analysis/fkrinner/yaml/include
  )
endif()


if( use_auto_diff )
  add_definitions( -DADOL_ON )
  set( INCLUDE_DIRS
       ${INCLUDE_DIRS}
      /nfs/hicran/home/fkrinner/private/include
  )
endif()

if( use_multiNest )
  set( INCLUDE_DIRS
      ${INCLUDE_DIRS}
      /nfs/hicran/project/compass/analysis/fkrinner/MultiNest-master/include
  )
endif()


include_directories(
     ${INCLUDE_DIRS}	
)

#set libraries
add_library(waveset               SHARED waveset.cxx            )
add_library(method                SHARED method.cxx             )
add_library(full_SDM              SHARED full_SDM.cxx           )
add_library(anchor_t              SHARED anchor_t.cxx           )
add_library(full_covariance       SHARED full_covariance.cxx    )
add_library(old_method            SHARED old_method.cxx         )
add_library(minimize              SHARED minimize.cxx           )


#set execuatbles
add_executable(fitmd.exe         fitmd.cc                       )
#add_executable(de_isobarring     test_de-isobar.cc              )
add_executable(standard_test.exe standard_test.cc               )
add_executable(amplitude         test_amplitude.cc              )

if ( use_multiNest )
  add_executable(multinest       multinest.cc                   )
endif()


#link own libraries
target_link_libraries(anchor_t   waveset                        )
#target_link_libraries(full_covariance waveset                  )
target_link_libraries(minimize   method                         )
target_link_libraries(minimize   anchor_t                       )
target_link_libraries(minimize   full_covariance                )
target_link_libraries(minimize   old_method                     )
target_link_libraries(minimize   full_SDM                       )
target_link_libraries(fitmd.exe  method                         )
target_link_libraries(fitmd.exe  old_method                     )
target_link_libraries(fitmd.exe  minimize                       )
target_link_libraries(standard_test.exe minimize                )
#target_link_libraries(de_isobarring anchor_t                    )

#link ROOT
target_link_libraries(minimize ${ROOT_LIBS} )

#link YAML
if ( use_YAML )
  target_link_libraries(fitmd.exe /nfs/hicran/project/compass/analysis/fkrinner/yaml/build/libyaml-cpp.so )
  target_link_libraries(standard_test.exe /nfs/hicran/project/compass/analysis/fkrinner/yaml/build/libyaml-cpp.so )
  target_link_libraries(waveset /nfs/hicran/project/compass/analysis/fkrinner/yaml/build/libyaml-cpp.so )
endif()

#link ADOL
if ( use_auto_diff )
  target_link_libraries(minimize  /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
  target_link_libraries(method    /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
  target_link_libraries(anchor_t  /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
  target_link_libraries(minimize  /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
  target_link_libraries(amplitude /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
endif()

if ( use_multiNest )
  target_link_libraries(multinest  /nfs/hicran/project/compass/analysis/fkrinner/MultiNest-master/lib/libmultinest.so )
  target_link_libraries(multinest minimize)
endif()
