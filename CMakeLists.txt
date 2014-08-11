cmake_minimum_required (VERSION 2.6)
set(CMAKE_CXX_FLAGS "-O3")
project (Mdep)

set(use_auto_diff ON) 

set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR} )
set( INCLUDE_DIRS
    ${PROJECT_BINARY_DIR} 
    ${PROJECT_SOURCE_DIR} 
    /nfs/hicran/project/compass/analysis/fkrinner/eigen-eigen-6b38706d90a9
    /nfs/mnemosyne/sys/slc6/sw/root/x86-64/5.34.18/root/include
    /nfs/hicran/project/compass/analysis/fkrinner/yaml-cpp-0.5.1/include
) 

if( use_auto_diff )
  add_definitions( -DADOL_ON )
  set( INCLUDE_DIRS
       ${INCLUDE_DIRS}
      /nfs/hicran/home/fkrinner/private/include
  )
endif()

include_directories(
     ${INCLUDE_DIRS}	
)

#set libraries
add_library(chi2 SHARED chi2.cxx)
add_library(chi2_2d SHARED chi2_2d.cxx)
add_library(anchor SHARED anchor.cxx)
add_library(anchor_t SHARED anchor_t.cxx)
add_library(minimize SHARED minimize.cxx)


#set execuatbles
add_executable(minimizers minimizers_6waves.cc)
add_executable(fitmd.exe fitmd.cc)

#link own libraries
target_link_libraries(minimize anchor_t)
target_link_libraries(anchor chi2)
target_link_libraries(anchor_t chi2_2d)
target_link_libraries(chi2_2d chi2)
target_link_libraries(minimizers minimize)
target_link_libraries(fitmd.exe minimize)

#link ROOT
target_link_libraries(minimize /nfs/mnemosyne/sys/slc6/sw/root/x86-64/5.34.18/root/lib/libMathCore.so.5.34 )
target_link_libraries(minimize /nfs/mnemosyne/sys/slc6/sw/root/x86-64/5.34.18/root/lib/libMathMore.so.5.34 )

#link YAML
target_link_libraries(fitmd.exe /nfs/hicran/project/compass/analysis/fkrinner/yaml-cpp-0.5.1/build/libyaml-cpp.a )

#link ADOL
if ( use_auto_diff )
  target_link_libraries(chi2 /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
  target_link_libraries(chi2_2d /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
  target_link_libraries(anchor_t /nfs/hicran/home/fkrinner/private/lib/libadolc.so )
endif()
