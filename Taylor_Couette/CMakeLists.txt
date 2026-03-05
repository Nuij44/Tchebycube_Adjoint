cmake_minimum_required(VERSION 3.10)
project(TchebyCube)
cmake_policy(SET CMP0074 NEW)

enable_language(Fortran)

set(case_name "optimisation")
set(DepDir $ENV{DepDir})

set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -r8 -check bounds -debug extended")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -r8")
set(CMAKE_BUILD_TYPE DEBUG)  # Spécifier le type de build


# Add m_fourier_tools.f90 as a source file for the executable
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)
set(BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/src)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# external libraries
#find_package(DECOMP2D REQUIRED)
#include_directories(${DECOMP2D_INCLUDE_DIR})

find_package(BLAS REQUIRED)

#find_package(FFTW3 REQUIRED)
#include_directories(${FFTW3_INCLUDE_DIR})

find_package(HDF5 COMPONENTS Fortran)
include_directories(${HDF5_INCLUDE_DIRS})

# Add subdirectory for Tcheb3 library if applicable
add_subdirectory(${SRC_DIR} ${BIN_DIR})

#set(COMMON_SOURCE_FILES m_udf_mms.f90)

set(EXECUTABLES
  DAL_TC.x:DAL_TC.f90
  visu_TC.x:visu_TC.f90
  lin_visu_TC.x:lin_visu_TC.f90
  test_dealiazing.x:test_dealiazing.f90
)

# Add source file paths
#set(COMMON_SOURCE_FILES m_udf_mms.f90)


# Boucle pour créer les exécutables
foreach(EXEC ${EXECUTABLES})
    string(REPLACE ":" ";" EXEC_SPLIT ${EXEC}) # Divise le nom et le fichier source
    list(GET EXEC_SPLIT 0 EXEC_NAME)           # Nom de l'exécutable
    list(GET EXEC_SPLIT 1 SOURCE_FILE)         # Fichier source spécifique

    # Créer l'exécutable
    add_executable(${EXEC_NAME} ${SOURCE_FILE} ${COMMON_SOURCE_FILES})

    # Ajout du répertoire contenant les .mod de 2decomp
    target_include_directories(${EXEC_NAME}
      PRIVATE
      ${DepDir}/2decomp-fft-2.0.4-icc20-mpi/include
      ${DepDir}/fftw-3.3.10-icc20-mpi/include
      ${HDF5_INCLUDE_DIRS}
    )
    
    # Lier les bibliothèques nécessaires
    target_link_libraries(${EXEC_NAME} PRIVATE Tcheb3 BLAS::BLAS HDF5::HDF5 ${DepDir}/2decomp-fft-2.0.4-icc20-mpi/lib64/libdecomp2d.a ${DepDir}/fftw-3.3.10-icc20-mpi/lib64/libfftw3.so)

    # Définir le répertoire de sortie
    set_target_properties(${EXEC_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endforeach()
