cmake_minimum_required(VERSION 3.10)
project(TchebyCube)
cmake_policy(SET CMP0074 NEW)
enable_language(Fortran)

set(case_name "optimisation")
set(DepDir $ENV{DepDir})

set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -r8 -check bounds -debug extended")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -r8")
set(CMAKE_BUILD_TYPE Release)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)
set(BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/src)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Chemins des dépendances
set(DECOMP2D_INCLUDE_DIR ${DepDir}/lib/intel/decomp2d/include)
set(DECOMP2D_LIB         ${DepDir}/lib/intel/decomp2d/lib64/libdecomp2d.a)
set(FFTW3_INCLUDE_DIR    ${DepDir}/lib/intel/fftw3/include)
set(FFTW3_LIB            ${DepDir}/lib/intel/fftw3/lib64/libfftw3.so)  # adapte lib ou lib64

# Includes globaux — s'appliquent à Tcheb3 ET aux exécutables
include_directories(
  ${DECOMP2D_INCLUDE_DIR}
  ${FFTW3_INCLUDE_DIR}
)

find_package(BLAS REQUIRED)
find_package(HDF5 COMPONENTS Fortran REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

find_package(MPI REQUIRED)
include_directories(${MPI_Fortran_INCLUDE_PATH})

add_subdirectory(${SRC_DIR} ${BIN_DIR})

set(EXECUTABLES
    kerswell.x:kerswell.f90
    Q1_kerswell.x:Q1_kerswell.f90
    post_kers.x:post_kers.f90
)

set(COMMON_SOURCE_FILES m_udf_mms.f90)

foreach(EXEC ${EXECUTABLES})
    string(REPLACE ":" ";" EXEC_SPLIT ${EXEC})
    list(GET EXEC_SPLIT 0 EXEC_NAME)
    list(GET EXEC_SPLIT 1 SOURCE_FILE)

    add_executable(${EXEC_NAME} ${SOURCE_FILE} ${COMMON_SOURCE_FILES})

    target_link_libraries(${EXEC_NAME} PRIVATE
        Tcheb3
        BLAS::BLAS
        HDF5::HDF5
        MPI::MPI_Fortran
        ${DECOMP2D_LIB}
        ${FFTW3_LIB}
    )

    set_target_properties(${EXEC_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endforeach()
