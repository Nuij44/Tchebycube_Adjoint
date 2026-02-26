
cmake_minimum_required(VERSION 3.10)
project(TchebyCube)

enable_language(Fortran)

set(case_name "check-poisson")

cmake_policy(SET CMP0074 NEW)

#add_compile_options(-fdefault-real-8 -O2 -ffree-line-length-256)


#set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -r8 -warn unused")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -r8 -check bounds")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -r8")
set(CMAKE_BUILD_TYPE DEBUG)  # Spécifier le type de build


# Add m_fourier_tools.f90 as a source file for the executable
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)
set(BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/src)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../cmake" ${CMAKE_MODULE_PATH})
message(STATUS "DBG ${CMAKE_SOURCE_DIR}/cmake  ")

# external libraries
find_package(FFTW3 REQUIRED)

message(STATUS "FFTW3_INCLUDE_DIR: ${FFTW3_INCLUDE_DIR}")
message(STATUS "FFTW3_LIBRARIES: ${FFTW3_LIBRARIES}")
include_directories(${FFTW3_INCLUDE_DIR})

find_package(DECOMP2D REQUIRED)
include_directories(${DECOMP2D_INCLUDE_DIRS})
message(STATUS "DECOMP2D_INCLUDE_DIRS: ${DECOMP2D_INCLUDE_DIRS}")
message(STATUS "DECOMP2D_LIBRARIES: ${DECOMP2D_LIBRARIES}")



find_package(BLAS REQUIRED)

find_package(HDF5 COMPONENTS Fortran)
include_directories(${HDF5_INCLUDE_DIRS})

# Add subdirectory for Tcheb3 library if applicable
add_subdirectory(${SRC_DIR} ${BIN_DIR})

set(EXECUTABLES
        DAL_TC.x:DAL_TC.f90
	visu_TC.x:visu_TC.f90
)


# Boucle pour créer les exécutables
foreach(EXEC ${EXECUTABLES})
    string(REPLACE ":" ";" EXEC_SPLIT ${EXEC}) # Divise le nom et le fichier source
    list(GET EXEC_SPLIT 0 EXEC_NAME)           # Nom de l'exécutable
    list(GET EXEC_SPLIT 1 SOURCE_FILE)         # Fichier source spécifique

    # Créer l'exécutable
    add_executable(${EXEC_NAME} ${SOURCE_FILE} ${COMMON_SOURCE_FILES})

    # Lier les bibliothèques nécessaires
    target_link_libraries(${EXEC_NAME} PRIVATE Tcheb3 BLAS::BLAS FFTW3::FFTW3 HDF5::HDF5 DECOMP2D::DECOMP2D)

    # Définir le répertoire de sortie
    set_target_properties(${EXEC_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endforeach()
