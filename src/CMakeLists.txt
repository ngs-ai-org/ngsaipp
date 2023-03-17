# compiler options
add_compile_options(-std=c++2a)
add_compile_options(-O3)
add_compile_options(-Wall)
add_compile_options(-Wextra)
add_compile_options(-Werror)
add_compile_options(-Wfatal-errors)
add_compile_options(-pedantic)
add_compile_options(-pg)

# include file locations
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories(${INCLUDE_DIRECTORIES})
include_directories(${Boost_INCLUDE_DIRS})

# linking to external lib
# Boost, pbcopper, pbbam, gtest are there
link_directories(${LINK_DIRECTORY})

# list of src files for lib
set(FILES_CPP_LIB 
        "algorithms/segmentation.cpp"
        "dna/dna_utility.cpp"
        "dna/SequenceEnumerator.cpp"
        "genome/constants.cpp"
        "genome/GenomeRegion.cpp"
        "genome/CpGRegion.cpp"
        "epigenetics/CcsKineticExtractor.cpp"
        "epigenetics/KmerData.cpp"
        "epigenetics/KmerMap.cpp"
        "epigenetics/KineticSignal.cpp"
        "epigenetics/KineticModel.cpp"
        "epigenetics/RawKineticModel.cpp"
        "epigenetics/NormalizedKineticModel.cpp"
        "epigenetics/PairWiseKineticModel.cpp"
        "epigenetics/PairWiseNormalizedKineticModel.cpp"
        "epigenetics/DiPositionKineticModel.cpp"
        "epigenetics/DiPositionNormalizedKineticModel.cpp"
        "epigenetics/KineticClassifier.cpp"
        "epigenetics/model_utility.cpp"
        "io/BedRecord.cpp"
        "io/BedReader.cpp"
        "parallel/ThreadPool.cpp"
        "utility/string_utility.cpp")
# list of src files for unittests
set(FILES_CPP_TESTS
        "unittests.cpp"
        "unittests/genome/constants_unittests.cpp"
        "unittests/genome/GenomeRegion_unittests.cpp"
        "unittests/genome/CpGRegion_unittests.cpp"
        "unittests/io/BedRecord_unittests.cpp"
        "unittests/io/BedReader_unittests.cpp"
        "unittests/dna/dna_utility_unittests.cpp"
        "unittests/dna/SequenceEnumerator_unittests.cpp"
        "unittests/epigenetics/utility.cpp"
        "unittests/epigenetics/CcsKineticExtractor_unittests.cpp"
        "unittests/epigenetics/KmerData_unittests.cpp"
        "unittests/epigenetics/KmerMap_unittests.cpp"
        "unittests/epigenetics/KineticSignal_unittests.cpp"
        "unittests/epigenetics/RawKineticModel_unittests.cpp"
        "unittests/epigenetics/NormalizedKineticModel_unittests.cpp"
        "unittests/epigenetics/PairWiseKineticModel_unittests.cpp"
        "unittests/epigenetics/PairWiseNormalizedKineticModel_unittests.cpp"
        "unittests/epigenetics/DiPositionKineticModel_unittests.cpp"
        "unittests/epigenetics/DiPositionNormalizedKineticModel_unittests.cpp")


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${ngsaipp_lib_SOURCE_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${ngsaipp_lib_SOURCE_DIR}/lib")


# compile libraries
## static library
add_library(lib_static STATIC ${FILES_CPP_LIB})
set_target_properties(lib_static PROPERTIES OUTPUT_NAME ngsaipp)
target_link_libraries(lib_static pthread
                                 pbbam
                                 pbcopper
                                 boost_serialization)
## shared library
add_library(lib_shared SHARED ${FILES_CPP_LIB})
set_target_properties(lib_shared PROPERTIES OUTPUT_NAME ngsaipp)
target_link_libraries(lib_shared pthread
                                 pbbam
                                 pbcopper
                                 boost_serialization)

# compile executables
## unittests
set(EXE_UNITTEST "ngsaipp_unittests")
add_executable(${EXE_UNITTEST} ${FILES_CPP_TESTS})
target_link_libraries(${EXE_UNITTEST} lib_static
                                      gtest
                                      pthread)
set_target_properties(${EXE_UNITTEST} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${ngsaipp_lib_SOURCE_DIR}/bin")


# install
install(TARGETS lib_static lib_shared
        LIBRARY DESTINATION ${INSTALL_LIB_DIRECTORY}
        ARCHIVE DESTINATION ${INSTALL_LIB_DIRECTORY})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${INSTALL_HEADER_DIRECTORY}
        FILES_MATCHING PATTERN "*.hpp")