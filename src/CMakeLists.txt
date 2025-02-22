set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -DLOG_USE_COLOR")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -DLOG_USE_COLOR")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

# Create executable
add_executable(gff_splice_sites_to_vcf gff_splice_sites_to_vcf.c)

find_package(Sanitizers)

# Link external libraries
target_link_libraries(gff_splice_sites_to_vcf ${HTSLIB_LIBRARIES})

# Add logging
add_library(Logging STATIC logging/log.c)
target_compile_definitions(Logging PRIVATE -DLOG_USE_COLOR)
target_include_directories(gff_splice_sites_to_vcf PRIVATE Logging)
target_link_libraries(gff_splice_sites_to_vcf Logging)

# Add BCFTools
add_library(Bcftools STATIC bcftools/regidx.c bcftools/gff.c)
target_link_libraries(Bcftools PRIVATE ${HTSLIB_LIBRARIES})
target_include_directories(gff_splice_sites_to_vcf PRIVATE Bcftools)
target_link_libraries(gff_splice_sites_to_vcf Bcftools)

