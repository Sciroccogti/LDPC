cmake_minimum_required(VERSION 3.10)
project(LDPC)

set(SOURCES LDPC/GF.cpp LDPC/GFLinAlg.cpp LDPC/NBLdpcBuilder.cpp
            LDPC/NBLdpcCodec.cpp LDPC/NBLdpcDecoder.cpp LDPC/testbed.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/Include
                                                   ${PROJECT_SOURCE_DIR}/LDPC)
