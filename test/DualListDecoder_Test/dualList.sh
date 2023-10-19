#!/bin/bash
echo "Google Test started"
# Define your remove command
REMOVE_COMMAND="rm dualList_Test"
$REMOVE_COMMAND
# Define your compile command
COMPILE_COMMAND="clang++ -g -o dualList_Test dualList_test.cpp ../../src/feedForwardTrellis.cpp ../../src/minHeap.cpp ../../src/dualListDecoder.cpp ../../src/viterbiCodec.cpp ../../src/listDecoder.cpp -I /usr/local/Cellar/googletest/1.14.0/include -I ../../src/ -lgtest -lgtest_main -pthread -std=c++14"
$COMPILE_COMMAND
# Define your test command
TEST_COMMAND="./dualList_Test"
$TEST_COMMAND