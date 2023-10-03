#!/bin/bash
echo "Google Test started"
# Define your remove command
REMOVE_COMMAND="rm duoCodec_test"
$REMOVE_COMMAND
# Define your compile command
COMPILE_COMMAND="clang++ -g -o duoCodec_test duoCodec_test.cpp ../../src/feedForwardTrellis.cpp ../../src/minHeap.cpp ../../src/viterbiCodec.cpp ../../src/listDecoder.cpp -I /usr/local/Cellar/googletest/1.14.0/include -I ../../src/ -lgtest -lgtest_main -pthread -std=c++14"
$COMPILE_COMMAND
# Define your test command
TEST_COMMAND="./duoCodec_test"
$TEST_COMMAND