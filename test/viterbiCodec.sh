#!/bin/bash
echo "Google Test started"
# Define your remove command
REMOVE_COMMAND="rm decoder_test"
$REMOVE_COMMAND
# Define your compile command
COMPILE_COMMAND="clang++ -g -o decoder_test viterbiCodec_test.cpp ../src/listDecoder.cpp ../src/feedForwardTrellis.cpp ../src/minHeap.cpp ../src/viterbiCodec.cpp -I /usr/local/Cellar/googletest/1.14.0/include -I ../src/ -lgtest -lgtest_main -pthread -std=c++14"
$COMPILE_COMMAND
# Define your test command
TEST_COMMAND="./decoder_test"
$TEST_COMMAND