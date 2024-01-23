#pragma once
#include <iostream>
#include <chrono>
#include <thread>

class Stopwatch {
private:
    std::chrono::steady_clock::time_point then{};
    std::chrono::microseconds accumulator{};

public:
    void tic();
    void toc();
    void reset();
    std::chrono::microseconds getElapsed() const;
};