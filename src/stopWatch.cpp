#include <iostream>
#include <chrono>
#include <thread>
#include "../include/stopWatch.h"


void Stopwatch::tic() {
    auto now = std::chrono::steady_clock::now();

    if (then != std::chrono::steady_clock::time_point{}) {
        std::cerr << "You have not ended your previous stopwatch." << std::endl;
        std::cout << "Restarting stopwatch now.." << std::endl;
    } else {
        //std::cout << "###" << std::endl;
    }
    then = now;
}

void Stopwatch::toc() {
    auto now = std::chrono::steady_clock::now();

    if (then != std::chrono::steady_clock::time_point{}) {
        auto delta = std::chrono::duration_cast<std::chrono::microseconds>(now - then);
        accumulator += delta;
    } else {
        std::cerr << "You have not started stopwatch" << std::endl;
    }
    then = std::chrono::steady_clock::time_point{};
}

void Stopwatch::reset() {
    then = std::chrono::steady_clock::time_point{};
    accumulator = std::chrono::microseconds(0);
}

std::chrono::microseconds Stopwatch::getElapsed() const {
    return accumulator;
}
