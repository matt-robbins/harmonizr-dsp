#include "Window.hpp"
#include "Util.hpp"
#include <math.h>
#include <iostream>
#include <vector>

Window::Window(winType type, int table_size) : type{type}, N{table_size} {
    table.resize(N+2*offset);
    float a0 = (type == Window::Hann) ? 0.5 : 0.54;
    float a1 = 1.0 - a0;
    
    for (int k = 0; k < N; k++){
        table [offset + k] = a0 - a1 * cosf (2 * M_PI * k / (N-1));
    }
}

float Window::operator[](int ix){
    return table[ix+offset];
}

void Window::print() {
    for (float v : table) {
        std::cout << v << std::endl;
    }
}

float Window::value(float f) {
    if (f < 0 || f > 1){
        return 0.0;
    }
    
    float wi = offset + (N-1) * f;
    int i = (int) wi;
    float w = cubic (table.begin() + i - 1, wi - i);    
    return w;
}

void Window::apply(float *data, int N) {
    for (int i = 0; i < N; i++) {
        data[i] *= Window::value((float)i / (N-1));
    }
}

// raw pointer implementation

void Window2::calc() {
    float a0 = (type == Window2::Hann) ? 0.5 : 0.54;
    float a1 = 1.0 - a0;
    
    for (int k = 0; k < N; k++){
        table [offset + k] = a0 - a1 * cosf (2 * M_PI * k / (N-1));
    }
}

Window2::Window2(winType type, int table_size=64) : type{type}, N{table_size} {
    table = new float[N+2*offset];
    calc();
}

Window2::Window2(const Window2& other) : type{other.type}, N{other.N} {
    table = new float[N+2*offset];
    calc();
}

Window2::~Window2() {
    //delete[] table;
}

float Window2::operator[](int ix){
    return table[ix+offset];
}

void Window2::print() {
    for (int i = 0; i < N; i++) {
        std::cout << table[i] << std::endl;
    }
}

float Window2::value(float f) {
    if (f < 0 || f > 1){
        return 0.0;
    }
    
    float wi = offset + (N-1) * f;
    int i = (int) floorf(wi);

    float w = cubic (table + i - 1, wi - i);    
    return w;
}

void Window2::apply(float *data, int N) {
    for (int i = 0; i < N; i++) {
        data[i] *= Window2::value((float)i / (N-1));
    }
}
