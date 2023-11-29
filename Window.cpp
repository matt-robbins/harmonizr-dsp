#include "Window.hpp"
#include "Util.hpp"
#include <math.h>
#include <iostream>

Window::Window(winType type, int table_size) : type{type}, N{table_size} {
    table = new float[N+6]();
    float a0 = 0.5;
    if (type == Window::Hamm) {
        a0 = 0.54;
    }
    
    float a1 = 1.0 - a0;
    
    for (int k = 0; k < N; k++){
        table [3 + k] = a0 - a1 * cosf (2 * M_PI * k / (N-1));
    }
}

Window::~Window() {
    delete(table);
}

float Window::operator[](int ix){
    return table[ix+3];
}

float Window::value(float f) {
    if (f < 0 || f > 1){
        return 0.0;
    }
    
    float wi = 3 + (N-1) * f;
    int i = (int) wi;
    float w = cubic (table + i - 1, wi - i);    
    return w;
}

void Window::apply(float *data, int N) {
    for (int i = 0; i < N; i++) {
        data[i] *= Window::value((float)i / (N-1));
    }
}

