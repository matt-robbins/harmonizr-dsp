#ifndef _window_
#define _window_
#include <vector>

class Window {
public:
    enum winType { Hann=0, Hamm=1 };

    Window(winType type, int table_size);
    ~Window() {};

    float operator[](int ix);
    float value(float position); // compute window value
    void apply(float * data, int N);
    void print();
    
private:
    winType type;
    int N;
    int offset = 2;
    std::vector<float> table;
};

class Window2 {
public:
    enum winType { Hann=0, Hamm=1 };

    Window2(winType type, int table_size);
    Window2(const Window2& other);
    ~Window2();

    float operator[](int ix);
    float value(float position); // compute window value
    void apply(float * data, int N);
    void print();
    
private:
    void calc();
    winType type;
    int N;
    int offset = 2;
    float* table;
};

#endif
