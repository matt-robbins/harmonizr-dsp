#ifndef _window_
#define _window_

class Window {
public:
    enum winType { Hann=0, Hamm=1 };
    Window(winType type, int table_size);
    ~Window();

    float operator[](int ix);
    float value(float position); // compute window value
    void apply(float * data, int N);
    
private:
    winType type;
    int N;
    float * table;
};

#endif
