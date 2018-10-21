#ifndef ASSIGN1CPP_COLOR_H
#define ASSIGN1CPP_COLOR_H

class Color {
public:
    float red;
    float green;
    float blue;
    Color() {
        red = 1.0f;
        green = 1.0f;
        blue = 1.0f;
    }

    Color(float red, float green, float blue): red(red), green(green), blue(blue) {}
};

#endif //ASSIGN1CPP_COLOR_H
