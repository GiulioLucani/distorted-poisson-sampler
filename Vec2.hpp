#pragma once

template <typename T>
class Vec2 {
public:
    T x;
    T y;
    Vec2(T _x = T(), T _y = T()) : x(_x), y(_y) {}
};
