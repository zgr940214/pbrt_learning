#pragma once
#include "math.hpp"

#include <algorithm>

namespace eric {
class Interval {
    public:
        Float low, high;        
    public:
        Interval(Float v): low(v), high(v) {};
        Interval(Float l, Float h): low(l), high(h) {};

        Interval _intersect(const Interval& i) const {
            return Interval(std::max(low, i.low), std::min(high, i.high));
        }

        bool _degenerated() const {
            return high < low;
        }

        Interval operator+(const Interval& i) {
            return Interval(low + i.low, high + i.high);
        }

        Interval operator-(const Interval& i) {
            return Interval(low - i.high, high - i.low);
        }

        Interval operator*(const Interval& i) {
            return Interval(std::min(std::min(low * i.low, high * i.low),
                    std::min(low * i.high, high * i.high)),
                    std::max(std::max(low * i.low, high * i.low),
                    std::max(low * i.high, high * i.high)));
        }

};

inline Interval Sin(const Interval& i) {// 
    if (i.low > 2 * M_PI) {
        Float nl = i.low - (int)(i.low / 2 * M_PI) * 2 * M_PI;
        Float nh = i.high - (int)(i.high / 2 * M_PI) * 2 * M_PI;
        return Sin(Interval(nl, nh));
    } else if (i.high < 0) {
        Float nl = i.low + ((int)((2 * M_PI - i.low) / 2 * M_PI) + 1) * 2 * M_PI;
        Float nh = i.high + ((int)((2 * M_PI - i.high) / 2 * M_PI) + 1) * 2 * M_PI;
        return Sin(Interval(nl, nh));
    }

    Float sinLow = sin(i.low);
    Float sinHigh = sin(i.high);
    if (sinLow > sinHigh) {
        std::swap(sinLow, sinHigh);
    }
    if (i.low < M_PI / 2 && i.high > M_PI / 2) {
        sinHigh = 1.f;
    }
    if (i.low < 3 / 2 * M_PI && i.high > 3 / 2 * M_PI) {
        sinLow = -1.f;
    }
    return Interval(sinLow, sinHigh);
}

inline Interval Cos(const Interval& i) {// interval i must be within (0, 2pi);
    if (i.low > 2 * M_PI) {
        Float nl = i.low - (int)(i.low / 2 * M_PI) * 2 * M_PI;
        Float nh = i.high - (int)(i.high / 2 * M_PI) * 2 * M_PI;
        return Cos(Interval(nl, nh));
    } else if (i.high < 0) {
        Float nl = i.low + ((int)((2 * M_PI - i.low) / 2 * M_PI) + 1) * 2 * M_PI;
        Float nh = i.high + ((int)((2 * M_PI - i.high) / 2 * M_PI) + 1) * 2 * M_PI;
        return Cos(Interval(nl, nh));
    }

    Float sinLow = cos(i.low);
    Float sinHigh = cos(i.high);
    if (sinLow > sinHigh) {
        std::swap(sinLow, sinHigh);
    }
    if (i.low < 0 && i.high > 0) {
        sinHigh = 1.f;
    }
    if (i.low <  M_PI && i.high > M_PI) {
        sinLow = -1.f;
    }
    return Interval(sinLow, sinHigh);
}

void IntervalFindZeros(Float c1, Float c2, Float c3, Float c4, Float c5, 
    Float theta, Interval tInterval, Float *zeros, int *zeroCount, int depth = 8) {
        Interval range = Interval(c1) +
            (Interval(c2) + Interval(c3) * tInterval) *
            Cos(Interval(2 * theta) * tInterval) +
            (Interval(c4) + Interval(c5) * tInterval) *
            Sin(Interval(2 * theta) * tInterval);
        if (range.low > 0. || range.high < 0. || range.low == range.high)
            return;
        else if (depth > 0) {
            Float mid = (tInterval.low + tInterval.high) * 0.5f;
            IntervalFindZeros(c1, c2, c3, c4, c5, theta,
            Interval(tInterval.low, mid), zeros, zeroCount, depth - 1);
            IntervalFindZeros(c1, c2, c3, c4, c5, theta,
            Interval(mid, tInterval.high), zeros, zeroCount, depth - 1);
        } else {
             // Use Newton's method to refine zero
            Float tNewton = (tInterval.low + tInterval.high) * 0.5f;
            for (int i = 0; i < 4; ++i) {
                Float fNewton =
                    c1 + (c2 + c3 * tNewton) * std::cos(2.f * theta * tNewton) +
                    (c4 + c5 * tNewton) * std::sin(2.f * theta * tNewton);
                Float fPrimeNewton = (c3 + 2 * (c4 + c5 * tNewton) * theta) *
                                        std::cos(2.f * tNewton * theta) +
                                    (c5 - 2 * (c2 + c3 * tNewton) * theta) *
                                        std::sin(2.f * tNewton * theta);
                if (fNewton == 0 || fPrimeNewton == 0) break;
                tNewton = tNewton - fNewton / fPrimeNewton;
            }
            if (tNewton >= tInterval.low - 1e-3f &&
                tNewton < tInterval.high + 1e-3f) {
                zeros[*zeroCount] = tNewton;
                (*zeroCount)++;
            }
        }
    };

} // namespace eric 
