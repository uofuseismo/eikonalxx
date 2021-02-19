#ifndef PRIVATE_TIMER_HPP
#define PRIVATE_TIMER_HPP
#include <chrono>
#include <ctime>
namespace
{
class Timer
{
public:
    /// C'tor
    Timer() :
        mStartTime(std::chrono::high_resolution_clock::now()),
        mEndTime(std::chrono::high_resolution_clock::now()),
        mDuration(0),
        mHaveDuration(false)
    {
    }
    /// Starts the timer
    void start()
    {
        mHaveDuration = false;
        mStartTime = std::chrono::high_resolution_clock::now();
        mEndTime = mStartTime;
    }
    /// Stops the clock
    void end()
    {
        mHaveDuration = false;
        mEndTime = std::chrono::high_resolution_clock::now();
    }
    /// Gets duration between start and end times in seconds
    double getDuration() const
    {
        if (!mHaveDuration)
        {
            auto timeSpan
                = std::chrono::duration_cast<std::chrono::microseconds>
                  (mEndTime - mStartTime);
            mDuration = timeSpan.count()*1.e-6;
            mHaveDuration = true;
        }
        return mDuration;
    }
private:
    std::chrono::high_resolution_clock::time_point mStartTime;
    std::chrono::high_resolution_clock::time_point mEndTime;
    mutable double mDuration;
    mutable bool mHaveDuration;
};
}
#endif
