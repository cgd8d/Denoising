/*
Create a thread-safe stopwatch -- multiple threads can use this watch concurrently without any sort of locking.
The resulting time will be the sum of the times accumulated by all of the component stopwatches.

Proper use:
void f() {
  static SafeStopwatch my_watch("my stopwatch name");
  SafeStopwatch::tag my_tag = my_watch.Start();
  ...
  my_watch.Stop(my_tag);
}

Upon program exit, the watch will print out its total accumulated time as it is destroyed.
Note that static local initialization is guaranteed thread-safe in C++11; it has also been
thread-safe in gcc for a long time.  Therefore, it is fine for our purposes.
*/
#ifndef SafeStopwatch_hpp
#define SafeStopwatch_hpp

#include <boost/atomic.hpp>
#include <boost/timer/timer.hpp>
#include <iostream>
#include <string>
#include <cassert>

class SafeStopwatch
{
 public:
  typedef boost::timer::cpu_timer tag;

  SafeStopwatch(std::string id) {
    // Constructor -- *NOT* guaranteed thread-safe (but C++11 protects you from this issue).
    assert(fNanoseconds.is_lock_free());
    fNanoseconds.store(0, boost::memory_order_relaxed);
    fID = id;
  }

  ~SafeStopwatch() {Print();}

  tag Start() const {
    // Get a "tag", which identifies the start time.
    // This is what will let the class compute a duration when "stop" is later called.
    // This function is thread-safe.
    return tag();
  }

  void Stop(tag& t) {
    // Increment the stopwatch time by the duration between when tag was created and now.
    // This function is thread-safe.
    t.stop();
    underlying_rep duration = t.elapsed().wall;
    fNanoseconds.fetch_add(duration, boost::memory_order_relaxed);
  }

  void Print() const {
    // Print the time of the stopwatch, terminated by an endline.
    // Any timers which have not been stopped will have no effect on the time printed.
    // Printing from multiple threads may yield garbage,
    // but technically this function is thread-safe.
    std::cout<<fID<<" has accumulated ";
    std::cout<<double(fNanoseconds.load(boost::memory_order_relaxed))/1e9;
    std::cout<<" sec."<<std::endl;
  }

 private:
  typedef boost::timer::nanosecond_type underlying_rep;
  typedef boost::atomic<underlying_rep> rep;
  rep fNanoseconds;
  std::string fID;
};
#endif
