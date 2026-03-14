/******************************************
Copyright (C) 2023 Authors of GANAK, see AUTHORS file

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#pragma once

#include <atomic>
#include <chrono>
#include <memory>
#include <thread>

class Timer {
    std::atomic<bool> finished = false;
    std::unique_ptr<std::thread> tp;
public:
    void set_timeout(std::atomic<bool>& appmc_timeout_fired, double delay) {
      tp = std::make_unique<std::thread>([this, &appmc_timeout_fired, delay]() {
          auto start = std::chrono::steady_clock::now();
          auto end = start + std::chrono::milliseconds(static_cast<int>(delay * 1000.0));

          // Sleep in smaller intervals to respond to stop() faster
          while(std::chrono::steady_clock::now() < end) {
              if (finished.load()) return;
              std::this_thread::sleep_for(std::chrono::milliseconds(100));
          }
          if (!finished.load()) {
            finished.store(true);
            appmc_timeout_fired = true;
            /* cout << "[appmc] Timeout fired, stopping ganak." << endl; */
          } else {
            /* cout << "[appmc] Timeout fired, but ganak already finished." << endl; */
          }
      });
    }

    void wait_all() {
        finished.store(true);
        if(tp && tp->joinable()) {
            tp->join();
            tp.reset();
        }
    }

    ~Timer() { wait_all(); }
};
