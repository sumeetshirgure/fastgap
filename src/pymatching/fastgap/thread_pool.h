// Copyright 2026 fastgap contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PYMATCHING2_FASTGAP_THREAD_POOL_H
#define PYMATCHING2_FASTGAP_THREAD_POOL_H

#include <atomic>
#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

#ifdef __linux__
#include <pthread.h>
#include <sched.h>
#endif

namespace pm {
namespace fastgap {

/// A persistent pool created once per GapDecoder (never per shot). `run(n, fn)` calls
/// fn(worker, task) for task in [0, n), with worker in [0, size()); the calling thread acts as
/// worker 0 and blocks until all tasks are done. Tasks are handed out dynamically.
class ThreadPool {
   public:
    explicit ThreadPool(size_t num_threads, bool pin = false) : num_threads_(num_threads < 1 ? 1 : num_threads) {
        if (pin)
            pin_current_thread(0);
        for (size_t w = 1; w < num_threads_; w++) {
            threads_.emplace_back([this, w, pin] {
                if (pin)
                    pin_current_thread(w);
                worker_loop(w);
            });
        }
    }
    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
            generation_++;
        }
        cv_.notify_all();
        for (auto& t : threads_)
            t.join();
    }
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    size_t size() const {
        return num_threads_;
    }

    void run(size_t num_tasks, const std::function<void(size_t, size_t)>& fn) {
        if (num_threads_ == 1 || num_tasks <= 1 || inside_job()) {
            // Serial fallback; also taken by a nested call from inside a task, since the pool is
            // not re-entrant. A nested caller keeps its own worker slot (0 here is only correct
            // for the outermost caller, so nested callers must not rely on the worker index).
            for (size_t t = 0; t < num_tasks; t++)
                fn(0, t);
            return;
        }
        {
            std::lock_guard<std::mutex> lock(mutex_);
            job_ = &fn;
            num_tasks_ = num_tasks;
            next_task_.store(0);
            pending_workers_ = num_threads_ - 1;
            error_ = nullptr;
            generation_++;
        }
        cv_.notify_all();
        drain(0);
        std::unique_lock<std::mutex> lock(mutex_);
        done_cv_.wait(lock, [this] { return pending_workers_ == 0; });
        job_ = nullptr;
        if (error_)
            std::rethrow_exception(error_);
    }

    /// Best effort: pins on Linux, no-op elsewhere (macOS has no hard affinity API).
    static bool pin_current_thread(size_t cpu) {
#ifdef __linux__
        cpu_set_t set;
        CPU_ZERO(&set);
        CPU_SET(cpu % std::max(1u, std::thread::hardware_concurrency()), &set);
        return pthread_setaffinity_np(pthread_self(), sizeof(set), &set) == 0;
#else
        (void)cpu;
        return false;
#endif
    }

   private:
    static bool& inside_job() {
        thread_local bool flag = false;
        return flag;
    }

    void drain(size_t worker) {
        while (true) {
            size_t t = next_task_.fetch_add(1);
            if (t >= num_tasks_)
                break;
            try {
                inside_job() = true;
                (*job_)(worker, t);
                inside_job() = false;
            } catch (...) {
                inside_job() = false;
                std::lock_guard<std::mutex> lock(mutex_);
                if (!error_)
                    error_ = std::current_exception();
                next_task_.store(num_tasks_);
            }
        }
    }

    void worker_loop(size_t worker) {
        uint64_t seen = 0;
        while (true) {
            {
                std::unique_lock<std::mutex> lock(mutex_);
                cv_.wait(lock, [&] { return generation_ != seen; });
                seen = generation_;
                if (stop_)
                    return;
            }
            drain(worker);
            {
                std::lock_guard<std::mutex> lock(mutex_);
                pending_workers_--;
            }
            done_cv_.notify_one();
        }
    }

    size_t num_threads_;
    std::vector<std::thread> threads_;
    std::mutex mutex_;
    std::condition_variable cv_;
    std::condition_variable done_cv_;
    const std::function<void(size_t, size_t)>* job_ = nullptr;
    size_t num_tasks_ = 0;
    std::atomic<size_t> next_task_{0};
    size_t pending_workers_ = 0;
    uint64_t generation_ = 0;
    bool stop_ = false;
    std::exception_ptr error_;
};

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_THREAD_POOL_H
