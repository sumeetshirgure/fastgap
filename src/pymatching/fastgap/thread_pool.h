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

#include <algorithm>
#include <atomic>
#include <chrono>
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
#ifdef __APPLE__
#include <pthread.h>
#include <pthread/qos.h>
#endif
#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
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

/// Hint to the CPU that the caller is spinning.
inline void cpu_relax() {
#if defined(__x86_64__) || defined(_M_X64)
    _mm_pause();
#elif defined(__aarch64__)
    asm volatile("yield" ::: "memory");
#endif
}

/// A team of threads that work together on *one* shot (intra-shot parallelism), as opposed to
/// ThreadPool, which hands out independent tasks. `run(active, fn)` calls fn(t) concurrently on
/// threads t in [0, active) (t = 0 is the caller) and returns once all of them have returned, so
/// fn may synchronise its threads with spin barriers.
///
/// Waking a sleeping thread costs several microseconds, more than the whole search on a typical
/// shot, so idle workers spin for `SPIN_BEFORE_SLEEP` after their last job before they block.
/// Back-to-back shots (a decode_batch, or a real-time stream) therefore never pay a wake-up.
class SpinTeam {
   public:
    static constexpr std::chrono::microseconds SPIN_BEFORE_SLEEP{2000};

    explicit SpinTeam(size_t num_threads, bool pin = false) : num_threads_(num_threads < 1 ? 1 : num_threads) {
        if (pin)
            ThreadPool::pin_current_thread(0);
        for (size_t t = 1; t < num_threads_; t++) {
            threads_.emplace_back([this, t, pin] {
                if (pin)
                    ThreadPool::pin_current_thread(t);
#ifdef __APPLE__
                // Ask for performance cores; no hard affinity exists on macOS.
                pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
                worker_loop(t);
            });
        }
    }
    ~SpinTeam() {
        stop_.store(true);
        job_seq_.fetch_add(1);
        {
            std::lock_guard<std::mutex> lock(mutex_);
        }
        cv_.notify_all();
        for (auto& t : threads_)
            t.join();
    }
    SpinTeam(const SpinTeam&) = delete;
    SpinTeam& operator=(const SpinTeam&) = delete;

    size_t size() const {
        return num_threads_;
    }

    /// fn must not throw.
    void run(size_t active, const std::function<void(size_t)>& fn) {
        active = std::min(std::max<size_t>(active, 1), num_threads_);
        if (active == 1) {
            fn(0);
            return;
        }
        job_ = &fn;
        active_ = active;
        remaining_.store(active - 1, std::memory_order_relaxed);
        job_seq_.fetch_add(1);  // seq_cst: pairs with the sleepers_ check in worker_loop
        if (sleepers_.load() > 0) {
            std::lock_guard<std::mutex> lock(mutex_);
            cv_.notify_all();
        }
        fn(0);
        while (remaining_.load(std::memory_order_acquire) != 0)
            cpu_relax();
    }

   private:
    void worker_loop(size_t t) {
        uint64_t seen = 0;
        while (true) {
            uint64_t now = wait_for_job(seen);
            seen = now;
            if (stop_.load())
                return;
            if (t < active_) {
                (*job_)(t);
                remaining_.fetch_sub(1, std::memory_order_release);
            }
        }
    }

    uint64_t wait_for_job(uint64_t seen) {
        auto deadline = std::chrono::steady_clock::now() + SPIN_BEFORE_SLEEP;
        while (true) {
            for (int i = 0; i < 1024; i++) {
                uint64_t s = job_seq_.load(std::memory_order_acquire);
                if (s != seen)
                    return s;
                cpu_relax();
            }
            if (std::chrono::steady_clock::now() >= deadline)
                break;
        }
        std::unique_lock<std::mutex> lock(mutex_);
        sleepers_.fetch_add(1);
        uint64_t s;
        while ((s = job_seq_.load()) == seen)
            cv_.wait(lock);
        sleepers_.fetch_sub(1);
        return s;
    }

    size_t num_threads_;
    std::vector<std::thread> threads_;
    const std::function<void(size_t)>* job_ = nullptr;
    size_t active_ = 0;
    std::atomic<uint64_t> job_seq_{0};
    std::atomic<size_t> remaining_{0};
    std::atomic<int> sleepers_{0};
    std::atomic<bool> stop_{false};
    std::mutex mutex_;
    std::condition_variable cv_;
};

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_THREAD_POOL_H
