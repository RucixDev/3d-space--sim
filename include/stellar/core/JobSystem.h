#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace stellar::core {

// A small, dependency-free thread pool / job system.
//
// Design notes:
//  - header-only: keeps the build simple for this repo.
//  - avoids nested-parallelism deadlocks by running one worker loop inline
//    in parallelFor() (so the calling thread always participates).
//  - intended for CPU-bound procedural generation, scanning, planning, etc.
//
// This is deliberately not a fully-featured task graph / work-stealing runtime.
// Keep it small, predictable, and easy to unit test.
class JobSystem {
public:
  // If threadCount == 0, uses a sensible default (hardware_concurrency(), fallback 4).
  explicit JobSystem(std::size_t threadCount = 0) {
    if (threadCount == 0) threadCount = defaultThreadCount();
    if (threadCount < 1) threadCount = 1;
    threadCount_ = threadCount;

    threads_.reserve(threadCount_);
    for (std::size_t i = 0; i < threadCount_; ++i) {
      threads_.emplace_back([this]() { workerLoop(); });
    }
  }

  ~JobSystem() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stopping_ = true;
    }
    cv_.notify_all();
    for (auto& t : threads_) {
      if (t.joinable()) t.join();
    }
  }

  JobSystem(const JobSystem&) = delete;
  JobSystem& operator=(const JobSystem&) = delete;
  JobSystem(JobSystem&&) = delete;
  JobSystem& operator=(JobSystem&&) = delete;

  std::size_t threadCount() const { return threadCount_; }

  static std::size_t defaultThreadCount() {
    const unsigned hc = std::thread::hardware_concurrency();
    if (hc > 0) return static_cast<std::size_t>(hc);
    return 4;
  }

  // Submit a job and get a future.
  template <class F, class... Args>
  auto submit(F&& f, Args&&... args)
    -> std::future<std::invoke_result_t<F, Args...>> {
    using R = std::invoke_result_t<F, Args...>;

    auto task = std::make_shared<std::packaged_task<R()>>(
      [func = std::forward<F>(f), tup = std::make_tuple(std::forward<Args>(args)...)]() mutable -> R {
        return std::apply(std::move(func), std::move(tup));
      });

    std::future<R> fut = task->get_future();

    {
      std::unique_lock<std::mutex> lock(mutex_);
      if (stopping_) {
        // Pool is stopping; run inline to avoid dropping work.
        lock.unlock();
        (*task)();
        return fut;
      }
      queue_.emplace_back([task]() { (*task)(); });
    }
    cv_.notify_one();
    return fut;
  }

  // Block until the queue is empty and no workers are executing tasks.
  void waitIdle() {
    std::unique_lock<std::mutex> lock(mutex_);
    idleCv_.wait(lock, [&]() { return queue_.empty() && active_ == 0; });
  }

  // Parallel for-loop over i in [0, count).
  //
  // grain:
  //   - 0: choose automatically
  //   - >0: process indices in blocks of `grain`
  template <class Fn>
  void parallelFor(std::size_t count, Fn&& fn, std::size_t grain = 0) {
    if (count == 0) return;

    // Fast path: single-thread or tiny work.
    if (threadCount_ <= 1 || count <= 1) {
      for (std::size_t i = 0; i < count; ++i) fn(i);
      return;
    }

    const std::size_t workers = threadCount_;
    const std::size_t desiredTasks = workers * 4;
    if (grain == 0) {
      grain = count / desiredTasks;
      if (grain < 1) grain = 1;
    }

    std::atomic<std::size_t> next{0};

    auto workLoop = [&]() {
      for (;;) {
        const std::size_t start = next.fetch_add(grain, std::memory_order_relaxed);
        if (start >= count) break;
        const std::size_t end = (start + grain < count) ? (start + grain) : count;
        for (std::size_t i = start; i < end; ++i) fn(i);
      }
    };

    // Avoid nested-parallelism deadlocks by ensuring the calling thread also works.
    // Launch (workers - 1) helper tasks and run one workLoop inline.
    std::vector<std::future<void>> futs;
    futs.reserve(workers - 1);
    for (std::size_t i = 0; i + 1 < workers; ++i) {
      futs.push_back(submit(workLoop));
    }

    workLoop();
    for (auto& f : futs) f.get();
  }

private:
  void workerLoop() {
    for (;;) {
      std::function<void()> task;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [&]() { return stopping_ || !queue_.empty(); });
        if (stopping_ && queue_.empty()) return;

        task = std::move(queue_.front());
        queue_.pop_front();
        ++active_;
      }

      task();

      {
        std::lock_guard<std::mutex> lock(mutex_);
        --active_;
        if (queue_.empty() && active_ == 0) idleCv_.notify_all();
      }
    }
  }

  std::size_t threadCount_{1};
  std::vector<std::thread> threads_{};

  std::mutex mutex_{};
  std::condition_variable cv_{};
  std::condition_variable idleCv_{};

  std::deque<std::function<void()>> queue_{};
  std::size_t active_{0};
  bool stopping_{false};
};

} // namespace stellar::core
