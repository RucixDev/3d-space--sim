#pragma once

#include <atomic>
#include <chrono>
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

  // True if the calling thread is one of this pool's worker threads.
  // Useful for avoiding oversubscription in nested parallel regions.
  bool isWorkerThread() const { return tlsCurrent_ == this; }

  static std::size_t defaultThreadCount() {
    const unsigned hc = std::thread::hardware_concurrency();
    if (hc > 0) return static_cast<std::size_t>(hc);
    return 4;
  }

  // Schedule a job without creating a future.
  //
  // This is a lower-overhead "fire-and-forget" path that is useful for internal
  // helpers like TaskGroup and parallelFor().
  //
  // NOTE: The callable must be safe to execute on any worker thread. If the pool
  // is stopping, the job is executed inline on the calling thread.
  template <class F>
  void enqueue(F&& f) {
    using FnT = std::decay_t<F>;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      if (stopping_) {
        lock.unlock();
        f();
        return;
      }

      if constexpr (std::is_copy_constructible_v<FnT>) {
        queue_.emplace_back(std::forward<F>(f));
      } else {
        // std::function requires copyable callables. For move-only tasks,
        // store the functor behind a shared_ptr.
        auto fn = std::make_shared<FnT>(std::forward<F>(f));
        queue_.emplace_back([fn]() { (*fn)(); });
      }
    }
    cv_.notify_one();
  }

  // Attempt to run exactly one queued task inline.
  // Returns true if a task was executed.
  //
  // This "help while waiting" primitive is important to avoid nested-parallelism
  // deadlocks when workers need to wait for other work in the same pool.
  bool tryRunOne() {
    std::function<void()> task;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      if (queue_.empty()) return false;
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

    return true;
  }

  // A small helper for "submit N jobs, then wait" patterns without futures.
  //
  // TaskGroup is safe to destroy without calling wait(): tasks will still run.
  // However, waiting is typically what you want for deterministic batching.
  class TaskGroup {
  public:
    explicit TaskGroup(JobSystem& js) : st_(std::make_shared<State>(&js)) {}

    TaskGroup(const TaskGroup&) = default;
    TaskGroup& operator=(const TaskGroup&) = default;

    // Submit a task into this group.
    template <class F>
    void run(F&& f) {
      auto st = st_;
      st->remaining.fetch_add(1, std::memory_order_relaxed);

      using FnT = std::decay_t<F>;
      if constexpr (std::is_copy_constructible_v<FnT>) {
        st->js->enqueue([st, func = std::forward<F>(f)]() mutable {
          func();
          if (st->remaining.fetch_sub(1, std::memory_order_acq_rel) == 1) {
            std::lock_guard<std::mutex> lock(st->mutex);
            st->cv.notify_all();
          }
        });
      } else {
        // std::function requires copyable callables; store move-only functors behind a shared_ptr.
        auto fn = std::make_shared<FnT>(std::forward<F>(f));
        st->js->enqueue([st, fn]() {
          (*fn)();
          if (st->remaining.fetch_sub(1, std::memory_order_acq_rel) == 1) {
            std::lock_guard<std::mutex> lock(st->mutex);
            st->cv.notify_all();
          }
        });
      }
    }

    // Wait for all tasks in the group to complete.
    void wait() {
      auto st = st_;
      while (st->remaining.load(std::memory_order_acquire) > 0) {
        // Help the pool make progress while we're waiting.
        if (!st->js->tryRunOne()) {
          std::unique_lock<std::mutex> lock(st->mutex);
          st->cv.wait_for(lock, std::chrono::milliseconds(1), [&]() {
            return st->remaining.load(std::memory_order_acquire) == 0;
          });
        }
      }
    }

    std::size_t pending() const {
      return (std::size_t)std::max(0, st_->remaining.load(std::memory_order_relaxed));
    }

  private:
    struct State {
      explicit State(JobSystem* jsIn) : js(jsIn) {}
      JobSystem* js{nullptr};
      std::atomic<int> remaining{0};
      std::mutex mutex{};
      std::condition_variable cv{};
    };

    std::shared_ptr<State> st_{};
  };

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

    // Keep the user callable by reference (matches old semantics) while still
    // allowing move-only functors (since we never copy it into tasks).
    auto* fnPtr = &fn;

    auto workLoop = [&, fnPtr]() {
      for (;;) {
        const std::size_t start = next.fetch_add(grain, std::memory_order_relaxed);
        if (start >= count) break;
        const std::size_t end = (start + grain < count) ? (start + grain) : count;
        for (std::size_t i = start; i < end; ++i) (*fnPtr)(i);
      }
    };

    // Avoid nested-parallelism deadlocks by ensuring the calling thread also works.
    // Launch (workers - 1) helper tasks and run one workLoop inline.
    //
    // Using TaskGroup avoids the allocations/futures overhead of submit().
    TaskGroup g(*this);
    for (std::size_t i = 0; i + 1 < workers; ++i) {
      g.run(workLoop);
    }

    workLoop();
    g.wait();
  }

private:
  void workerLoop() {
    tlsCurrent_ = this;
    for (;;) {
      std::function<void()> task;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [&]() { return stopping_ || !queue_.empty(); });
        if (stopping_ && queue_.empty()) break;

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
    tlsCurrent_ = nullptr;
  }

  std::size_t threadCount_{1};
  std::vector<std::thread> threads_{};

  std::mutex mutex_{};
  std::condition_variable cv_{};
  std::condition_variable idleCv_{};

  std::deque<std::function<void()>> queue_{};
  std::size_t active_{0};
  bool stopping_{false};

  inline static thread_local JobSystem* tlsCurrent_{nullptr};
};

} // namespace stellar::core
