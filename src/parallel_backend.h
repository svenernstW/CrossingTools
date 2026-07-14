#pragma once

#include <functional>

#ifdef _OPENMP

  // ===== OpenMP path =====
  #include <omp.h>
  #include <atomic>
  #include <algorithm>

  inline std::atomic<int>& ct_thread_cap() {
    static std::atomic<int> cap{1};
    return cap;
  }

  inline void ct_set_threads(int n) {
    if (n < 1) n = 1;

    // CrossingTools-specific setting only.
    ct_thread_cap().store(n, std::memory_order_relaxed);
  }

  template <typename F>
  inline void ct_parallel_for(int begin, int end, F fn) {
    int n_threads =
      ct_thread_cap().load(std::memory_order_relaxed);

    // Never request more threads than loop iterations.
    n_threads = std::min(n_threads, std::max(1, end - begin));

    // Remember the user's existing OpenMP dynamic setting.
    const int old_dynamic = omp_get_dynamic();

    // Allow OpenMP to use fewer than n_threads when appropriate.
    omp_set_dynamic(1);

    #pragma omp parallel for schedule(static) num_threads(n_threads)
    for (int i = begin; i < end; ++i) {
      fn(i);
    }

    // Restore the previous setting immediately after this loop.
    omp_set_dynamic(old_dynamic);
  }

#else
  // ===== RcppParallel (oneTBB) path =====
  #include <RcppParallel.h>

  // oneTBB controls
  #include <tbb/global_control.h>
  #include <tbb/task_arena.h>
  #include <atomic>
  #include <memory>

  // Store the requested thread cap; 0 means "no cap".
  inline std::atomic<int>& ct_thread_cap() {
    static std::atomic<int> cap{0};
    return cap;
  }

  inline void ct_set_threads(int n) {
    // n <= 0 -> let TBB decide; n > 0 -> cap parallelism to n
    ct_thread_cap().store(
      n > 0 ? n : 0,
      std::memory_order_relaxed
    );
  }

  template <typename F>
  inline void ct_parallel_for(int begin, int end, F fn) {
    // Apply a temporary global cap if requested.
    std::unique_ptr<tbb::global_control> gc;

    int cap =
      ct_thread_cap().load(std::memory_order_relaxed);

    if (cap > 0) {
      gc = std::make_unique<tbb::global_control>(
        tbb::global_control::max_allowed_parallelism,
        cap
      );
    }

    struct Worker : RcppParallel::Worker {
      int b, e;
      F f;

      Worker(int _b, int _e, F _f)
        : b(_b), e(_e), f(_f) {}

      void operator()(std::size_t i, std::size_t j) {
        for (std::size_t k = i; k < j; ++k) {
          f(static_cast<int>(k));
        }
      }
    } w(begin, end, fn);

    RcppParallel::parallelFor(begin, end, w);

    // gc is destroyed here, releasing the cap.
  }

#endif