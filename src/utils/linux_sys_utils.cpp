#include <unistd.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>
#include <chrono>
#include <functional>
#include <thread>

namespace linux_sys_utils {
    long perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                         int cpu, int group_fd, unsigned long flags) {
        int ret;

        ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,
                      group_fd, flags);
        return ret;
    }

    uint64_t timing(std::function<void()> fn) {
        const auto start = std::chrono::high_resolution_clock::now();
        fn();
        const auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count();
    }
}

