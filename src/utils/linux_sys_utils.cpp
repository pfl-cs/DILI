#include <unistd.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>

namespace linux_sys_utils {
    long perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                         int cpu, int group_fd, unsigned long flags) {
        int ret;

        ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,
                      group_fd, flags);
        return ret;
    }
}