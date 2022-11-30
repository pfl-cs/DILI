#include <unistd.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>

#ifndef UTILS_LINUX_SYS_UTILS_H
#define UTILS_LINUX_SYS_UTILS_H

namespace linux_sys_utils {
    long perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                         int cpu, int group_fd, unsigned long flags);
}
#endif //UTILS_LINUX_SYS_UTILS_H
