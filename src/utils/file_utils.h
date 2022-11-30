#include <string>
#include <algorithm>
#include <regex>
#include <fstream>
#include <vector>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

#ifndef UTILS_FILE_UTILS_H
#define UTILS_FILE_UTILS_H

namespace file_utils {
    void readLines(const std::string &filepath, std::vector<std::string> &lines);
    void writeLines(const std::string &filepath, std::vector<std::string> &lines);

    std::string path_join(const std::string &_dir, const std::string &_filename);

    bool file_exists(const std::string &path_str);

    int path_status(const std::string &path_str);

    bool detect_and_create_dir(const std::string &_dir);
}

#endif //UTILS_FILE_UTILS_H
