#include <string>
#include <algorithm>
#include <regex>
#include <fstream>
#include <vector>
#include <iostream>
#include <filesystem>
#include "file_utils.h"
namespace fs = std::filesystem;

namespace file_utils {
    void readLines(const std::string &filepath, std::vector<std::string> &lines) {
        std::ifstream fin(filepath);
        std::string line;
        while (getline(fin, line)) {
            lines.emplace_back(line);
        }
        fin.close();
    }

    void writeLines(const std::string &filepath, std::vector<std::string> &lines) {
        std::ofstream fout(filepath);
        for (const std::string &line: lines) {
            fout << line << std::endl;
        }
        fout << std::flush;
        fout.close();
    }

    std::string path_join(const std::string &_dir, const std::string &_filename) {
        fs::path dir = _dir;
        auto filepath = dir / _filename;
        return filepath.string();
    }

//    std::string path_join(const fs::path &dir, const std::string &_filename) {
//        auto filepath = dir / _filename;
//        return filepath.string();
//    }

//    bool file_exists(const char *path) {
//        FILE *fp = NULL;
//
//        if (NULL == (fp = fopen(path, "rb"))) {
//            return false;
//        }
//
//        fclose(fp);
//        return true;
//    }

    bool file_exists(const std::string &path_str) {
        fs::path _path = path_str;
        return fs::exists(_path);
    }

//    int path_status(const string &path_name) {
//        struct stat info;
//        if (stat(path_name.c_str(), &info) != 0) {
//            return 0; // path not exists
//        }
//        else if (info.st_mode & S_IFDIR) {
//            return 1; // path is a directory
//        }
//        else {
//            return 2; // path is not a directory
//        }
//    }

    int path_status(const std::string &path_str) {
        fs::path _path = path_str;
        if (!fs::exists(_path)) {
            return 0; // path not exists
        }
        if (fs::is_directory(_path)) {
            return 1; // path is a directory
        } else {
            return 2; // path is not a directory
        }
    }

    bool detect_and_create_dir(const std::string &_dir) {
        fs::path dir = _dir;
        if (!fs::exists(dir)) {
            return fs::create_directories(dir);
        } else {
            return true;
        }
    }

//    bool detect_and_create_dir(const fs::path &dir) {
//        if (!fs::exists(dir)) {
//            return fs::create_directory(dir);
//        } else {
//            return true;
//        }
//    }
}

