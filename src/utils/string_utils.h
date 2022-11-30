//
// Created by Pengfei Li on 11/18/22.
//
#include <string>
#include <algorithm>
#include <regex>
#include <cassert>
#ifndef UTILS_STRING_UTILS_H
#define UTILS_STRING_UTILS_H

namespace string_utils {
    // trim from start (in place)
    static inline void ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
    }

    // trim from end (in place)
    static inline void rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), s.end());
    }

    // trim from both ends (in place)
    static inline void trim(std::string &s) {
        ltrim(s);
        rtrim(s);
    }

    // trim from start (copying)
    static inline std::string ltrim_copy(std::string s) {
        ltrim(s);
        return s;
    }

    // trim from end (copying)
    static inline std::string rtrim_copy(std::string s) {
        rtrim(s);
        return s;
    }

    // trim from both ends (copying)
    static inline std::string trim_copy(std::string s) {
        trim(s);
        return s;
    }

    // tolower
    static inline void lowercase(std::string &s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); }
        );
    }

    static inline std::string lowercase_copy(std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); }
        );
        return s;
    }

    static inline std::string replace_substr(const std::string &s, const std::string &a, const std::string &b) {
        return std::regex_replace(s, std::regex(a), b);
    }

    static void split(const std::string &str, const std::string &delim, std::vector<std::string> &items) {
        size_t prev = 0, pos = 0;
        do {
            pos = str.find(delim, prev);
            if (pos == str.npos) {
                pos = str.size();
            }
            std::string item = str.substr(prev, pos-prev);
            if (!item.empty()) {
                items.emplace_back(item);
            }
            prev = pos + delim.size();
        } while (pos < str.size() && prev < str.size());
    }

    static std::string join(std::vector<std::string> &items, const std::string &delim) {
        assert (items.size() > 0);
        std::string s = items[0];
        for (auto i = 1; i < items.size(); ++i) {
            s += delim + items[i];
        }
        return s;
    }

}

#endif //UTILS_STRING_UTILS_H
