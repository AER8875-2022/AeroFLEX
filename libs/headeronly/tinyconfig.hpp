/*
MIT License

Copyright (c) 2023 Samuel Ayala

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#include <unordered_map>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <functional>

namespace tiny {

class config {
    public:
        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> config;
        std::unordered_map<std::string, std::vector<std::unordered_map<std::string, std::string>>> config_vec;
        std::vector<std::string> sections;

        template<typename T> T get(const std::string &section, const std::string &key);
        template<typename T> T get(const std::string &section, const std::string &key, const T default_value);
        template<typename T> T get_i(const std::string &section, const std::string &key, const int idx);
        template<typename T> T get_i(const std::string &section, const std::string &key, const int idx, const T default_value);
        bool has(const std::string &section, const std::string &key);
        bool has_i(const std::string &section, const std::string &key, const int idx);
        int how_many(const std::string &section);

        bool read(const std::string& ifilename);
        bool write(const std::string& ofilename);

    private:
        void read_file(std::ifstream &ifile);
        void write_file(std::ofstream& ofile);
};

template<typename T>
inline T convert_value(std::string s) {
    if constexpr (std::is_same<T, std::string>()) {
        return s;
    } else {
        T value;
        std::istringstream(s) >> value;
        return value;
    }
}

template<>
inline bool convert_value(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    static const std::unordered_map<std::string, bool> s2b{
        {"1", true},  {"true", true},   {"yes", true}, {"on", true},
        {"0", false}, {"false", false}, {"no", false}, {"off", false},
    };
    auto const val = s2b.find(s);
    if (val == s2b.end()) {
        throw std::runtime_error("'" + s + "' is not a valid boolean value.");
    }
    return val->second;
}

constexpr auto err_missing_section = [](const std::string &s) {throw std::runtime_error("Section [" + s + "] not found");};
constexpr auto err_missing_setting = [](const std::string &s) {throw std::runtime_error("Setting (" + s + ") not found");};
constexpr auto err_missing_config = [](const std::string &s) {throw std::runtime_error("Config file: " + s + " not found");};
constexpr auto err_missing_closing = [](const std::string &s) {throw std::runtime_error("Closing character is missing on line: " + s);};

constexpr auto err_failed_read = [](const std::string &s) {throw std::runtime_error("Failed to read :" + s);};
constexpr auto err_failed_write = [](const std::string &s) {throw std::runtime_error("Failed to write :" + s);};

template<typename T>
inline auto get_map_elem(std::unordered_map<std::string, T> &map, const std::string &key, const std::function<void(const std::string&)> &&err) {
    auto elem = map.find(key);
    if (elem == map.end()) err(key);
    return elem;
}

template<typename T>
inline auto get_map_elem_default_value(std::unordered_map<std::string, T> &map, const std::string &key, const std::string default_value) {
    auto elem = map.find(key);
    if (elem == map.end()) elem = map.insert({key, default_value}).first;
    return elem;
}

template<typename T>
inline T config::get(const std::string &section, const std::string &key) {
    const auto section_elem = get_map_elem(config, section, err_missing_section);
    const auto value_elem = get_map_elem(section_elem->second, key, err_missing_setting);
    return convert_value<T>(value_elem->second);
}

template<typename T>
inline T config::get(const std::string &section, const std::string &key, const T default_value) {
    const auto section_elem = get_map_elem(config, section, err_missing_section);
    const auto value_elem = get_map_elem_default_value(section_elem->second, key, std::to_string(default_value));
    return convert_value<T>(value_elem->second);
}

template<typename T>
inline T config::get_i(const std::string &section, const std::string &key, const int idx) {
    const auto section_elem = get_map_elem(config_vec, section, err_missing_section);
    const auto value_elem = get_map_elem(section_elem->second[idx], key, err_missing_setting);
    return convert_value<T>(value_elem->second);
}

template<typename T>
inline T config::get_i(const std::string &section, const std::string &key, const int idx, const T default_value) {
    const auto section_elem = get_map_elem(config_vec, section, err_missing_section);
    const auto value_elem = get_map_elem_default_value(section_elem->second[idx], key, std::to_string(default_value));
    return convert_value<T>(value_elem->second);
}

inline bool config::has(const std::string &section, const std::string &key) {
    const auto section_elem = get_map_elem(config, section, err_missing_section);
    const auto value_elem = section_elem->second.find(key);
    return value_elem == section_elem->second.end() ? false : true;
}

inline bool config::has_i(const std::string &section, const std::string &key, const int idx) {
    const auto section_elem = get_map_elem(config_vec, section, err_missing_section);
    const auto value_elem = section_elem->second[idx].find(key);
    return value_elem == section_elem->second[idx].end() ? false : true;
}

inline int config::how_many(const std::string &section) {
    return static_cast<int>(get_map_elem(config_vec, section, err_missing_section)->second.size());
}

inline bool config::read(const std::string& ifilename) {
    std::filesystem::path path(ifilename);
    if (!std::filesystem::exists(path)) err_missing_config(ifilename);

    std::ifstream f(path);
	if (f.is_open()) {
        read_file(f);
        f.close();
	} else return false;
    return !f.bad();
}

inline bool config::write(const std::string& ofilename) {
    std::filesystem::path path(ofilename);
    //if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    std::ofstream f(path);
	if (f.is_open()) {
        write_file(f);
        f.close();
	} else return false;
    return !f.bad();
}

inline static std::string clean_str(std::string str) {
    str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
    return str;
};

inline static void extract(std::unordered_map<std::string, std::string> &map, const std::string& line, const size_t sep, const size_t end, const int line_nb) {
    if (sep == std::string::npos) return;
    if (end == std::string::npos) err_missing_closing(std::to_string(line_nb));
    map[clean_str(line.substr(0, sep))] = clean_str(line.substr(sep + 1, end - sep - 1));
}

inline void config::read_file(std::ifstream &f) {
    std::string line;
    std::string section = "_";
    const std::string whitespaces (" \t\f\v\n\r");
    int line_nb = 0;
    while (std::getline(f, line)) {
        line_nb++;
        if (auto first_bracket = line.find("["); first_bracket != std::string::npos) {
            auto second_bracket = line.find("]");
            if (second_bracket == std::string::npos) err_missing_closing(std::to_string(line_nb));
            section = clean_str(line.substr(first_bracket + 1, second_bracket - first_bracket - 1));
            config[section];
            sections.push_back(section);
        } else if (auto first_abracket = line.find("<"); first_abracket != std::string::npos) {
            std::unordered_map<std::string, std::string> entry;
            if (first_abracket == line.find_last_not_of(whitespaces)) {
                while (std::getline(f, line) && line.find(">") == std::string::npos) {
                    line_nb++;
                    extract(entry, line, line.find("="), line.find(","), line_nb);
                }
            } else if (line.find(">") != std::string::npos) {
                line = line.substr(first_abracket + 1);
                auto comma = line.find(",");
                while (comma != std::string::npos) {
                    extract(entry, line, line.find("="), comma, line_nb);
                    line = line.substr(comma + 1, std::string::npos);
                    comma = line.find(",");
                }
                extract(entry, line, line.find("="), line.find(">"), line_nb);
            }
            config_vec[section].push_back(std::move(entry));
        } else {
            extract(config[section], line, line.find("="), line.length(), line_nb);
        }
    }
}

inline void config::write_file(std::ofstream &f) {
    for (const auto &section : sections) {
        f << "[" << section << "] \n";
        for (const auto & [setting, value] : config[section]) {
            f << setting << " = " << value << "\n";
        }
        for (const auto &entry : config_vec[section]) {
            if (entry.size() > 1) {
                f << "<\n";
                for (const auto& [setting, value] : entry) {
                    f << "\t" << setting << " = " << value << ",\n";
                }
                f << ">\n";
            } else if (entry.size() == 1){
                f << "<" << entry.begin()->first << " = " << entry.begin()->second << ">\n";
            }
        }
        f << "\n";
    }
}
} // namespace tiny