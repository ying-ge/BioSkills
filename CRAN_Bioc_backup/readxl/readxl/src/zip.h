#pragma once

#include "rapidxml/rapidxml.h"
#include <string>

std::string zip_buffer(const std::string& zip_path, const std::string& file_path);
bool zip_has_file(const std::string& zip_path, const std::string& file_path);
