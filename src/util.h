/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#pragma once 

#include <vector>
#include <string>
#include <unordered_map>

#include <sys/types.h>
#include <sys/stat.h>

/******************************************************************************/

mode_t stat_file(const std::string &path);

std::vector<std::string> split(const std::string &s, char delim);

std::string rc(const std::string &s);

/******************************************************************************/

double tau(double edit_error, int kmer_size);

int relaxed_jaccard_estimate(int s, int kmer_size, std::unordered_map<int, int> &mm);

