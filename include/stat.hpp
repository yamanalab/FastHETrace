#pragma once

#include "csv2.hpp"

template <class T>
T ComputeAvgDropOne(const std::vector<T>& v) {
  std::size_t e = v.size();
  T sum = 0.0;
  for (std::size_t i = 1; i < e; ++i) {
    sum += v[i];
  }
  return sum / e;
}

struct StatDropFirst {
  StatDropFirst(std::vector<double>& v) { Setup(v); }

  void Setup(const std::vector<double>& v) {
    raw_datas = v;
    // Drop the first stat
    raw_datas.erase(raw_datas.begin());
    std::size_t n = raw_datas.size();
    auto sum = std::accumulate(raw_datas.begin(), raw_datas.end(), 0);

    mean = sum / n;

    stdev = 0;
    for (const auto& value : raw_datas) {
      stdev += (value - mean) * (value - mean);
    }
    stdev /= n;
    stdev = std::sqrt(stdev);

    minimum = *std::min_element(raw_datas.begin(), raw_datas.end());
    maximum = *std::max_element(raw_datas.begin(), raw_datas.end());
  }

  StatDropFirst(const std::vector<double>& v,
                const std::vector<double> oth_info) {
    Setup(v);
    other_info = oth_info;
  }

  void Write(std::ostream& os) {
    os << mean << "," << stdev << "," << minimum << "," << maximum << ",";

    if (other_info.size() == 0) {
      os << std::endl;
      return;
    }

    for (size_t i = 0, e = other_info.size(); i < e - 1; ++i)
      os << other_info[i] << ",";

    os << other_info.back() << std::endl;
  }

  /**
   * Write the stat as an csv file, including the first trial.
   */
  void WriteCSV(std::ofstream& ofs) {
    csv2::Writer<csv2::delimiter<','>> writer(ofs);
    std::vector<std::vector<std::string>> datas;
    datas.resize(raw_datas.size());
    std::transform(
        raw_datas.begin(), raw_datas.end(), datas.begin(),
        [](double c) { return std::vector<std::string>{std::to_string(c)}; });

    writer.write_row(vector<std::string>{"time"});
    writer.write_rows(datas);
    ofs.close();
  }

  double mean;
  double stdev;
  double minimum;
  double maximum;
  vector<double> other_info;
  vector<double> raw_datas;
};
