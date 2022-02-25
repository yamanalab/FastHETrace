#pragma once

#include "palisade.h"
#include <cereal/cereal.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>

// CEREAL_REGISTER_TYPE(fastertracetype::InvEvks);
// CEREAL_REGISTER_TYPE(fastertracetype::AutoIndicesPar);

namespace fastertracetype {
class InvEvks {
 public:
  using T = map<usint, std::pair<std::vector<DCRTPoly>, std::vector<DCRTPoly>>>;
  InvEvks() = default;
  InvEvks(T v) : data_(v){};
  ~InvEvks() = default;

  T& GetData() { return data_; }

 private:
  T data_;

  friend class cereal::access;

  template <class Archive>
  void save(Archive& ar, std::uint32_t const version) const {
    ar(::cereal::make_nvp("data", data_));
  }

  template <class Archive>
  void load(Archive& ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }
    ar(::cereal::make_nvp("data", data_));
  }

  std::string SerializedObjectName() const { return "InvEvks"; }
  static uint32_t SerializedVersion() { return 1; }
};

class AutoIndicesPar {
 public:
  using T = vector<vector<usint>>;
  AutoIndicesPar() = default;
  AutoIndicesPar(T v) : data_(v){};
  ~AutoIndicesPar() = default;

  T& GetData() { return data_; }

 private:
  T data_;

  friend class cereal::access;

  template <class Archive>
  void save(Archive& ar, std::uint32_t const version) const {
    ar(::cereal::make_nvp("data", data_));
  }

  template <class Archive>
  void load(Archive& ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }
    ar(::cereal::make_nvp("data", data_));
  }

  std::string SerializedObjectName() const { return "AutoIndicesPar"; }
  static uint32_t SerializedVersion() { return 1; }
};
}  // namespace fastertracetype
