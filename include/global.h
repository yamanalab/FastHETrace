// 2021 Kotaro Inoue

#pragma once

#include <filesystem>

namespace fastertracetype {
const std::filesystem::path CONFIG_FILE_PATH = "config";
const std::filesystem::path CONTEXT_FILE_NAME = "context.txt";
const std::filesystem::path PUBKEY_FILE_NAME = "public.key";
const std::filesystem::path SECKEY_FILE_NAME = "secret.key";
const std::filesystem::path EVALMULTKEY_FILE_NAME = "evalmult.key";
const std::filesystem::path EVALSUMKEY_FILE_NAME = "evalsum.key";
}  // namespace fastertracetype
