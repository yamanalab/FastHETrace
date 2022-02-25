// Copyright 2021 Yu Ishimaki

// Define PROFILE to enable TIC-TOC timing measurements
#define PROFILE

#include <cassert>
#include <fstream>

#include <omp.h>       // NOLINT
#include "palisade.h"  // NOLINT
#include "utils/serialize-binary.h"
#include "ciphertext-ser.h"     // NOLINT
#include "cryptocontext-ser.h"  // NOLINT
#include "scheme/ckks/ckks-ser.h"
#include "pubkeylp-ser.h"  // NOLINT

#include "argparse.hpp"  // NOLINT

using namespace std;  // NOLINT

using namespace lbcrypto;  // NOLINT
using PolyType = DCRTPoly;
using Ctxt = Ciphertext<PolyType>;
using EvkAut = shared_ptr<map<usint, LPEvalKey<PolyType>>>;
using FHEContext = CryptoContext<PolyType>;

#include "stat.hpp"
#include "mem_usage.h"  // NOLINT
#include "trace_keygen.hpp"
#include "automo_opt.hpp"
#include "hoist.hpp"
#include "global.h"  // NOLINT
#include "types.hpp"

double EvalTraceUnrollHKMult(
    const size_t M, const CryptoContext<DCRTPoly> cc,
    const LPKeyPair<DCRTPoly>& keys, const uint32_t L, const size_t num_unroll,
    const size_t num_mod_reduce,
    const vector<vector<usint>>& auto_indices_parallel,
    map<usint, std::pair<std::vector<DCRTPoly>, std::vector<DCRTPoly>>>&
        inv_evks) {
  // Input
  vector<complex<double>> x = {0, 0, 0, 0, 0, 0, 0, 1};
  Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x);

  cout << "Input x: " << ptxt << endl;
  auto c = cc->Encrypt(keys.publicKey, ptxt);
  vector<complex<double>> x1(cc->GetRingDimension() / 2, 1);
  auto c1 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(x1));
  for (size_t i = 0; i < num_mod_reduce; ++i) {
    c = cc->EvalMult(c, c1);
    c = cc->Rescale(c);
  }
  cout << "# Dropped = " << c->GetLevel()
       << " # MODs Available = " << (L + 1) - c->GetLevel() << endl;

  TimeVar t;
  TIC(t);

  for (size_t numitr = 0; numitr < num_unroll; ++numitr)
    c += EvalHoistedAutomorph(c, auto_indices_parallel[numitr], inv_evks);
  double timing = TOC(t);

  Plaintext result;
  cc->Decrypt(keys.secretKey, c, &result);
  result->SetLength(5);
  cout << "\t Result = " << result << endl;
  cout << "Timing = " << timing << " ms...";

  auto m = cc->GetCyclotomicOrder();

  size_t num_ans = 0;
  int ans = 1;
  size_t num_expected_ans = 1UL << M;
  if (M == size_t(log2(m / 2))) {
    ans = 2;
    num_expected_ans = (1UL << (M - 1));
  }

  // See all the slot values
  for (size_t i = 0; i < m / 4; ++i) {
    if (round(result->GetCKKSPackedValue()[i].real()) == ans) num_ans += 1;
  }

  if (num_ans == num_expected_ans)
    cout << " correct :) #1 = " << num_ans << " Expected = " << num_expected_ans
         << endl;
  else
    cout << "      wrong #1 = " << num_ans << " Expected = " << num_expected_ans
         << endl;

  return timing;
}

double EvalTraceUnrollSingleTh(
    const size_t M, const CryptoContext<DCRTPoly> cc,
    const LPKeyPair<DCRTPoly>& keys, const uint32_t L, const size_t num_unroll,
    const size_t num_mod_reduce,
    const vector<vector<usint>>& auto_indices_parallel,
    map<usint, std::pair<std::vector<DCRTPoly>, std::vector<DCRTPoly>>>&
        inv_evks) {
  // Input
  vector<complex<double>> x = {0, 0, 0, 0, 0, 0, 0, 1};
  Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x);

  cout << "Input x: " << ptxt << endl;
  auto c = cc->Encrypt(keys.publicKey, ptxt);
  vector<complex<double>> x1(cc->GetRingDimension() / 2, 1);
  auto c1 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(x1));
  for (size_t i = 0; i < num_mod_reduce; ++i) {
    c = cc->EvalMult(c, c1);
    c = cc->Rescale(c);
  }
  cout << "# Dropped = " << c->GetLevel()
       << " # MODs Available = " << (L + 1) - c->GetLevel() << endl;

  // First, we perform 7 regular (non-hoisted) rotations
  // and measure the runtime.
  TimeVar t;
  TIC(t);

  for (size_t numitr = 0; numitr < num_unroll; ++numitr)
    c += EvalHoistedAutomorphHKSingleTh(c, auto_indices_parallel[numitr],
                                        inv_evks);

  double timing = TOC(t);

  Plaintext result;
  cc->Decrypt(keys.secretKey, c, &result);
  result->SetLength(5);
  cout << "\t Result = " << result << endl;
  cout << "Timing = " << timing << " ms...";

  auto m = cc->GetCyclotomicOrder();

  size_t num_ans = 0;
  int ans = 1;
  size_t num_expected_ans = 1UL << M;
  if (M == size_t(log2(m / 2))) {
    ans = 2;
    num_expected_ans = (1UL << (M - 1));
  }

  // See all the slot values
  for (size_t i = 0; i < m / 4; ++i) {
    if (round(result->GetCKKSPackedValue()[i].real()) == ans) num_ans += 1;
  }

  if (num_ans == num_expected_ans)
    cout << " correct :) #1 = " << num_ans << " Expected = " << num_expected_ans
         << endl;
  else
    cout << "      wrong #1 = " << num_ans << " Expected = " << num_expected_ans
         << endl;

  return timing;
}

template <class T>
void ShowVec(const std::vector<T>& v) {
  std::size_t e = v.size();
  std::cout << "[";
  for (std::size_t i = 0, e = v.size(); i < e - 1; ++i) {
    std::cout << v[i] << ",";
  }
  std::cout << v[e - 1] << "]" << std::endl;
}

int main(int argc, char const* argv[]) {
#ifdef USE_FAST_AUTOMORPH
  cout << "Fast Automorph" << endl;
#endif

  int max_num_thread = omp_get_max_threads();

  // try to unroll evenly distributed
  // In case of evalsum, logN-1/#Unroll
  argparse::ArgumentParser parser("bench_unroll");

  parser.add_argument("logN").help("logN").required().scan<'i', int>();
  parser.add_argument("M")
      .help("iterations for rotations-and-sums")
      .required()
      .scan<'i', int>();
  parser.add_argument("h").help("# of unrolling").required().scan<'i', int>();
  parser.add_argument("L")
      .help("maximal level (# of RNS moduli for this level = L + 1)")
      .required()
      .scan<'i', uint32_t>();
  parser.add_argument("ell")
      .help(
          "level to examine (# of RNS moduli for this level = \ell + 1 where "
          "\ell <= L)")
      .required()
      .scan<'i', uint32_t>();
  parser.add_argument("delta")
      .help("scaling factor for the CKKS encoding")
      .required()
      .scan<'i', uint32_t>();
  parser.add_argument("dnum")
      .help(
          "maximal number of digits for key-switching (an integer that divides "
          "L+1)")
      .required()
      .scan<'i', uint32_t>();
  parser.add_argument("T")
      .help("the number of experiments to perform")
      .required()
      .scan<'i', int>();
  parser.add_argument("-t", "--thread")
      .help("# of threads used for calculation")
      .default_value(max_num_thread)
      .scan<'i', int>();
  parser.add_argument("-c", "--config-path")
      .help("Path to the directory which contains keys and crypto context");

  try {
    parser.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << parser;
    std::exit(1);
  }

  auto logN = parser.get<int>("logN");
  auto M = parser.get<int>("M");
  auto num_unroll = parser.get<int>("h");
  auto depth = parser.get<uint32_t>("L");
  auto ell = parser.get<uint32_t>("ell");
  auto scale = parser.get<uint32_t>("delta");
  auto dnum = parser.get<uint32_t>("dnum");
  auto num_exp = parser.get<int>("T");
  auto num_thread = parser.is_used("-t") ? parser.get<int>("-t") : 1;
  auto config_store_path =
      parser.is_used("-c")
          ? static_cast<filesystem::path>(parser.get<std::string>("-c"))
          : fastertracetype::CONFIG_FILE_PATH;

  assert(logN >= M);
  assert(M >= num_unroll);
  assert(depth >= ell);

  // config/<per-parameter directory>/*.key
  if (filesystem::exists(config_store_path) &&
          !filesystem::is_directory(config_store_path) ||
      config_store_path.has_extension()) {
    std::cerr << config_store_path
              << " is not a directory. Please specify a valid path."
              << std::endl;
    std::exit(1);
  }
  filesystem::create_directories(config_store_path);

  std::cout << "Config path is set to " << config_store_path << std::endl;

  omp_set_num_threads(num_thread);

#pragma omp parallel
  {
#pragma omp single
    { std::cout << "# of threads: " << omp_get_num_threads() << std::endl; }
  }

  size_t num_mod_reduce = depth - ell;

  SecurityLevel securityLevel = HEStd_128_classic;
  RescalingTechnique rsTech = APPROXRESCALE;
  KeySwitchTechnique ksTech = HYBRID;

  // Create or Load crypto context
  CryptoContext<DCRTPoly> cc;
  LPKeyPair<DCRTPoly> keys;

  std::string config_dirname =
      "logN_"s + std::to_string(logN) + "_L_"s + std::to_string(depth) +
      "_delta_"s + std::to_string(scale) + "_dnum_"s + std::to_string(dnum);
  auto config_path = config_store_path / config_dirname;

  filesystem::create_directories(config_path);

  auto context_file_path = config_path / fastertracetype::CONTEXT_FILE_NAME;
  auto pubkey_file_path = config_path / fastertracetype::PUBKEY_FILE_NAME;
  auto seckey_file_path = config_path / fastertracetype::SECKEY_FILE_NAME;
  auto evalmultkey_file_path =
      config_path / fastertracetype::EVALMULTKEY_FILE_NAME;

  if (filesystem::exists(pubkey_file_path) &&
      filesystem::exists(seckey_file_path) &&
      filesystem::exists(evalmultkey_file_path)) {
#ifdef DEBUG
    std::cout << "Found context file." << std::endl;
#endif
    cc->ClearEvalMultKeys();
    cc->ClearEvalAutomorphismKeys();
    cc->ClearEvalSumKeys();
    lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
    if (!Serial::DeserializeFromFile(context_file_path, cc, SerType::BINARY)) {
      std::cerr << "Failed to load crypto context file" << std::endl;
      std::exit(1);
    }
    if (!Serial::DeserializeFromFile(pubkey_file_path, keys.publicKey,
                                     SerType::BINARY)) {
      std::cerr << "Failed to load public key" << std::endl;
      std::exit(1);
    }
    if (!Serial::DeserializeFromFile(seckey_file_path, keys.secretKey,
                                     SerType::BINARY)) {
      std::cerr << "Failed to load secret key" << std::endl;
      std::exit(1);
    }
    {
      std::ifstream mult_key_file(evalmultkey_file_path,
                                  std::ios::out | std::ios::binary);
      if (mult_key_file.is_open()) {
        if (!cc->DeserializeEvalMultKey(mult_key_file, SerType::BINARY)) {
          std::cerr << "Failed to load eval mult keys" << std::endl;
          std::exit(1);
        }
        std::cout << "EvalMult/relinearization keys have been loaded"
                  << std::endl;
      } else {
        std::cerr << "Error deserializing EvalMult keys" << std::endl;
        std::exit(1);
      }
    }
  } else {
    cc = CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
        depth, scale, 1 << (logN - 1), securityLevel, 1 << logN, rsTech, ksTech,
        dnum);
    cc->Enable(ENCRYPTION);
    cc->Enable(SHE);
    cc->Enable(LEVELEDSHE);
    // Serialize cryptocontext
    if (!Serial::SerializeToFile(context_file_path, cc, SerType::BINARY)) {
      std::cerr << "Error writing serialization of the crypto context"
                << std::endl;
      std::exit(1);
    }

    // KeyGen
    keys = cc->KeyGen();
    if (!Serial::SerializeToFile(pubkey_file_path, keys.publicKey,
                                 SerType::BINARY)) {
      std::cerr << "Error writing serialization of the public key" << std::endl;
      std::exit(1);
    }
    if (!Serial::SerializeToFile(seckey_file_path, keys.secretKey,
                                 SerType::BINARY)) {
      std::cerr << "Error writing serialization of the secret key" << std::endl;
      std::exit(1);
    }

    // Relinealization key gen
    cc->EvalMultKeyGen(keys.secretKey);
    {
      std::ofstream mult_key_file(evalmultkey_file_path,
                                  std::ios::out | std::ios::binary);
      if (mult_key_file.is_open()) {
        if (!cc->SerializeEvalMultKey(mult_key_file, SerType::BINARY)) {
          std::cerr << "Error writing eval mult keys" << std::endl;
          std::exit(1);
        }
        std::cout << "EvalMult/ relinearization keys have been serialized"
                  << std::endl;
      } else {
        std::cerr << "Error serializing EvalMult keys" << std::endl;
        std::exit(1);
      }
    }
  }

  // evk
  fastertracetype::InvEvks inv_evks;
  fastertracetype::AutoIndicesPar auto_indices_par;

  std::string rotate_key_file_name = "rotate_"s + std::to_string(M) + "_"s +
                                     std::to_string(num_unroll) + ".key";
  std::string inv_evks_file_name = "inv_evks_"s + std::to_string(M) + "_"s +
                                   std::to_string(num_unroll) + ".bin";
  std::string auto_indices_par_file_name = "auto_indices_par_"s +
                                           std::to_string(M) + "_"s +
                                           std::to_string(num_unroll) + ".bin";

  auto rotate_key_path = config_path / rotate_key_file_name;
  auto inv_evks_path = config_path / inv_evks_file_name;
  auto auto_indices_par_path = config_path / auto_indices_par_file_name;

  if (filesystem::exists(rotate_key_path)) {
    std::ifstream rot_key(rotate_key_path, std::ios::in | std::ios::binary);
    if (!rot_key.is_open()) {
      std::cerr << "Cannot read serialization of rotate key" << std::endl;
      std::exit(1);
    }
    if (!cc->DeserializeEvalAutomorphismKey(rot_key, SerType::BINARY)) {
      std::cerr << "Failed to load rotation key." << std::endl;
      std::exit(1);
    }
    if (!Serial::DeserializeFromFile(inv_evks_path, inv_evks,
                                     SerType::BINARY)) {
      std::cerr << "Failed to load inv_evks." << std::endl;
      std::exit(1);
    }
    if (!Serial::DeserializeFromFile(auto_indices_par_path, auto_indices_par,
                                     SerType::BINARY)) {
      std::cerr << "Failed to load auto_indices_par." << std::endl;
      std::exit(1);
    }
  } else {
    // ctxt
    vector<complex<double>> x = {0};
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x);
    auto cipher = cc->Encrypt(keys.publicKey, ptxt);

    // Generate rotation key (which used in unrolled EvalSum;summaization all
    // vector element)
    EvalSumKeyGenUnroll(cc, keys, M, num_unroll, auto_indices_par.GetData());
    std::ofstream rot_key(rotate_key_path, std::ios::out | std::ios::binary);
    rot_key.exceptions(std::ios::failbit);
    if (rot_key.is_open()) {
      if (!cc->SerializeEvalAutomorphismKey(rot_key, SerType::BINARY)) {
        std::cerr << "Error writing rotation keys" << std::endl;
        std::exit(1);
      }
      std::cout << "Rotation keys have been serialized" << std::endl;
    } else {
      std::cerr << "Error serializing Rotation keys" << std::endl;
      std::exit(1);
    }

    // This step should be performed at lower-level
    for (auto vec : auto_indices_par.GetData())
      for (auto id : vec) PreComputeEvk(cc, id, cipher, inv_evks.GetData());

    // save inv_evks and auto_indices_par
    if (!Serial::SerializeToFile(inv_evks_path, inv_evks, SerType::BINARY)) {
      std::cerr << "Error serializing inv_evks" << std::endl;
      std::exit(1);
    }
    if (!Serial::SerializeToFile(auto_indices_par_path, auto_indices_par,
                                 SerType::BINARY)) {
      std::cerr << "Error serializing auto_indices_par" << std::endl;
      std::exit(1);
    }

#ifdef DEBUG
    int numkeys = 0;
    for (auto key_vec : auto_indices_par.GetData()) {
      std::cout << "[" << key_vec.size() << ",";
      numkeys += key_vec.size();
    }
    std::cout << "]";
    std::cout << "# keys =" << numkeys << std::endl;
#endif
  }

  vector<double> timings;
  if (num_thread == 1) {
    for (auto i = 0; i < num_exp; ++i)
      timings.push_back(EvalTraceUnrollSingleTh(
          M, cc, keys, depth, num_unroll, num_mod_reduce,
          auto_indices_par.GetData(), inv_evks.GetData()));
  } else {
    for (auto i = 0; i < num_exp; ++i)
      timings.push_back(EvalTraceUnrollHKMult(
          M, cc, keys, depth, num_unroll, num_mod_reduce,
          auto_indices_par.GetData(), inv_evks.GetData()));
  }

  int num_special_mod = (depth + 1) / dnum;

  string outprefix = "";
  std::string outfile =
      "/tmp/result/TraceRuntime/" + outprefix + "dynamic_dnum_unroll_HK_logN" +
      to_string(logN) + "_M" + to_string(M) + "_L" + to_string(depth) + "_ell" +
      to_string(ell) + "_numunroll" + to_string(num_unroll) + "_scale" +
      to_string(scale) + "_numth" + to_string(num_thread) + "_k" +
      to_string(num_special_mod) + "_numexp" + to_string(num_exp) + ".csv";

  StatDropFirst s_unroll(timings,
                         {static_cast<double>(getPeakRSS()) / (1 << 20)});
  s_unroll.Write(cout);
  ofstream ofs(outfile);
  s_unroll.WriteCSV(ofs);
  ofs.close();

  return 0;
}
