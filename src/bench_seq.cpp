// Define PROFILE to enable TIC-TOC timing measurements
#define PROFILE

#include "palisade.h"

#include <cassert>
#include <fstream>
#include <omp.h>

using namespace lbcrypto;
using PolyType = DCRTPoly;
using Ctxt = Ciphertext<PolyType>;
using EvkAut = shared_ptr<map<usint, LPEvalKey<PolyType>>>;
using FHEContext = CryptoContext<PolyType>;

#include "stat.hpp"
#include "mem_usage.h"
#include "trace_keygen.hpp"
#include "automo_opt.hpp"
#include "argparse.hpp"  // NOLINT

Ciphertext<DCRTPoly> EvalAutoKSHybrid(ConstCiphertext<DCRTPoly> cipherText,
                                      const usint autoIndex) {
  Ciphertext<DCRTPoly> newCiphertext = cipherText->CloneEmpty();

  // Retrieve the automorphism key that corresponds to the auto index.
  auto ek = cipherText->GetCryptoContext()
                ->GetEvalAutomorphismKeyMap(cipherText->GetKeyTag())
                .find(autoIndex)
                ->second;

  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(
          ek->GetCryptoParameters());

  LPEvalKeyRelin<DCRTPoly> evalKey =
      std::static_pointer_cast<LPEvalKeyRelinImpl<DCRTPoly>>(ek);

  std::vector<DCRTPoly> c;
// = cipherText->GetElements();
#ifdef USE_FAST_AUTOMORPH
  vector<usint> perm;
  GenAutomorphTable(cipherText->GetElements()[0].GetRingDimension(), autoIndex,
                    perm);
  c.push_back(std::move(cipherText->GetElements()[0].Permute(perm)));
  c.push_back(std::move(cipherText->GetElements()[1].Permute(perm)));
#else
  c.push_back(
      std::move(cipherText->GetElements()[0].AutomorphismTransform(autoIndex)));
  c.push_back(
      std::move(cipherText->GetElements()[1].AutomorphismTransform(autoIndex)));
#endif

  const std::vector<DCRTPoly>& b = evalKey->GetBVector();
  const std::vector<DCRTPoly>& a = evalKey->GetAVector();

  DCRTPoly ct0;
  DCRTPoly ct1;

  const shared_ptr<typename DCRTPoly::Params> paramsQ = c[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP =
      cryptoParamsLWE->GetAuxElementParams();

  size_t cipherTowers = c[0].GetParams()->GetParams().size();
  size_t towersToSkip =
      cryptoParamsLWE->GetElementParams()->GetParams().size() - cipherTowers;

  DCRTPoly cTmp;
  cTmp = c[1].Clone();

  uint32_t l = cipherTowers - 1;
  uint32_t alpha = cryptoParamsLWE->GetNumberOfTowersPerDigit();
  uint32_t beta =
      ceil(((double)(l + 1)) /
           alpha);  // The number of digits of the current ciphertext
  // uint32_t digits = cryptoParamsLWE->GetNumberOfDigits();
  if (beta > cryptoParamsLWE->GetNumberOfQPartitions())
    beta = cryptoParamsLWE->GetNumberOfQPartitions();

  vector<DCRTPoly> digitsCTmp(beta);

  // Digit decomposition
  // Zero-padding and split
  uint32_t numTowersLastDigit =
      cryptoParamsLWE->GetQPartition(beta - 1)->GetParams().size();
  for (uint32_t j = 0; j < beta; j++) {
    if (j == beta - 1) {
      auto part = cryptoParamsLWE->GetQPartition(j);
      part->GetParams();

      numTowersLastDigit = cipherTowers - alpha * j;

      vector<NativeInteger> moduli(numTowersLastDigit);
      vector<NativeInteger> roots(numTowersLastDigit);

      for (uint32_t i = 0; i < numTowersLastDigit; i++) {
        moduli[i] = part->GetParams()[i]->GetModulus();
        roots[i] = part->GetParams()[i]->GetRootOfUnity();
      }

      auto params = DCRTPoly::Params(part->GetCyclotomicOrder(), moduli, roots,
                                     {}, {}, 0);

      digitsCTmp[j] =
          DCRTPoly(std::make_shared<typename DCRTPoly::Params>(params),
                   EVALUATION, true);

    } else
      digitsCTmp[j] =
          DCRTPoly(cryptoParamsLWE->GetQPartition(j), Format::EVALUATION, true);

    uint32_t iters = (j == beta - 1) ? numTowersLastDigit : alpha;
    for (uint32_t i = 0; i < iters; i++) {
      if (j * alpha + i <= l) {
        auto tmp = cTmp.GetElementAtIndex(j * alpha + i);
        digitsCTmp[j].SetElementAtIndex(i, tmp);
      }
    }
  }
  // RNS decompose
  for (uint32_t j = 0; j < beta; j++) {
    for (uint32_t i = 0; i < alpha; i++) {
      if (j * alpha + i <= l) {
        auto tmp = digitsCTmp[j].GetElementAtIndex(i).Times(
            cryptoParamsLWE->GetQHatInvModqTable()[j][j * alpha + i]);
        digitsCTmp[j].SetElementAtIndex(i, tmp);
      }
    }
  }

  vector<DCRTPoly> pPartExtC(digitsCTmp.size());
  vector<DCRTPoly> expandedC(digitsCTmp.size());
  for (uint32_t j = 0; j < digitsCTmp.size(); j++) {
    auto tmpDigit = digitsCTmp[j].Clone();

    tmpDigit.SetFormat(Format::COEFFICIENT);

    const shared_ptr<typename DCRTPoly::Params> params =
        cryptoParamsLWE->GetComplementaryPartition(cipherTowers - 1, j);

    pPartExtC[j] = tmpDigit.ApproxSwitchCRTBasis(
        cryptoParamsLWE->GetQPartition(j), params,  // paramsP,
        cryptoParamsLWE->GetPartitionQHatInvModQTable(
            j)[digitsCTmp[j].GetNumOfElements() - 1],
        cryptoParamsLWE->GetPartitionQHatInvModQPreconTable(
            j)[digitsCTmp[j].GetNumOfElements() - 1],
        cryptoParamsLWE->GetPartitionQHatModPTable(cipherTowers - 1)[j],
        cryptoParamsLWE->GetPartitionPrecon(cipherTowers - 1)[j]);

    pPartExtC[j].SetFormat(Format::EVALUATION);

    expandedC[j] =
        DCRTPoly(cTmp.GetExtendedCRTBasis(paramsP), Format::EVALUATION, true);

    for (usint i = 0; i < cipherTowers; i++) {
      if (i / alpha == j)
        expandedC[j].SetElementAtIndex(
            i, digitsCTmp[j].GetElementAtIndex(i % alpha));
      else {
        if (i / alpha < j) {
          expandedC[j].SetElementAtIndex(i, pPartExtC[j].GetElementAtIndex(i));
        } else {
          expandedC[j].SetElementAtIndex(
              i, pPartExtC[j].GetElementAtIndex(i - alpha));
        }
      }
    }

    for (usint i = 0; i < paramsP->GetParams().size(); i++) {
      expandedC[j].SetElementAtIndex(
          i + cipherTowers,
          pPartExtC[j].GetElementAtIndex(i + params->GetParams().size() -
                                         paramsP->GetParams().size()));
    }
  }

  DCRTPoly cTilda0(cTmp.GetExtendedCRTBasis(paramsP), Format::EVALUATION, true);
  DCRTPoly cTilda1(cTmp.GetExtendedCRTBasis(paramsP), Format::EVALUATION, true);

#if USE_FAST_IP
  // This only works for ell = L
  uint32_t j2 = expandedC.size() - 1;
  for (uint32_t j = 0; j < j2; j++) {
    cTilda0 += expandedC[j] * b[j];
    cTilda1 += expandedC[j] * a[j];
  }
  for (usint i = 0; i < expandedC[j2].GetNumOfElements(); i++) {
    usint idx = (i < cipherTowers) ? i : i + towersToSkip;
    cTilda0.SetElementAtIndex(
        i, cTilda0.GetElementAtIndex(i) + expandedC[j2].GetElementAtIndex(i) *
                                              b[j2].GetElementAtIndex(idx));
    cTilda1.SetElementAtIndex(
        i, cTilda1.GetElementAtIndex(i) + expandedC[j2].GetElementAtIndex(i) *
                                              a[j2].GetElementAtIndex(idx));
  }

#else

  for (uint32_t j = 0; j < digitsCTmp.size(); j++) {
    for (usint i = 0; i < expandedC[j].GetNumOfElements(); i++) {
      // The following skips the switch key elements that are missing from the
      // ciphertext
      usint idx = (i < cipherTowers) ? i : i + towersToSkip;
      cTilda0.SetElementAtIndex(
          i, cTilda0.GetElementAtIndex(i) + expandedC[j].GetElementAtIndex(i) *
                                                b[j].GetElementAtIndex(idx));
      cTilda1.SetElementAtIndex(
          i, cTilda1.GetElementAtIndex(i) + expandedC[j].GetElementAtIndex(i) *
                                                a[j].GetElementAtIndex(idx));
    }
  }
#endif

  cTilda0.SetFormat(Format::COEFFICIENT);
  cTilda1.SetFormat(Format::COEFFICIENT);

  DCRTPoly cHat0 = cTilda0.ApproxModDown(
      paramsQ, paramsP, cryptoParamsLWE->GetPInvModQTable(),
      cryptoParamsLWE->GetPInvModQPreconTable(),
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  DCRTPoly cHat1 = cTilda1.ApproxModDown(
      paramsQ, paramsP, cryptoParamsLWE->GetPInvModQTable(),
      cryptoParamsLWE->GetPInvModQPreconTable(),
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  cHat0.SetFormat(Format::EVALUATION);
  cHat1.SetFormat(Format::EVALUATION);

  ct0 = c[0] + cHat0;
  ct1 = cHat1;

  newCiphertext->SetElements({ct0, ct1});
  newCiphertext->SetDepth(cipherText->GetDepth());
  newCiphertext->SetScalingFactor(cipherText->GetScalingFactor());
  newCiphertext->SetLevel(cipherText->GetLevel());

  return newCiphertext;
}

void CtxtStatus(const Ciphertext<DCRTPoly>& c, const string&& msg) {
  std::cout << "\t" << msg;
  std::cout << ": depth " << c->GetDepth() << std::endl;

  std::cout << "\t" << msg;
  std::cout << ": level " << c->GetLevel() << std::endl;

  std::cout << "\t" << msg;
  std::cout << ": ScaleBit " << std::log2(c->GetScalingFactor()) << std::endl;

  std::cout << "\t" << msg;
  std::cout << ": Modulus (bit) "
            << log2(c->GetElements()[0].GetWorkingModulus().ConvertToDouble())
            << std::endl;
}

double TraceSeqHybridOnline(const CryptoContext<DCRTPoly> cc,
                            const LPKeyPair<DCRTPoly>& keys,
                            const uint32_t num_drop, const uint32_t L,
                            const vector<usint>& auto_indices_seq) {
  // Input
  vector<complex<double>> x = {0, 0, 0, 0, 0, 0, 0, 1};
  Plaintext ptxt = cc->MakeCKKSPackedPlaintext(x);

  std::cout << "Input x: " << ptxt << std::endl;
  auto c = cc->Encrypt(keys.publicKey, ptxt);
  vector<complex<double>> x1(cc->GetRingDimension() / 2, 1);
  auto c1 = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(x1));
  for (uint32_t i = 0; i < num_drop; ++i) {
    c = cc->EvalMult(c, c1);
    c = cc->Rescale(c);
  }
  std::cout << "# Mods (Dropped = " << c->GetLevel()
            << ", Available = " << (L + 1) - c->GetLevel() << ")" << std::endl;

  TimeVar t;
  TIC(t);

  long num_seq = auto_indices_seq.size();
#ifdef USE_FAST_AUTOMORPH
  for (long i = 0; i < num_seq; ++i)
    c += EvalAutoKSHybrid(c, auto_indices_seq[i]);
#else
  for (long i = 0; i < num_seq; ++i)
    c += cc->EvalAutomorphism(c, auto_indices_seq[i],
                              cc->GetEvalAutomorphismKeyMap(c->GetKeyTag()));
#endif
  double timing = TOC(t);

  std::cout << "Timing: " << timing << " ms" << std::endl;

  Plaintext result;
  cc->Decrypt(keys.secretKey, c, &result);
  result->SetLength(8);
  std::cout << "\tResult dec = " << result << std::endl;

  auto m = cc->GetCyclotomicOrder();

  size_t num_ans = 0;
  int ans = 1;
  size_t num_expected_ans = 1UL << num_seq;
  if (num_seq == int(log2(m / 2))) {
    ans = 2;
    num_expected_ans = (1UL << (num_seq - 1));
  }

  // See all the slot values
  for (size_t i = 0; i < m / 4; ++i) {
    if (round(result->GetCKKSPackedValue()[i].real()) == ans) num_ans += 1;
  }

  if (num_ans == num_expected_ans)
    std::cout << " correct :) #" << ans << " = " << num_ans
              << " Expected = " << num_expected_ans << std::endl;
  else
    std::cout << "      wrong #" << ans << " = " << num_ans
              << " Expected = " << num_expected_ans << std::endl;

  return timing;
}

int main(int argc, char* argv[]) {
#ifdef USE_FAST_AUTOMORPH
  std::cout << "Fast Automorph" << std::endl;
#endif

  int max_num_thread = omp_get_max_threads();

  argparse::ArgumentParser parser("bench_seq");

  parser.add_argument("logN").help("logN").required().scan<'i', int>();
  parser.add_argument("M")
      .help("iterations for rotations-and-sums")
      .required()
      .scan<'i', int>();
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

  try {
    parser.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << parser;
    std::exit(1);
  }

  auto logN = parser.get<int>("logN");
  auto M = parser.get<int>("M");
  auto depth = parser.get<uint32_t>("L");
  auto ell = parser.get<uint32_t>("ell");
  auto scale = parser.get<uint32_t>("delta");
  auto dnum = parser.get<uint32_t>("dnum");
  auto num_exp = parser.get<int>("T");
  auto num_thread = parser.is_used("-t") ? parser.get<int>("-t") : 1;

  assert(logN >= M);
  assert(depth >= ell);

  omp_set_num_threads(num_thread);

  uint32_t num_drop = depth - ell;

  CryptoContext<DCRTPoly> cc;
  LPKeyPair<DCRTPoly> keys;
  vector<usint> auto_indices_seq;
  TraceSeqHKOffline(logN, M, depth, scale, dnum, cc, keys, auto_indices_seq,
                    false);
  cc->EvalMultKeyGen(keys.secretKey);

  vector<double> timings;
  for (auto i = 0; i < num_exp; ++i)
    timings.push_back(
        TraceSeqHybridOnline(cc, keys, num_drop, depth, auto_indices_seq));

  int num_special_mod = (depth + 1) / dnum;
  ofstream ofs("/tmp/result/TraceRuntime/trace_seq_dynamic_dnum_HK_logN" +
               to_string(logN) + "_M" + to_string(M) + "_L" + to_string(depth) +
               "_ell" + to_string(ell) + "_scale" + to_string(scale) +
               "_numth" + to_string(num_thread) + "_k" +
               to_string(num_special_mod) + "_numexp" + to_string(num_exp) +
               ".csv");

  StatDropFirst trace(timings, {(double)(getPeakRSS()) / (1 << 20)});
  trace.Write(std::cout);
  trace.WriteCSV(ofs);
  ofs.close();

  return 0;
}
