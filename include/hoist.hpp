// Copyright 2021 Yu Ishimaki

#pragma once

#include <vector>
#include <memory>
#include <utility>

Ciphertext<DCRTPoly> EvalHoistedAutomorphHKSingleTh(
    ConstCiphertext<DCRTPoly> ciphertext,  // original ciphertext to  rotate
    const vector<usint>& autoids,
    map<usint, std::pair<std::vector<DCRTPoly>, std::vector<DCRTPoly>>>&
        inv_evks) {
  Ciphertext<DCRTPoly> newCiphertext = ciphertext->CloneEmpty();

  const shared_ptr<vector<DCRTPoly>> expandedCiphertext =
      ciphertext->GetCryptoContext()->EvalFastRotationPrecompute(ciphertext);

  // 1. Retrieve the automorphism key that corresponds to the auto index.
  std::vector<DCRTPoly> c(2);

  c[0] = ciphertext->GetElements()[0].Clone();
  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(
          ciphertext->GetCryptoParameters());

  const shared_ptr<typename DCRTPoly::Params> paramsQ =
      ciphertext->GetElements()[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP =
      cryptoParamsLWE->GetAuxElementParams();
  size_t cipherTowers = paramsQ->GetParams().size();
  size_t towersToSkip =
      cryptoParamsLWE->GetElementParams()->GetParams().size() - cipherTowers;

  size_t n = autoids.size();
  DCRTPoly c_sum(paramsQ, Format::EVALUATION, true);
  DCRTPoly cTilda0_sum((*expandedCiphertext)[0].GetParams(), Format::EVALUATION,
                       true);
  DCRTPoly cTilda1_sum((*expandedCiphertext)[0].GetParams(), Format::EVALUATION,
                       true);
  for (size_t ii = 0; ii < n; ++ii) {
    const usint autoIndex = autoids[ii];

    auto evalKey = ciphertext->GetCryptoContext()
                       ->GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag())
                       .find(autoIndex)
                       ->second;
#ifdef GET_RAM
    const std::vector<DCRTPoly>& b = evalKey->GetBVector();
    const std::vector<DCRTPoly>& a = evalKey->GetAVector();
#else
    const std::vector<DCRTPoly>& b = inv_evks[autoIndex].first;
    const std::vector<DCRTPoly>& a = inv_evks[autoIndex].second;
#endif

    DCRTPoly ct0;
    DCRTPoly ct1;

    // 1. Retrieve the automorphism key that corresponds to the auto index.

    // Inner Prod with Evk
    // Now, we are preparing the resulting storage mod PQ (filled with 0)
    DCRTPoly cTilda0((*expandedCiphertext)[0].GetParams(), Format::EVALUATION,
                     true);
    DCRTPoly cTilda1((*expandedCiphertext)[0].GetParams(), Format::EVALUATION,
                     true);

    for (uint32_t j = 0; j < expandedCiphertext->size(); j++) {
      for (usint i = 0; i < (*expandedCiphertext)[j].GetNumOfElements(); i++) {
        usint idx = (i < cipherTowers) ? i : i + towersToSkip;
        cTilda0.SetElementAtIndex(
            i, b[j].GetElementAtIndex(idx).FMA(
                   (*expandedCiphertext)[j].GetElementAtIndex(i),
                   cTilda0.GetElementAtIndex(i)));
        cTilda1.SetElementAtIndex(
            i, a[j].GetElementAtIndex(idx).FMA(
                   (*expandedCiphertext)[j].GetElementAtIndex(i),
                   cTilda1.GetElementAtIndex(i)));
      }
    }

#ifdef TRACE_DEBUG
    // We see that the last k RNS components are relative to P
    uint32_t num_digits = expandedCiphertext->size();
    cout << " k                             = " << paramsP->GetParams().size()
         << endl;
    cout << " cipherTows (ell+1)            = " << cipherTowers << endl;
    cout << " k + ell + 1                   = "
         << cipherTowers + paramsP->GetParams().size() << endl;
    cout << " TowsToSkip (#Drop = L - ell ) = " << towersToSkip << endl;
    cout << " # Digtit                      = " << num_digits << endl;
    cout << " # RSN on EVK                  = "
         << towersToSkip + cipherTowers + paramsP->GetParams().size() << endl;
    for (uint32_t j = 0; j < num_digits; j++) {
      cout << "\tDigit No: " << j << endl;
      usint num_rns_on_digit = (*expandedCiphertext)[j].GetNumOfElements();
      cout << "\t\t #RNS Components on digit: " << num_rns_on_digit << endl;
      for (usint i = 0; i < num_rns_on_digit; i++) {
        usint idx = (i < cipherTowers) ? i : i + towersToSkip;
        cout << "\t\t EVk RNS-id to refer     : " << idx << endl;
      }
    }
#endif

#ifdef USE_FAST_AUTOMORPH
    vector<usint> perm;
    GenAutomorphTable(cTilda0.GetRingDimension(), autoIndex, perm);
    c_sum += std::move(ciphertext->GetElements()[0].Permute(perm));
    cTilda0_sum += std::move(cTilda0.Permute(perm));
    cTilda1_sum += std::move(cTilda1.Permute(perm));
#else
    c_sum += std::move(
        ciphertext->GetElements()[0].AutomorphismTransform(autoIndex));
    cTilda0_sum += std::move(cTilda0.AutomorphismTransform(autoIndex));
    cTilda1_sum += std::move(cTilda1.AutomorphismTransform(autoIndex));
#endif
  }

  cTilda0_sum.SetFormat(Format::COEFFICIENT);
  cTilda1_sum.SetFormat(Format::COEFFICIENT);

  DCRTPoly cHat0 = cTilda0_sum.ApproxModDown(
      paramsQ, paramsP,
      cryptoParamsLWE->GetPInvModQTable(),        // P^{-1} mod q_i
      cryptoParamsLWE->GetPInvModQPreconTable(),  // Barrett Const to multiply
                                                  // with P^{-1} mod q_i
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  DCRTPoly cHat1 = cTilda1_sum.ApproxModDown(
      paramsQ, paramsP, cryptoParamsLWE->GetPInvModQTable(),
      cryptoParamsLWE->GetPInvModQPreconTable(),
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  cHat0.SetFormat(Format::EVALUATION);
  cHat1.SetFormat(Format::EVALUATION);

  DCRTPoly ct0;
  DCRTPoly ct1;
  ct0 = c_sum + cHat0;
  ct1 = cHat1;

  newCiphertext->SetElements({ct0, ct1});
  newCiphertext->SetDepth(ciphertext->GetDepth());
  newCiphertext->SetLevel(ciphertext->GetLevel());
  newCiphertext->SetScalingFactor(ciphertext->GetScalingFactor());

  return newCiphertext;
}

Ciphertext<DCRTPoly> EvalHoistedAutomorph(
    ConstCiphertext<DCRTPoly> ciphertext,  // original ciphertext to  rotate
    const vector<usint>& autoids,
    map<usint, std::pair<std::vector<DCRTPoly>, std::vector<DCRTPoly>>>&
        inv_evks) {
  Ciphertext<DCRTPoly> newCiphertext = ciphertext->CloneEmpty();

  const shared_ptr<vector<DCRTPoly>> expandedCiphertext =
      ciphertext->GetCryptoContext()->EvalFastRotationPrecompute(ciphertext);

  // 1. Retrieve the automorphism key that corresponds to the auto index.
  std::vector<DCRTPoly> c(2);

  c[0] = ciphertext->GetElements()[0].Clone();
  const shared_ptr<LPCryptoParametersCKKS<DCRTPoly>> cryptoParamsLWE =
      std::dynamic_pointer_cast<LPCryptoParametersCKKS<DCRTPoly>>(
          ciphertext->GetCryptoParameters());

  const shared_ptr<typename DCRTPoly::Params> paramsQ =
      ciphertext->GetElements()[0].GetParams();
  const shared_ptr<typename DCRTPoly::Params> paramsP =
      cryptoParamsLWE->GetAuxElementParams();
  const shared_ptr<typename DCRTPoly::Params> paramsPQ =
      (*expandedCiphertext)[0].GetParams();
  size_t cipherTowers = paramsQ->GetParams().size();
  size_t towersToSkip =
      cryptoParamsLWE->GetElementParams()->GetParams().size() - cipherTowers;

  size_t n = autoids.size();
  vector<DCRTPoly> c_vec(n);
  vector<DCRTPoly> cTilda0_vec(n);
  vector<DCRTPoly> cTilda1_vec(n);
#pragma omp parallel for
  for (size_t ii = 0; ii < n; ++ii) {
    const usint autoIndex = autoids[ii];

    auto evalKey = ciphertext->GetCryptoContext()
                       ->GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag())
                       .find(autoIndex)
                       ->second;
#ifdef GET_RAM
    const std::vector<DCRTPoly>& b = evalKey->GetBVector();
    const std::vector<DCRTPoly>& a = evalKey->GetAVector();
#else
    const std::vector<DCRTPoly>& b = inv_evks[autoIndex].first;
    const std::vector<DCRTPoly>& a = inv_evks[autoIndex].second;
#endif

    DCRTPoly ct0;
    DCRTPoly ct1;

    // 1. Retrieve the automorphism key that corresponds to the auto index.

    // Inner Prod with Evk
    // Now, we are preparing the resulting storage mod PQ (filled with 0)
    DCRTPoly cTilda0(paramsPQ, Format::EVALUATION, true);
    DCRTPoly cTilda1(paramsPQ, Format::EVALUATION, true);

    for (uint32_t j = 0; j < expandedCiphertext->size(); j++) {
      for (usint i = 0; i < (*expandedCiphertext)[j].GetNumOfElements(); i++) {
        usint idx = (i < cipherTowers) ? i : i + towersToSkip;
        cTilda0.SetElementAtIndex(
            i, b[j].GetElementAtIndex(idx).FMA(
                   (*expandedCiphertext)[j].GetElementAtIndex(i),
                   cTilda0.GetElementAtIndex(i)));
        cTilda1.SetElementAtIndex(
            i, a[j].GetElementAtIndex(idx).FMA(
                   (*expandedCiphertext)[j].GetElementAtIndex(i),
                   cTilda1.GetElementAtIndex(i)));
      }
    }

#ifdef USE_FAST_AUTOMORPH
    vector<usint> perm;
    GenAutomorphTable(cTilda0.GetRingDimension(), autoIndex, perm);
    c_vec[ii] = std::move(ciphertext->GetElements()[0].Permute(perm));
    cTilda0_vec[ii] = std::move(cTilda0.Permute(perm));
    cTilda1_vec[ii] = std::move(cTilda1.Permute(perm));
#else
    c_vec[ii] = std::move(
        ciphertext->GetElements()[0].AutomorphismTransform(autoIndex));
    cTilda0_vec[ii] = std::move(cTilda0.AutomorphismTransform(autoIndex));
    cTilda1_vec[ii] = std::move(cTilda1.AutomorphismTransform(autoIndex));
#endif
  }
  DCRTPoly c_sum(paramsQ, Format::EVALUATION, true);
  DCRTPoly cTilda0_sum(paramsPQ, Format::EVALUATION, true);
  DCRTPoly cTilda1_sum(paramsPQ, Format::EVALUATION, true);

#if 0
#pragma omp parallel for reduction(+ : c_sum)
  for (size_t ii = 0; ii < n; ++ii)
     c_sum += c_vec[ii];

#pragma omp parallel for reduction(+ : cTilda0_sum, cTilda1_sum)
  for (size_t ii = 0; ii < n; ++ii) {
     cTilda0_sum += cTilda0_vec[ii];
     cTilda1_sum += cTilda1_vec[ii];
  }
#else
  for (size_t ii = 0; ii < n; ++ii) {
    c_sum += c_vec[ii];
    cTilda0_sum += cTilda0_vec[ii];
    cTilda1_sum += cTilda1_vec[ii];
  }
#endif

  cTilda0_sum.SetFormat(Format::COEFFICIENT);
  cTilda1_sum.SetFormat(Format::COEFFICIENT);

  DCRTPoly cHat0 = cTilda0_sum.ApproxModDown(
      paramsQ, paramsP,
      cryptoParamsLWE->GetPInvModQTable(),        // P^{-1} mod q_i
      cryptoParamsLWE->GetPInvModQPreconTable(),  // Barrett Const to multiply
                                                  // with P^{-1} mod q_i
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  DCRTPoly cHat1 = cTilda1_sum.ApproxModDown(
      paramsQ, paramsP, cryptoParamsLWE->GetPInvModQTable(),
      cryptoParamsLWE->GetPInvModQPreconTable(),
      cryptoParamsLWE->GetPHatInvModPTable(),
      cryptoParamsLWE->GetPHatInvModPPreconTable(),
      cryptoParamsLWE->GetPHatModQTable(),
      cryptoParamsLWE->GetModBarretPreconQTable());

  cHat0.SetFormat(Format::EVALUATION);
  cHat1.SetFormat(Format::EVALUATION);

  DCRTPoly ct0;
  DCRTPoly ct1;
  ct0 = c_sum + cHat0;
  ct1 = cHat1;

  newCiphertext->SetElements({ct0, ct1});
  newCiphertext->SetDepth(ciphertext->GetDepth());
  newCiphertext->SetLevel(ciphertext->GetLevel());
  newCiphertext->SetScalingFactor(ciphertext->GetScalingFactor());

  return newCiphertext;
}
