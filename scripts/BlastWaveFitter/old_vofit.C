double BW_V0_pT(double pT, double mass, double Tkin, double betaT, double n_flow, double rho2, double Rx, double Ry, double A, double gIntegralA, double gIntegralB, double sT2,
                double sB2) {
  double Np = A * BW_Spectrum(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry);
  double dNdT = BW_dSpectrum_dT(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A);
  double dNdB = BW_dSpectrum_dBeta(pT, mass, Tkin, betaT, n_flow, rho2, Rx, Ry, A);
  double num = gIntegralA * dNdT * sT2 + gIntegralB * dNdB * sB2;
  double den = Np * std::sqrt(gIntegralA * gIntegralA * sT2 + gIntegralB * gIntegralB * sB2);
  return num / (den);
}
