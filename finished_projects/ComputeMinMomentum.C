// ComputeMinMomentum.C

#include <iostream>
#include <cmath>


const Double_t e_charge           = 1.602e-19;   // Coulomb
const Double_t GeV_c_to_kg_m_s    = 5.344e-19;   // [kg·m/s] per (GeV/c)

//=== Compute momentum p [GeV/c] required to achieve curvature radius R [cm] in B-field ===
Double_t MomentumFromRadius(Double_t R_cm, Double_t B) {
    Double_t R_m = R_cm * 1e-2;
    // p_SI [kg·m/s] = q * B * R
    Double_t p_SI = e_charge * B * R_m;
    // Convert to GeV/c: p_GeV = p_SI / GeV_c_to_kg_m_s
    return p_SI / GeV_c_to_kg_m_s;
}

void ComputeMinMomentum() {
    // Detector radii in cm
    const Double_t rTPC = 85.0;
    const Double_t rTOF = 370.0;
    const Double_t Bfield = 0.5;  // Tesla

    Double_t pMinTPC = MomentumFromRadius(rTPC, Bfield);
    Double_t pMinTOF = MomentumFromRadius(rTOF, Bfield);

    std::cout << "Minimum momentum to reach TPC (R = "
              << rTPC << " cm): " << pMinTPC << " GeV/c" << std::endl;
    std::cout << "Minimum momentum to reach TOF (R = "
              << rTOF << " cm): " << pMinTOF << " GeV/c" << std::endl;
}