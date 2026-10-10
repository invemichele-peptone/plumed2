#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <random>
#include <string>

// OPES with the exact neighbor list must give the same bias and forces as OPES without it
int main() {
  std::ofstream output("output");
  bool all_passed=true;
  for(const std::string mode : {
        "OPES_METAD","OPES_METAD_EXPLORE"
      }) {
    //adaptive sigma: kernels of different widths, many of them outside the neighbor list
    const std::string opes=" ARG=d1,d2 PACE=5 BARRIER=20 TEMP=300 CUTOFF_SCHEME=PR";
    PLMD::Plumed p[2];
    int precision=8, natoms=3;
    double timestep=0.001, kbt=2.4943387854;
    for(unsigned n=0; n<2; n++) {
      p[n].cmd("setRealPrecision",&precision);
      p[n].cmd("setNatoms",&natoms);
      p[n].cmd("setTimestep",&timestep);
      p[n].cmd("setKbT",&kbt);
      p[n].cmd("init");
      p[n].cmd("readInputLine","d1: DISTANCE ATOMS=1,2 NOPBC");
      p[n].cmd("readInputLine","d2: DISTANCE ATOMS=1,3 NOPBC");
    }
    p[0].cmd("readInputLine",("b: "+mode+opes+" NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL FILE="+mode+"-nlist.kernels").c_str());
    p[1].cmd("readInputLine",("b: "+mode+opes+" FILE="+mode+"-ref.kernels").c_str());
    std::mt19937 rng(12345);
    auto uniform=[&]() {
      return rng()/4294967296.0-0.5;
    };
    std::array<double,9> positions{{0,0,0,1,0,0,0,1,0}};
    std::array<double,9> box{{3,0,0,0,3,0,0,0,3}};
    std::array<double,3> masses{{1,1,1}};
    double max_dbias=0, max_dforce=0;
    for(int step=0; step<8000; ++step) {
      for(unsigned j=3; j<9; j++) {
        positions[j]+=0.03*uniform();
      }
      double bias[2];
      std::array<double,9> forces[2];
      for(unsigned n=0; n<2; n++) {
        std::array<double,9> virial{};
        forces[n].fill(0);
        p[n].cmd("setStep",&step);
        p[n].cmd("setPositions",positions.data());
        p[n].cmd("setMasses",masses.data());
        p[n].cmd("setBox",box.data());
        p[n].cmd("setForces",forces[n].data());
        p[n].cmd("setVirial",virial.data());
        p[n].cmd("calc");
        p[n].cmd("getBias",&bias[n]);
      }
      max_dbias=std::max(max_dbias,std::abs(bias[0]-bias[1])/(1+std::abs(bias[1])));
      for(unsigned j=0; j<9; j++) {
        max_dforce=std::max(max_dforce,std::abs(forces[0][j]-forces[1][j])/(1+std::abs(forces[1][j])));
      }
    }
    const bool passed=max_dbias<1e-9 && max_dforce<1e-9;
    all_passed=all_passed && passed;
    output<<mode<<" exact neighbor list "<<(passed ? "PASS" : "FAIL")<<"\n";
  }
  return all_passed ? 0 : 1;
}
