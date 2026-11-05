#ifndef maskedFaceIds_H
#define maskedFaceIds_H 1

#include "platform.hpp"
#include "mesh.h"
#include <array>
#include <map>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <sstream>

static std::tuple<dlong, occa::memory, dlong, occa::memory, dlong, occa::memory>
maskedFaceIds(mesh_t *mesh,
              dlong mNlocal,
              int nFields,
              dlong offset,
              const std::vector<int>& EToB,
              int dirichletBcTypeId,
              const std::string &ellipticName = "")
{
  dlong Nmasked;
  occa::memory o_maskIds;
  dlong NmaskedLocal;
  occa::memory o_maskIdsLocal;
  dlong NmaskedGlobal;
  occa::memory o_maskIdsGlobal;

  const int Nlocal = (nFields == 1) ? mNlocal : nFields * offset;
  const int largeNumber = 1 << 20;

  std::vector<int> mapB(Nlocal);
 
  // assumes mapB is identical for all nFields
  for (int fld = 0; fld < nFields; fld++) {

    for (dlong e = 0; e < mesh->Nelements; e++) {
      for (int n = 0; n < mesh->Np; n++) {
        mapB[n + e * mesh->Np + fld * offset] = largeNumber;
      }
      for (int f = 0; f < mesh->Nfaces; f++) {
        const int fOffset = fld * mesh->Nelements * mesh->Nfaces;
        int bc = EToB[f + e * mesh->Nfaces + fOffset];
        if (bc > 0) {
          for (int n = 0; n < mesh->Nfp; n++) {
            int fid = mesh->faceNodes[n + f * mesh->Nfp];
            mapB[fid + e * mesh->Np + fld * offset] =
                std::min(bc, mapB[fid + e * mesh->Np + fld * offset]); // DIRICHLET wins
          }
        }
      }
    }

  }

  // STAGE4B: fixed-particle interior nodal Dirichlet constraints.
  // Only active if the Slurm job supplies NEKRS_PHYSALIS_CAGE_FILE.
  // The file is generated from the validated 3D P6 stage2 cage.
  // Values are set by UDF callback BEFORE the elliptic solves.
  if (const char *cagePath = std::getenv("NEKRS_PHYSALIS_CAGE_FILE")) {
    const auto keyName = upperCase(ellipticName);
    const bool pressure = keyName == "FLUID PRESSURE";
    const bool velocity = keyName == "FLUID VELOCITY";
    if (pressure || velocity) {
      using Coord = std::array<long long,3>;
      static std::map<Coord, int> cage;
      static std::string loadedPath;
      if (loadedPath != std::string(cagePath)) {
        cage.clear();
        std::ifstream in(cagePath);
        if (!in) throw std::runtime_error(std::string("Stage4B missing cage mask file: ")+cagePath);
        std::string line;
        while (std::getline(in,line)) {
          if (line.empty() || line[0]=='#') continue;
          std::istringstream s(line);
          char kind=0; double x,y,z;
          if(!(s>>kind>>x>>y>>z) || !(kind=='V'||kind=='B'))
            throw std::runtime_error("Stage4B malformed cage mask line");
          cage[{std::llround(x*1e9),std::llround(y*1e9),std::llround(z*1e9)}] =
              kind == 'B' ? 3 : 1;
        }
        if (cage.empty())
          throw std::runtime_error("Stage5C cage mask file has no constrained coordinates");
        loadedPath = cagePath;
      }
      const auto xyz = mesh->xyzHost();
      const auto &x = std::get<0>(xyz);
      const auto &y = std::get<1>(xyz);
      const auto &z = std::get<2>(xyz);
      dlong localExtra=0;
      for (dlong n=0; n<mesh->Nlocal; ++n) {
        auto it=cage.find({std::llround(double(x[n])*1e9),
                           std::llround(double(y[n])*1e9),
                           std::llround(double(z[n])*1e9)});
        if (it == cage.end()) continue;
        const int bits = it->second;
        if (!(velocity ? (bits & 1) : (bits & 2))) continue;
        for (int fld=0; fld<nFields; ++fld)
          mapB[n + fld * offset] =
            std::min(mapB[n + fld * offset],dirichletBcTypeId);
        localExtra++;
      }
      dlong globalExtra=0;
      MPI_Allreduce(&localExtra,&globalExtra,1,MPI_LONG_LONG,MPI_SUM,
                    platform->comm.mpiComm());
      if (!globalExtra)
        throw std::runtime_error("Stage4B cage mask has zero matching GLL nodes");
      if (platform->comm.mpiRank()==0)
        printf("[STAGE4B CORE] %s: element-local cage matches=%lld\n",
               keyName.c_str(),static_cast<long long>(globalExtra));
    }
  }

  ogsGatherScatterMany(mapB.data(),
                       nFields,
                       offset,
                       ogsInt,
                       ogsMin, // DIRICHLET wins
                       mesh->ogs);

  Nmasked = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong n = 0; n < mesh->Nlocal; n++) {
      if (mapB[n + fld * offset] == largeNumber) {
        mapB[n + fld * offset] = 0;
      } else if (mapB[n + fld * offset] == dirichletBcTypeId) {
        Nmasked++;
      }
    }
  }
  std::vector<dlong> maskIds(Nmasked);

  Nmasked = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong n = 0; n < mesh->Nlocal; n++) {
      if (mapB[n + fld * offset] == dirichletBcTypeId) {
        maskIds[Nmasked++] = n + fld * offset;
      }
    }
  }
  if (Nmasked) {
    o_maskIds = platform->device.malloc<dlong>(Nmasked);
    o_maskIds.copyFrom(maskIds.data());
  }

  NmaskedLocal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong el = 0; el < mesh->NlocalGatherElements; ++el) {
      const dlong elemOffset = mesh->localGatherElementList[el] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == dirichletBcTypeId) {
          NmaskedLocal++;
        }
      }
    }
  }

  std::vector<dlong> localMaskIds(NmaskedLocal);
  NmaskedLocal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong el = 0; el < mesh->NlocalGatherElements; ++el) {
      const dlong elemOffset = mesh->localGatherElementList[el] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == dirichletBcTypeId) {
          localMaskIds[NmaskedLocal++] = n + fld * offset;
        }
      }
    }
  }
  if (NmaskedLocal) {
    o_maskIdsLocal = platform->device.malloc<dlong>(NmaskedLocal);
    o_maskIdsLocal.copyFrom(localMaskIds.data());
  }

  NmaskedGlobal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong eg = 0; eg < mesh->NglobalGatherElements; ++eg) {
      const dlong elemOffset = mesh->globalGatherElementList[eg] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == dirichletBcTypeId) {
          NmaskedGlobal++;
        }
      }
    }
  }
  std::vector<dlong> globalMaskIds(NmaskedGlobal);
  NmaskedGlobal = 0;
  for (int fld = 0; fld < nFields; fld++) {
    for (dlong eg = 0; eg < mesh->NglobalGatherElements; ++eg) {
      const dlong elemOffset = mesh->globalGatherElementList[eg] * mesh->Np;
      for (dlong qp = 0; qp < mesh->Np; qp++) {
        const dlong n = elemOffset + qp;
        if (mapB[n + fld * offset] == dirichletBcTypeId) {
          globalMaskIds[NmaskedGlobal++] = n + fld * offset;
        }
      }
    }
  }
  if (NmaskedGlobal) {
    o_maskIdsGlobal = platform->device.malloc<dlong>(NmaskedGlobal);
    o_maskIdsGlobal.copyFrom(globalMaskIds.data());
  }


  return {Nmasked, o_maskIds, NmaskedLocal, o_maskIdsLocal, NmaskedGlobal, o_maskIdsGlobal};  
}

#endif
