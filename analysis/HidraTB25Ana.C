#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <iostream>
#include <array>
#include <stdint.h>
#include <string>
#include <fstream>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <sys/file.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "HidraGeo.h"
#include <map>
#include <TRandom3.h>
#include <cctype>

using json = nlohmann::json;

// -----------------------------------------------------------------------------
// Global constants and parameters
// -----------------------------------------------------------------------------

const double MinVisiblePheS = 0.0;
const double MinVisiblePheC = 0.0;

const unsigned int grouping = 8;
const unsigned int EventDisplayEvery = 100;
const unsigned int FersMultiplicity = 2;

const double chi = 0.38;
const double sciPheGeV = 119.001; // tb24
const double cerPheGeV = 29.4;  // tb24
//const double sciPheGeV = 178.501;
//const double cerPheGeV = 43;

const double NofSipmCells_sci = 7772;
const double NofSipmCells_cer = 3443;
const double sci_pde = 0.22;
const double cer_pde = 0.38;

const double elcont = 1.005;
const double picont = 1.028;


enum class NoiseCorrelationMode {
  UncorrelatedByChannel,
  CorrelatedWithinFers
};

enum class NoiseDistributionMode {
  Gaussian,
  LogNormal
};

/////////////////////////////////////////////////////////
// Paremeters for changing simulation behaviour
const bool ApplySaturation = true;
const bool ApplyTbNoise = false;
const bool AddPedestalToEnergyContribution = false;
const bool ApplyPedestalSubtraction = false;
const double ScalePedestalSubtractionFactorS = 1.0;  // 1.09 for correlated
const double ScalePedestalSubtractionFactorC = 1.0; // 1.055 for correlated
//const double ScalePedestalSubtractionFactorS = 1.1;  // 1.09 for correlated
//const double ScalePedestalSubtractionFactorC = 1.055; // 1.055 for correlated
const double FersThresholdScaleFactor = 1.;
const NoiseCorrelationMode NoiseMode = NoiseCorrelationMode::UncorrelatedByChannel; // use uncorrelated noise
//const NoiseCorrelationMode NoiseMode = NoiseCorrelationMode::CorrelatedWithinFers; // use correlated noise within FERS

//const NoiseDistributionMode NoiseDistribution = NoiseDistributionMode::LogNormal;
const NoiseDistributionMode NoiseDistribution = NoiseDistributionMode::Gaussian;
const bool printSmearingLog = false;
const bool printFersLog = false;

///////////////////////////////////////////////////////////

const unsigned int ChannelsPerFers = 64;
const unsigned int TotalFers = 16;
const unsigned int TotalCalibChannels = ChannelsPerFers * TotalFers;






using NoiseGeneratorBank = std::vector<TRandom3>;

// -----------------------------------------------------------------------------
// Utility helpers
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// 3) Add deterministic RNG builders
//    Note: ROOT treats seed=0 specially, so channel seed uses idx+1.
// -----------------------------------------------------------------------------

unsigned int MakeChannelNoiseSeed(unsigned int calibIndex){return calibIndex + 1U;}
unsigned int MakeFersNoiseSeed(unsigned int fersId){return fersId;}

NoiseGeneratorBank BuildChannelNoiseGenerators()
{
  NoiseGeneratorBank generators(TotalCalibChannels);
  for (unsigned int calibIndex = 0; calibIndex < TotalCalibChannels; ++calibIndex) {
    generators[calibIndex].SetSeed(MakeChannelNoiseSeed(calibIndex));
  }
  return generators;
}

NoiseGeneratorBank BuildFersNoiseGenerators()
{
  NoiseGeneratorBank generators(TotalFers + 1U);
  for (unsigned int fersId = 1; fersId <= TotalFers; ++fersId) {
    generators[fersId].SetSeed(MakeFersNoiseSeed(fersId));
  }
  return generators;
}


double apply_sipm_saturation_from_pe(double pe,
                                     double phePerGeV,
                                     double nCells)
{
  if (pe <= 0.0 || phePerGeV <= 0.0 || nCells <= 0.0) {
    return 0.0;
  }

  const double firedCells = -nCells * std::expm1(-pe / nCells);
  return firedCells / phePerGeV;
}


// FERS activation values (pedestal subtracted)
std::map<unsigned int, double> fers_to_thr_map = {
    { 1, 0.17 },
    { 2, 0.20 },
    { 3, 0.30 },
    { 4, 0.28 },
    { 5, 0.105},
    { 6, 0.095},
    { 7, 0.08 },
    { 8, 0.08 },
    { 9, 0.05 },
    {10, 0.07 },
    {11, 0.05 },
    {12, 0.025},
    {13, 0.06 },
    {14, 0.05 },
    {15, 0.12 },
    {16, 0.06 }
};

/*
// FERS activation values ( non-pedestal subtracted)
std::map<unsigned int, double> fers_to_thr_map = {
{1, 0.2225},
{2, 0.24},
{3, 0.3508},
{4, 0.3508},
{5, 0.1641},
{6, 0.1525},
{7, 0.1405},
{8, 0.0942},
{9, 0.0592},
{10, 0.0592},
{11, 0.0358},
{12, 0.0825},
{13, 0.0825},
{14, 0.0825},
{15, 0.1413},
{16, 0.1413}
};*/

struct PedestalHgEntry {
  double medianAdc = 0.0;
  double rmsAdc = 0.0;
};

using PedestalHgVector = std::vector<PedestalHgEntry>;
using AdcToGeVVector = std::vector<double>;

struct SipmMapEntry {
  int calibIndex = -1;
  int boardID = -1;
  std::string type;
  int row = -1;
  int column = -1;
  double x = 0.0;
  double y = 0.0;
  int fersId = -1;
  int ch = -1;
  std::string module_name;
  double x_local = 0.0;
  double y_local = 0.0;

  double pedestalMedianHgAdc = 0.0;
  double pedestalRmsHgAdc = 0.0;
  double pedestalMedianHgGeV = 0.0;
  double pedestalRmsHgGeV = 0.0;
};

struct EventDisplayHit {
  uint64_t key = 0;
  double x = 0.0;
  double y = 0.0;
  double rawPhe = 0.0;
  double signal = 0.0;
};

struct GaussianFitSummary {
  double mean = std::numeric_limits<double>::quiet_NaN();
  double meanErr = std::numeric_limits<double>::quiet_NaN();
  double rms = std::numeric_limits<double>::quiet_NaN();
  double rmsErr = std::numeric_limits<double>::quiet_NaN();
};

struct HistMomentSummary {
  double mean = std::numeric_limits<double>::quiet_NaN();
  double meanErr = std::numeric_limits<double>::quiet_NaN();
  double rms = std::numeric_limits<double>::quiet_NaN();
  double rmsErr = std::numeric_limits<double>::quiet_NaN();
};

using SipmLookup = std::unordered_map<std::string, SipmMapEntry>;
using FersKey = uint64_t;

// Constrct a unique key for each FERS+channel combination for easier lookup and accumulation
struct ChannelKey {
  FersKey fersKey = 0;
  int ch = -1;

  bool operator==(const ChannelKey& other) const
  {
    return fersKey == other.fersKey && ch == other.ch;
  }
};

struct ChannelKeyHash {
  std::size_t operator()(const ChannelKey& key) const noexcept
  {
    const std::size_t h1 = std::hash<FersKey>{}(key.fersKey);
    const std::size_t h2 = std::hash<int>{}(key.ch);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6U) + (h1 >> 2U));
  }
};

struct ChannelAccumulatedSignal {
  FersKey fersKey = 0;
  int fersId = -1;
  int ch = -1;
  int calibIndex = -1;
  bool isSci = false;
  double threshold = 0.0;
  double thresholdSignal = 0.0;
  double outputSignal = 0.0;
  double pedestalMean = 0.0;
  double noiseSigma = 0.0;
  double smearedThresholdSignal = 0.0;
  double smearedOutputSignal = 0.0;
};


using ChannelSignalMap =
  std::unordered_map<ChannelKey, ChannelAccumulatedSignal, ChannelKeyHash>;

std::string make_sipm_key(const std::string& tower,
                          const std::string& type,
                          int row,
                          int column)
{
  std::ostringstream os;
  os << tower << '|' << type << '|' << row << '|' << column;
  return os.str();
}

FersKey make_fers_key(int boardID, int fersId)
{
  return (static_cast<FersKey>(static_cast<uint32_t>(boardID)) << 32) |
         static_cast<uint32_t>(fersId);
}

int board_id_from_fers_key(FersKey key)
{
  return static_cast<int>(static_cast<uint32_t>(key >> 32));
}

int grouped_column(int colID_original, unsigned int groupingValue)
{
  return static_cast<int>(colID_original / groupingValue);
}

double GetFersThreshold(int fersId)
{
  if (fersId < 0) {
    throw std::runtime_error("Invalid FERS id: " + std::to_string(fersId));
  }

  const auto it = fers_to_thr_map.find(static_cast<unsigned int>(fersId));
  if (it == fers_to_thr_map.end()) {
    throw std::runtime_error("Missing threshold for FERS id: " + std::to_string(fersId));
  }
  
  //return it->second;
  return it->second * FersThresholdScaleFactor; // try scaling to consider pedestal non-subtracted values
}

// -----------------------------------------------------------------------------
// Calibration loaders
// -----------------------------------------------------------------------------

PedestalHgVector LoadHgPedestalsFromJson(const std::string& jsonPath)
{
  std::ifstream in(jsonPath);
  if (!in) {
    throw std::runtime_error("Cannot open pedestal JSON: " + jsonPath);
  }

  json j;
  in >> j;

  if (!j.is_object()) {
    throw std::runtime_error("Pedestal JSON must be an object: " + jsonPath);
  }

  std::size_t maxIndex = 0;
  for (auto it = j.begin(); it != j.end(); ++it) {
    const int idx = std::stoi(it.key());
    if (idx < 0) {
      throw std::runtime_error("Negative pedestal index in: " + jsonPath);
    }
    maxIndex = std::max(maxIndex, static_cast<std::size_t>(idx));
  }

  PedestalHgVector values(maxIndex + 1);
  std::vector<bool> found(maxIndex + 1, false);

  for (auto it = j.begin(); it != j.end(); ++it) {
    const int idx = std::stoi(it.key());
    const json& node = it.value();

    values[static_cast<std::size_t>(idx)].medianAdc = node.value("median_HG", 0.0);
    values[static_cast<std::size_t>(idx)].rmsAdc = node.value("iqr_eff_HG", 0.0);
    found[static_cast<std::size_t>(idx)] = true;
  }

  for (std::size_t i = 0; i < found.size(); ++i) {
    if (!found[i]) {
      std::ostringstream os;
      os << "Missing HG pedestal entry for channel " << i
         << " in " << jsonPath;
      throw std::runtime_error(os.str());
    }
  }

  return values;
}

AdcToGeVVector LoadAdcToGeVFromJson(const std::string& jsonPath)
{
  std::ifstream in(jsonPath);
  if (!in) {
    throw std::runtime_error("Cannot open ADC->GeV JSON: " + jsonPath);
  }

  json j;
  in >> j;

  if (!j.is_array()) {
    throw std::runtime_error("ADC->GeV JSON must be an array: " + jsonPath);
  }

  AdcToGeVVector values;
  values.reserve(j.size());

  for (const auto& item : j) {
    values.push_back(item.get<double>());
  }

  if (values.empty()) {
    throw std::runtime_error("ADC->GeV JSON is empty: " + jsonPath);
  }

  return values;
}

// not used in current version
std::size_t BuildCalibrationIndex(const SipmMapEntry& entry)
{
  if (entry.fersId <= 0 || entry.fersId > static_cast<int>(TotalFers)) {
    std::ostringstream os;
    os << "Invalid FERS id " << entry.fersId;
    throw std::runtime_error(os.str());
  }

  if (entry.ch < 0 || entry.ch >= static_cast<int>(ChannelsPerFers)) {
    std::ostringstream os;
    os << "Invalid channel " << entry.ch
       << " for FERS " << entry.fersId;
    throw std::runtime_error(os.str());
  }

  const std::size_t idx =
    (static_cast<std::size_t>(entry.fersId) - 1U) * ChannelsPerFers +
    static_cast<std::size_t>(entry.ch);

  if (idx >= TotalCalibChannels) {
    std::ostringstream os;
    os << "Calibration index out of range: " << idx
       << " from FERS " << entry.fersId
       << ", ch " << entry.ch;
    throw std::runtime_error(os.str());
  }

  return idx;
}

void AttachHgPedestalNoiseToSipmMap(SipmLookup& lookup,
                                    const PedestalHgVector& pedestalsHg,
                                    const AdcToGeVVector& adcToGeV)
{
  if (pedestalsHg.size() != TotalCalibChannels) {
    std::ostringstream os;
    os << "Expected " << TotalCalibChannels
       << " HG pedestal entries, got " << pedestalsHg.size();
    throw std::runtime_error(os.str());
  }

  if (adcToGeV.size() != TotalCalibChannels) {
    std::ostringstream os;
    os << "Expected " << TotalCalibChannels
       << " ADC->GeV entries, got " << adcToGeV.size();
    throw std::runtime_error(os.str());
  }

  for (auto& [key, entry] : lookup) {
    if (entry.calibIndex < 0 || entry.calibIndex >= static_cast<int>(TotalCalibChannels)) {
      std::ostringstream os;
      os << "Invalid calibIndex " << entry.calibIndex
         << " for map entry " << key;
      throw std::runtime_error(os.str());
    }

    const std::size_t idx = static_cast<std::size_t>(entry.calibIndex);

    entry.pedestalMedianHgAdc = pedestalsHg[idx].medianAdc;
    entry.pedestalRmsHgAdc = pedestalsHg[idx].rmsAdc;
    entry.pedestalMedianHgGeV = pedestalsHg[idx].medianAdc * adcToGeV[idx];
    entry.pedestalRmsHgGeV = pedestalsHg[idx].rmsAdc * adcToGeV[idx];
  }
}



// -----------------------------------------------------------------------------
// Signal accumulation in each channel and smearing
// Apply saturation to a single signal path for both threshold and output consistency
// -----------------------------------------------------------------------------
void AccumulateChannelSignal(ChannelSignalMap& channelSignals,
                             const SipmMapEntry& info,
                             double signalContribution)
{
  if (info.boardID < 0 || info.fersId < 0 || info.ch < 0 || info.calibIndex < 0) {
    return;
  }

  const FersKey fersKey = make_fers_key(info.boardID, info.fersId);
  const ChannelKey channelKey{fersKey, info.ch};

  auto [it, inserted] = channelSignals.emplace(
    channelKey,
    ChannelAccumulatedSignal{
      fersKey,
      info.fersId,
      info.ch,
      info.calibIndex,
      info.type == "S",
      GetFersThreshold(info.fersId),
      0.0,
      0.0,
      info.pedestalMedianHgGeV,
      info.pedestalRmsHgGeV,
      0.0,
      0.0
    }
  );

  if (!inserted && it->second.calibIndex != info.calibIndex) {
    std::ostringstream os;
    os << "Inconsistent calibIndex for accumulated channel: existing="
       << it->second.calibIndex << ", new=" << info.calibIndex
       << " for FERS " << info.fersId << ", ch " << info.ch;
    throw std::runtime_error(os.str());
  }

  // Use single signal for both threshold and output to ensure consistency
  it->second.thresholdSignal += signalContribution;
  it->second.outputSignal += signalContribution;
}

double DrawUnitGaussianNoiseForChannel(const ChannelAccumulatedSignal& channel,
                                       NoiseGeneratorBank& channelNoiseGenerators)
{
  if (channel.calibIndex < 0 ||
      channel.calibIndex >= static_cast<int>(channelNoiseGenerators.size())) {
    std::ostringstream os;
    os << "Invalid calibIndex " << channel.calibIndex
       << " while drawing channel noise";
    throw std::runtime_error(os.str());
  }

  return channelNoiseGenerators[static_cast<std::size_t>(channel.calibIndex)].Gaus(0.0, 1.0);
}

double DrawUnitGaussianNoiseForFers(
  const ChannelAccumulatedSignal& channel,
  NoiseGeneratorBank& fersNoiseGenerators,
  std::unordered_map<FersKey, double>& fersEventNoiseCache)
{
  const auto it = fersEventNoiseCache.find(channel.fersKey);
  if (it != fersEventNoiseCache.end()) {
    return it->second;
  }

  if (channel.fersId <= 0 ||
      channel.fersId >= static_cast<int>(fersNoiseGenerators.size()))
      {
        throw std::runtime_error(
          "Invalid FERS id " + std::to_string(channel.fersId) +
          " while drawing correlated FERS noise"
        ); 
      }


  const double sharedUnitNoise =
    fersNoiseGenerators[static_cast<std::size_t>(channel.fersId)].Gaus(0.0, 1.0);

  fersEventNoiseCache[channel.fersKey] = sharedUnitNoise;
  return sharedUnitNoise;
}


struct LogNormalParameters {
  double mu = 0.0;
  double sigma = 0.0;
};

LogNormalParameters GetLogNormalParametersFromMeanRms(double mean, double rms)
{
  if (mean <= 0.0 || rms <= 0.0) {
    return {0.0, 0.0};
  }

  const double variance = rms * rms;

  LogNormalParameters pars;
  pars.sigma = std::sqrt(std::log(1.0 + variance / (mean * mean)));
  pars.mu = std::log(mean) - 0.5 * pars.sigma * pars.sigma;

  return pars;
}

double DrawLogNormalWithMeanRms(TRandom3& rng, double mean, double rms)
{
  if (mean <= 0.0) {
    return 0.0;
  }

  if (rms <= 0.0) {
    return mean;
  }

  const LogNormalParameters pars =
    GetLogNormalParametersFromMeanRms(mean, rms);

  return std::exp(rng.Gaus(pars.mu, pars.sigma));
}



// -----------------------------------------------------------------------------
//Smearing: noise first, activation decision first, no pedestal subtraction here
void SmearAndCountChannelsOverThreshold(
  ChannelSignalMap& channelSignals,
  std::unordered_map<FersKey, unsigned int>& fersChannelsOverThreshold,
  NoiseCorrelationMode noiseMode,
  NoiseGeneratorBank& channelNoiseGenerators,
  NoiseGeneratorBank& fersNoiseGenerators,
  std::unordered_map<FersKey, double>& fersEventNoiseCache) 
{

  for (auto& [channelKey, channel] : channelSignals) {

    double pedestalContribution = 0.0;

    if (AddPedestalToEnergyContribution) {

      if (!ApplyTbNoise || channel.noiseSigma <= 0.0) {
        pedestalContribution = channel.pedestalMean;
      }

      else if (NoiseDistribution == NoiseDistributionMode::Gaussian) {
        double unitNoise = 0.0;

        if (noiseMode == NoiseCorrelationMode::UncorrelatedByChannel) {
          unitNoise = DrawUnitGaussianNoiseForChannel(
            channel,
            channelNoiseGenerators
          );
        } else if (noiseMode == NoiseCorrelationMode::CorrelatedWithinFers) {
          unitNoise = DrawUnitGaussianNoiseForFers(
            channel,
            fersNoiseGenerators,
            fersEventNoiseCache
          );
        }

        pedestalContribution =
          channel.pedestalMean + unitNoise * channel.noiseSigma;
      }

      else if (NoiseDistribution == NoiseDistributionMode::LogNormal) {
        if (noiseMode == NoiseCorrelationMode::UncorrelatedByChannel) {
          pedestalContribution = DrawLogNormalWithMeanRms(
            channelNoiseGenerators[static_cast<std::size_t>(channel.calibIndex)],
            channel.pedestalMean,
            channel.noiseSigma
          );
        } else if (noiseMode == NoiseCorrelationMode::CorrelatedWithinFers) {
          pedestalContribution = DrawLogNormalWithMeanRms(
            fersNoiseGenerators[static_cast<std::size_t>(channel.fersId)],
            channel.pedestalMean,
            channel.noiseSigma
          );
        }
      }
    }

    channel.smearedThresholdSignal =
      channel.thresholdSignal + pedestalContribution;

    channel.smearedOutputSignal =
      channel.outputSignal + pedestalContribution;





    if (channel.smearedThresholdSignal > channel.threshold) {
      ++fersChannelsOverThreshold[channel.fersKey];
    }




  }
}




bool IsActivatedFers(const std::unordered_map<FersKey, unsigned int>& fersChannelsOverThreshold,
                     FersKey key)
{
  const auto it = fersChannelsOverThreshold.find(key);
  return it != fersChannelsOverThreshold.end() &&
         it->second >= FersMultiplicity;
}

double SumActivatedSmearedChannelOutput(
  const ChannelSignalMap& channelSignals,
  const std::unordered_map<FersKey, unsigned int>& fersChannelsOverThreshold)
{
  double total = 0.0;

  for (const auto& [channelKey, channel] : channelSignals) {
    if (!IsActivatedFers(fersChannelsOverThreshold, channel.fersKey)) {
      continue;
    }

    total += channel.smearedOutputSignal;
  }

  return total;
}

unsigned int CountActivatedFers(
  const std::unordered_map<FersKey, unsigned int>& fersChannelsOverThreshold)
{
  unsigned int activated = 0;

  for (const auto& [fersKey, nChannels] : fersChannelsOverThreshold) {
    if (nChannels >= FersMultiplicity) {
      ++activated;
    }
  }

  return activated;
}


// -----------------------------------------------------------------------------
// 4) subtract pedestal only for activated FERS
void SubtractPedestalMedianFromActivatedChannels(
  ChannelSignalMap& channelSignals,
  const std::unordered_map<FersKey, unsigned int>& fersChannelsOverThreshold)
{
  for (auto& [channelKey, channel] : channelSignals) {
    if (!IsActivatedFers(fersChannelsOverThreshold, channel.fersKey)) {
      continue;
    }

    if(ApplyPedestalSubtraction)
    {
      const double scale = channel.isSci
        ? ScalePedestalSubtractionFactorS
        : ScalePedestalSubtractionFactorC;

      //channel.smearedOutputSignal -= channel.pedestalMean;
      if(channel.outputSignal>0 and printSmearingLog){
      std::cout << "Subtracting pedestal median for FERS " << channel.fersId
                << ", ch " << channel.ch
                << ", true signal: " << channel.outputSignal
                << "\t smeared signal: " << channel.smearedOutputSignal
                << "\t pedestal median: " << channel.pedestalMean
                << "\t scale: " << scale << std::endl;}
      channel.smearedOutputSignal -= channel.pedestalMean * scale;
    }
    // Allow negative values: downward fluctuations naturally subtract from total energy
  }
}

void FillActivatedChannelOutputArrays(
  const ChannelSignalMap& channelSignals,
  const std::unordered_map<FersKey, unsigned int>& fersChannelsOverThreshold,
  std::vector<double>& activatedChannelOutput,
  std::vector<double>& activatedChannelOutputTrue,
  std::vector<double>& activatedChannelOutputSmeared)
{
  if (activatedChannelOutput.size() != activatedChannelOutputTrue.size() ||
      activatedChannelOutput.size() != activatedChannelOutputSmeared.size()) {
    throw std::runtime_error("Activated channel output arrays must have the same size");
  }

  for (const auto& [channelKey, channel] : channelSignals) {
    if (channel.calibIndex < 0 ||
        channel.calibIndex >= static_cast<int>(activatedChannelOutput.size())) {
      std::ostringstream os;
      os << "Invalid calibIndex " << channel.calibIndex
         << " while filling activated channel output arrays";
      throw std::runtime_error(os.str());
    }

    if (!IsActivatedFers(fersChannelsOverThreshold, channel.fersKey)) {
      continue;
    }

    const std::size_t idx = static_cast<std::size_t>(channel.calibIndex);

    activatedChannelOutputTrue[idx] = channel.outputSignal;
    activatedChannelOutputSmeared[idx] = channel.smearedOutputSignal;

    // keep the legacy branch identical to the smeared output
    activatedChannelOutput[idx] = channel.smearedOutputSignal;
  }
}




void PrintPerEventFersChannelCounts(
  unsigned int eventIndex,
  const ChannelSignalMap& sciChannelSignals,
  const ChannelSignalMap& cerChannelSignals,
  const std::unordered_map<FersKey, unsigned int>& fersChannelsOverThreshold)
{
  std::map<std::pair<int, int>, FersKey> fersSeen;

  auto collectFers = [&](const ChannelSignalMap& channelSignals) {
    for (const auto& [channelKey, channel] : channelSignals) {
      const int boardID = board_id_from_fers_key(channel.fersKey);
      fersSeen[{boardID, channel.fersId}] = channel.fersKey;
    }
  };

  collectFers(sciChannelSignals);
  collectFers(cerChannelSignals);

  if (printFersLog) {
    std::cout << "Event " << eventIndex << ":\n";
    for (const auto& [id, key] : fersSeen) {
      const int boardID = id.first;
      const int fersId = id.second;

      const auto it = fersChannelsOverThreshold.find(key);
      const unsigned int nChannelsOverThreshold =
        (it != fersChannelsOverThreshold.end()) ? it->second : 0U;

      std::cout
        << "  board " << boardID
        << ", FERS " << fersId
        << " -> channels above threshold = " << nChannelsOverThreshold;

      if (nChannelsOverThreshold >= FersMultiplicity) {
        std::cout << " [ON]";
      } else {
        std::cout << " [OFF]";
      }

      std::cout << '\n';
    }
  }
}




// -----------------------------------------------------------------------------
// SiPM map and event display helpers
// -----------------------------------------------------------------------------
SipmLookup LoadSipmMap(const std::string& jsonPath)
{
  std::ifstream in(jsonPath);
  if (!in) {
    throw std::runtime_error("Cannot open SiPM map JSON: " + jsonPath);
  }

  json j;
  in >> j;

  SipmLookup lookup;
  for (auto it = j.begin(); it != j.end(); ++it) {
    const json& node = it.value();

    SipmMapEntry entry;
    entry.calibIndex  = std::stoi(it.key());
    entry.boardID     = node.value("boardID", -1);
    entry.type        = node.value("type", "");
    entry.row         = node.value("row", -1);
    entry.column      = node.value("column", -1);
    entry.x           = node.value("x", 0.0);
    entry.y           = node.value("y", 0.0);
    entry.fersId      = node.value("fersId", -1);
    entry.ch          = node.value("ch", -1);
    entry.module_name = node.value("module_name", "");
    entry.x_local     = node.value("x_local", 0.0);
    entry.y_local     = node.value("y_local", 0.0);

    lookup[make_sipm_key(entry.module_name, entry.type, entry.row, entry.column)] = entry;
  }

  return lookup;
}

const SipmMapEntry* FindSipmInfo(const SipmLookup& lookup,
                                 const std::string& tower,
                                 const std::string& type,
                                 int row,
                                 int groupedCol,
                                 unsigned int /*groupingValue*/)
{
  const int jsonColumn = groupedCol;
  const auto key = make_sipm_key(tower, type, row, jsonColumn);
  const auto it = lookup.find(key);
  return (it == lookup.end()) ? nullptr : &it->second;
}

std::string make_event_display_name(const std::string& prefix, unsigned int entry)
{
  std::ostringstream os;
  os << prefix << "_evt" << std::setw(6) << std::setfill('0') << entry;
  return os.str();
}

TH2F* CreateEventDisplayHist(const std::string& name,
                             const std::string& title,
                             unsigned int groupingValue)
{
  const int xBins = static_cast<int>(NofSiPMTowersX * NofFiberscolumn / groupingValue);
  const int yBins = NofSiPMTowersY * NofFibersrow;

  return new TH2F(name.c_str(), title.c_str(),
                  xBins,
                  -NofSiPMTowersX * moduleX / 2.0, NofSiPMTowersX * moduleX / 2.0,
                  yBins,
                  -NofSiPMTowersY * moduleY / 2.0, NofSiPMTowersY * moduleY / 2.0);
}

void FillActivatedEventDisplay(TH2F* hist,
                               const std::vector<EventDisplayHit>& hits,
                               const std::unordered_map<FersKey, unsigned int>& fersChannelsOverThreshold,
                               double minVisiblePhe)
{
  if (!hist) {
    return;
  }

  for (const auto& hit : hits) {
    if (IsActivatedFers(fersChannelsOverThreshold, hit.key) &&
        hit.rawPhe >= minVisiblePhe) {
      hist->Fill(hit.x, hit.y, hit.signal);
    }
  }
}

void GetSiPMcoordinate(int TowID,
                       int rowID,
                       int colID_original,
                       double& SiPM_X,
                       double& SiPM_Y,
                       std::string fiber,
                       unsigned int groupingValue)
{
  double TowerOffsetY = -((NoModulesSiPM - 1) * moduleY) / 2 + TowID * moduleY;
  unsigned int channel = static_cast<unsigned int>(colID_original / groupingValue);
  double colID = (static_cast<double>(groupingValue) - 1) / 2 + channel * groupingValue;

  if (fiber == "S") {
    SiPM_X = +moduleX / 2 - tuberadius - (tuberadius * 2) * colID;
    SiPM_Y = -moduleY / 2 + tuberadius + (sq3 * tuberadius) * rowID + tuberadius * (2. * sq3m1 - 1.);
  }

  if (fiber == "C") {
    SiPM_X = +moduleX / 2 + tuberadius - (tuberadius * 2) * colID;
    SiPM_Y = -moduleY / 2 + tuberadius + (sq3 * tuberadius) * rowID + tuberadius * (2. * sq3m1 - 1.);
  }

  SiPM_Y = TowerOffsetY + SiPM_Y;
}

// -----------------------------------------------------------------------------
// Analysis/output helpers
// -----------------------------------------------------------------------------

std::string simTower_to_tbTower(std::string simTower)
{
  return std::to_string(306 + std::stoi(simTower));
}

GaussianFitSummary FitGaussianSummary(TH1* hist)
{
  GaussianFitSummary out;
  if (!hist || hist->GetEntries() < 10) {
    return out;
  }

  hist->Fit("gaus", "Q0");
  TF1* fit = hist->GetFunction("gaus");
  if (!fit) {
    return out;
  }

  out.mean = fit->GetParameter(1);
  out.meanErr = fit->GetParError(1);
  out.rms = fit->GetParameter(2);
  out.rmsErr = fit->GetParError(2);
  return out;
}

HistMomentSummary GetHistMomentSummary(const TH1* hist)
{
  HistMomentSummary out;
  if (!hist || hist->GetEntries() < 1) {
    return out;
  }

  out.mean = hist->GetMean();
  out.meanErr = hist->GetMeanError();
  out.rms = hist->GetRMS();
  out.rmsErr = hist->GetRMSError();
  return out;
}

std::string CsvEscape(const std::string& value)
{
  std::string escaped = "\"";
  for (char c : value) {
    if (c == '"') {
      escaped += "\"\"";
    } else {
      escaped += c;
    }
  }
  escaped += "\"";
  return escaped;
}

void AppendSummaryCsvLocked(const std::string& csvPath,
                            const std::string& inputFile,
                            double truthEnergy,
                            int nentries,
                            const GaussianFitSummary& sFit,
                            const GaussianFitSummary& cFit,
                            const GaussianFitSummary& combFit,
                            const GaussianFitSummary& combChiFit,
                            const HistMomentSummary& sciX,
                            const HistMomentSummary& sciY,
                            const HistMomentSummary& cerX,
                            const HistMomentSummary& cerY,
                            const GaussianFitSummary& sFitFersOn,
                            const GaussianFitSummary& cFitFersOn)
{
  const int fd = open(csvPath.c_str(), O_CREAT | O_WRONLY | O_APPEND, 0644);
  if (fd == -1) {
    throw std::runtime_error("Cannot open summary CSV: " + csvPath);
  }

  if (flock(fd, LOCK_EX) != 0) {
    close(fd);
    throw std::runtime_error("Cannot lock summary CSV: " + csvPath);
  }

  struct stat st;
  if (fstat(fd, &st) != 0) {
    flock(fd, LOCK_UN);
    close(fd);
    throw std::runtime_error("Cannot stat summary CSV: " + csvPath);
  }

  std::ostringstream out;
  out << std::setprecision(17);

  if (st.st_size == 0) {
    out
      << "input_file,truth_energy,nentries,"
      << "s_fit_mean,s_fit_mean_err,s_fit_rms,s_fit_rms_err,"
      << "c_fit_mean,c_fit_mean_err,c_fit_rms,c_fit_rms_err,"
      << "comb_fit_mean,comb_fit_mean_err,comb_fit_rms,comb_fit_rms_err,"
      << "combchi_fit_mean,combchi_fit_mean_err,combchi_fit_rms,combchi_fit_rms_err,"
      << "res_sci_x_mean,res_sci_x_mean_err,res_sci_x_rms,res_sci_x_rms_err,"
      << "res_sci_y_mean,res_sci_y_mean_err,res_sci_y_rms,res_sci_y_rms_err,"
      << "res_cer_x_mean,res_cer_x_mean_err,res_cer_x_rms,res_cer_x_rms_err,"
      << "res_cer_y_mean,res_cer_y_mean_err,res_cer_y_rms,res_cer_y_rms_err,"
      << "s_fit_fers_on_mean,s_fit_fers_on_mean_err,s_fit_fers_on_rms,s_fit_fers_on_rms_err,"
      << "c_fit_fers_on_mean,c_fit_fers_on_mean_err,c_fit_fers_on_rms,c_fit_fers_on_rms_err\n";
  }

  out
    << CsvEscape(inputFile) << ','
    << truthEnergy << ','
    << nentries << ','
    << sFit.mean << ',' << sFit.meanErr << ',' << sFit.rms << ',' << sFit.rmsErr << ','
    << cFit.mean << ',' << cFit.meanErr << ',' << cFit.rms << ',' << cFit.rmsErr << ','
    << combFit.mean << ',' << combFit.meanErr << ',' << combFit.rms << ',' << combFit.rmsErr << ','
    << combChiFit.mean << ',' << combChiFit.meanErr << ',' << combChiFit.rms << ',' << combChiFit.rmsErr << ','
    << sciX.mean << ',' << sciX.meanErr << ',' << sciX.rms << ',' << sciX.rmsErr << ','
    << sciY.mean << ',' << sciY.meanErr << ',' << sciY.rms << ',' << sciY.rmsErr << ','
    << cerX.mean << ',' << cerX.meanErr << ',' << cerX.rms << ',' << cerX.rmsErr << ','
    << cerY.mean << ',' << cerY.meanErr << ',' << cerY.rms << ',' << cerY.rmsErr << ','
    << sFitFersOn.mean << ',' << sFitFersOn.meanErr << ',' << sFitFersOn.rms << ',' << sFitFersOn.rmsErr << ','
    << cFitFersOn.mean << ',' << cFitFersOn.meanErr << ',' << cFitFersOn.rms << ',' << cFitFersOn.rmsErr << '\n';

  const std::string text = out.str();
  const ssize_t written = write(fd, text.c_str(), text.size());
  if (written != static_cast<ssize_t>(text.size())) {
    flock(fd, LOCK_UN);
    close(fd);
    throw std::runtime_error("Cannot write summary CSV: " + csvPath);
  }

  flock(fd, LOCK_UN);
  close(fd);
}

// -----------------------------------------------------------------------------
// Main analysis
// -----------------------------------------------------------------------------

void HidraTB25Ana(double energy, const std::string& input)
{
  // read SiPM map and calibration files
  const std::string sipmMapPath =
    "/home/apareti/HidraSim2025/TBDataPreparation/2025_SPS/MapAndCalibration/sipm_map.json";

  const std::string pedestalJsonPath =
    "/home/apareti/HidraSim2025/TBDataPreparation/2025_SPS/MapAndCalibration/SiPM_pedestals_v1.json";

  const std::string adcToGeVJsonPath =
    "/home/apareti/HidraSim2025/TBDataPreparation/2025_SPS/MapAndCalibration/SiPM_ADCtoGeV_v1.json";

  // Prepare lookup map for SiPMs, with attached pedestal and noise info
  // key: "tower|type|row|column", value: map entry with FERS/channel info and pedestal/noise
  SipmLookup sipmLookup;
  try {
    sipmLookup = LoadSipmMap(sipmMapPath);

    const PedestalHgVector pedestalsHg = LoadHgPedestalsFromJson(pedestalJsonPath);
    const AdcToGeVVector adcToGeV = LoadAdcToGeVFromJson(adcToGeVJsonPath);

    AttachHgPedestalNoiseToSipmMap(sipmLookup, pedestalsHg, adcToGeV);

    std::cout << "Loaded " << sipmLookup.size()
              << " SiPM map entries from " << sipmMapPath << std::endl;
    std::cout << "Loaded " << pedestalsHg.size()
              << " HG pedestal entries from " << pedestalJsonPath << std::endl;
    std::cout << "Loaded " << adcToGeV.size()
              << " ADC->GeV entries from " << adcToGeVJsonPath << std::endl;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return;
  }

  // read input file and prepare output file
  const std::string infile = "../build/" + input;
  std::cout << "Using file: " << infile << std::endl;

  TFile* simfile = TFile::Open(infile.c_str(), "READ");
  if (!simfile || simfile->IsZombie()) {
    std::cerr << "Cannot open input file " << infile << std::endl;
    return;
  }

  TTree* simtree = static_cast<TTree*>(simfile->Get("DREMTubesout"));
  if (!simtree) {
    std::cerr << "Cannot find TTree DREMTubesout in " << infile << std::endl;
    simfile->Close();
    delete simfile;
    return;
  }

  std::ostringstream os;
  os << energy;
  const std::string enstr = os.str();
  const std::string outfile = "hidra" + enstr + ".root";
  TFile f(outfile.c_str(), "RECREATE");

  TDirectory* eventDisplayDir = f.mkdir("EventDisplays");
  f.cd();

  TTree* activatedChannelsTree =
    new TTree("ActivatedChannels", "Per-event activated channel output");

  unsigned int outEvent = 0;
  std::vector<double> activatedChannelOutput(TotalCalibChannels, 0.0);
  std::vector<double> activatedChannelOutputTrue(TotalCalibChannels, 0.0);
  std::vector<double> activatedChannelOutputSmeared(TotalCalibChannels, 0.0);

  activatedChannelsTree->Branch("event", &outEvent);
  activatedChannelsTree->Branch("activatedChannelOutput", &activatedChannelOutput);
  activatedChannelsTree->Branch("activatedChannelOutputTrue", &activatedChannelOutputTrue);
  activatedChannelsTree->Branch("activatedChannelOutputSmeared", &activatedChannelOutputSmeared);

  // ---------------------------------------
  // create the RNG banks once
  // Prepare random generators for noise smearing:
  // one per calibration channel for uncorrelated noise, one per FERS for correlated noise
  // ----------------------------------------
  NoiseGeneratorBank channelNoiseGenerators = BuildChannelNoiseGenerators();
  NoiseGeneratorBank fersNoiseGenerators = BuildFersNoiseGenerators();



  // Load module mapping to fill the event display maps
  int modcol[NofmodulesX * NofmodulesY];
  int modrow[NofmodulesX * NofmodulesY];
  for (int i = 0; i < NofmodulesX * NofmodulesY; i++) {
    int row = i / NofmodulesX;
    int col = i % NofmodulesX;
    int imod = modflag[i];
    if (imod >= 0) {
      modcol[imod] = col;
      modrow[imod] = row;
    }
  }

  // ----------------------------------------
  // Prepare output histograms
  // ----------------------------------------
  const double bmin = energy - 0.4 * std::sqrt(energy) * 10.;
  const double bmax = energy + 0.4 * std::sqrt(energy) * 10.;

  auto sciene = new TH1F("sciene", "sciene", 100, bmin, bmax);
  auto cerene = new TH1F("cerene", "cerene", 100, bmin, bmax);
  auto totene = new TH1F("totene", "totene", 100, bmin, bmax);
  auto totenec = new TH1F("totenec", "totenec", 100, bmin, bmax);
  auto totdep = new TH1F("totdep", "totdep", 100, 0., bmax);
  auto leakene = new TH1F("leakene", "leakene", 100, 0., 0.1);
  auto chidist = new TH1F("chidist", "chidist", 100, 0., 1.);
  auto mapcalo  = new TH2F("mapcalo", "mapcalo",
                           NofmodulesX, 0., NofmodulesX,
                           NofmodulesY, 0., NofmodulesY);
  auto SipmMapS = new TH2F("SipmMapS", "SipmS; Col; Row",
                           NofSiPMTowersX * NofFiberscolumn, 0, NofSiPMTowersX * NofFiberscolumn,
                           NofSiPMTowersY * NofFibersrow, 0, NofSiPMTowersY * NofFibersrow);
  auto SipmMapC = new TH2F("SipmMapC", "SipmC; Col; Row",
                           NofSiPMTowersX * NofFiberscolumn, 0, NofSiPMTowersX * NofFiberscolumn,
                           NofSiPMTowersY * NofFibersrow, 0, NofSiPMTowersY * NofFibersrow);

  auto SciSiPMCoordinates = new TH2F("SciSiPMCoordinates", "Sci SiPM Coordinates; X [mm]; Y[mm]",
                                     NofSiPMTowersX * NofFiberscolumn / grouping,
                                     -NofSiPMTowersX * moduleX / 2, NofSiPMTowersX * moduleX / 2,
                                     NofSiPMTowersY * NofFibersrow,
                                     -NofSiPMTowersY * moduleY / 2, NofSiPMTowersY * moduleY / 2);
  auto CerSiPMCoordinates = new TH2F("CerSiPMCoordinates", "Cer SiPM Coordinates; X [mm]; Y[mm]",
                                     NofSiPMTowersX * NofFiberscolumn / grouping,
                                     -NofSiPMTowersX * moduleX / 2, NofSiPMTowersX * moduleX / 2,
                                     NofSiPMTowersY * NofFibersrow,
                                     -NofSiPMTowersY * moduleY / 2, NofSiPMTowersY * moduleY / 2);

  auto ResidualHistSciX = new TH1F("ResidualHistSciX", "Residual Sci X; Residual [mm]; Entries",
                                   100, -10, 10);
  auto ResidualHistSciY = new TH1F("ResidualHistSciY", "Residual Sci Y; Residual [mm]; Entries",
                                   100, -10, 10);
  auto ResidualHistCerX = new TH1F("ResidualHistCerX", "Residual Cer X; Residual [mm]; Entries",
                                   100, -10, 10);
  auto ResidualHistCerY = new TH1F("ResidualHistCerY", "Residual Cer Y; Residual [mm]; Entries",
                                   100, -10, 10);

  auto SignalSfersID = new TH1F("SignalSfersID", "Signal S fibers FERS ID; FERS; Signal", 16, 0.5, 16.5);
  auto SignalCfersID = new TH1F("SignalCfersID", "Signal C fibers FERS ID; FERS; Signal", 16, 0.5, 16.5);

  auto scieneFersOn = new TH1F("scieneFersOn", "S energy for Activated FERS", 100, 0., bmax);
  auto cereneFersOn = new TH1F("cereneFersOn", "C energy for Activated FERS", 100, 0., bmax);

  auto h2_SvsC = new TH2F("SvsC", "S vs C energy (all fibres); S/E; C/E", 100, -0.5, 1.5, 100, 0., 1.5);
  auto h2_SvsCFersOn = new TH2F("SvsCFersOn", "S vs C energy for Activated FERS; S/E; C/E", 100, 0., 1.5, 100, 0., 1.5);

  // Count number of activated channels per FERS, to check the multiplicity distribution and the effect of different thresholds
  auto ActivatedChannelsPerFers = new TH2F("ActivatedChannelsPerFers", "Activated channels per FERS; FERS ID; Channels above threshold",16, 0.5, 16.5, 65, -0.5, 64.5);


  const int nentries = simtree->GetEntries();
  std::cout << "Entries " << nentries << std::endl;

  int pdg; simtree->SetBranchAddress("PrimaryPDGID", &pdg);
  double venergy; simtree->SetBranchAddress("PrimaryParticleEnergy", &venergy);
  double lenergy; simtree->SetBranchAddress("EscapedEnergyl", &lenergy);
  double denergy; simtree->SetBranchAddress("EscapedEnergyd", &denergy);
  double edep; simtree->SetBranchAddress("EnergyDepCalo", &edep);
  double Stot; simtree->SetBranchAddress("NofPMTScinDet", &Stot);
  double Ctot; simtree->SetBranchAddress("NofPMTCherDet", &Ctot);
  double PSdep; simtree->SetBranchAddress("PSEnergy", &PSdep);
  double beamX; simtree->SetBranchAddress("PrimaryX", &beamX);
  double beamY; simtree->SetBranchAddress("PrimaryY", &beamY);
  std::vector<double>* TowerE = nullptr; simtree->SetBranchAddress("VecTowerE", &TowerE);
  std::vector<double>* SPMT = nullptr; simtree->SetBranchAddress("VecSPMT", &SPMT);
  std::vector<double>* CPMT = nullptr; simtree->SetBranchAddress("VecCPMT", &CPMT);
  std::vector<double>* SSiPM = nullptr; simtree->SetBranchAddress("VectorSignals", &SSiPM);
  std::vector<double>* CSiPM = nullptr; simtree->SetBranchAddress("VectorSignalsCher", &CSiPM);

  for (unsigned int i = 0; i < static_cast<unsigned int>(nentries); i++) {
    simtree->GetEntry(i);

    // store one event display every N events
    const bool writeEventDisplay =
      (EventDisplayEvery > 0) && (((i + 1) % EventDisplayEvery) == 0);

    double ecalo = energy - lenergy / 1000.0;
    double totsci = 0.;
    double totcer = 0.;
    double tottow = 0.;

    double barX_sci = 0.;
    double barY_sci = 0.;
    double barX_cer = 0.;
    double barY_cer = 0.;

    outEvent = i;
    // prepare the per-event data structures for the smearing and activation steps
    std::fill(activatedChannelOutput.begin(), activatedChannelOutput.end(), 0.0);
    std::fill(activatedChannelOutputTrue.begin(), activatedChannelOutputTrue.end(), 0.0);
    std::fill(activatedChannelOutputSmeared.begin(), activatedChannelOutputSmeared.end(), 0.0);
    std::unordered_map<FersKey, unsigned int> fersChannelsOverThreshold;

    std::vector<EventDisplayHit> eventDisplaySciHits;
    std::vector<EventDisplayHit> eventDisplayCerHits;

    ChannelSignalMap sciChannelSignals;
    ChannelSignalMap cerChannelSignals;

    sciChannelSignals.reserve(SSiPM->size());
    cerChannelSignals.reserve(CSiPM->size());

    if (writeEventDisplay) {
      eventDisplaySciHits.reserve(SSiPM->size());
      eventDisplayCerHits.reserve(CSiPM->size());
    }

    // Fill energy in PMT towers
    for (unsigned int j = 0; j < SPMT->size(); j++) {
      totsci += SPMT->at(j) / sciPheGeV;
      totcer += CPMT->at(j) / cerPheGeV;
      tottow += TowerE->at(j);
      mapcalo->Fill(modcol[j], modrow[j], TowerE->at(j) / 1000. / nentries);
    }

    double sciPosWeight = 0.0;
    double cerPosWeight = 0.0;

    // Prepare shared random map to have correlated noise between S/C channels
    std::unordered_map<FersKey, double> fersEventNoiseCache;


    // Start looping on SiPM signals
    for (unsigned int n = 0; n < SSiPM->size(); n++) {
      const double rawPhe = SSiPM->at(n);        // number of photoelectrons
      const double content = rawPhe / sciPheGeV; // converted to GeV

      totsci += content;

      const unsigned int towID  = static_cast<unsigned int>(n / (NofFiberscolumn * NofFibersrow / 2));
      const unsigned int SiPMID = n % (NofFiberscolumn * NofFibersrow / 2);
      //const unsigned int colID  = NofFiberscolumn - static_cast<unsigned int>(SiPMID / (NofFibersrow / 2));
      const unsigned int colID  = NofFiberscolumn -1U - static_cast<unsigned int>(SiPMID / (NofFibersrow / 2));
      const unsigned int rowID  = 2 * static_cast<unsigned int>(SiPMID % (NofFibersrow / 2));

      double SiPM_X = 0.0;
      double SiPM_Y = 0.0;

      GetSiPMcoordinate(towID, rowID, colID, SiPM_X, SiPM_Y, "S", grouping); // Currently not used, JSON map is read instead

      const std::string tbTower = simTower_to_tbTower(std::to_string(towID)); // convert from sim tower ID to test beam tower ID
      const int groupedCol = grouped_column(static_cast<int>(colID), grouping); // single fibre to channel mapping

      SipmMapS->Fill(colID, towID * NofFibersrow + rowID, content);

      // now associate the SiPM signal to the corresponding FERS and channel
      const SipmMapEntry* info =
        FindSipmInfo(sipmLookup, tbTower, "S", static_cast<int>(rowID), groupedCol, grouping);

      if (info) {
        SciSiPMCoordinates->Fill(info->x, info->y, content); // add fibre content to mapped channel signal
        SignalSfersID->Fill(info->fersId, content / nentries);
        barX_sci += info->x * content;
        barY_sci += info->y * content;
        sciPosWeight += content;

        // accumulate the signal in the corresponding channel for smearing and activation steps
        if (info->boardID >= 0 && info->fersId >= 0) {
          // create the key to identify the FERS this SiPM belongs to
          const FersKey key = make_fers_key(info->boardID, info->fersId); 
          
          // Apply saturation consistently for both threshold and output energy
          const double saturatedSignal =
            ApplySaturation
              ? apply_sipm_saturation_from_pe(rawPhe, sciPheGeV, NofSipmCells_sci)
              : content;

          // accumulate the signal for this channel
          // reads threshold values for this channel from the map
          // and updates the signal and FERS->channel count maps
          AccumulateChannelSignal(
            sciChannelSignals,
            *info,
            saturatedSignal
          );

          if (writeEventDisplay) {
            eventDisplaySciHits.push_back({key, info->x, info->y, rawPhe, content});
          }
        }
      }
    }

    for (unsigned int n = 0; n < CSiPM->size(); n++) {
      const double rawPhe = CSiPM->at(n);
      const double content = rawPhe / cerPheGeV;

      totcer += content;

      const unsigned int towID  = static_cast<unsigned int>(n / (NofFiberscolumn * NofFibersrow / 2));
      const unsigned int SiPMID = n % (NofFiberscolumn * NofFibersrow / 2);
      //const unsigned int colID  = NofFiberscolumn - static_cast<unsigned int>(SiPMID / (NofFibersrow / 2));
      const unsigned int colID  = NofFiberscolumn -1U - static_cast<unsigned int>(SiPMID / (NofFibersrow / 2));
      const unsigned int rowID  = 2 * static_cast<unsigned int>(SiPMID % (NofFibersrow / 2)) + 1;

      SipmMapC->Fill(colID, towID * NofFibersrow + rowID, content);

      const std::string tbTower = simTower_to_tbTower(std::to_string(towID));
      const int groupedCol = grouped_column(static_cast<int>(colID), grouping);

      const SipmMapEntry* info =
        FindSipmInfo(sipmLookup, tbTower, "C", static_cast<int>(rowID), groupedCol, grouping);

      if (info) {
        CerSiPMCoordinates->Fill(info->x, info->y, content);
        SignalCfersID->Fill(info->fersId, content / nentries);
        barX_cer += info->x * content;
        barY_cer += info->y * content;
        cerPosWeight += content;

        // accumulate the signal for this channel
        if (info->boardID >= 0 && info->fersId >= 0) {
          const FersKey key = make_fers_key(info->boardID, info->fersId);
          
          // Apply saturation consistently for both threshold and output energy
          const double saturatedSignal =
            ApplySaturation
              ? apply_sipm_saturation_from_pe(rawPhe, cerPheGeV, NofSipmCells_cer)
              : content;

          AccumulateChannelSignal(
            cerChannelSignals,
            *info,
            saturatedSignal
          );

          if (writeEventDisplay) {
            eventDisplayCerHits.push_back({key, info->x, info->y, rawPhe, content});
          }
        }
      }
    } // end C Sipm loop


    //------------------------------
    // FERS ACTIVATION LOGIC
    //------------------------------

    // Now apply noise smearing and thresholding to the accumulated channel signals
    // and fill the activated channel output array for this event
    // Count number of channels over threshold per FERS
    SmearAndCountChannelsOverThreshold( // S channel
      sciChannelSignals,
      fersChannelsOverThreshold,
      NoiseMode,
      channelNoiseGenerators,
      fersNoiseGenerators,
      fersEventNoiseCache 
    );

    SmearAndCountChannelsOverThreshold( // C channel
      cerChannelSignals,
      fersChannelsOverThreshold,
      NoiseMode,
      channelNoiseGenerators,
      fersNoiseGenerators,
      fersEventNoiseCache 
    );

    for (unsigned int fersId = 1; fersId <= TotalFers; ++fersId) {
      unsigned int nOverThreshold = 0;

      for (const auto& [fersKey, count] : fersChannelsOverThreshold) {
        if (static_cast<unsigned int>(static_cast<uint32_t>(fersKey)) == fersId) {
          nOverThreshold += count;
        }
      }

      ActivatedChannelsPerFers->Fill(fersId, nOverThreshold);
    }


    // For activated FERS channels, subtract the pedestal median
    SubtractPedestalMedianFromActivatedChannels(
      sciChannelSignals,
      fersChannelsOverThreshold
    );
    SubtractPedestalMedianFromActivatedChannels(
      cerChannelSignals,
      fersChannelsOverThreshold
    );

    FillActivatedChannelOutputArrays(
      sciChannelSignals,
      fersChannelsOverThreshold,
      activatedChannelOutput,
      activatedChannelOutputTrue,
      activatedChannelOutputSmeared
    );

    FillActivatedChannelOutputArrays(
      cerChannelSignals,
      fersChannelsOverThreshold,
      activatedChannelOutput,
      activatedChannelOutputTrue,
      activatedChannelOutputSmeared
    );

    activatedChannelsTree->Fill();


    PrintPerEventFersChannelCounts(
      i,
      sciChannelSignals,
      cerChannelSignals,
      fersChannelsOverThreshold
    );

    const double totsciFersOn =
      SumActivatedSmearedChannelOutput(sciChannelSignals, fersChannelsOverThreshold);

    const double totcerFersOn =
      SumActivatedSmearedChannelOutput(cerChannelSignals, fersChannelsOverThreshold);

    if (sciPosWeight > 0.) {
      ResidualHistSciX->Fill((barX_sci / sciPosWeight) - beamX);
      ResidualHistSciY->Fill((barY_sci / sciPosWeight) - beamY);
    }

    if (cerPosWeight > 0.) {
      ResidualHistCerX->Fill((barX_cer / cerPosWeight) - beamX);
      ResidualHistCerY->Fill((barY_cer / cerPosWeight) - beamY);
    }

    if (writeEventDisplay && eventDisplayDir) {
      eventDisplayDir->cd();

      std::ostringstream sciTitle;
      sciTitle << "Sci event display, entry " << i
               << " (grouped channels, activated FERS only); X [mm]; Y [mm]";
      TH2F* sciEventDisplay = CreateEventDisplayHist(
        make_event_display_name("EventDisplaySci", i),
        sciTitle.str(),
        grouping
      );

      std::ostringstream cerTitle;
      cerTitle << "Cer event display, entry " << i
               << " (grouped channels, activated FERS only); X [mm]; Y [mm]";
      TH2F* cerEventDisplay = CreateEventDisplayHist(
        make_event_display_name("EventDisplayCer", i),
        cerTitle.str(),
        grouping
      );

      FillActivatedEventDisplay(
        sciEventDisplay,
        eventDisplaySciHits,
        fersChannelsOverThreshold,
        MinVisiblePheS
      );
      FillActivatedEventDisplay(
        cerEventDisplay,
        eventDisplayCerHits,
        fersChannelsOverThreshold,
        MinVisiblePheC
      );

      sciEventDisplay->Write();
      cerEventDisplay->Write();

      delete sciEventDisplay;
      delete cerEventDisplay;

      f.cd();
    }



    scieneFersOn->Fill(totsciFersOn);
    cereneFersOn->Fill(totcerFersOn);

    sciene->Fill(totsci);
    cerene->Fill(totcer);
    totene->Fill(0.5 * (totsci + totcer));
    totenec->Fill(picont * (totsci - chi * totcer) / (1 - chi));
    totdep->Fill(tottow / 1000.);
    leakene->Fill(lenergy / 1000. / energy);
    chidist->Fill((totsci - ecalo) / (totcer - ecalo));


    h2_SvsC->Fill(totsci / energy, totcer / energy);
    h2_SvsCFersOn->Fill(totsciFersOn / energy, totcerFersOn / energy);


    //break;
  }

  const GaussianFitSummary sFit = FitGaussianSummary(sciene);
  const GaussianFitSummary cFit = FitGaussianSummary(cerene);
  const GaussianFitSummary combFit = FitGaussianSummary(totene);
  const GaussianFitSummary combChiFit = FitGaussianSummary(totenec);

  const HistMomentSummary resSciX = GetHistMomentSummary(ResidualHistSciX);
  const HistMomentSummary resSciY = GetHistMomentSummary(ResidualHistSciY);
  const HistMomentSummary resCerX = GetHistMomentSummary(ResidualHistCerX);
  const HistMomentSummary resCerY = GetHistMomentSummary(ResidualHistCerY);

  const GaussianFitSummary sFitFersOn = FitGaussianSummary(scieneFersOn);
  const GaussianFitSummary cFitFersOn = FitGaussianSummary(cereneFersOn);

  try {
    AppendSummaryCsvLocked("hidra_summary.csv",
                           input,
                           energy,
                           nentries,
                           sFit,
                           cFit,
                           combFit,
                           combChiFit,
                           resSciX,
                           resSciY,
                           resCerX,
                           resCerY,
                           sFitFersOn,
                           cFitFersOn);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }

  f.Write();
  f.Close();

  simfile->Close();
  delete simfile;
}
