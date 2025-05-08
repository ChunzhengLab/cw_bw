#include <getopt.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "AnalyzerCVE.h"
#include "Config.h"
#include "EventFactory.h"

static void PrintUsage() {
  std::cout << "Usage: bwgen [options]\n"
            << "  -h, --help               Show this help message\n"
            << "  -c, --config <file>      Configuration YAML file (default: "
               "../configs/default.yaml)\n"
            << "  -o, --output-file <file> Output ROOT file name (default: "
               "result.root)\n"
            << "  -C, --centrality <int>    Centrality % (one of 15,25,35,45,55)\n"
            << "  -n, --events <int>       Number of events to simulate (overrides config)\n"
            << "  -r, --ratio-proton-lambda <float>    Override ratioProtonLambda\n"
            << "  -i, --ratio-proton-inclusive <float> Override ratioProtonInclusive\n"
            << "  -f, --frac-lbc <float>               Override fracLBC\n"
            << "  -R, --rho2scale <float>  Scale all rho2_p and rho2_L by given factor\n"
            << "  -B, --betascale <float>  Scale all betaT by given factor\n";
}

int main(int argc, char** argv) {
  // 1. 解析命令行参数
  std::string configPath = "../configs/default.yaml";
  std::string outputFile = "result.root";
  int centralityArg = -1;
  int nEventsArg = -1;
  float ratioProtonLambdaArg = -1.0f;
  float ratioProtonInclusiveArg = -1.0f;
  float fracLBCArg = -1.0f;
  float rho2scaleArg = -1.0f;
  float betascaleArg = -1.0f;
  const struct option longOpts[] = {{"help", no_argument, nullptr, 'h'},
                                    {"config", required_argument, nullptr, 'c'},
                                    {"output-file", required_argument, nullptr, 'o'},
                                    {"centrality", required_argument, nullptr, 'C'},
                                    {"events", required_argument, nullptr, 'n'},
                                    {"ratio-proton-lambda", required_argument, nullptr, 'r'},
                                    {"ratio-proton-inclusive", required_argument, nullptr, 'i'},
                                    {"frac-lbc", required_argument, nullptr, 'f'},
                                    {"rho2scale", required_argument, nullptr, 'R'},
                                    {"betascale", required_argument, nullptr, 'B'},
                                    {nullptr, 0, nullptr, 0}};
  int opt;
  while ((opt = getopt_long(argc, argv, "hc:o:C:n:r:i:f:R:B:", longOpts, nullptr)) != -1) {
    switch (opt) {
      case 'h':
        PrintUsage();
        return 0;
      case 'c':
        configPath = optarg;
        break;
      case 'o':
        outputFile = optarg;
        break;
      case 'C':
        centralityArg = std::stoi(optarg);
        break;
      case 'n':
        nEventsArg = std::stoi(optarg);
        break;
      case 'r':
        ratioProtonLambdaArg = std::stof(optarg);
        break;
      case 'i':
        ratioProtonInclusiveArg = std::stof(optarg);
        break;
      case 'f':
        fracLBCArg = std::stof(optarg);
        break;
      case 'R':
        rho2scaleArg = std::stof(optarg);
        break;
      case 'B':
        betascaleArg = std::stof(optarg);
        break;
      default:
        PrintUsage();
        return 1;
    }
  }

  if (centralityArg != 15 && centralityArg != 25 && centralityArg != 35 && centralityArg != 45 && centralityArg != 55) {
    std::cerr << "Error: --centrality must be one of {15, 25, 35, 45, 55}\n";
    return 1;
  }

  // 2. 加载配置
  Config cfg = Config::Load(configPath);

  if (ratioProtonLambdaArg > 0.0f) {
    std::cout << "Overriding ratioProtonLambda from " << cfg.ratioProtonLambda << " to " << ratioProtonLambdaArg
              << std::endl;
    cfg.ratioProtonLambda = ratioProtonLambdaArg;
  }
  if (ratioProtonInclusiveArg > 0.0f) {
    std::cout << "Overriding ratioProtonInclusive from " << cfg.ratioProtonInclusive << " to "
              << ratioProtonInclusiveArg << std::endl;
    cfg.ratioProtonInclusive = ratioProtonInclusiveArg;
  }
  if (rho2scaleArg > 0.0f) {
    std::cout << "Scaling rho2_p and rho2_L by factor " << rho2scaleArg << std::endl;
    for (auto& v : cfg.rho2_p) v *= rho2scaleArg;
    for (auto& v : cfg.rho2_L) v *= rho2scaleArg;
  }
  if (betascaleArg > 0.0f) {
    std::cout << "Scaling betaT by factor " << betascaleArg << std::endl;
    for (auto& v : cfg.betaT) v *= betascaleArg;
  }
  size_t idx = 0;
  if (centralityArg == 25) idx = 1;
  else if (centralityArg == 35) idx = 2;
  else if (centralityArg == 45) idx = 3;
  else if (centralityArg == 55) idx = 4;
  if (fracLBCArg > 0.0f) {
    std::cout << "Overriding fracLBC[" << idx << "] from " << cfg.fracLBC[idx]
              << " to " << fracLBCArg << std::endl;
    cfg.fracLBC[idx] = fracLBCArg;
  }
  // Generate output file name based on centrality, fracLBC and optional scales
  {
    std::ostringstream ofname;
    ofname << "results_cent" << centralityArg << "_fLBC" << cfg.fracLBC[idx];
    if (rho2scaleArg > 0.0f) ofname << "_rho2scale" << rho2scaleArg;
    if (betascaleArg > 0.0f) ofname << "_betascale" << betascaleArg;
    ofname << ".root";
    outputFile = ofname.str();
  }

  // 如果指定了事件数，则覆盖配置文件中的值
  if (nEventsArg > 0) {
    std::cout << "Overriding number of events from " << cfg.nEvents << " to " << nEventsArg << std::endl;
    cfg.nEvents = nEventsArg;
  }

  // Startup banner
  std::cout << "==================================================================================\n";
  std::cout << "     Chunzheng Wang's Blast Wave Model with Local Baryon Number Conservation \n";
  std::cout << "             Author: Chunzheng Wang (chunzheng.wang@icloud.com) \n";
  std::cout << "==================================================================================\n";

  // === Simulation Parameters ===
  std::cout << ">>>Config file:            " << configPath << "\n";
  std::cout << ">>>Output ROOT file:       " << outputFile << "\n";
  std::cout << ">>>Centrality:             " << centralityArg << "\n";
  std::cout << ">>>Number of events:       " << cfg.nEvents << "\n";
  std::cout << ">>>ratioProtonLambda:      " << cfg.ratioProtonLambda << "\n";
  std::cout << ">>>ratioProtonInclusive:   " << cfg.ratioProtonInclusive << "\n";
  std::cout << ">>>fracLBC:                " << cfg.fracLBC[idx] << "\n";
  std::cout << ">>>rho2scale:              " << (rho2scaleArg > 0.0f ? rho2scaleArg : 1.0f) << "\n";
  std::cout << ">>>betascale:              " << (betascaleArg > 0.0f ? betascaleArg : 1.0f) << "\n";
  std::cout << "==================================================================================\n";

  // 3. 实例化事件生成器和分析器
  EventFactory factory(cfg);
  AnalyzerCVE analyzer;
  analyzer.Init();

  // Progress bar
  auto printProgress = [&](unsigned int current) {
    const unsigned int width = 80;
    unsigned int pos = static_cast<unsigned int>(width * current / cfg.nEvents);
    std::cout << "\r[";
    for (unsigned int i = 0; i < width; ++i) { std::cout << (i < pos ? '=' : ' '); }
    std::cout << "] " << std::setw(3) << (100 * current / cfg.nEvents) << "% (" << current << "/" << cfg.nEvents << ")"
              << std::flush;
    if (current == cfg.nEvents) std::cout << std::endl;
  };

  // 4. 随机数引擎，用随机设备进行初始化
  std::random_device rd;
  std::mt19937 rng(rd());
  // std::uniform_real_distribution<float> centralityDist(0.0f, 100.0f);

  // 5. 事件循环
  for (unsigned int i = 0; i < static_cast<unsigned int>(cfg.nEvents); ++i) {
    float centrality = centralityArg;
    Event evt = factory.GetEvent(centrality);
    analyzer.Process(evt);
    printProgress(i + 1);
  }

  // 6. 写出所有直方图到 ROOT 文件
  analyzer.Write(outputFile);

  std::cout << "Simulation completed. Output saved to " << outputFile << std::endl;
  return 0;
}
