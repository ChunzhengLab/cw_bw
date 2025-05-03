#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <cstdlib>  // for std::stof
#include <fstream>

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
            << "  -f, --frac-lbc <float>               Override fracLBC\n";
}

int main(int argc, char **argv) {
  // 1. 解析命令行参数
  std::string configPath = "../configs/default.yaml";
  std::string outputFile = "result.root";
  int centralityArg = -1;
  int nEventsArg = -1;
  float ratioProtonLambdaArg = -1.0f;
  float ratioProtonInclusiveArg = -1.0f;
  float fracLBCArg = -1.0f;
  const struct option longOpts[] = {
      {"help", no_argument, nullptr, 'h'},
      {"config", required_argument, nullptr, 'c'},
      {"output-file", required_argument, nullptr, 'o'},
      {"centrality", required_argument, nullptr, 'C'},
      {"events", required_argument, nullptr, 'n'},
      {"ratio-proton-lambda", required_argument, nullptr, 'r'},
      {"ratio-proton-inclusive", required_argument, nullptr, 'i'},
      {"frac-lbc", required_argument, nullptr, 'f'},
      {nullptr, 0, nullptr, 0}};
  int opt;
  while ((opt = getopt_long(argc, argv, "hc:o:C:n:r:i:f:", longOpts, nullptr)) != -1) {
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
    default:
      PrintUsage();
      return 1;
    }
  }

  if (centralityArg != 15 && centralityArg != 25 && centralityArg != 35 &&
      centralityArg != 45 && centralityArg != 55) {
    std::cerr << "Error: --centrality must be one of {15, 25, 35, 45, 55}\n";
    return 1;
  }

  // 2. 加载配置
  Config cfg = Config::Load(configPath);

  if (ratioProtonLambdaArg > 0.0f) {
    std::cout << "Overriding ratioProtonLambda from " << cfg.ratioProtonLambda
              << " to " << ratioProtonLambdaArg << std::endl;
    cfg.ratioProtonLambda = ratioProtonLambdaArg;
  }
  if (ratioProtonInclusiveArg > 0.0f) {
    std::cout << "Overriding ratioProtonInclusive from " << cfg.ratioProtonInclusive
              << " to " << ratioProtonInclusiveArg << std::endl;
    cfg.ratioProtonInclusive = ratioProtonInclusiveArg;
  }
  if (fracLBCArg > 0.0f) {
    std::cout << "Overriding fracLBC from " << cfg.fracLBC
              << " to " << fracLBCArg << std::endl;
    cfg.fracLBC = fracLBCArg;
  }
  // Write parameters to par.log
  {
    std::ofstream parLog("par.log");
    parLog << "configPath=" << configPath << "\n";
    parLog << "outputFile=" << outputFile << "\n";
    parLog << "centrality=" << centralityArg << "\n";
    parLog << "nEvents=" << cfg.nEvents << "\n";
    parLog << "ratioProtonLambda=" << cfg.ratioProtonLambda << "\n";
    parLog << "ratioProtonInclusive=" << cfg.ratioProtonInclusive << "\n";
    parLog << "fracLBC=" << cfg.fracLBC << "\n";
    parLog.close();
  }

  // 如果指定了事件数，则覆盖配置文件中的值
  if (nEventsArg > 0) {
    std::cout << "Overriding number of events from " << cfg.nEvents
              << " to " << nEventsArg << std::endl;
    cfg.nEvents = nEventsArg;
  }

  // === Simulation Parameters ===
  std::cout << "============================================================\n";
  std::cout << "bwgen run parameters:\n";
  std::cout << "  Config file:            " << configPath << "\n";
  std::cout << "  Output ROOT file:       " << outputFile << "\n";
  std::cout << "  Centrality:             " << centralityArg << "\n";
  std::cout << "  Number of events:       " << cfg.nEvents << "\n";
  std::cout << "  ratioProtonLambda:      " << cfg.ratioProtonLambda << "\n";
  std::cout << "  ratioProtonInclusive:   " << cfg.ratioProtonInclusive << "\n";
  std::cout << "  fracLBC:                " << cfg.fracLBC << "\n";
  std::cout << "============================================================\n";

  // 3. 实例化事件生成器和分析器
  EventFactory factory(cfg);
  AnalyzerCVE analyzer;
  analyzer.Init();

  // Progress bar
  auto printProgress = [&](unsigned int current) {
    const unsigned int width = 50;
    unsigned int pos = static_cast<unsigned int>(width * current / cfg.nEvents);
    std::cout << "\r[";
    for (unsigned int i = 0; i < width; ++i) {
      std::cout << (i < pos ? '=' : ' ');
    }
    std::cout << "] " << std::setw(3) << (100 * current / cfg.nEvents)
              << "% (" << current << "/" << cfg.nEvents << ")" << std::flush;
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

  std::cout << "Simulation completed. Output saved to " << outputFile
            << std::endl;
  return 0;
}
