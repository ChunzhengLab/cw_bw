#include <getopt.h>
#include <iostream>
#include <string>
#include <random>

#include "Config.h"
#include "EventFactory.h"
#include "AnalyzerCVE.h"

static void PrintUsage() {
    std::cout << "Usage: bwgen [options]\n"
              << "  -h, --help               Show this help message\n"
              << "  -c, --config <file>      Configuration YAML file (default: configs/default.yaml)\n"
              << "  -o, --output-file <file> Output ROOT file name (default: result.root)\n";
}

int main(int argc, char** argv) {
    // 1. 解析命令行参数
    std::string configPath = "configs/default.yaml";
    std::string outputFile = "result.root";
    const struct option longOpts[] = {
        {"help",        no_argument,       nullptr, 'h'},
        {"config",      required_argument, nullptr, 'c'},
        {"output-file", required_argument, nullptr, 'o'},
        {nullptr,       0,                 nullptr,  0 }
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "hc:o:", longOpts, nullptr)) != -1) {
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
            default:
                PrintUsage();
                return 1;
        }
    }

    // 2. 加载配置
    Config cfg = Config::Load(configPath);

    // 3. 实例化事件生成器和分析器
    EventFactory  factory(cfg);
    AnalyzerCVE   analyzer;
    analyzer.Init();

    // 4. 随机数引擎，用随机设备进行初始化
    std::random_device rd;
    std::mt19937       rng(rd());
    std::uniform_real_distribution<float> centralityDist(0.0f, 100.0f);

    // 5. 事件循环
    for (unsigned int i = 0; i < static_cast<unsigned int>(cfg.nEvents); ++i) {
        float centrality = centralityDist(rng);
        Event evt = factory.GetEvent(centrality);
        analyzer.Process(evt);
    }

    // 6. 写出所有直方图到 ROOT 文件
    analyzer.Write(outputFile);

    std::cout << "Simulation completed. Output saved to " << outputFile << std::endl;
    return 0;
}
