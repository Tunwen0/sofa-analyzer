***向下滚动页面以浏览中文版***

# sofa-analyzer

**sofa-analyzer** is a lightweight, robust Python utility designed to parse and analyze **SOFA (Spatially Oriented Format for Acoustics)** files. It specifically focuses on extracting, converting, and  visualizing measurement positions (Source Positions), providing a clear  breakdown of Head-Related Transfer Function (HRTF) or Head-Related  Impulse Response (HRIR) measurement points located within the frontal  hemisphere of the subject.

------

### Introduction

In spatial audio and acoustics research, understanding the distribution of measurement points in a SOFA file is crucial. **sofa-analyzer** reads standard `.sofa` files, validates the data dimensions, and outputs a human-readable list of measurement points. It is particularly optimized to display points  within the **frontal 90-degree range** ($-90^\circ \le \text{Azimuth} \le +90^\circ$), grouped by elevation.

### Features

- **Format Support:** Fully supports the SOFA standard (Spatially Oriented Format for Acoustics).
- **Coordinate Conversion:** Automatically detects the  coordinate system (Spherical or Cartesian) and converts Cartesian  coordinates to Spherical (Azimuth, Elevation, Radius) for consistent  analysis.
- **Robust Dimension Handling:** Enhanced stability for handling multi-receiver data structures (Dimensions: `M, R, C`), ensuring compatibility with complex SOFA files containing multiple ears or receivers.
- **Numerical Stability:** Implements safe trigonometric calculations to prevent floating-point errors (e.g., `arcsin` domain errors).
- **Smart Filtering:** Automatically filters and displays points relevant to the frontal hemisphere.
- **Formatted Output:** Sorts and groups data by elevation and azimuth with an intuitive display format (0° centered, alternating left/right).

### Prerequisites & Installation

Ensure you have Python 3.6+ installed. You will need to install the required dependencies before running the script.

1. **Clone or Download** this repository/script.
2. **Install Dependencies** using `pip`:

bash

```bash
pip install pysofaconventions numpy netCDF4
```

### Usage

1. Open your terminal or command prompt.
2. Navigate to the directory containing `sofa-analyzer.py`.
3. Run the script:

bash

```bash
python sofa-analyzer.py
```

1. Follow the on-screen prompts:
   - Enter the absolute path to your `.sofa` file (e.g., `C:\AudioData\HRTF_Subject1.sofa`).
   - The tool will analyze the file and print the distribution of measurement points.

#### Example Output

text

```text
============================================================
       sofa-analyzer v1.1.0 - SOFA File Analysis Tool
============================================================

Please enter the storage path of your sofa file:
>>> ./data/HRTF_Subject_001.sofa

Loading and analyzing SOFA file...

============================================================

Import Successful. The file contains 1200 measurement points.

Among them, there are 850 points within the frontal 90-degree range
(Distributed across 15 elevations):

  1、Elevation -40°：Azimuth 0°，+15°，-15°，+30°，-30°
  2、Elevation 0°：Azimuth 0°，+5°，-5°，+10°，-10°，+15°，-15°
  ...
  
  Remaining points not listed require algorithmic conversion.

============================================================
```

### Dependencies

This project relies on the following open-source libraries:

- **[pysofaconventions](https://github.com/andresperezlopez/pysofaconventions):** For reading and parsing SOFA files according to the AES69 standard.
- **[numpy](https://numpy.org/):** For efficient numerical operations and matrix manipulations.
- **[netCDF4](https://github.com/Unidata/netcdf4-python):** The underlying I/O library used by SOFA conventions.

### License

This project is licensed under the **MIT License**.

------

## 项目简介

**sofa-analyzer** 是一个 Python 小工具，专为解析和分析 **SOFA (Spatially Oriented Format for Acoustics)** 文件而设计。它主要用于提取、转换和可视化声源测量位置（Source Positions），能够清晰地列出受试者头部前方半球内的 HRTF（头相关传输函数）或 HRIR（头相关脉冲响应）测量点位分布情况。

### 功能特性

- **标准支持：** 完美支持 SOFA 声学数据存储标准。
- **坐标转换：** 自动识别坐标系类型（球坐标或笛卡尔坐标），并将笛卡尔坐标自动转换为球坐标（方位角、俯仰角、距离），以便于统一分析。
- **维度健壮性：** 针对多接收器数据结构（维度：`M, R, C`）进行了增强处理，确保在处理包含双耳或多接收器的复杂 SOFA 文件时程序不会崩溃。
- **数值稳定性：** 优化了三角函数计算逻辑，有效防止因浮点数精度问题导致的计算错误（如 `arcsin` 域错误）。
- **智能筛选：** 自动筛选并展示人头前方 90 度范围内的有效测量点。
- **格式化输出：** 按俯仰角分组，并对方位角进行智能排序（0° 居首，左右交替排列），输出结果清晰直观。

### 安装与配置

请确保您的环境中已安装 Python 3.6 或更高版本。在使用本脚本之前，需要安装以下依赖库。

1. **下载** 本项目代码。
2. 使用 `pip` **安装依赖**：

bash

```bash
pip install pysofaconventions numpy netCDF4
```

### 使用方法

1. 打开终端（Terminal）或命令提示符（CMD）。
2. 进入 `sofa-analyzer.py` 所在的目录。
3. 运行脚本：

bash

```bash
python sofa-analyzer.py
```

1. 按照屏幕提示操作：
   - 输入 `.sofa` 文件的绝对路径（例如：`C:\Users\Admin\Desktop\kemar.sofa`）。
   - 程序将自动分析文件并打印点位分布报告。

#### 使用示例

text

```text
============================================================
       sofa-analyzer v1.1.0 - SOFA文件点位分析工具
============================================================

请输入您的sofa文件的储存路径，然后按回车：
（例如：C:\Users\Tunwen\Desktop\kemar1.sofa）

>>> D:\AudioData\HRTF_Subject_001.sofa

正在加载和分析SOFA文件...

============================================================

导入成功。该文件共记载了 1200 个测量点位。

其中，人头前方90度范围内共有 850 个点位
（分布在 15 个俯仰角上）：

  一、俯仰角-40°：左右方向0°，+15°，-15°，+30°，-30°
  二、俯仰角0°：左右方向0°，+5°，-5°，+10°，-10°，+15°，-15°
  ...
  
  其余未列出的点位需要经算法换算得到。

============================================================
```

### 依赖库

本项目基于以下开源库构建：

- **[pysofaconventions](https://github.com/andresperezlopez/pysofaconventions):** 用于遵循 AES69 标准读取和解析 SOFA 文件。
- **[numpy](https://numpy.org/):** 用于高效的数值计算和矩阵操作。
- **[netCDF4](https://github.com/Unidata/netcdf4-python):** SOFA 格式底层的 I/O 依赖库。

### 开源协议


本项目采用 **MIT 许可证** (MIT License)。
