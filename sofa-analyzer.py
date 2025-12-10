#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sofa-analyzer
分析SOFA文件并列出其中直接记载的测量点位

功能：读取SOFA文件，显示人头前方90度范围内的所有测量点位

依赖库：
    - pysofaconventions
    - numpy
    - netCDF4
    
版本：1.1.0
作者：tunwen0
"""

import os
import sys
import warnings
import numpy as np
from collections import defaultdict

# 抑制警告
warnings.filterwarnings('ignore')

# 程序版本
VERSION = "1.1.0"

# ============================================================
# 库导入检查
# ============================================================
try:
    from pysofaconventions import SOFAFile
except ImportError:
    print("=" * 60)
    print("错误：未找到 pysofaconventions 库")
    print("请使用以下命令安装：")
    print("    pip install pysofaconventions")
    print("=" * 60)
    sys.exit(1)

try:
    import netCDF4
except ImportError:
    print("=" * 60)
    print("错误：未找到 netCDF4 库")
    print("请使用以下命令安装：")
    print("    pip install netCDF4")
    print("=" * 60)
    sys.exit(1)


def clear_screen():
    """清屏函数"""
    os.system('cls' if os.name == 'nt' else 'clear')


def print_header():
    """打印程序头部"""
    print("=" * 60)
    print(f"       sofa-analyzer v{VERSION} - SOFA文件点位分析工具")
    print("=" * 60)
    print()


def load_sofa_file(filepath):
    """
    加载SOFA文件
    
    参数：
        filepath: SOFA文件路径
    
    返回：
        SOFAFile对象 或 None（如果失败）
    """
    try:
        sofa = SOFAFile(filepath, 'r')
        return sofa
    except Exception as e:
        print(f"错误：无法打开SOFA文件")
        print(f"详细信息：{e}")
        return None


def get_source_positions(sofa):
    """
    获取所有测量点的位置
    
    返回：
        positions: numpy数组 (M, 3)，[azimuth, elevation, distance]
        coord_type: 坐标系类型
    """
    # 获取坐标系信息
    try:
        pos_units, coord_type = sofa.getSourcePositionInfo()
        coord_type = coord_type if coord_type else 'spherical'
    except:
        coord_type = 'spherical'
    
    # 获取位置数据
    positions = sofa.getSourcePositionValues()
    
    # ============================================================
    # 维度检查与处理 (改进部分)
    # ============================================================
    # SOFA文件可能包含多个接收器 (M, R, C)，例如 (测量点数, 2耳, 3坐标)
    # 如果是这种情况，我们通常取第一个接收器(R=0)的数据作为声源位置参考
    if positions.ndim == 3:
        positions = positions[:, 0, :]
    
    # 确保是二维数组 (M, 3)
    if positions.ndim == 1:
        positions = positions.reshape(1, -1)
    
    # ============================================================
    # 坐标转换逻辑
    # ============================================================
    # 如果是笛卡尔坐标，转换为球面坐标
    if 'cartesian' in coord_type.lower():
        spherical_positions = []
        for pos in positions:
            x, y, z = pos[0], pos[1], pos[2]
            distance = np.sqrt(x**2 + y**2 + z**2)
            
            if distance < 1e-10:
                # 距离极小，视为原点
                spherical_positions.append([0.0, 0.0, 0.0])
            else:
                # 使用 np.clip 确保数值在 [-1, 1] 之间，防止浮点误差导致 arcsin 报错
                # (这也是之前提到的数值稳定性优化)
                sin_val = np.clip(z / distance, -1.0, 1.0)
                elevation = np.degrees(np.arcsin(sin_val))
                azimuth = np.degrees(np.arctan2(y, x))
                spherical_positions.append([azimuth, elevation, distance])
        positions = np.array(spherical_positions)
    
    return positions, coord_type


def normalize_azimuth(azimuth):
    """
    将方位角归一化到 [-180, 180] 范围
    """
    while azimuth > 180:
        azimuth -= 360
    while azimuth < -180:
        azimuth += 360
    return azimuth


def format_azimuth(azimuth):
    """
    格式化方位角显示
    SOFA约定：正值向左，负值向右
    """
    azimuth = round(azimuth, 1)
    # 去除不必要的小数位
    if azimuth == int(azimuth):
        azimuth = int(azimuth)
    
    if azimuth == 0:
        return "0°"
    elif azimuth > 0:
        return f"+{azimuth}°"
    else:
        return f"{azimuth}°"


def format_elevation(elevation):
    """
    格式化俯仰角显示
    """
    elevation = round(elevation, 1)
    if elevation == int(elevation):
        elevation = int(elevation)
    
    if elevation == 0:
        return "0°"
    elif elevation > 0:
        return f"+{elevation}°"
    else:
        return f"{elevation}°"


def analyze_front_positions(positions):
    """
    分析人头前方90度范围内的点位
    
    前方90度范围定义：
    - 方位角绝对值 <= 90°，即 -90° <= azimuth <= +90°
    
    返回：
        grouped_positions: 按俯仰角分组的字典
            key: 俯仰角
            value: 方位角集合
    """
    grouped = defaultdict(set)
    
    for pos in positions:
        azimuth = normalize_azimuth(pos[0])
        elevation = pos[1]
        
        # 筛选前方90度范围（方位角绝对值 <= 90°）
        if -90 <= azimuth <= 90:
            # 四舍五入到1位小数
            az_rounded = round(azimuth, 1)
            el_rounded = round(elevation, 1)
            grouped[el_rounded].add(az_rounded)
    
    return grouped


def sort_azimuths(azimuths):
    """
    对方位角列表进行排序
    排序规则：0°排最前面，然后正负交替排列
    例如：0°, +10°, -10°, +20°, -20°, ...
    """
    azimuths = list(azimuths)
    
    # 分离0°、正值和负值
    zeros = [a for a in azimuths if a == 0]
    positives = sorted([a for a in azimuths if a > 0])
    negatives = sorted([a for a in azimuths if a < 0], key=lambda x: abs(x))
    
    # 交替排列
    result = zeros.copy()
    
    max_len = max(len(positives), len(negatives)) if positives or negatives else 0
    for i in range(max_len):
        if i < len(positives):
            result.append(positives[i])
        if i < len(negatives):
            result.append(negatives[i])
    
    return result


def sort_elevations(elevations):
    """
    对俯仰角列表进行排序
    排序规则：0°排最前面，然后正负交替排列
    """
    elevations = list(elevations)
    
    # 分离0°、正值和负值
    zeros = [e for e in elevations if e == 0]
    positives = sorted([e for e in elevations if e > 0])
    negatives = sorted([e for e in elevations if e < 0], key=lambda x: abs(x))
    
    # 交替排列
    result = zeros.copy()
    
    max_len = max(len(positives), len(negatives)) if positives or negatives else 0
    for i in range(max_len):
        if i < len(positives):
            result.append(positives[i])
        if i < len(negatives):
            result.append(negatives[i])
    
    return result


def get_chinese_number(n):
    """
    获取中文序号（1-99）
    """
    chinese_digits = ["零", "一", "二", "三", "四", "五", "六", "七", "八", "九"]
    
    if n <= 0:
        return str(n)
    elif n <= 10:
        return ["一", "二", "三", "四", "五", "六", "七", "八", "九", "十"][n - 1]
    elif n < 20:
        return "十" + (chinese_digits[n - 10] if n > 10 else "")
    elif n < 100:
        tens = n // 10
        ones = n % 10
        result = chinese_digits[tens] + "十"
        if ones > 0:
            result += chinese_digits[ones]
        return result
    else:
        return str(n)


def display_positions(grouped_positions):
    """
    显示分组后的点位信息
    """
    if not grouped_positions:
        print("  该SOFA文件在前方90度范围内没有记载任何点位。")
        return
    
    # 获取排序后的俯仰角列表
    elevations = sort_elevations(list(grouped_positions.keys()))
    
    print()
    for i, elevation in enumerate(elevations):
        azimuths = grouped_positions[elevation]
        sorted_azimuths = sort_azimuths(azimuths)
        
        # 格式化方位角列表
        az_strs = [format_azimuth(az) for az in sorted_azimuths]
        az_display = "，".join(az_strs)
        
        # 获取中文序号
        num_str = get_chinese_number(i + 1)
        
        print(f"  {num_str}、俯仰角{format_elevation(elevation)}：左右方向{az_display}")
    
    print()
    print("  其余未列出的点位需要经算法换算得到。")


def main():
    """主程序入口"""
    
    while True:
        # 清屏并打印头部
        clear_screen()
        print_header()
        
        # ============================================================
        # 第一步：输入SOFA文件路径
        # ============================================================
        print("请输入您的sofa文件的储存路径，然后按回车：")
        print("（例如：C:\\Users\\Tunwen\\Desktop\\kemar1.sofa）")
        print()
        
        sofa_path = input(">>> ").strip()
        
        # 去除可能的引号
        sofa_path = sofa_path.strip('"').strip("'")
        
        # 检查输入是否为空
        if not sofa_path:
            print()
            print("错误：请输入有效的文件路径")
            print()
            input("按回车键继续...")
            continue
        
        # 检查文件是否存在
        if not os.path.exists(sofa_path):
            print()
            print(f"错误：文件不存在")
            print(f"路径：{sofa_path}")
            print()
            input("按回车键继续...")
            continue
        
        # ============================================================
        # 第二步：加载和分析SOFA文件
        # ============================================================
        print()
        print("正在加载和分析SOFA文件...")
        
        sofa = load_sofa_file(sofa_path)
        if sofa is None:
            print()
            input("按回车键继续...")
            continue
        
        try:
            # 获取测量点位置
            positions, coord_type = get_source_positions(sofa)
            
            # 获取总测量点数
            total_points = len(positions)
            
            # 分析前方90度范围内的点位
            grouped_positions = analyze_front_positions(positions)
            
            # 计算前方点位数
            front_points = sum(len(azs) for azs in grouped_positions.values())
            
            # 计算俯仰角数量
            elevation_count = len(grouped_positions)
            
            # 关闭SOFA文件
            try:
                sofa.close()
            except:
                pass
            
            # 显示结果
            print()
            print("=" * 60)
            print()
            print(f"导入成功。该文件共记载了 {total_points} 个测量点位。")
            print()
            print(f"其中，人头前方90度范围内共有 {front_points} 个点位")
            print(f"（分布在 {elevation_count} 个俯仰角上）：")
            
            display_positions(grouped_positions)
            
            print()
            print("=" * 60)
            
        except Exception as e:
            print(f"错误：分析SOFA文件时出错 - {e}")
            print()
            # 打印详细错误堆栈，方便调试
            import traceback
            traceback.print_exc()
            
            try:
                sofa.close()
            except:
                pass
            input("按回车键继续...")
            continue
        
        # ============================================================
        # 第三步：询问是否继续
        # ============================================================
        print()
        print("分析完毕，继续分析其他文件请按1，退出程序请按2：")
        print()
        
        try:
            choice = input(">>> ").strip()
            
            if choice == "2":
                print()
                print("感谢使用 sofa-analyzer，再见！")
                print()
                break
            # 其他任何输入都继续循环
            
        except (KeyboardInterrupt, EOFError):
            print()
            print("程序已退出。")
            break


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n程序被用户中断，已退出。")
        sys.exit(0)
    except Exception as e:
        print(f"\n程序发生未预期的错误：{e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)