#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import argparse
from collections import defaultdict
import sys
import os

class EpitopeAnalyzer:
    def __init__(self):
        self.software_predictions = {}
        self.position_counts = defaultdict(int)
        self.shared_regions = []
        self.protein_sequence = ""
        self.predefined_peptides = []
        
    def set_protein_sequence(self, sequence):
        """设置蛋白质序列，用于提取肽段序列"""
        self.protein_sequence = sequence.strip()
        
    def get_peptide_sequence(self, start, end):
        """根据起始和结束位置提取肽段序列"""
        if not self.protein_sequence:
            return "序列未提供"
        
        # 将位置转为0-based索引
        start_idx = start - 1  # 转为0-based索引
        end_idx = end  # end位置是闭区间，需要转为开区间
        
        if start_idx < 0 or end_idx > len(self.protein_sequence):
            return "位置超出序列范围"
            
        return self.protein_sequence[start_idx:end_idx]
    
    def add_predefined_peptide(self, name, sequence, start, end):
        """添加一个预定义的肽段"""
        length = end - start + 1
        self.predefined_peptides.append({
            'name': name,
            'sequence': sequence,
            'start': start,
            'end': end,
            'length': length
        })
        print(f"已添加预定义肽段: {name} ({start}-{end})")
    
    def load_peptides_from_file(self, file_path):
        """从文件中加载肽段数据"""
        if not os.path.exists(file_path):
            print(f"错误: 文件 '{file_path}' 不存在")
            return False
            
        # 清空当前肽段列表
        self.predefined_peptides = []
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                
            count = 0
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):  # 跳过空行和注释行
                    continue
                    
                # 支持CSV格式和空格分隔格式
                if ',' in line:
                    parts = [part.strip() for part in line.split(',')]
                else:
                    parts = [part.strip() for part in line.split()]
                    
                if len(parts) >= 4:
                    name = parts[0]
                    sequence = parts[1]
                    
                    try:
                        start = int(parts[2])
                        end = int(parts[3])
                        self.add_predefined_peptide(name, sequence, start, end)
                        count += 1
                    except ValueError:
                        print(f"警告: 行 '{line}' 包含无效的位置数据，已跳过")
                else:
                    print(f"警告: 行 '{line}' 格式不正确，已跳过")
            
            print(f"从文件 '{file_path}' 成功加载了 {count} 个肽段")
            return count > 0
        except Exception as e:
            print(f"读取肽段文件出错: {e}")
            return False
        
    def create_example_peptides_file(self, file_path="peptides_example.txt"):
        """创建一个肽段数据文件示例"""
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write("# 肽段数据文件示例\n")
                f.write("# 格式: 名称,序列,起始位置,结束位置\n")
                f.write("# 或者: 名称 序列 起始位置 结束位置\n\n")
                f.write("Peptide 1,EQLEQKLRDVEETKAKAE,23,40\n")
                f.write("Peptide 2,TLLQKKYTNLENEFDQVNEK,44,63\n")
                f.write("Peptide 3,LEVSEKRVTEAED,71,83\n")
                f.write("Peptide 4,RRIQLLLEDDLERT,90,102\n")
                f.write("Peptide 5,DATKTADESERGRKYLES,115,132\n")
                f.write("Peptide 6,EKQVKDAKYVAEEADRKYDE,145,164\n")
                f.write("Peptide 7,EVDLERSETRLEA,173,185\n")
                f.write("Peptide 8,LQNAVDQASQREDSYEE,207,223\n")
                f.write("Peptide 9,KDAENRAAEAERVVNKLQ,233,250\n")
                f.write("Peptide 10,ELLAEKEKYKAISDELDQ,259,276\n")
            
            print(f"已创建肽段数据文件示例: {file_path}")
            return True
        except Exception as e:
            print(f"创建肽段数据文件示例出错: {e}")
            return False
        
    def parse_range(self, range_str):
        """解析类似'10~40'或'10-40'的范围格式为起始和结束位置"""
        match = re.search(r'(\d+)[~\-](\d+)', range_str)
        if match:
            start, end = int(match.group(1)), int(match.group(2))
            return (start, end)
        return None
    
    def add_software_prediction(self, software_name, ranges_text):
        """添加一个软件的预测结果"""
        ranges = []
        # 支持多种分隔符：逗号、分号、空格
        for part in re.split(r'[,;\s]+', ranges_text):
            if not part.strip():
                continue
            range_tuple = self.parse_range(part)
            if range_tuple:
                ranges.append(range_tuple)
        
        self.software_predictions[software_name] = ranges
        print(f"已添加 {software_name} 的 {len(ranges)} 个预测区间")
        
        # 更新每个位置的计数
        for start, end in ranges:
            for pos in range(start, end + 1):
                self.position_counts[pos] += 1
    
    def find_shared_regions(self, min_software_count=2):
        """查找至少被指定数量软件预测的连续区域"""
        positions = sorted(self.position_counts.keys())
        if not positions:
            print("未找到任何预测区间，请至少添加一个软件的预测结果")
            return []
        
        shared_regions = []
        start = None
        
        for i in range(min(positions), max(positions) + 2):  # +2 to handle the last region
            if self.position_counts[i] >= min_software_count:
                if start is None:
                    start = i
            elif start is not None:
                end = i - 1
                shared_regions.append((start, end))
                start = None
        
        self.shared_regions = shared_regions
        
        if not shared_regions:
            print(f"未找到至少被{min_software_count}个软件共同预测的区域")
            
        return shared_regions
    
    def print_results(self, use_predefined=False):
        """打印分析结果（按照表4的格式）"""
        if use_predefined and self.predefined_peptides:
            peptides = self.predefined_peptides
        elif not use_predefined and self.shared_regions:
            # 将共享区域转换为肽段格式
            peptides = []
            for i, (start, end) in enumerate(self.shared_regions, 1):
                peptide_name = f"Peptide {i}"
                peptide_sequence = self.get_peptide_sequence(start, end)
                length = end - start + 1
                peptides.append({
                    'name': peptide_name,
                    'sequence': peptide_sequence,
                    'start': start,
                    'end': end,
                    'length': length
                })
        else:
            print("没有可显示的肽段数据")
            return
        
        print("\n表4 联合生物信息学方法预测的线性表位")
        print("Table 4 Linear epitopes predicted by combined bioinformatics methods")
        print("=" * 80)
        print(f"{'名称':<12} {'序列':<40} {'位置':<15} {'长度':<8}")
        print("-" * 80)
        
        for peptide in peptides:
            location = f"{peptide['start']}~{peptide['end']}"
            print(f"{peptide['name']:<12} {peptide['sequence']:<40} {location:<15} {peptide['length']:<8}")
    
    def save_results(self, filename="linear_epitopes.txt", use_predefined=False):
        """保存结果到文件（按照表4的格式）"""
        if use_predefined and self.predefined_peptides:
            peptides = self.predefined_peptides
        elif not use_predefined and self.shared_regions:
            # 将共享区域转换为肽段格式
            peptides = []
            for i, (start, end) in enumerate(self.shared_regions, 1):
                peptide_name = f"Peptide {i}"
                peptide_sequence = self.get_peptide_sequence(start, end)
                length = end - start + 1
                peptides.append({
                    'name': peptide_name,
                    'sequence': peptide_sequence,
                    'start': start,
                    'end': end,
                    'length': length
                })
        else:
            print("没有可保存的肽段数据，不生成结果文件")
            return
        
        with open(filename, 'w', encoding='utf-8') as f:
            f.write("表4 联合生物信息学方法预测的线性表位\n")
            f.write("Table 4 Linear epitopes predicted by combined bioinformatics methods\n")
            f.write("=" * 80 + "\n")
            f.write(f"{'名称':<12} {'序列':<40} {'位置':<15} {'长度':<8}\n")
            f.write("-" * 80 + "\n")
            
            for peptide in peptides:
                location = f"{peptide['start']}~{peptide['end']}"
                f.write(f"{peptide['name']:<12} {peptide['sequence']:<40} {location:<15} {peptide['length']:<8}\n")
            
            if not use_predefined:
                f.write("\n\n软件预测详情:\n")
                for software, ranges in self.software_predictions.items():
                    f.write(f"{software}: ")
                    range_strs = [f"{start}-{end}" for start, end in ranges]
                    f.write(", ".join(range_strs))
                    f.write("\n")
        
        print(f"结果已保存到 {filename}")
        
    def analyze_overlapping_regions(self):
        """分析不同软件预测区间的重叠情况"""
        if len(self.software_predictions) < 2:
            print("需要至少两个软件的预测数据才能分析重叠情况")
            return
            
        print("\n不同软件预测区间的重叠情况:")
        print("=" * 60)
        
        software_names = list(self.software_predictions.keys())
        
        for i in range(len(software_names)):
            for j in range(i+1, len(software_names)):
                name1 = software_names[i]
                name2 = software_names[j]
                
                print(f"{name1} 和 {name2} 的重叠区域:")
                
                has_overlap = False
                for r1_start, r1_end in self.software_predictions[name1]:
                    for r2_start, r2_end in self.software_predictions[name2]:
                        # 检查两个区间是否重叠
                        if max(r1_start, r2_start) <= min(r1_end, r2_end):
                            overlap_start = max(r1_start, r2_start)
                            overlap_end = min(r1_end, r2_end)
                            print(f"  {overlap_start}-{overlap_end} (长度: {overlap_end - overlap_start + 1})")
                            has_overlap = True
                
                if not has_overlap:
                    print("  没有重叠区域")
                print("-" * 40)

def print_usage():
    """打印使用说明"""
    print("\n表位预测分析工具使用说明")
    print("=" * 60)
    print("此工具需要通过参数或交互模式提供数据，以下是常用选项：")
    print("\n1. 交互模式（推荐新用户使用）:")
    print("   python epitope_analyzer.py --interactive")
    print("\n2. 使用肽段数据文件:")
    print("   python epitope_analyzer.py --peptides-file peptides.txt")
    print("\n3. 创建肽段数据文件示例:")
    print("   python epitope_analyzer.py --create-example")
    print("\n4. 结合蛋白质序列和肽段数据文件:")
    print("   python epitope_analyzer.py --peptides-file peptides.txt --sequence \"PROTEIN_SEQUENCE...\"")
    print("\n5. 分析软件预测区间的重叠情况:")
    print("   python epitope_analyzer.py --interactive --analyze-overlap")
    print("\n详细使用说明请参考README.md文件")
    print("=" * 60)

def main():
    parser = argparse.ArgumentParser(description='分析多个软件的表位预测结果，找出共同预测的区域')
    parser.add_argument('--interactive', action='store_true', help='使用交互模式输入数据')
    parser.add_argument('--analyze-overlap', action='store_true', help='分析不同软件预测区间的重叠情况')
    parser.add_argument('--sequence', type=str, help='蛋白质序列，用于提取肽段序列')
    parser.add_argument('--sequence-file', type=str, help='包含蛋白质序列的文件')
    parser.add_argument('--peptides-file', type=str, help='包含预定义肽段数据的文件')
    parser.add_argument('--create-example', action='store_true', help='创建肽段数据文件示例')
    args = parser.parse_args()
    
    # 检查是否提供了任何参数
    if len(sys.argv) == 1:
        print_usage()
        return
    
    analyzer = EpitopeAnalyzer()
    
    # 创建肽段数据文件示例
    if args.create_example:
        analyzer.create_example_peptides_file()
        return
    
    # 处理蛋白质序列
    if args.sequence:
        analyzer.set_protein_sequence(args.sequence)
    elif args.sequence_file:
        try:
            with open(args.sequence_file, 'r') as f:
                sequence = f.read()
                analyzer.set_protein_sequence(sequence)
        except Exception as e:
            print(f"读取序列文件出错: {e}")
    
    # 使用肽段数据文件
    if args.peptides_file:
        if analyzer.load_peptides_from_file(args.peptides_file):
            analyzer.print_results(use_predefined=True)
            analyzer.save_results(use_predefined=True)
            return
    
    if args.interactive:
        print("=== 表位预测分析工具 ===")
        
        # 询问用户是否使用预定义肽段
        use_predefined = input("是否使用预定义肽段数据? (y/n): ")
        if use_predefined.lower() == 'y':
            peptides_source = input("从文件加载(1)还是创建示例文件(2)? (1/2): ")
            
            if peptides_source == '1':
                file_path = input("请输入肽段数据文件路径: ")
                if analyzer.load_peptides_from_file(file_path):
                    # 询问用户是否提供蛋白质序列
                    if not analyzer.protein_sequence:
                        seq_answer = input("是否提供蛋白质序列? (y/n): ")
                        if seq_answer.lower() == 'y':
                            seq_input_type = input("通过键盘输入(1)还是文件读取(2)? (1/2): ")
                            if seq_input_type == '1':
                                sequence = input("请输入蛋白质序列: ")
                                analyzer.set_protein_sequence(sequence)
                            elif seq_input_type == '2':
                                file_path = input("请输入序列文件路径: ")
                                try:
                                    with open(file_path, 'r') as f:
                                        sequence = f.read()
                                        analyzer.set_protein_sequence(sequence)
                                except Exception as e:
                                    print(f"读取序列文件出错: {e}")
                    
                    analyzer.print_results(use_predefined=True)
                    analyzer.save_results(use_predefined=True)
                    return
            elif peptides_source == '2':
                example_file = input("请输入示例文件名 (默认: peptides_example.txt): ") or "peptides_example.txt"
                analyzer.create_example_peptides_file(example_file)
                print(f"请编辑 {example_file} 文件后再次运行程序")
                return
        
        # 询问用户是否提供蛋白质序列
        if not analyzer.protein_sequence:
            seq_answer = input("是否提供蛋白质序列? (y/n): ")
            if seq_answer.lower() == 'y':
                seq_input_type = input("通过键盘输入(1)还是文件读取(2)? (1/2): ")
                if seq_input_type == '1':
                    sequence = input("请输入蛋白质序列: ")
                    analyzer.set_protein_sequence(sequence)
                elif seq_input_type == '2':
                    file_path = input("请输入序列文件路径: ")
                    try:
                        with open(file_path, 'r') as f:
                            sequence = f.read()
                            analyzer.set_protein_sequence(sequence)
                    except Exception as e:
                        print(f"读取序列文件出错: {e}")
        
        print("\n请输入各软件的预测区间，格式为: 起始位置-结束位置,起始位置-结束位置...")
        print("每个软件输入完成后按回车，输入'done'结束输入")
        
        while True:
            software_name = input("\n请输入软件名称 (输入'done'结束): ")
            if software_name.lower() == 'done':
                break
                
            ranges_text = input(f"请输入 {software_name} 的预测区间: ")
            analyzer.add_software_prediction(software_name, ranges_text)
    else:
        # 不再使用默认数据，必须通过参数或交互模式提供输入
        if not args.peptides_file and not args.interactive and len(analyzer.software_predictions) == 0:
            print("错误: 未提供任何软件预测数据。请使用--interactive参数进入交互模式，或指定--peptides-file参数。")
            print_usage()
            return
    
    if args.analyze_overlap:
        analyzer.analyze_overlapping_regions()
    
    # 仅当有软件预测数据时才执行以下操作
    if analyzer.software_predictions:
        shared_regions = analyzer.find_shared_regions(min_software_count=3)
        
        if shared_regions:
            analyzer.print_results(use_predefined=False)
            analyzer.save_results(use_predefined=False)

if __name__ == "__main__":
    main() 